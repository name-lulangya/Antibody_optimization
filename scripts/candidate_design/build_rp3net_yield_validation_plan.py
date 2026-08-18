#!/usr/bin/env python3
"""Build the fixed 47-sequence RP3Net–BL21 yield validation plan."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.rp3net_yield import build_rp3net_validation_inputs  # noqa: E402


NAMES = {
    "samples": "rp3net_validation_samples.csv",
    "fasta": "rp3net_validation_sequences.fasta",
    "contract": "rp3net_yield_validation_contract.json",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression-records", type=Path, required=True)
    parser.add_argument("--numbering-review", type=Path, required=True)
    parser.add_argument("--numbering-positions", type=Path, required=True)
    parser.add_argument("--allowed-use-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--check_only", action="store_true")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    sources = [
        args.expression_records.resolve(strict=True),
        args.numbering_review.resolve(strict=True),
        args.numbering_positions.resolve(strict=True),
        args.allowed_use_manifest.resolve(strict=True),
    ]
    allowed = _json(sources[3])
    if allowed["gates"]["cross_assay_pooling_gate"] != "pass" or allowed["counts"]["samples"] != 47:
        raise ValueError("Expression audit does not release the fixed 47-sample exploratory analysis")
    samples = build_rp3net_validation_inputs(_csv(sources[0]), _csv(sources[1]), _csv(sources[2]))["sample_rows"]
    if args.check_only:
        print(json.dumps({"status": "pass", "samples": len(samples)}, ensure_ascii=False))
        return 0

    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite RP3Net validation plan")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    contract = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated,
        "expression_system": "BL21_E_coli",
        "software": {"name": "RP3Net", "version": "0.0.2", "model": "RP3Esm2_650m"},
        "remote_environment": "/data/software/env/luly25/rp3net",
        "entry_point": "/data/software/env/luly25/rp3net/bin/rp3",
        "checkpoint": {
            "path": "/homes/Tianlab/luly25/software/RP3Net/checkpoints/rp3net_v0.1_d.ckpt",
            "sha256": "443743bd031689aaf17dc6f7c22c5da3d23cf87b38e10341f114b27d651e6d2b",
        },
        "primary_score": "predicted_probability_of_recombinant_small_scale_expression_in_E_coli",
        "score_direction": "higher_is_better",
        "continuous_analysis": "31_LTT_WCC_individual_approximate_plus_16_LLJ_ordinal_censored",
        "classification": {
            "numeric_samples": 31,
            "outcome": "high_yield_at_or_above_matching_provider_median_fitted_in_outer_training_fold",
            "score_threshold": "maximize_training_MCC_then_balanced_accuracy_then_higher_threshold",
            "outer_schemes": ["leave_one_out", "leave_one_cluster_out"],
            "reported_metrics": ["ROC-AUC", "PR-AUC", "MCC", "balanced_accuracy", "sensitivity", "specificity", "threshold_stability"],
        },
        "llj_use": "group_ordinal_or_censored_only_no_individual_point_estimates",
        "release": "ready_for_remote_rp3net_scoring",
    }
    with tempfile.TemporaryDirectory(prefix=".rp3net-plan-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(path).name for key, path in final.items()}
        _write_csv(staged["samples"], samples)
        staged["fasta"].write_text("".join(f'>{row["sample_uid"]}\n{row["sequence_raw"]}\n' for row in samples), encoding="utf-8", newline="\n")
        _write_json(staged["contract"], contract)
        _write_json(staged["summary"], {
            "schema_version": 1, "status": "pass", "generated_at": generated,
            "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter() - started, 6),
            "sample_count": 47, "numeric_count": 31, "llj_ordinal_censored_count": 16,
            "outputs": {key: str(value) for key, value in final.items() if key != "summary"},
        })
        replace_staged_files({staged[key]: final[key] for key in staged}, project_root=ROOT, protected_source_paths=valid.source_paths)
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
