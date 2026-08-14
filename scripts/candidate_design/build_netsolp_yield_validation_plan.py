#!/usr/bin/env python3
"""Build the fixed 47-sequence NetSolP–BL21 yield validation plan."""

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
from antibody_optimization.nanobert_yield import BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED, PERMUTATION_REPLICATES  # noqa: E402
from antibody_optimization.netsolp_yield import (  # noqa: E402
    SEQUENCE_CLUSTER_IDENTITY,
    build_netsolp_validation_inputs,
)


NAMES = {
    "samples": "netsolp_validation_samples.csv",
    "fasta": "netsolp_validation_sequences.fasta",
    "contract": "netsolp_yield_validation_contract.json",
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
    result = build_netsolp_validation_inputs(_csv(sources[0]), _csv(sources[1]), _csv(sources[2]))
    samples = result["sample_rows"]
    if args.check_only:
        print(
            json.dumps(
                {
                    "status": "pass",
                    "samples": len(samples),
                    "sequence_clusters_90": len({row["sequence_cluster_90"] for row in samples}),
                },
                ensure_ascii=False,
            )
        )
        return 0

    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite NetSolP validation plan")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    contract = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated,
        "expression_system": "BL21_E_coli",
        "software": "NetSolP-1.0_official_full_distribution",
        "remote_environment": "/data/software/env/luly25/netsolp",
        "remote_working_directory": "/homes/Tianlab/luly25/software/netsolp",
        "remote_entry_point": "predict.py",
        "model_type": "Distilled",
        "prediction_type": "SU",
        "primary_score": "predicted_usability",
        "secondary_score": "predicted_solubility",
        "primary_numeric_analysis": "31_LTT_WCC_individual_approximate",
        "llj_use": "16_group_ordinal_or_censored_observations_only",
        "sequence_cluster_definition": {
            "algorithm": "single_linkage_on_1_minus_normalized_levenshtein_distance",
            "identity_threshold": SEQUENCE_CLUSTER_IDENTITY,
        },
        "bootstrap": {"replicates": BOOTSTRAP_REPLICATES, "seed": BOOTSTRAP_SEED},
        "permutation": {"replicates": PERMUTATION_REPLICATES, "seed": BOOTSTRAP_SEED + 1},
        "high_capacity_model_training": False,
        "evidence_levels": ["weak_ranking_evidence", "compatibility_filter_only", "no_supported_use"],
        "release": "ready_for_remote_netsolp_scoring",
    }
    with tempfile.TemporaryDirectory(prefix=".netsolp-plan-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(path).name for key, path in final.items()}
        _write_csv(staged["samples"], samples)
        _write_fasta(staged["fasta"], samples)
        _write_json(staged["contract"], contract)
        _write_json(
            staged["summary"],
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated,
                "python": platform.python_version(),
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "sample_count": 47,
                "numeric_count": 31,
                "llj_ordinal_count": 16,
                "sequence_cluster_count_90": len({row["sequence_cluster_90"] for row in samples}),
                "outputs": {key: str(value) for key, value in final.items() if key != "summary"},
            },
        )
        replace_staged_files(
            {staged[key]: final[key] for key in staged},
            project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    text = "".join(f'>{row["sample_uid"]}\n{row["sequence_raw"]}\n' for row in rows)
    path.write_text(text, encoding="utf-8", newline="\n")


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
