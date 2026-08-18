#!/usr/bin/env python3
"""Analyze RP3Net scores against collaborator-reported BL21 yield."""

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
from antibody_optimization.rp3net_yield import analyze_rp3net_associations  # noqa: E402
from antibody_optimization.rp3net_yield_plot import render_rp3net_yield_figure  # noqa: E402


NAMES = {
    "samples": "rp3net_yield_sample_evidence.csv",
    "metrics": "rp3net_yield_associations.csv",
    "classification": "rp3net_yield_classification.csv",
    "predictions": "rp3net_yield_classification_predictions.csv",
    "gate": "rp3net_yield_validation_gate.json",
    "png": "rp3net_yield_validation.png",
    "svg": "rp3net_yield_validation.svg",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--score-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan = args.plan_dir.resolve(strict=True)
    scores = args.score_dir.resolve(strict=True)
    sources = [
        plan / "rp3net_validation_samples.csv", plan / "rp3net_yield_validation_contract.json",
        scores / "rp3net_sample_scores.csv", scores / "rp3net_model_run.json",
    ]
    contract, model_run = _json(sources[1]), _json(sources[3])
    if contract.get("status") != "pass" or model_run.get("status") != "pass":
        raise ValueError("RP3Net plan/model-run status mismatch")
    if model_run.get("checkpoint_sha256") != contract["checkpoint"]["sha256"]:
        raise ValueError("RP3Net plan/model-run checkpoint mismatch")
    result = analyze_rp3net_associations(_csv(sources[0]), _csv(sources[2]))
    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite RP3Net validation outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    release = {
        "weak_ranking_evidence": "ready_for_weak_rp3net_ranking_use",
        "compatibility_filter_only": "rp3net_compatibility_filter_only",
        "no_supported_use": "rp3net_not_supported_for_candidate_use",
    }[result["evidence_level"]]
    gate = {
        "schema_version": 1, "status": "pass", "generated_at": generated,
        "gate_name": "nb252_rp3net_bl21_reported_yield_validation",
        "sample_count": 47, "numeric_individual_count": 31, "llj_ordinal_censored_count": 16,
        "primary_feature": "rp3net_expression_probability", "evidence_level": result["evidence_level"],
        "decision_reasons": result["decision_reasons"], "continuous_statistics": result["primary"],
        "classification_statistics": result["classification_rows"], "high_capacity_model_trained": False,
        "nb252_mutant_expression_prediction_validated": False, "release": release,
        "interpretation": "Model association with reported BL21 yield; not measured yield or a universal expression cutoff.",
    }
    with tempfile.TemporaryDirectory(prefix=".rp3net-analysis-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(path).name for key, path in final.items()}
        _write_csv(staged["samples"], result["sample_rows"])
        _write_csv(staged["metrics"], result["metric_rows"])
        _write_csv(staged["classification"], result["classification_rows"])
        _write_csv(staged["predictions"], result["classification_prediction_rows"])
        _write_json(staged["gate"], gate)
        render_rp3net_yield_figure(result["sample_rows"], result["metric_rows"], result["classification_rows"], result["classification_prediction_rows"], png_path=staged["png"], svg_path=staged["svg"])
        _write_json(staged["summary"], {
            "schema_version": 1, "status": "pass", "generated_at": generated,
            "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter() - started, 6),
            "evidence_level": result["evidence_level"], "counts": {"samples": 47, "numeric": 31, "llj_ordinal_censored": 16},
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
