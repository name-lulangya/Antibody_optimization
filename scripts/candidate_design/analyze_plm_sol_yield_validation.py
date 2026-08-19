#!/usr/bin/env python3
"""Analyze fixed PLM_Sol scores against collaborator-reported BL21 yield."""

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
from antibody_optimization.plm_sol_yield import analyze_plm_sol_associations  # noqa: E402
from antibody_optimization.plm_sol_yield_plot import (  # noqa: E402
    render_plm_sol_fixed5_figure,
    render_plm_sol_yield_figure,
)


NAMES = {
    "samples": "plm_sol_yield_sample_evidence.csv",
    "metrics": "plm_sol_yield_associations.csv",
    "classification": "plm_sol_yield_classification.csv",
    "predictions": "plm_sol_yield_classification_predictions.csv",
    "comparison": "plm_sol_predictor_comparison.csv",
    "fixed_metrics": "plm_sol_fixed5mg_metrics.csv",
    "fixed_predictions": "plm_sol_fixed5mg_predictions.csv",
    "gate": "plm_sol_yield_validation_gate.json",
    "png": "plm_sol_yield_validation.png",
    "svg": "plm_sol_yield_validation.svg",
    "fixed_png": "plm_sol_fixed5mg_display.png",
    "fixed_svg": "plm_sol_fixed5mg_display.svg",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--score-dir", type=Path, required=True)
    parser.add_argument("--netsolp-samples", type=Path, required=True)
    parser.add_argument("--rp3net-samples", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan = args.plan_dir.resolve(strict=True)
    scores = args.score_dir.resolve(strict=True)
    sources = [
        plan / "plm_sol_validation_samples.csv",
        plan / "plm_sol_yield_validation_contract.json",
        scores / "plm_sol_sample_scores.csv",
        scores / "plm_sol_model_run.json",
        args.netsolp_samples.resolve(strict=True),
        args.rp3net_samples.resolve(strict=True),
    ]
    contract, model_run = _json(sources[1]), _json(sources[3])
    if contract.get("status") != "pass" or model_run.get("status") != "pass":
        raise ValueError("PLM_Sol plan/model-run status mismatch")
    result = analyze_plm_sol_associations(_csv(sources[0]), _csv(sources[2]), _csv(sources[4]), _csv(sources[5]))
    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite PLM_Sol validation outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    release = {
        "weak_ranking_evidence": "ready_for_weak_plm_sol_ranking_use",
        "compatibility_filter_only": "plm_sol_compatibility_filter_only",
        "no_supported_use": "plm_sol_not_supported_for_candidate_use",
    }[result["evidence_level"]]
    gate = {
        "schema_version": 1, "status": "pass", "generated_at": generated,
        "gate_name": "nb252_plm_sol_bl21_reported_yield_validation",
        "sample_count": 47, "numeric_individual_count": 31, "llj_ordinal_censored_count": 16,
        "primary_feature": "plm_sol_solubility_score", "evidence_level": result["evidence_level"],
        "decision_reasons": result["decision_reasons"], "continuous_statistics": result["primary"],
        "classification_statistics": result["classification_rows"],
        "independent_cluster_cv_increment_over_netsolp_s": result["independent_cluster_cv_increment"],
        "fixed_5mg_use": "display_only_not_a_predictor_gate_or_candidate_filter",
        "model_trained_or_finetuned_on_project_data": False,
        "nb252_mutant_expression_prediction_validated": False,
        "release": release,
        "interpretation": "PLM_Sol association with reported BL21 yield; not measured solubility, yield, or a universal cutoff.",
    }
    with tempfile.TemporaryDirectory(prefix=".plm-sol-analysis-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(path).name for key, path in final.items()}
        for key, rows in (
            ("samples", result["sample_rows"]), ("metrics", result["metric_rows"]),
            ("classification", result["classification_rows"]), ("predictions", result["classification_prediction_rows"]),
            ("comparison", result["comparison_rows"]), ("fixed_metrics", result["fixed5_metric_rows"]),
            ("fixed_predictions", result["fixed5_prediction_rows"]),
        ):
            _write_csv(staged[key], rows)
        _write_json(staged["gate"], gate)
        render_plm_sol_yield_figure(
            result["sample_rows"], result["metric_rows"], result["classification_rows"], result["comparison_rows"],
            png_path=staged["png"], svg_path=staged["svg"],
        )
        render_plm_sol_fixed5_figure(
            result["sample_rows"], result["fixed5_metric_rows"], png_path=staged["fixed_png"], svg_path=staged["fixed_svg"],
        )
        _write_json(staged["summary"], {
            "schema_version": 1, "status": "pass", "generated_at": generated,
            "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter() - started, 6),
            "evidence_level": result["evidence_level"], "counts": {"samples": 47, "numeric": 31, "llj_ordinal_censored": 16},
            "outputs": {key: str(value) for key, value in final.items() if key != "summary"},
        })
        replace_staged_files({staged[key]: final[key] for key in staged}, project_root=ROOT, protected_source_paths=valid.source_paths)
    return 0


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
