#!/usr/bin/env python3
"""Analyze fixed NanoMelt predictions against collaborator-reported BL21 yield."""

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
from antibody_optimization.nanomelt_yield import analyze_nanomelt_associations  # noqa: E402
from antibody_optimization.nanomelt_yield_plot import render_nanomelt_yield_figure  # noqa: E402


NAMES = {
    "samples": "nanomelt_yield_sample_evidence.csv",
    "metrics": "nanomelt_yield_associations.csv",
    "cv": "nanomelt_yield_cv_comparison.csv",
    "influence": "nanomelt_yield_leave_one_out.csv",
    "gate": "nanomelt_yield_validation_gate.json",
    "png": "nanomelt_yield_validation.png",
    "svg": "nanomelt_yield_validation.svg",
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
        plan / "nanomelt_validation_samples.csv",
        plan / "nanomelt_yield_validation_contract.json",
        scores / "nanomelt_sample_scores.csv",
        scores / "nanomelt_model_run.json",
    ]
    contract, model_run = _json(sources[1]), _json(sources[3])
    if contract.get("schema_version") != 1 or contract.get("status") != "pass":
        raise ValueError("NanoMelt plan is not released")
    if model_run.get("schema_version") != 1 or model_run.get("status") != "pass" or model_run.get("sample_count") != 47:
        raise ValueError("NanoMelt model run is incomplete")
    result = analyze_nanomelt_associations(_csv(sources[0]), _csv(sources[2]))

    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite NanoMelt validation outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    release = {
        "weak_ranking_evidence": "nanomelt_ready_for_weak_candidate_ranking_use",
        "compatibility_filter_only": "nanomelt_stability_compatibility_only",
        "no_supported_use": "nanomelt_not_supported_for_yield_use",
    }[result["evidence_level"]]
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_nanomelt_predicted_tm_bl21_reported_yield_validation",
        "status": "pass",
        "generated_at": generated,
        "coverage": {"planned": 47, "scored": 47, "numeric": 31, "llj_ordinal_censored": 16},
        "primary_feature": "nanomelt_predicted_apparent_tm_c",
        "expected_direction": "higher_predicted_apparent_tm_higher_reported_yield",
        "primary_statistics": result["primary"],
        "evidence_level": result["evidence_level"],
        "decision_reasons": result["decision_reasons"],
        "high_capacity_model_trained": False,
        "experimental_tm_available": False,
        "nb252_expression_prediction_validated": False,
        "release": release,
        "interpretation": "Association of NanoMelt-predicted apparent Tm for tool-aligned VHH domains with collaborator-reported BL21 yield; not measured Tm, causal evidence, or an mg/L predictor.",
    }
    with tempfile.TemporaryDirectory(prefix=".nanomelt-analysis-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(value).name for key, value in final.items()}
        _write_csv(staged["samples"], result["sample_rows"])
        _write_csv(staged["metrics"], result["metric_rows"])
        _write_csv(staged["cv"], result["cv_rows"])
        _write_csv(staged["influence"], result["leave_one_out_rows"])
        _write_json(staged["gate"], gate)
        render_nanomelt_yield_figure(
            result["sample_rows"],
            result["primary"],
            result["cv_rows"],
            result["leave_one_out_rows"],
            png_path=staged["png"],
            svg_path=staged["svg"],
        )
        _write_json(
            staged["summary"],
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated,
                "python": platform.python_version(),
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "evidence_level": result["evidence_level"],
                "coverage": gate["coverage"],
                "outputs": {key: str(value) for key, value in final.items() if key != "summary"},
            },
        )
        replace_staged_files(
            {staged[key]: final[key] for key in staged},
            project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


def _csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
