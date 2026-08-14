#!/usr/bin/env python3
"""Analyze fixed TNP scores against collaborator-reported BL21 yield."""

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
from antibody_optimization.tnp_yield import analyze_tnp_associations  # noqa: E402
from antibody_optimization.tnp_yield_plot import render_tnp_yield_figure  # noqa: E402


NAMES = {"samples": "tnp_yield_sample_evidence.csv", "metrics": "tnp_yield_associations.csv", "cv": "tnp_yield_cv_comparison.csv", "gate": "tnp_yield_validation_gate.json", "png": "tnp_yield_validation.png", "svg": "tnp_yield_validation.svg"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--score-dir", type=Path, required=True)
    parser.add_argument("--netsolp-evidence", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan, scores = args.plan_dir.resolve(strict=True), args.score_dir.resolve(strict=True)
    sources = [plan / "tnp_validation_samples.csv", plan / "tnp_yield_validation_contract.json", scores / "tnp_sample_scores.csv", scores / "tnp_model_run.json", args.netsolp_evidence.resolve(strict=True)]
    if _json(sources[1]).get("status") != "pass" or _json(sources[3]).get("status") != "pass": raise ValueError("TNP plan/model run status mismatch")
    result = analyze_tnp_associations(_csv(sources[0]), _csv(sources[2]), _csv(sources[4]))
    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths): raise FileExistsError("Refusing to overwrite TNP validation outputs")
    for path in valid.target_paths: path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    release = {"weak_ranking_evidence": "ready_for_weak_tnp_psh_yield_ranking_use", "compatibility_filter_only": "tnp_developability_compatibility_only", "no_supported_yield_use": "tnp_not_supported_for_yield_use"}[result["evidence_level"]]
    gate = {"schema_version": 1, "gate_name": "nb252_tnp_bl21_reported_yield_validation", "status": "pass", "generated_at": generated, "coverage": result["coverage"], "primary_feature": "tnp_psh", "expected_direction": "higher_PSH_lower_yield", "evidence_level": result["evidence_level"], "decision_reasons": result["decision_reasons"], "primary_statistics": result["primary"], "high_capacity_model_trained": False, "nb252_expression_prediction_validated": False, "release": release, "interpretation": "Association with collaborator-reported BL21 yield; TNP remains a developability-risk profiler and is not a measured expression or mg/L predictor."}
    with tempfile.TemporaryDirectory(prefix=".tnp-analysis-", dir=ROOT) as temp:
        stage = Path(temp); staged = {key: stage / Path(value).name for key, value in final.items()}
        _write_csv(staged["samples"], result["sample_rows"]); _write_csv(staged["metrics"], result["metric_rows"]); _write_csv(staged["cv"], result["cv_rows"]); _write_json(staged["gate"], gate)
        render_tnp_yield_figure(result["sample_rows"], result["metric_rows"], result["cv_rows"], png_path=staged["png"], svg_path=staged["svg"])
        _write_json(staged["summary"], {"schema_version": 1, "status": "pass", "generated_at": generated, "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter() - started, 6), "evidence_level": result["evidence_level"], "coverage": result["coverage"], "outputs": {key: str(value) for key, value in final.items() if key != "summary"}})
        replace_staged_files({staged[key]: final[key] for key in staged}, project_root=ROOT, protected_source_paths=valid.source_paths)
    return 0


def _csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))
def _json(path): return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
def _write_json(path, value): path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__": raise SystemExit(main())
