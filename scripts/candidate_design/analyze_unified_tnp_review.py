#!/usr/bin/env python3
"""Analyze the fixed unified TNP review and render compact artifacts."""

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
from antibody_optimization.unified_tnp_plot import render_unified_tnp_review  # noqa: E402
from antibody_optimization.unified_tnp_review import (  # noqa: E402
    SCORE_COUNT,
    analyze_unified_tnp_scores,
)


NAMES = {
    "evidence": "unified_tnp_candidate_evidence.csv",
    "summary_table": "unified_tnp_review_summary.csv",
    "gate": "unified_tnp_review_gate.json",
    "png": "unified_tnp_review.png",
    "svg": "unified_tnp_review.svg",
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
        plan / "unified_tnp_samples.csv",
        plan / "unified_tnp_review_contract.json",
        plan / "unified_tnp_review_plan_gate.json",
        scores / "tnp_sample_scores.csv",
        scores / "tnp_model_run.json",
    ]
    for path in sources:
        path.resolve(strict=True)
    contract = _json(sources[1])
    plan_gate = _json(sources[2])
    model_run = _json(sources[4])
    if (
        contract.get("release") != "ready_for_remote_single_process_unified_tnp_review"
        or plan_gate.get("status") != "pass"
        or model_run.get("status") != "pass"
        or int(model_run.get("pass_count", 0)) != SCORE_COUNT
    ):
        raise ValueError("Unified TNP plan or model-run gate is not passed")
    evidence, summary_table, gate = analyze_unified_tnp_scores(
        _csv(sources[0]), _csv(sources[3])
    )
    gate["generated_at"] = generated

    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite unified TNP review outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "run_summary"), valid.target_paths, strict=True))
    with tempfile.TemporaryDirectory(prefix=".unified-tnp-analysis-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / path.name for key, path in final.items()}
        _write_csv(staged["evidence"], evidence)
        _write_csv(staged["summary_table"], summary_table)
        _write_json(staged["gate"], gate)
        render_unified_tnp_review(evidence, png_path=staged["png"], svg_path=staged["svg"])
        _write_json(staged["run_summary"], {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated,
            "python": platform.python_version(),
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "planned_count": SCORE_COUNT,
            "candidate_count": len(evidence),
            "source_counts": gate["source_counts"],
            "flag_regression_count": gate["flag_regression_count"],
            "new_red_flag_candidate_count": gate["new_red_flag_candidate_count"],
            "yield_prediction_performed": False,
            "candidate_selection_performed": False,
            "outputs": {key: str(path) for key, path in final.items() if key != "run_summary"},
        })
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


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
