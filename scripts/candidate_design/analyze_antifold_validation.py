#!/usr/bin/env python3
"""Analyze the three-view AntiFold validation for the released affinity core."""

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

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.antifold_plot import render_antifold_validation  # noqa: E402
from antibody_optimization.antifold_validation import (  # noqa: E402
    build_candidate_evidence,
    normalize_antifold_rows,
    validate_result_gate,
)
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402


OUTPUTS = {
    "evidence": "antifold_candidate_view_evidence.csv",
    "summary": "antifold_candidate_summary.csv",
    "plot_data": "antifold_validation_plot_data.csv",
    "gate": "antifold_validation_gate.json",
    "png": "antifold_validation.png",
    "svg": "antifold_validation.svg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--score-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan_dir = args.plan_dir.resolve(strict=True)
    score_dir = args.score_dir.resolve(strict=True)
    output_dir = args.output_dir.absolute()
    run_summary = args.run_summary.absolute()
    if output_dir.exists() or run_summary.exists():
        raise FileExistsError("Refusing to overwrite AntiFold validation result outputs")
    plan_gate = _json(plan_dir / "antifold_validation_plan_gate.json")
    model_run = _json(score_dir / "antifold_model_run.json")
    if plan_gate.get("status") != "pass" or model_run.get("status") != "pass":
        raise ValueError("AntiFold plan or remote model run did not pass")
    candidates = _csv(plan_dir / "antifold_candidate_panel.csv")
    views = _csv(plan_dir / "antifold_structure_views.csv")
    indexed = {
        view["view_id"]: normalize_antifold_rows(
            _csv(score_dir / f"{view['view_id']}.csv"),
            view_id=view["view_id"],
            vhh_chain=view["vhh_chain"],
        )
        for view in views
    }
    evidence, summaries = build_candidate_evidence(candidates, indexed)
    gate_facts = validate_result_gate(evidence, summaries)
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_antifold_minimal_validation",
        "generated_at": generated_at,
        **gate_facts,
        "all_view_direction_concordant_count": sum(_bool(row["all_view_directions_concordant"]) for row in summaries),
        "experimental_context_direction_change_count": sum(_bool(row["experimental_context_direction_change"]) for row in summaries),
        "experimental_vs_af3_direction_change_count": sum(_bool(row["experimental_vs_af3_direction_change"]) for row in summaries),
        "interpretation": "AntiFold structure-conditioned compatibility only; no affinity, stability, expression, yield, or experimental-effect threshold was applied.",
    }
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".antifold-analysis-", dir=PROJECT_ROOT) as temp_name:
        staging = Path(temp_name)
        staged = {key: staging / name for key, name in OUTPUTS.items()}
        _write_csv(staged["evidence"], evidence)
        _write_csv(staged["summary"], summaries)
        _write_csv(staged["plot_data"], summaries)
        _write_json(staged["gate"], gate)
        render_antifold_validation(summaries, png_path=staged["png"], svg_path=staged["svg"])
        summary_stage = staging / "run_summary.json"
        _write_json(summary_stage, {
            "schema_version": 1, "status": gate["status"], "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6), "python": platform.python_version(),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
            "counts": {
                "candidates": len(summaries), "views": len(views), "evidence_rows": len(evidence),
                "all_view_direction_concordant": gate["all_view_direction_concordant_count"],
            },
            "release": gate["release"],
            "outputs": {key: str(output_dir / name) for key, name in OUTPUTS.items()},
        })
        final_pairs = {staged[key]: output_dir / name for key, name in OUTPUTS.items()}
        final_pairs[summary_stage] = run_summary
        sources = [
            plan_dir / "antifold_validation_plan.json", plan_dir / "antifold_candidate_panel.csv",
            plan_dir / "antifold_structure_views.csv", score_dir / "antifold_model_run.json",
            *[score_dir / f"{view['view_id']}.csv" for view in views],
        ]
        validated = validate_file_paths(
            project_root=PROJECT_ROOT, source_paths=sources, target_paths=list(final_pairs.values())
        )
        for path in validated.target_paths:
            path.parent.mkdir(parents=True, exist_ok=True)
        replace_staged_files(final_pairs, project_root=PROJECT_ROOT, protected_source_paths=validated.source_paths)
    return 0 if gate["status"] == "pass" else 2


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


def _bool(value: object) -> bool:
    return value is True or str(value).lower() == "true"


if __name__ == "__main__":
    raise SystemExit(main())
