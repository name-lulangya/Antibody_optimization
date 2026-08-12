#!/usr/bin/env python3
"""Validate and summarize all 1000 Flex ddG production tasks without filtering."""

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

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.flex_ddg import TASK_METRIC_FIELDS  # noqa: E402
from antibody_optimization.flex_ddg_production import summarize_production_results  # noqa: E402
from antibody_optimization.flex_ddg_production_plot import render_flex_ddg_production_figure  # noqa: E402


OUTPUT_NAMES = {
    "tasks": "flex_ddg_production_task_metrics.csv",
    "candidates": "flex_ddg_production_candidate_summary.csv",
    "gate": "flex_ddg_production_gate.json",
    "png": "flex_ddg_production_qc.png",
    "svg": "flex_ddg_production_qc.svg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--task-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan_dir = _project_directory(args.plan_dir)
    task_root = _project_directory(args.task_root)
    manifest_path = plan_dir / "flex_ddg_production_manifest.csv"
    plan_path = plan_dir / "flex_ddg_production_plan.json"
    manifest = _load_csv(manifest_path)
    task_paths = [task_root / str(row["task_id"]) / "task_result.json" for row in manifest]
    task_results = [_load_json(path) for path in task_paths]
    summary = summarize_production_results(manifest_rows=manifest, task_results=task_results)
    candidate_fields = list(summary["candidate_rows"][0])
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_flex_ddg_production",
        "status": summary["gate_status"],
        "generated_at": generated_at,
        "task_count": len(summary["task_rows"]),
        "candidate_count": len(summary["candidate_rows"]),
        "samples_per_candidate": 20,
        "scope": "all_tier_1_and_2_plus_two_selected_tier_3",
        "candidate_selection_performed": False,
        "release": "ready_for_ensemble_filter_design",
        "interpretation": "Complete model-specific ensemble evidence; not measured affinity.",
    }
    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    targets = [output_dir / name for name in OUTPUT_NAMES.values()] + [run_summary]
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[manifest_path, plan_path, *task_paths],
        target_paths=targets,
    )
    if any(path.exists() for path in validated.target_paths):
        raise FileExistsError("Refusing to overwrite production summary outputs")
    for path in validated.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*OUTPUT_NAMES, "run_summary"), validated.target_paths, strict=True))
    with tempfile.TemporaryDirectory(prefix=".flex-ddg-production-summary-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / Path(path).name for key, path in final.items()}
        _write_csv(staged["tasks"], summary["task_rows"], TASK_METRIC_FIELDS)
        _write_csv(staged["candidates"], summary["candidate_rows"], candidate_fields)
        _write_json(staged["gate"], gate)
        render_flex_ddg_production_figure(
            summary["candidate_rows"], png_path=staged["png"], svg_path=staged["svg"]
        )
        _write_json(
            staged["run_summary"],
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "task_count": len(summary["task_rows"]),
                "candidate_count": len(summary["candidate_rows"]),
                "candidate_selection_performed": False,
                "outputs": {key: str(path) for key, path in final.items() if key != "run_summary"},
            },
        )
        replace_staged_files(
            {staged[key]: final[key] for key in staged},
            project_root=PROJECT_ROOT,
            protected_source_paths=validated.source_paths,
        )
    return 0


def _project_directory(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    if resolved.is_symlink() or not resolved.is_dir():
        raise ValueError(f"Expected regular project directory: {resolved}")
    return resolved


def _load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_csv(path: Path, rows, fields) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
