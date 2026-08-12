#!/usr/bin/env python3
"""Merge eight Flex ddG pilot tasks and project full-scope wall time."""

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

from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.flex_ddg import (  # noqa: E402
    TASK_METRIC_FIELDS,
    summarize_pilot_results,
)
from antibody_optimization.flex_ddg_plot import render_flex_ddg_pilot_figure  # noqa: E402


OUTPUT_NAMES = {
    "tasks": "flex_ddg_pilot_task_metrics.csv",
    "timing": "flex_ddg_pilot_timing_summary.csv",
    "projections": "flex_ddg_scope_projections.csv",
    "gate": "flex_ddg_pilot_gate.json",
    "figure_png": "flex_ddg_pilot_timing.png",
    "figure_svg": "flex_ddg_pilot_timing.svg",
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
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    plan_dir = _project_directory(args.plan_dir)
    task_root = _project_directory(args.task_root)
    manifest_path = plan_dir / "flex_ddg_pilot_manifest.csv"
    plan_path = plan_dir / "flex_ddg_pilot_plan.json"
    plan = _load_json(plan_path)
    if plan.get("status") != "pass" or int(plan.get("task_count", 0)) != 8:
        raise ValueError("Flex ddG pilot plan is not released")
    manifest = _load_csv(manifest_path)
    task_paths = [
        task_root / f"task_{index:02d}" / "task_result.json" for index in range(8)
    ]
    task_results = [_load_json(path) for path in task_paths]
    result = summarize_pilot_results(manifest_rows=manifest, task_results=task_results)
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_flex_ddg_timing_pilot",
        "status": result["gate_status"],
        "generated_at": generated_at,
        "task_count": len(result["task_metric_rows"]),
        "status_counts": result["status_counts"],
        "candidate_selection_performed": False,
        "tier_3_scope_decision_performed": False,
        "scope_decision_release": (
            "ready_for_tier_scope_decision"
            if result["gate_status"] == "pass"
            else "blocked"
        ),
        "interpretation": (
            "Timing and protocol feasibility only; two samples per candidate "
            "are insufficient for affinity ranking."
        ),
    }
    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[manifest_path, plan_path, *task_paths],
        target_paths=[
            *[output_dir / name for name in OUTPUT_NAMES.values()],
            run_summary,
        ],
    )
    final_paths = dict(zip(OUTPUT_NAMES, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*final_paths.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n"
            + "\n".join(map(str, existing))
        )
    for path in [*final_paths.values(), run_summary]:
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".flex-ddg-summary-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = staging / "run_summary.json"
        _write_csv(staged["tasks"], result["task_metric_rows"], TASK_METRIC_FIELDS)
        _write_csv(
            staged["timing"],
            result["timing_summary_rows"],
            list(result["timing_summary_rows"][0]),
        )
        _write_csv(
            staged["projections"],
            result["projection_rows"],
            list(result["projection_rows"][0]),
        )
        _write_json(staged["gate"], gate)
        render_flex_ddg_pilot_figure(
            task_rows=result["task_metric_rows"],
            projection_rows=result["projection_rows"],
            png_path=staged["figure_png"],
            svg_path=staged["figure_svg"],
        )
        _write_json(
            staged_summary,
            {
                "schema_version": 1,
                "status": result["gate_status"],
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    *sys.argv[1:],
                ],
                "task_count": len(result["task_metric_rows"]),
                "status_counts": result["status_counts"],
                "candidate_selection_performed": False,
                "tier_3_scope_decision_performed": False,
                "outputs": {key: str(path) for key, path in final_paths.items()},
            },
        )
        replace_staged_files(
            {
                **{staged[key]: final_paths[key] for key in OUTPUT_NAMES},
                staged_summary: run_summary,
            },
            project_root=PROJECT_ROOT,
            protected_source_paths=validated.source_paths,
        )
    return 0 if result["gate_status"] == "pass" else 2


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
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
