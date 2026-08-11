#!/usr/bin/env python3
"""Merge all 12 unfiltered affinity scan shards into one complete result."""

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

from antibody_optimization.affinity_full_scan import (  # noqa: E402
    PLOT_FIELDS,
    SHARD_COUNT,
    merge_full_scan_shards,
)
from antibody_optimization.affinity_full_scan_plot import (  # noqa: E402
    render_full_scan_figure,
)
from antibody_optimization.affinity_scoring import (  # noqa: E402
    PAIRED_FIELDS,
    SUMMARY_FIELDS,
    WT_CONTROL_FIELDS,
)
from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)


OUTPUT_NAMES = {
    "wt_controls": "wt_replicate_metrics.csv",
    "paired": "candidate_replicate_metrics.csv",
    "summary": "candidate_summary.csv",
    "plot_data": "full_scan_plot_data.csv",
    "gate": "full_scan_merge_gate.json",
    "figure_png": "affinity_full_scan_qc.png",
    "figure_svg": "affinity_full_scan_qc.svg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--scan-plan-dir", type=Path, required=True)
    parser.add_argument("--shard-root", type=Path, required=True)
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
    candidate_dir = _project_directory(args.candidate_dir)
    plan_dir = _project_directory(args.scan_plan_dir)
    shard_root = _project_directory(args.shard_root)
    candidate_path = candidate_dir / "affinity_single_mutants.csv"
    assignment_path = plan_dir / "full_scan_shard_manifest.csv"
    plan_path = plan_dir / "full_scan_plan.json"
    candidates = _load_csv(candidate_path)
    assignments = _load_csv(assignment_path)
    plan = _load_json(plan_path)
    if (
        plan.get("status") != "pass"
        or plan.get("shard_count") != SHARD_COUNT
        or plan.get("maximum_array_concurrency") != 4
        or plan.get("candidate_filtering_applied") is not False
    ):
        raise ValueError("Full-scan plan is not released")

    shards = []
    shard_sources: list[Path] = []
    for index in range(SHARD_COUNT):
        shard_id = f"shard_{index:02d}"
        shard_dir = shard_root / shard_id
        paths = {
            "gate": shard_dir / "affinity_scoring_shard_gate.json",
            "wt_rows": shard_dir / "wt_replicate_metrics.csv",
            "paired_rows": shard_dir / "candidate_replicate_metrics.csv",
            "summary_rows": shard_dir / "candidate_summary.csv",
            "run_summary": shard_dir / "shard_run_summary.json",
        }
        shard_sources.extend(paths.values())
        shards.append(
            {
                "shard_id": shard_id,
                "gate": _load_json(paths["gate"]),
                "wt_rows": _load_csv(paths["wt_rows"]),
                "paired_rows": _load_csv(paths["paired_rows"]),
                "summary_rows": _load_csv(paths["summary_rows"]),
                "run_summary": _load_json(paths["run_summary"]),
            }
        )
    merged = merge_full_scan_shards(
        candidates=candidates,
        assignments=assignments,
        shards=shards,
    )
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_affinity_pyrosetta_full_scan_merge",
        "status": "pass",
        "generated_at": generated_at,
        **merged["counts"],
        "aggregate_shard_elapsed_seconds": merged[
            "aggregate_shard_elapsed_seconds"
        ],
        "candidate_filtering_applied": False,
        "full_scan_contract": "score_all_declared_candidates_then_filter_once",
        "merged_candidate_set_complete": True,
        "merged_replicate_keys_complete": True,
        "cross_shard_wt_controls_consistent_within_tolerance": True,
        "cross_shard_wt_numeric_tolerance": 1e-6,
        "unified_filtering_release": "ready_for_post_scan_filter_implementation",
        "score_semantics": "mutant_minus_paired_WT_Rosetta_ranking_signal",
    }

    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[candidate_path, assignment_path, plan_path, *shard_sources],
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
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )
    for path in [*final_paths.values(), run_summary]:
        path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".affinity-full-merge-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = staging / "run_summary.json"
        _write_csv(staged["wt_controls"], merged["wt_rows"], WT_CONTROL_FIELDS)
        _write_csv(staged["paired"], merged["paired_rows"], PAIRED_FIELDS)
        _write_csv(staged["summary"], merged["summary_rows"], SUMMARY_FIELDS)
        _write_csv(staged["plot_data"], merged["plot_rows"], PLOT_FIELDS)
        _write_json(staged["gate"], gate)
        render_full_scan_figure(
            rows=merged["plot_rows"],
            png_path=staged["figure_png"],
            svg_path=staged["figure_svg"],
        )
        _write_json(
            staged_summary,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "counts": merged["counts"],
                "candidate_filtering_applied": False,
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


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
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
