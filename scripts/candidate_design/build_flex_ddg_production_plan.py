#!/usr/bin/env python3
"""Build the fixed 50-candidate by 20-sample Flex ddG production plan."""

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
    BACKRUB_NEIGHBORHOOD_ANGSTROM,
    BACKRUB_SEGMENT_LENGTH_RANGE,
    BACKRUB_TEMPERATURE,
    BACKRUB_TRIALS,
    PILOT_CANDIDATES,
)
from antibody_optimization.flex_ddg_production import (  # noqa: E402
    DEFAULT_ARRAY_CHUNK_SIZE,
    DEFAULT_ARRAY_CONCURRENCY,
    PRODUCTION_CANDIDATE_COUNT,
    PRODUCTION_MANIFEST_FIELDS,
    PRODUCTION_SAMPLES_PER_CANDIDATE,
    PRODUCTION_TIER3_CANDIDATES,
    PRODUCTION_TASK_COUNT,
    build_production_manifest,
)


OUTPUT_NAMES = {
    "manifest": "flex_ddg_production_manifest.csv",
    "plan": "flex_ddg_production_plan.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tier-dir", type=Path, required=True)
    parser.add_argument("--pilot-result-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--check_only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    tier_dir = _project_directory(args.tier_dir)
    pilot_dir = _project_directory(args.pilot_result_dir)
    tier_path = tier_dir / "affinity_candidate_tiers.csv"
    tier_gate_path = tier_dir / "affinity_post_scan_gate.json"
    pilot_review_path = pilot_dir / "flex_ddg_pilot_scientific_review.json"
    tier_gate = _load_json(tier_gate_path)
    pilot_review = _load_json(pilot_review_path)
    if (
        tier_gate.get("status") != "pass"
        or tier_gate.get("candidate_selection_performed") is not False
        or pilot_review.get("execution", {}).get("status") != "pass"
        or pilot_review.get("integrity_review", {}).get("status") != "pass"
    ):
        raise ValueError("Tier or Flex ddG pilot evidence does not release production")

    tasks = build_production_manifest(_load_csv(tier_path))
    tier123_p90_at_8 = float(
        pilot_review["scope_projections_concurrency_8_overhead_factor_1_2"]
        ["tier_1_2_3_87_candidates_20_samples"]["p90_wall_hours"]
    )
    production_p90_at_8 = tier123_p90_at_8 * PRODUCTION_CANDIDATE_COUNT / 87
    plan = {
        "schema_version": 1,
        "plan_name": "nb252_flex_ddg_production",
        "status": "pass",
        "generated_at": generated_at,
        "purpose": "tier_1_2_plus_selected_tier_3_ensemble_affinity_review",
        "candidate_selection_performed": False,
        "tier_3_scope_decision_performed": True,
        "tier_3_scope_decision_source": "user_revision_for_30_sequence_delivery_2026-08-12",
        "scope": "all_tier_1_and_2_plus_two_selected_tier_3",
        "selected_tier_3_candidate_ids": list(PRODUCTION_TIER3_CANDIDATES),
        "candidate_count": PRODUCTION_CANDIDATE_COUNT,
        "samples_per_candidate": PRODUCTION_SAMPLES_PER_CANDIDATE,
        "task_count": PRODUCTION_TASK_COUNT,
        "precheck_candidate_ids": [
            item[0] for item in PILOT_CANDIDATES[:3]
        ] + list(PRODUCTION_TIER3_CANDIDATES),
        "backrub": {
            "trials": BACKRUB_TRIALS,
            "temperature_kT": BACKRUB_TEMPERATURE,
            "mutation_neighborhood_angstrom": BACKRUB_NEIGHBORHOOD_ANGSTROM,
            "segment_length_residues": BACKRUB_SEGMENT_LENGTH_RANGE,
            "sequence_during_sampling": "WT",
            "branching": "same_sampled_backbone_to_independent_WT_and_mutant_branches",
        },
        "score_function": "ref2015",
        "protocol_identity": "project_ref2015_paired_backbone_flex_ddg_production",
        "starting_structure_state": "selected_wt_prepared",
        "scheduler": {
            "default_array_concurrency": DEFAULT_ARRAY_CONCURRENCY,
            "concurrency_changes_scientific_task_identity": False,
            "default_array_chunk_size": DEFAULT_ARRAY_CHUNK_SIZE,
            "logs_directory": "logs/flex_ddg_production",
        },
        "resume": {
            "unit": "complete_task_output_set",
            "completed_files_required": 6,
            "identity_and_safety_validation_required": True,
            "partial_or_conflicting_outputs": "block_without_overwrite",
        },
        "runtime_projection": {
            "source": "flex_ddg_pilot_scientific_review.json",
            "overhead_factor": 1.2,
            "p90_wall_hours_at_concurrency_8": production_p90_at_8,
            "p90_wall_hours_at_default_concurrency_12": production_p90_at_8 * 8 / 12,
        },
        "interpretation": (
            "All 48 Tier-1/2 and two position-diverse Tier-3 single mutants are evaluated before ensemble "
            "filtering; model-specific Rosetta scores are not measured affinity."
        ),
    }
    if args.check_only:
        print(json.dumps({"status": "pass", "task_count": len(tasks)}, sort_keys=True))
        return 0

    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[tier_path, tier_gate_path, pilot_review_path],
        target_paths=[
            output_dir / OUTPUT_NAMES["manifest"],
            output_dir / OUTPUT_NAMES["plan"],
            run_summary,
        ],
    )
    if any(path.exists() for path in validated.target_paths):
        raise FileExistsError("Refusing to overwrite production plan outputs")
    for path in validated.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final_manifest, final_plan, final_summary = validated.target_paths
    with tempfile.TemporaryDirectory(prefix=".flex-ddg-production-plan-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged_manifest = staging / OUTPUT_NAMES["manifest"]
        staged_plan = staging / OUTPUT_NAMES["plan"]
        staged_summary = staging / "run_summary.json"
        _write_csv(staged_manifest, tasks, PRODUCTION_MANIFEST_FIELDS)
        _write_json(staged_plan, plan)
        _write_json(
            staged_summary,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "candidate_count": PRODUCTION_CANDIDATE_COUNT,
                "task_count": PRODUCTION_TASK_COUNT,
                "outputs": {"manifest": str(final_manifest), "plan": str(final_plan)},
            },
        )
        replace_staged_files(
            {
                staged_manifest: final_manifest,
                staged_plan: final_plan,
                staged_summary: final_summary,
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
