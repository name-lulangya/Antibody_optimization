#!/usr/bin/env python3
"""Build the fixed eight-task production-parameter Flex ddG timing pilot."""

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
    MANIFEST_FIELDS,
    MAXIMUM_ARRAY_CONCURRENCY,
    PILOT_CANDIDATES,
    PILOT_SAMPLES_PER_CANDIDATE,
    build_pilot_manifest,
)


OUTPUT_NAMES = {
    "manifest": "flex_ddg_pilot_manifest.csv",
    "plan": "flex_ddg_pilot_plan.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tier-dir", type=Path, required=True)
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
    tier_path = tier_dir / "affinity_candidate_tiers.csv"
    gate_path = tier_dir / "affinity_post_scan_gate.json"
    gate = _load_json(gate_path)
    if (
        gate.get("status") != "pass"
        or gate.get("release") != "ready_for_strict_review_protocol_design"
        or gate.get("candidate_selection_performed") is not False
    ):
        raise ValueError("Post-scan tier gate does not release strict-review planning")
    tasks = build_pilot_manifest(_load_csv(tier_path))
    plan = {
        "schema_version": 1,
        "plan_name": "nb252_flex_ddg_timing_pilot",
        "status": "pass",
        "generated_at": generated_at,
        "purpose": "production_parameter_runtime_and_protocol_feasibility_only",
        "candidate_selection_performed": False,
        "tier_3_scope_decision_performed": False,
        "task_count": len(tasks),
        "candidate_count": len(PILOT_CANDIDATES),
        "samples_per_candidate": PILOT_SAMPLES_PER_CANDIDATE,
        "maximum_array_concurrency": MAXIMUM_ARRAY_CONCURRENCY,
        "backrub": {
            "trials": BACKRUB_TRIALS,
            "temperature_kT": BACKRUB_TEMPERATURE,
            "mutation_neighborhood_angstrom": BACKRUB_NEIGHBORHOOD_ANGSTROM,
            "segment_length_residues": BACKRUB_SEGMENT_LENGTH_RANGE,
            "sequence_during_sampling": "WT",
            "branching": "same_sampled_backbone_to_independent_WT_and_mutant_branches",
        },
        "score_function": "ref2015",
        "protocol_identity": "project_ref2015_paired_backbone_flex_ddg_timing_pilot",
        "starting_structure_state": "selected_wt_prepared",
        "published_protocol_basis": {
            "paper_doi": "10.1371/journal.pone.0182479",
            "rosetta_tutorial": "https://github.com/Kortemme-Lab/flex_ddG_tutorial",
        },
        "project_adaptations": [
            "use project ref2015 instead of legacy talaris-era score weights",
            "reuse the released prepared WT complex and project interface metrics",
            "one final accepted 35000-trial backbone per independent sample",
        ],
        "interpretation": (
            "Timing and protocol feasibility only; two samples per candidate are "
            "insufficient for scientific candidate ranking."
        ),
    }
    if args.check_only:
        print(json.dumps({"status": "pass", "task_count": len(tasks)}, sort_keys=True))
        return 0

    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[tier_path, gate_path],
        target_paths=[
            output_dir / OUTPUT_NAMES["manifest"],
            output_dir / OUTPUT_NAMES["plan"],
            run_summary,
        ],
    )
    final_manifest, final_plan, run_summary = validated.target_paths
    existing = [path for path in validated.target_paths if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )
    for path in validated.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".flex-ddg-plan-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged_manifest = staging / OUTPUT_NAMES["manifest"]
        staged_plan = staging / OUTPUT_NAMES["plan"]
        staged_summary = staging / "run_summary.json"
        _write_csv(staged_manifest, tasks, MANIFEST_FIELDS)
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
                "task_count": len(tasks),
                "candidate_count": len(PILOT_CANDIDATES),
                "outputs": {"manifest": str(final_manifest), "plan": str(final_plan)},
            },
        )
        replace_staged_files(
            {
                staged_manifest: final_manifest,
                staged_plan: final_plan,
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
