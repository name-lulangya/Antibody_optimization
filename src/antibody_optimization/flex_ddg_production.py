"""Production Flex ddG task planning, recovery, and result aggregation.

The scientific task identity is independent of Slurm scheduling.  A fixed
50-candidate by 20-sample manifest is built once; later submissions may vary
their array throttle without changing candidates, samples, or seeds.  Resume
logic recognizes only complete, identity-matched, structurally safe task
outputs and never deletes or overwrites ambiguous partial results.
"""

from __future__ import annotations

import json
import math
import statistics
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

from .flex_ddg import TASK_METRIC_FIELDS


PRODUCTION_TIERS = ("tier_1", "tier_2", "tier_3")
PRODUCTION_TIER_COUNTS = {"tier_1": 18, "tier_2": 30, "tier_3": 2}
PRODUCTION_TIER3_CANDIDATES = (
    "Nb252_aff_seq033_D33N",
    "Nb252_aff_seq115_Y115F",
)
PRODUCTION_CANDIDATE_COUNT = 50
PRODUCTION_SAMPLES_PER_CANDIDATE = 20
PRODUCTION_TASK_COUNT = 1000
PRODUCTION_BASE_SEED = 8123000
DEFAULT_ARRAY_CONCURRENCY = 12
DEFAULT_ARRAY_CHUNK_SIZE = 900
TASK_OUTPUT_NAMES = (
    "task_result.json",
    "energy_terms.csv",
    "contact_qc.csv",
    "backrub_backbone.pdb",
    "wt_final.pdb",
    "mutant_final.pdb",
)
PRODUCTION_MANIFEST_FIELDS = [
    "task_index",
    "task_id",
    "candidate_id",
    "tier",
    "sample_index",
    "seed",
    "sequence_index_1based",
    "wt_residue",
    "mutant_residue",
    "experimental_auth_asym_id",
    "experimental_auth_seq_id",
    "experimental_insertion_code",
    "mutation_reported_label",
    "mutation_numbering_label",
    "mutation_source_auth_label",
    "risk_flags",
]


class FlexDdgProductionError(ValueError):
    """Raised when a production task, recovery, or aggregation contract fails."""


def build_production_manifest(
    tier_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Build the deterministic Tier-1/2/3 by 20-sample production manifest."""

    if len(tier_rows) != 456 or len({str(row["candidate_id"]) for row in tier_rows}) != 456:
        raise FlexDdgProductionError("Tier table must contain 456 unique candidates")
    selected = [
        row
        for row in tier_rows
        if str(row["tier"]) in ("tier_1", "tier_2")
        or str(row["candidate_id"]) in PRODUCTION_TIER3_CANDIDATES
    ]
    counts = Counter(str(row["tier"]) for row in selected)
    if counts != Counter(PRODUCTION_TIER_COUNTS):
        raise FlexDdgProductionError(f"Unexpected Tier-1/2/3 counts: {dict(counts)}")
    selected_tier3 = {
        str(row["candidate_id"]) for row in selected if str(row["tier"]) == "tier_3"
    }
    if selected_tier3 != set(PRODUCTION_TIER3_CANDIDATES):
        raise FlexDdgProductionError("Unexpected selected Tier-3 candidates")
    if any(str(row["hard_validity_status"]) != "pass" for row in selected):
        raise FlexDdgProductionError("Production candidates must all pass hard validity")

    selected.sort(
        key=lambda row: (
            PRODUCTION_TIERS.index(str(row["tier"])),
            int(row["sequence_index_1based"]),
            str(row["mutant_residue"]),
            str(row["candidate_id"]),
        )
    )
    tasks: list[dict[str, object]] = []
    for candidate_number, candidate in enumerate(selected):
        chain_id, auth_seq_id = _source_auth_identity(
            str(candidate["mutation_source_auth_label"])
        )
        for sample_index in range(1, PRODUCTION_SAMPLES_PER_CANDIDATE + 1):
            task_index = candidate_number * PRODUCTION_SAMPLES_PER_CANDIDATE + sample_index - 1
            tasks.append(
                {
                    "task_index": task_index,
                    "task_id": f"task_{task_index:04d}",
                    "candidate_id": candidate["candidate_id"],
                    "tier": candidate["tier"],
                    "sample_index": sample_index,
                    "seed": PRODUCTION_BASE_SEED + task_index + 1,
                    "sequence_index_1based": int(candidate["sequence_index_1based"]),
                    "wt_residue": candidate["wt_residue"],
                    "mutant_residue": candidate["mutant_residue"],
                    "experimental_auth_asym_id": chain_id,
                    "experimental_auth_seq_id": auth_seq_id,
                    "experimental_insertion_code": "",
                    "mutation_reported_label": candidate["mutation_reported_label"],
                    "mutation_numbering_label": candidate["mutation_numbering_label"],
                    "mutation_source_auth_label": candidate["mutation_source_auth_label"],
                    "risk_flags": candidate["risk_flags"],
                }
            )
    validate_production_manifest(tasks)
    return tasks


def validate_production_manifest(rows: Sequence[Mapping[str, object]]) -> None:
    """Validate the complete immutable production task identity."""

    if len(rows) != PRODUCTION_TASK_COUNT:
        raise FlexDdgProductionError(
            f"Production manifest must contain {PRODUCTION_TASK_COUNT} tasks"
        )
    indices = [int(row["task_index"]) for row in rows]
    ids = [str(row["task_id"]) for row in rows]
    seeds = [int(row["seed"]) for row in rows]
    if indices != list(range(PRODUCTION_TASK_COUNT)):
        raise FlexDdgProductionError("Production task indices must be contiguous")
    if ids != [f"task_{index:04d}" for index in indices] or len(set(ids)) != len(ids):
        raise FlexDdgProductionError("Production task IDs are invalid")
    if len(set(seeds)) != len(seeds) or min(seeds) <= 0:
        raise FlexDdgProductionError("Production seeds must be positive and unique")
    coverage = Counter(
        (str(row["candidate_id"]), int(row["sample_index"])) for row in rows
    )
    candidates = {str(row["candidate_id"]) for row in rows}
    expected = Counter(
        (candidate_id, sample_index)
        for candidate_id in candidates
        for sample_index in range(1, PRODUCTION_SAMPLES_PER_CANDIDATE + 1)
    )
    if len(candidates) != PRODUCTION_CANDIDATE_COUNT or coverage != expected:
        raise FlexDdgProductionError("Production candidate/sample coverage is incomplete")
    if Counter(str(row["tier"]) for row in rows) != Counter(
        {
            tier: count * PRODUCTION_SAMPLES_PER_CANDIDATE
            for tier, count in PRODUCTION_TIER_COUNTS.items()
        }
    ):
        raise FlexDdgProductionError("Production task tier counts are invalid")


def assess_task_outputs(
    *, manifest_rows: Sequence[Mapping[str, object]], task_root: Path
) -> dict[str, object]:
    """Classify complete tasks and safe-to-submit missing tasks without mutation."""

    validate_production_manifest(manifest_rows)
    completed: list[int] = []
    pending: list[int] = []
    invalid: list[dict[str, object]] = []
    for task in manifest_rows:
        index = int(task["task_index"])
        task_dir = task_root / str(task["task_id"])
        existing = [name for name in TASK_OUTPUT_NAMES if (task_dir / name).exists()]
        if not existing:
            pending.append(index)
            continue
        if set(existing) != set(TASK_OUTPUT_NAMES):
            invalid.append(
                {"task_index": index, "task_id": task["task_id"], "reason": "partial_output_set"}
            )
            continue
        if any((task_dir / name).stat().st_size == 0 for name in TASK_OUTPUT_NAMES):
            invalid.append(
                {
                    "task_index": index,
                    "task_id": task["task_id"],
                    "reason": "empty_completed_output",
                }
            )
            continue
        try:
            result = json.loads((task_dir / "task_result.json").read_text(encoding="utf-8"))
            _validate_completed_result(task, result)
        except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
            invalid.append(
                {
                    "task_index": index,
                    "task_id": task["task_id"],
                    "reason": f"invalid_completed_result:{exc}",
                }
            )
            continue
        completed.append(index)
    return {
        "completed_task_indices": completed,
        "pending_task_indices": pending,
        "invalid_tasks": invalid,
        "completed_count": len(completed),
        "pending_count": len(pending),
        "invalid_count": len(invalid),
    }


def chunk_task_indices(indices: Sequence[int], chunk_size: int) -> list[list[int]]:
    """Split pending global task indices into scheduler-safe sequential batches."""

    if chunk_size <= 0 or len(set(indices)) != len(indices):
        raise FlexDdgProductionError("Pending indices or chunk size are invalid")
    ordered = sorted(int(value) for value in indices)
    if any(value < 0 or value >= PRODUCTION_TASK_COUNT for value in ordered):
        raise FlexDdgProductionError("Pending task index is outside the manifest")
    return [ordered[start : start + chunk_size] for start in range(0, len(ordered), chunk_size)]


def summarize_production_results(
    *, manifest_rows: Sequence[Mapping[str, object]], task_results: Sequence[Mapping[str, object]]
) -> dict[str, object]:
    """Validate all 1000 task results and aggregate 20-sample candidate evidence."""

    validate_production_manifest(manifest_rows)
    by_id = {str(row["task_id"]): row for row in task_results}
    if len(by_id) != len(task_results) or set(by_id) != {
        str(row["task_id"]) for row in manifest_rows
    }:
        raise FlexDdgProductionError("Production task results are not unique and complete")
    task_rows = []
    for task in manifest_rows:
        result = by_id[str(task["task_id"])]
        _validate_completed_result(task, result)
        row = {field: result[field] for field in TASK_METRIC_FIELDS}
        task_rows.append(row)

    candidate_rows = []
    for candidate_id in dict.fromkeys(str(row["candidate_id"]) for row in manifest_rows):
        rows = [row for row in task_rows if str(row["candidate_id"]) == candidate_id]
        candidate_manifest = next(
            row for row in manifest_rows if str(row["candidate_id"]) == candidate_id
        )
        if len(rows) != PRODUCTION_SAMPLES_PER_CANDIDATE:
            raise FlexDdgProductionError(f"Incomplete candidate samples for {candidate_id}")
        first = rows[0]
        dg = [float(row["delta_dG_separated"]) for row in rows]
        cross = [float(row["delta_cross_interface_energy"]) for row in rows]
        rep = [float(row["delta_interface_fa_rep"]) for row in rows]
        candidate_rows.append(
            {
                "candidate_id": candidate_id,
                "tier": first["tier"],
                "sequence_index_1based": int(
                    candidate_manifest["sequence_index_1based"]
                ),
                "wt_residue": candidate_manifest["wt_residue"],
                "mutant_residue": candidate_manifest["mutant_residue"],
                "mutation_reported_label": candidate_manifest[
                    "mutation_reported_label"
                ],
                "mutation_numbering_label": candidate_manifest[
                    "mutation_numbering_label"
                ],
                "mutation_source_auth_label": candidate_manifest[
                    "mutation_source_auth_label"
                ],
                "risk_flags": candidate_manifest["risk_flags"],
                "sample_count": len(rows),
                "delta_dG_separated_median": statistics.median(dg),
                "delta_dG_separated_mad": _mad(dg),
                "negative_delta_dG_count": sum(value < 0 for value in dg),
                "delta_cross_interface_energy_median": statistics.median(cross),
                "delta_cross_interface_energy_mad": _mad(cross),
                "negative_delta_cross_interface_count": sum(value < 0 for value in cross),
                "delta_interface_fa_rep_median": statistics.median(rep),
                "minimum_vhh_contact_retention": min(
                    float(row["candidate_vs_paired_wt_vhh_contact_retention"])
                    for row in rows
                ),
                "minimum_receptor_epitope_retention": min(
                    float(row["candidate_vs_paired_wt_receptor_epitope_retention"])
                    for row in rows
                ),
                "candidate_selection_performed": False,
            }
        )
    return {"task_rows": task_rows, "candidate_rows": candidate_rows, "gate_status": "pass"}


def _validate_completed_result(task: Mapping[str, object], result: Mapping[str, object]) -> None:
    if result.get("run_kind") != "production":
        raise FlexDdgProductionError(
            f"Task result is not a production result: {task['task_id']}"
        )
    if result.get("tier_3_scope_decision_performed") is not True:
        raise FlexDdgProductionError(
            f"Task result lacks the Tier-3 scope decision: {task['task_id']}"
        )
    for field in ("task_index", "task_id", "candidate_id", "tier", "sample_index", "seed"):
        if str(result[field]) != str(task[field]):
            raise FlexDdgProductionError(f"Task result identity mismatch: {task['task_id']} {field}")
    if str(result["status"]) != "pass":
        raise FlexDdgProductionError(f"Task is not pass: {task['task_id']}")
    for field in TASK_METRIC_FIELDS:
        if field not in result:
            raise FlexDdgProductionError(f"Task result lacks {field}: {task['task_id']}")
    for field in (
        "wt_mapping_pass",
        "wt_breaks_pass",
        "wt_disulfide_pass",
        "mutant_mapping_pass",
        "mutant_breaks_pass",
        "mutant_disulfide_pass",
    ):
        if not _as_bool(result[field]):
            raise FlexDdgProductionError(f"Task safety failed: {task['task_id']} {field}")
    for field in (
        "total_elapsed_seconds",
        "delta_dG_separated",
        "delta_cross_interface_energy",
        "delta_interface_fa_rep",
    ):
        if not math.isfinite(float(result[field])):
            raise FlexDdgProductionError(f"Task metric is non-finite: {task['task_id']} {field}")


def _source_auth_identity(label: str) -> tuple[str, int]:
    parts = label.split()
    try:
        chain_id, auth_seq_id = parts[2], int(parts[4])
    except (IndexError, ValueError) as exc:
        raise FlexDdgProductionError(f"Cannot parse source-auth label {label!r}") from exc
    if len(chain_id) != 1:
        raise FlexDdgProductionError(f"Invalid chain in source-auth label {label!r}")
    return chain_id, auth_seq_id


def _mad(values: Sequence[float]) -> float:
    median = statistics.median(values)
    return statistics.median(abs(value - median) for value in values)


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1"}
