"""Contracts and timing summaries for the Nb252 Flex ddG pilot.

The pilot benchmarks four declared candidates with two production-parameter
samples each.  This module is PyRosetta-free: it validates tier provenance,
builds deterministic task manifests, validates task result handoff, and
projects wall time for candidate scopes.  It never selects candidates or
interprets Rosetta scores as measured affinity.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter
from typing import Mapping, Sequence


PILOT_CANDIDATES = (
    ("Nb252_aff_seq045_R45T", "tier_1", "routine_low_risk_baseline"),
    ("Nb252_aff_seq102_Y102W", "tier_1", "source_sidechain_incomplete_large_aromatic"),
    ("Nb252_aff_seq105_E105W", "tier_2", "large_aromatic_fa_rep_increase"),
    ("Nb252_aff_seq114_S114I", "tier_3", "conditional_receptor_epitope_contact_loss"),
)
PILOT_SAMPLES_PER_CANDIDATE = 2
PILOT_TASK_COUNT = len(PILOT_CANDIDATES) * PILOT_SAMPLES_PER_CANDIDATE
BACKRUB_TRIALS = 35_000
BACKRUB_TEMPERATURE = 1.2
BACKRUB_NEIGHBORHOOD_ANGSTROM = 8.0
BACKRUB_SEGMENT_LENGTH_RANGE = "3-12"
BASE_SEED = 8122000
MAXIMUM_ARRAY_CONCURRENCY = 8

MANIFEST_FIELDS = [
    "task_index",
    "task_id",
    "candidate_id",
    "tier",
    "pilot_role",
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

TASK_METRIC_FIELDS = [
    "task_index",
    "task_id",
    "candidate_id",
    "tier",
    "sample_index",
    "seed",
    "status",
    "total_elapsed_seconds",
    "initial_minimization_seconds",
    "backrub_seconds",
    "wt_branch_seconds",
    "mutant_branch_seconds",
    "measurement_seconds",
    "peak_rss_mb",
    "output_size_bytes",
    "backrub_trials",
    "backrub_neighborhood_residue_count",
    "delta_dG_separated",
    "delta_cross_interface_energy",
    "delta_interface_fa_rep",
    "candidate_vs_paired_wt_vhh_contact_retention",
    "candidate_vs_paired_wt_receptor_epitope_retention",
    "wt_mapping_pass",
    "wt_breaks_pass",
    "wt_disulfide_pass",
    "mutant_mapping_pass",
    "mutant_breaks_pass",
    "mutant_disulfide_pass",
]


class FlexDdgError(ValueError):
    """Raised when a Flex ddG pilot contract is invalid or incomplete."""


def build_pilot_manifest(
    tier_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Build the fixed four-candidate by two-sample pilot task manifest."""

    by_id = {str(row["candidate_id"]): row for row in tier_rows}
    if len(by_id) != len(tier_rows) or len(by_id) != 456:
        raise FlexDdgError("Tier table must contain 456 unique candidates")
    tasks: list[dict[str, object]] = []
    for candidate_number, (candidate_id, expected_tier, role) in enumerate(
        PILOT_CANDIDATES
    ):
        if candidate_id not in by_id:
            raise FlexDdgError(f"Missing declared pilot candidate {candidate_id}")
        candidate = by_id[candidate_id]
        if str(candidate["tier"]) != expected_tier:
            raise FlexDdgError(
                f"Pilot candidate {candidate_id} is not expected {expected_tier}"
            )
        if str(candidate["hard_validity_status"]) != "pass":
            raise FlexDdgError(f"Pilot candidate {candidate_id} failed hard validity")
        source_label = str(candidate["mutation_source_auth_label"])
        chain_id, auth_seq_id = _source_auth_identity(source_label)
        for sample_index in range(1, PILOT_SAMPLES_PER_CANDIDATE + 1):
            task_index = candidate_number * PILOT_SAMPLES_PER_CANDIDATE + sample_index - 1
            seed = BASE_SEED + task_index + 1
            tasks.append(
                {
                    "task_index": task_index,
                    "task_id": f"task_{task_index:02d}",
                    "candidate_id": candidate_id,
                    "tier": expected_tier,
                    "pilot_role": role,
                    "sample_index": sample_index,
                    "seed": seed,
                    "sequence_index_1based": int(candidate["sequence_index_1based"]),
                    "wt_residue": candidate["wt_residue"],
                    "mutant_residue": candidate["mutant_residue"],
                    "experimental_auth_asym_id": chain_id,
                    "experimental_auth_seq_id": auth_seq_id,
                    "experimental_insertion_code": "",
                    "mutation_reported_label": candidate["mutation_reported_label"],
                    "mutation_numbering_label": candidate["mutation_numbering_label"],
                    "mutation_source_auth_label": source_label,
                    "risk_flags": candidate["risk_flags"],
                }
            )
    _validate_manifest(tasks)
    return tasks


def summarize_pilot_results(
    *,
    manifest_rows: Sequence[Mapping[str, object]],
    task_results: Sequence[Mapping[str, object]],
    overhead_factor: float = 1.2,
    concurrency: int = MAXIMUM_ARRAY_CONCURRENCY,
) -> dict[str, object]:
    """Validate eight task results and project full-scope wall times."""

    _validate_manifest(manifest_rows)
    if overhead_factor < 1.0 or concurrency <= 0:
        raise FlexDdgError("Timing projection parameters are invalid")
    by_task = {str(item["task_id"]): item for item in task_results}
    if len(by_task) != len(task_results):
        raise FlexDdgError("Task results contain duplicate task IDs")
    expected_ids = {str(row["task_id"]) for row in manifest_rows}
    if set(by_task) != expected_ids:
        missing = sorted(expected_ids - set(by_task))
        extra = sorted(set(by_task) - expected_ids)
        raise FlexDdgError(f"Pilot task results are incomplete; missing={missing}, extra={extra}")

    metric_rows: list[dict[str, object]] = []
    for manifest in sorted(manifest_rows, key=lambda row: int(row["task_index"])):
        result = by_task[str(manifest["task_id"])]
        for field in ("task_index", "task_id", "candidate_id", "sample_index", "seed"):
            if str(result[field]) != str(manifest[field]):
                raise FlexDdgError(
                    f"Task result identity mismatch for {manifest['task_id']} field {field}"
                )
        row = {field: result[field] for field in TASK_METRIC_FIELDS}
        for field in (
            "total_elapsed_seconds",
            "initial_minimization_seconds",
            "backrub_seconds",
            "wt_branch_seconds",
            "mutant_branch_seconds",
            "measurement_seconds",
            "peak_rss_mb",
            "delta_dG_separated",
            "delta_cross_interface_energy",
            "delta_interface_fa_rep",
            "candidate_vs_paired_wt_vhh_contact_retention",
            "candidate_vs_paired_wt_receptor_epitope_retention",
        ):
            if not math.isfinite(float(row[field])):
                raise FlexDdgError(f"Non-finite task metric {field} for {manifest['task_id']}")
        safety_fields = (
            "wt_mapping_pass",
            "wt_breaks_pass",
            "wt_disulfide_pass",
            "mutant_mapping_pass",
            "mutant_breaks_pass",
            "mutant_disulfide_pass",
        )
        if str(row["status"]) == "pass" and not all(
            _as_bool(row[field]) for field in safety_fields
        ):
            raise FlexDdgError(
                f"Passing task has failed structural safety for {manifest['task_id']}"
            )
        metric_rows.append(row)

    elapsed = [float(row["total_elapsed_seconds"]) for row in metric_rows]
    timing_summary_rows = []
    for field in (
        "total_elapsed_seconds",
        "initial_minimization_seconds",
        "backrub_seconds",
        "wt_branch_seconds",
        "mutant_branch_seconds",
        "measurement_seconds",
        "peak_rss_mb",
    ):
        values = sorted(float(row[field]) for row in metric_rows)
        timing_summary_rows.append(
            {
                "metric": field,
                "minimum": min(values),
                "median": statistics.median(values),
                "p90": _quantile(values, 0.90),
                "maximum": max(values),
            }
        )
    projections = []
    for scope, candidate_count in (
        ("tier_1_2", 48),
        ("tier_1_2_3", 87),
    ):
        task_count = candidate_count * 20
        projections.append(
            {
                "scope": scope,
                "candidate_count": candidate_count,
                "samples_per_candidate": 20,
                "task_count": task_count,
                "concurrency": concurrency,
                "overhead_factor": overhead_factor,
                "projected_wall_hours_from_median": (
                    task_count * statistics.median(elapsed) / concurrency / 3600 * overhead_factor
                ),
                "projected_wall_hours_from_p90": (
                    task_count * _quantile(sorted(elapsed), 0.90)
                    / concurrency
                    / 3600
                    * overhead_factor
                ),
                "projected_wall_hours_from_maximum": (
                    task_count * max(elapsed) / concurrency / 3600 * overhead_factor
                ),
            }
        )
    status_counts = Counter(str(row["status"]) for row in metric_rows)
    return {
        "task_metric_rows": metric_rows,
        "timing_summary_rows": timing_summary_rows,
        "projection_rows": projections,
        "status_counts": dict(status_counts),
        "gate_status": "pass" if status_counts == {"pass": PILOT_TASK_COUNT} else "blocked",
    }


def _validate_manifest(rows: Sequence[Mapping[str, object]]) -> None:
    if len(rows) != PILOT_TASK_COUNT:
        raise FlexDdgError(f"Pilot manifest must contain {PILOT_TASK_COUNT} tasks")
    indices = [int(row["task_index"]) for row in rows]
    ids = [str(row["task_id"]) for row in rows]
    seeds = [int(row["seed"]) for row in rows]
    if indices != list(range(PILOT_TASK_COUNT)) or len(set(ids)) != len(ids):
        raise FlexDdgError("Pilot task indices or IDs are not deterministic and unique")
    if len(set(seeds)) != len(seeds) or any(seed <= 0 for seed in seeds):
        raise FlexDdgError("Pilot seeds must be positive and unique")
    observed = Counter((str(row["candidate_id"]), int(row["sample_index"])) for row in rows)
    expected = Counter(
        (candidate_id, sample_index)
        for candidate_id, _, _ in PILOT_CANDIDATES
        for sample_index in range(1, PILOT_SAMPLES_PER_CANDIDATE + 1)
    )
    if observed != expected:
        raise FlexDdgError("Pilot candidate/sample coverage is incorrect")


def _source_auth_identity(label: str) -> tuple[str, int]:
    # Fixed manifest grammar produced by build_affinity_single_mutants.py.
    parts = label.split()
    try:
        chain_id = parts[2]
        auth_seq_id = int(parts[4])
    except (IndexError, ValueError) as exc:
        raise FlexDdgError(f"Cannot parse source-auth mutation label {label!r}") from exc
    if len(chain_id) != 1:
        raise FlexDdgError(f"Invalid source-auth chain in {label!r}")
    return chain_id, auth_seq_id


def _quantile(sorted_values: Sequence[float], probability: float) -> float:
    if not sorted_values:
        raise FlexDdgError("Cannot calculate a quantile from an empty sequence")
    position = (len(sorted_values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(sorted_values[lower])
    fraction = position - lower
    return float(sorted_values[lower]) * (1.0 - fraction) + float(
        sorted_values[upper]
    ) * fraction


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if str(value).strip().lower() in {"true", "1"}:
        return True
    if str(value).strip().lower() in {"false", "0"}:
        return False
    raise FlexDdgError(f"Invalid boolean value {value!r}")
