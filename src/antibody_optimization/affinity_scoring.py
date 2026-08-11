"""Pure contracts for paired PyRosetta affinity-candidate scoring.

The runtime produces one WT control and one mutant result for each replicate.
This module validates those paired records, calculates mutant-minus-WT
differences, and summarizes candidates.  Rosetta energy units are retained as
ranking signals and are never converted to measured affinity.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter
from typing import Mapping, Sequence


class AffinityScoringError(ValueError):
    """Raised when paired scoring records violate the released contract."""


METRIC_FIELDS = (
    "total_score",
    "dG_separated",
    "cross_interface_energy",
    "interface_fa_atr",
    "interface_fa_rep",
    "vhh_contact_count",
    "receptor_epitope_count",
    "vhh_contact_retention",
    "receptor_epitope_retention",
    "interface_ca_rmsd",
    "minimum_interchain_distance",
)

CONTACT_SET_FIELDS = (
    "vhh_contact_auth_positions",
    "receptor_contact_auth_positions",
)

PAIRED_FIELDS = [
    "candidate_id",
    "mutation_reported_label",
    "mutation_numbering_label",
    "mutation_source_auth_label",
    "sequence_index_1based",
    "wt_residue",
    "mutant_residue",
    "region",
    "prepared_contact_sensitive",
    "replicate",
    "seed",
    "wt_control_id",
    *[f"mutant_{field}" for field in METRIC_FIELDS],
    *[f"mutant_{field}" for field in CONTACT_SET_FIELDS],
    "paired_wt_vhh_contact_count",
    "paired_wt_receptor_epitope_count",
    "candidate_vs_paired_wt_vhh_contact_retention",
    "candidate_vs_paired_wt_receptor_epitope_retention",
    "delta_dG_separated",
    "delta_cross_interface_energy",
    "delta_interface_fa_rep",
    "mutant_mapping_pass",
    "mutant_breaks_pass",
    "mutant_disulfide_pass",
    "mutant_finite_metrics",
    "mutant_runtime_valid",
    "status",
]

WT_CONTROL_FIELDS = [
    "wt_control_id",
    "replicate",
    "seed",
    *METRIC_FIELDS,
    *CONTACT_SET_FIELDS,
    "mapping_pass",
    "breaks_pass",
    "disulfide_pass",
    "finite_metrics",
    "status",
]

SUMMARY_FIELDS = [
    "candidate_id",
    "mutation_reported_label",
    "mutation_numbering_label",
    "mutation_source_auth_label",
    "sequence_index_1based",
    "wt_residue",
    "mutant_residue",
    "region",
    "prepared_contact_sensitive",
    "replicate_count",
    "runtime_valid_replicate_count",
    "delta_dG_separated_median",
    "delta_dG_separated_mad",
    "delta_cross_interface_energy_median",
    "delta_interface_fa_rep_median",
    "minimum_vhh_contact_retention",
    "minimum_receptor_epitope_retention",
    "minimum_candidate_vs_paired_wt_vhh_contact_retention",
    "minimum_candidate_vs_paired_wt_receptor_epitope_retention",
    "maximum_interface_ca_rmsd",
    "status",
    "selection_status",
    "interpretation",
]


def build_paired_row(
    candidate: Mapping[str, object],
    *,
    replicate: int,
    seed: int,
    wt_metrics: Mapping[str, object],
    mutant_metrics: Mapping[str, object],
) -> dict[str, object]:
    """Create one unfiltered paired mutant-minus-WT scoring record."""

    if replicate < 1 or seed <= 0:
        raise AffinityScoringError("Replicate and seed must be positive")
    for label, metrics in (("WT", wt_metrics), ("mutant", mutant_metrics)):
        missing = [field for field in METRIC_FIELDS if field not in metrics]
        if missing:
            raise AffinityScoringError(f"{label} metrics missing fields: {missing}")
        if not all(math.isfinite(float(metrics[field])) for field in METRIC_FIELDS):
            raise AffinityScoringError(f"{label} metrics contain non-finite values")
        missing_contacts = [field for field in CONTACT_SET_FIELDS if field not in metrics]
        if missing_contacts:
            raise AffinityScoringError(
                f"{label} metrics missing contact sets: {missing_contacts}"
            )
    wt_vhh_contacts = _contact_set(wt_metrics["vhh_contact_auth_positions"])
    wt_receptor_contacts = _contact_set(
        wt_metrics["receptor_contact_auth_positions"]
    )
    mutant_vhh_contacts = _contact_set(
        mutant_metrics["vhh_contact_auth_positions"]
    )
    mutant_receptor_contacts = _contact_set(
        mutant_metrics["receptor_contact_auth_positions"]
    )
    runtime_valid = all(
        bool(mutant_metrics.get(field))
        for field in ("mapping_pass", "breaks_pass", "disulfide_pass", "finite_metrics")
    )
    row = {
        key: candidate[key]
        for key in (
            "candidate_id",
            "mutation_reported_label",
            "mutation_numbering_label",
            "mutation_source_auth_label",
            "sequence_index_1based",
            "wt_residue",
            "mutant_residue",
            "region",
            "prepared_contact_sensitive",
        )
    }
    row.update(
        {
            "replicate": replicate,
            "seed": seed,
            "wt_control_id": f"Nb252_WT_rep{replicate:02d}_seed{seed}",
        }
    )
    row.update({f"mutant_{field}": mutant_metrics[field] for field in METRIC_FIELDS})
    row.update(
        {
            f"mutant_{field}": _serialize_contact_set(mutant_metrics[field])
            for field in CONTACT_SET_FIELDS
        }
    )
    row.update(
        {
            "paired_wt_vhh_contact_count": len(wt_vhh_contacts),
            "paired_wt_receptor_epitope_count": len(wt_receptor_contacts),
            "candidate_vs_paired_wt_vhh_contact_retention": _set_retention(
                wt_vhh_contacts, mutant_vhh_contacts
            ),
            "candidate_vs_paired_wt_receptor_epitope_retention": _set_retention(
                wt_receptor_contacts, mutant_receptor_contacts
            ),
            "delta_dG_separated": float(mutant_metrics["dG_separated"])
            - float(wt_metrics["dG_separated"]),
            "delta_cross_interface_energy": float(
                mutant_metrics["cross_interface_energy"]
            )
            - float(wt_metrics["cross_interface_energy"]),
            "delta_interface_fa_rep": float(mutant_metrics["interface_fa_rep"])
            - float(wt_metrics["interface_fa_rep"]),
            "mutant_mapping_pass": bool(mutant_metrics.get("mapping_pass")),
            "mutant_breaks_pass": bool(mutant_metrics.get("breaks_pass")),
            "mutant_disulfide_pass": bool(mutant_metrics.get("disulfide_pass")),
            "mutant_finite_metrics": bool(mutant_metrics.get("finite_metrics")),
            "mutant_runtime_valid": runtime_valid,
            "status": "pass" if runtime_valid else "blocked",
        }
    )
    return row


def summarize_paired_rows(
    rows: Sequence[Mapping[str, object]], *, expected_replicates: int
) -> list[dict[str, object]]:
    """Summarize every scored candidate without applying scientific selection."""

    if expected_replicates < 3:
        raise AffinityScoringError("At least three paired replicates are required")
    grouped: dict[str, list[Mapping[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["candidate_id"]), []).append(row)
    summaries = []
    for candidate_id, selected in sorted(grouped.items()):
        if len(selected) != expected_replicates:
            raise AffinityScoringError(
                f"{candidate_id} has {len(selected)} rather than {expected_replicates} replicates"
            )
        replicate_ids = [int(row["replicate"]) for row in selected]
        if len(set(replicate_ids)) != expected_replicates:
            raise AffinityScoringError(f"Duplicate replicate for {candidate_id}")
        first = selected[0]
        delta_dg = [float(row["delta_dG_separated"]) for row in selected]
        delta_cross = [float(row["delta_cross_interface_energy"]) for row in selected]
        delta_rep = [float(row["delta_interface_fa_rep"]) for row in selected]
        median_dg = statistics.median(delta_dg)
        runtime_valid = all(bool(row["mutant_runtime_valid"]) for row in selected)
        summaries.append(
            {
                **{
                    key: first[key]
                    for key in (
                        "candidate_id",
                        "mutation_reported_label",
                        "mutation_numbering_label",
                        "mutation_source_auth_label",
                        "sequence_index_1based",
                        "wt_residue",
                        "mutant_residue",
                        "region",
                        "prepared_contact_sensitive",
                    )
                },
                "replicate_count": len(selected),
                "runtime_valid_replicate_count": sum(
                    bool(row["mutant_runtime_valid"]) for row in selected
                ),
                "delta_dG_separated_median": median_dg,
                "delta_dG_separated_mad": statistics.median(
                    abs(value - median_dg) for value in delta_dg
                ),
                "delta_cross_interface_energy_median": statistics.median(delta_cross),
                "delta_interface_fa_rep_median": statistics.median(delta_rep),
                "minimum_vhh_contact_retention": min(
                    float(row["mutant_vhh_contact_retention"]) for row in selected
                ),
                "minimum_receptor_epitope_retention": min(
                    float(row["mutant_receptor_epitope_retention"]) for row in selected
                ),
                "minimum_candidate_vs_paired_wt_vhh_contact_retention": min(
                    float(row["candidate_vs_paired_wt_vhh_contact_retention"])
                    for row in selected
                ),
                "minimum_candidate_vs_paired_wt_receptor_epitope_retention": min(
                    float(
                        row[
                            "candidate_vs_paired_wt_receptor_epitope_retention"
                        ]
                    )
                    for row in selected
                ),
                "maximum_interface_ca_rmsd": max(
                    float(row["mutant_interface_ca_rmsd"]) for row in selected
                ),
                "status": "pass" if runtime_valid else "blocked",
                "selection_status": "not_applied_scan_stage",
                "interpretation": (
                    "unfiltered_relative_scoring_result"
                    if runtime_valid
                    else "runtime_failure"
                ),
            }
        )
    return summaries


def build_scoring_run_gate(
    *,
    wt_controls: Sequence[Mapping[str, object]],
    paired_rows: Sequence[Mapping[str, object]],
    summaries: Sequence[Mapping[str, object]],
    expected_candidate_count: int,
    expected_replicates: int,
    run_kind: str,
    shard_id: str | None = None,
) -> dict[str, object]:
    """Validate one unfiltered pilot or full-scan shard execution."""

    if run_kind not in {"pilot", "full_scan_shard"}:
        raise AffinityScoringError(f"Unsupported run kind: {run_kind}")
    if run_kind == "full_scan_shard" and not shard_id:
        raise AffinityScoringError("Full-scan shard runs require a shard ID")
    if run_kind == "pilot" and shard_id is not None:
        raise AffinityScoringError("Pilot runs must not declare a shard ID")

    expected_rows = expected_candidate_count * expected_replicates
    blockers = []
    if len(paired_rows) != expected_rows:
        blockers.append("paired_row_count_mismatch")
    if len(summaries) != expected_candidate_count:
        blockers.append("candidate_count_mismatch")
    if any(not bool(row["mutant_runtime_valid"]) for row in paired_rows):
        blockers.append("one_or_more_mutant_runtime_failures")
    if len(wt_controls) != expected_replicates or len(
        {int(row["replicate"]) for row in wt_controls}
    ) != expected_replicates:
        blockers.append("wt_control_count_mismatch")
    wt_control_failures = sum(
        not (
            bool(row["mapping_pass"])
            and bool(row["breaks_pass"])
            and bool(row["disulfide_pass"])
            and bool(row["finite_metrics"])
            and float(row["dG_separated"]) < 0
            and float(row["cross_interface_energy"]) < 0
        )
        for row in wt_controls
    )
    if wt_control_failures:
        blockers.append("paired_wt_control_failure")
    return {
        "schema_version": 2,
        "gate_name": (
            "nb252_affinity_pyrosetta_pilot_v2"
            if run_kind == "pilot"
            else "nb252_affinity_pyrosetta_full_scan_shard"
        ),
        "run_kind": run_kind,
        "shard_id": shard_id,
        "status": "pass" if not blockers else "blocked",
        "full_affinity_scan_release": (
            "pass"
            if not blockers and run_kind == "pilot"
            else "pending_complete_merge"
            if not blockers
            else "blocked"
        ),
        "blockers": blockers,
        "candidate_count": len(summaries),
        "wt_control_count": len(wt_controls),
        "replicate_count_per_candidate": expected_replicates,
        "paired_row_count": len(paired_rows),
        "candidate_status_counts": dict(Counter(str(row["status"]) for row in summaries)),
        "candidate_filtering_applied": False,
        "full_scan_contract": "score_all_declared_candidates_then_filter_once",
        "structure_metrics_are_nonblocking_qc": True,
        "candidate_energy_direction_is_not_a_workflow_gate": True,
        "score_semantics": "mutant_minus_paired_WT_Rosetta_ranking_signal",
    }


def build_pilot_gate(
    *,
    wt_controls: Sequence[Mapping[str, object]],
    paired_rows: Sequence[Mapping[str, object]],
    summaries: Sequence[Mapping[str, object]],
    expected_candidate_count: int,
    expected_replicates: int,
) -> dict[str, object]:
    """Backward-compatible pilot wrapper around the generic run gate."""

    return build_scoring_run_gate(
        wt_controls=wt_controls,
        paired_rows=paired_rows,
        summaries=summaries,
        expected_candidate_count=expected_candidate_count,
        expected_replicates=expected_replicates,
        run_kind="pilot",
    )


def build_wt_control_row(
    *, replicate: int, seed: int, metrics: Mapping[str, object]
) -> dict[str, object]:
    """Record one shared WT control exactly once for one replicate seed."""

    missing = [field for field in METRIC_FIELDS if field not in metrics]
    if missing:
        raise AffinityScoringError(f"WT metrics missing fields: {missing}")
    if not all(math.isfinite(float(metrics[field])) for field in METRIC_FIELDS):
        raise AffinityScoringError("WT metrics contain non-finite values")
    missing_contacts = [field for field in CONTACT_SET_FIELDS if field not in metrics]
    if missing_contacts:
        raise AffinityScoringError(f"WT metrics missing contact sets: {missing_contacts}")
    return {
        "wt_control_id": f"Nb252_WT_rep{replicate:02d}_seed{seed}",
        "replicate": replicate,
        "seed": seed,
        **{field: metrics[field] for field in METRIC_FIELDS},
        **{
            field: _serialize_contact_set(metrics[field])
            for field in CONTACT_SET_FIELDS
        },
        **{
            field: bool(metrics.get(field))
            for field in (
                "mapping_pass",
                "breaks_pass",
                "disulfide_pass",
                "finite_metrics",
            )
        },
        "status": str(metrics.get("status", "blocked")),
    }


def _contact_set(value: object) -> set[int]:
    if isinstance(value, str):
        return {int(item) for item in value.split(";") if item}
    if not isinstance(value, (set, frozenset, list, tuple)):
        raise AffinityScoringError(f"Invalid contact-set value: {value!r}")
    return {int(item) for item in value}


def _serialize_contact_set(value: object) -> str:
    return ";".join(str(item) for item in sorted(_contact_set(value)))


def _set_retention(reference: set[int], observed: set[int]) -> float:
    if not reference:
        raise AffinityScoringError("Paired WT contact set must not be empty")
    return len(reference & observed) / len(reference)
