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
    "delta_dG_separated",
    "delta_cross_interface_energy",
    "delta_interface_fa_rep",
    "mutant_mapping_pass",
    "mutant_breaks_pass",
    "mutant_disulfide_pass",
    "mutant_finite_metrics",
    "mutant_runtime_valid",
    "mutant_structure_gate",
    "status",
]

WT_CONTROL_FIELDS = [
    "wt_control_id",
    "replicate",
    "seed",
    *METRIC_FIELDS,
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
    "passing_replicate_count",
    "delta_dG_separated_median",
    "delta_dG_separated_mad",
    "delta_cross_interface_energy_median",
    "delta_interface_fa_rep_median",
    "minimum_vhh_contact_retention",
    "minimum_receptor_epitope_retention",
    "maximum_interface_ca_rmsd",
    "status",
    "interpretation",
]


def build_paired_row(
    candidate: Mapping[str, object],
    *,
    replicate: int,
    seed: int,
    wt_metrics: Mapping[str, object],
    mutant_metrics: Mapping[str, object],
    minimum_vhh_contact_retention: float = 0.80,
    minimum_receptor_epitope_retention: float = 0.90,
    maximum_interface_ca_rmsd: float = 0.50,
) -> dict[str, object]:
    """Create one paired mutant-minus-WT scoring record."""

    if replicate < 1 or seed <= 0:
        raise AffinityScoringError("Replicate and seed must be positive")
    for label, metrics in (("WT", wt_metrics), ("mutant", mutant_metrics)):
        missing = [field for field in METRIC_FIELDS if field not in metrics]
        if missing:
            raise AffinityScoringError(f"{label} metrics missing fields: {missing}")
        if not all(math.isfinite(float(metrics[field])) for field in METRIC_FIELDS):
            raise AffinityScoringError(f"{label} metrics contain non-finite values")
    runtime_valid = all(
        bool(mutant_metrics.get(field))
        for field in ("mapping_pass", "breaks_pass", "disulfide_pass", "finite_metrics")
    )
    structure_gate = runtime_valid and (
        float(mutant_metrics["vhh_contact_retention"])
        >= minimum_vhh_contact_retention
        and float(mutant_metrics["receptor_epitope_retention"])
        >= minimum_receptor_epitope_retention
        and float(mutant_metrics["interface_ca_rmsd"])
        <= maximum_interface_ca_rmsd
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
            "mutant_structure_gate": structure_gate,
            "status": "pass" if runtime_valid else "blocked",
        }
    )
    return row


def summarize_paired_rows(
    rows: Sequence[Mapping[str, object]], *, expected_replicates: int
) -> list[dict[str, object]]:
    """Summarize each candidate without treating unfavorable energy as run failure."""

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
        all_structure_pass = all(bool(row["mutant_structure_gate"]) for row in selected)
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
                "passing_replicate_count": sum(
                    bool(row["mutant_structure_gate"]) for row in selected
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
                "maximum_interface_ca_rmsd": max(
                    float(row["mutant_interface_ca_rmsd"]) for row in selected
                ),
                "status": "pass" if all_structure_pass else "blocked",
                "interpretation": (
                    "favorable_relative_signal"
                    if all_structure_pass and median_dg < 0
                    else "unfavorable_or_neutral_relative_signal"
                    if all_structure_pass
                    else "runtime_failure"
                    if not runtime_valid
                    else "structurally_not_evaluable"
                ),
            }
        )
    return summaries


def build_pilot_gate(
    *,
    wt_controls: Sequence[Mapping[str, object]],
    paired_rows: Sequence[Mapping[str, object]],
    summaries: Sequence[Mapping[str, object]],
    expected_candidate_count: int,
    expected_replicates: int,
) -> dict[str, object]:
    """Release full scoring only when the pilot route itself is evaluable."""

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
        "schema_version": 1,
        "gate_name": "nb252_affinity_pyrosetta_pilot",
        "status": "pass" if not blockers else "blocked",
        "full_affinity_scan_release": "pass" if not blockers else "blocked",
        "blockers": blockers,
        "candidate_count": len(summaries),
        "wt_control_count": len(wt_controls),
        "replicate_count_per_candidate": expected_replicates,
        "paired_row_count": len(paired_rows),
        "candidate_status_counts": dict(Counter(str(row["status"]) for row in summaries)),
        "structurally_rejected_candidates_do_not_block_workflow": True,
        "candidate_energy_direction_is_not_a_workflow_gate": True,
        "score_semantics": "mutant_minus_paired_WT_Rosetta_ranking_signal",
    }


def build_wt_control_row(
    *, replicate: int, seed: int, metrics: Mapping[str, object]
) -> dict[str, object]:
    """Record one shared WT control exactly once for one replicate seed."""

    missing = [field for field in METRIC_FIELDS if field not in metrics]
    if missing:
        raise AffinityScoringError(f"WT metrics missing fields: {missing}")
    if not all(math.isfinite(float(metrics[field])) for field in METRIC_FIELDS):
        raise AffinityScoringError("WT metrics contain non-finite values")
    return {
        "wt_control_id": f"Nb252_WT_rep{replicate:02d}_seed{seed}",
        "replicate": replicate,
        "seed": seed,
        **{field: metrics[field] for field in METRIC_FIELDS},
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
