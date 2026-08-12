"""Uniform post-scan interpretation for the complete Nb252 affinity scan.

This module validates the already merged 456-candidate result, assigns one
evidence tier to every candidate, records structural/reproducibility risks,
and calculates within-tier Pareto fronts.  It does not run PyRosetta, combine
mutations, or convert Rosetta scores into measured affinity.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Mapping, Sequence


EXPECTED_CANDIDATES = 456
EXPECTED_REPLICATES = 3
TIERS = ("tier_1", "tier_2", "tier_3", "tier_4", "tier_5")


class AffinityPostScanError(ValueError):
    """Raised when post-scan inputs or tier assignments violate the contract."""


def tier_affinity_candidates(
    *,
    candidates: Sequence[Mapping[str, object]],
    summaries: Sequence[Mapping[str, object]],
    replicates: Sequence[Mapping[str, object]],
    merge_gate: Mapping[str, object],
    scientific_review: Mapping[str, object],
    critical_residue_sets: Mapping[str, object],
    calibration_gate: Mapping[str, object],
) -> dict[str, object]:
    """Validate one complete scan and assign deterministic evidence tiers.

    Tiers use paired-WT score directions, contact retention, and fa_rep.  Risk
    flags remain explicit annotations and do not silently override a tier.
    Pareto fronts are calculated separately inside each valid tier.
    """

    _validate_releases(
        merge_gate=merge_gate,
        scientific_review=scientific_review,
        critical_residue_sets=critical_residue_sets,
        calibration_gate=calibration_gate,
    )
    candidate_by_id = _unique_rows(candidates, "candidate manifest")
    summary_by_id = _unique_rows(summaries, "candidate summaries")
    if len(candidate_by_id) != EXPECTED_CANDIDATES:
        raise AffinityPostScanError(
            f"Expected {EXPECTED_CANDIDATES} candidates, observed {len(candidate_by_id)}"
        )
    if set(candidate_by_id) != set(summary_by_id):
        raise AffinityPostScanError("Candidate manifest and summaries disagree")

    replicate_groups: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    replicate_keys: set[tuple[str, int, int]] = set()
    for row in replicates:
        candidate_id = str(row["candidate_id"])
        if candidate_id not in candidate_by_id:
            raise AffinityPostScanError(f"Unknown replicate candidate {candidate_id}")
        key = (candidate_id, int(row["replicate"]), int(row["seed"]))
        if key in replicate_keys:
            raise AffinityPostScanError(f"Duplicate replicate key {key}")
        replicate_keys.add(key)
        replicate_groups[candidate_id].append(row)
    if len(replicate_keys) != EXPECTED_CANDIDATES * EXPECTED_REPLICATES:
        raise AffinityPostScanError("Replicate result count is incomplete")

    incomplete_source_positions = set(
        int(value)
        for value in calibration_gate["raw_import_metrics"][
            "source_incomplete_vhh_interface_positions"
        ]
    )
    maximum_rmsd = float(
        calibration_gate["thresholds"]["maximum_interface_ca_rmsd_angstrom"]
    )
    high_mad_limit = float(calibration_gate["thresholds"]["maximum_dg_mad_reu"])

    output_rows: list[dict[str, object]] = []
    for candidate_id, candidate in candidate_by_id.items():
        summary = summary_by_id[candidate_id]
        selected = sorted(
            replicate_groups[candidate_id], key=lambda row: int(row["replicate"])
        )
        if len(selected) != EXPECTED_REPLICATES or {
            int(row["replicate"]) for row in selected
        } != {1, 2, 3}:
            raise AffinityPostScanError(f"Incomplete replicates for {candidate_id}")
        _validate_identity(candidate, summary, selected)

        hard_blockers = _hard_blockers(selected, summary, maximum_rmsd)
        dg_values = [float(row["delta_dG_separated"]) for row in selected]
        cross_values = [float(row["delta_cross_interface_energy"]) for row in selected]
        dg_negative_count = sum(value < 0.0 for value in dg_values)
        cross_negative_count = sum(value < 0.0 for value in cross_values)
        dg_median = float(summary["delta_dG_separated_median"])
        cross_median = float(summary["delta_cross_interface_energy_median"])
        fa_rep_median = float(summary["delta_interface_fa_rep_median"])
        paired_vhh = float(
            summary["minimum_candidate_vs_paired_wt_vhh_contact_retention"]
        )
        paired_receptor = float(
            summary["minimum_candidate_vs_paired_wt_receptor_epitope_retention"]
        )
        position = int(candidate["sequence_index_1based"])
        prepared_sensitive = _as_bool(candidate["prepared_contact_sensitive"])

        risk_flags: list[str] = []
        if prepared_sensitive:
            risk_flags.append("prepared_contact_sensitive")
        if position in incomplete_source_positions:
            risk_flags.append("source_sidechain_incomplete_position")
        if paired_receptor < 1.0:
            risk_flags.append("receptor_epitope_contact_loss")
        if paired_vhh < 1.0:
            risk_flags.append("vhh_contact_reorganization")
        if fa_rep_median > 0.0:
            risk_flags.append("fa_rep_increase")
        if dg_negative_count not in {0, EXPECTED_REPLICATES} or (
            cross_negative_count not in {0, EXPECTED_REPLICATES}
        ):
            risk_flags.append("replicate_direction_flip")
        if (dg_median < 0.0) != (cross_median < 0.0):
            risk_flags.append("energy_metric_disagreement")
        if float(summary["delta_dG_separated_mad"]) > high_mad_limit:
            risk_flags.append("dg_mad_above_calibration_limit")

        tier, tier_reason = _assign_tier(
            hard_blockers=hard_blockers,
            dg_negative_count=dg_negative_count,
            cross_negative_count=cross_negative_count,
            dg_median=dg_median,
            cross_median=cross_median,
            paired_vhh=paired_vhh,
            paired_receptor=paired_receptor,
            fa_rep_median=fa_rep_median,
        )
        row = {
            "candidate_id": candidate_id,
            "mutation_reported_label": candidate["mutation_reported_label"],
            "mutation_numbering_label": candidate["mutation_numbering_label"],
            "mutation_source_auth_label": candidate["mutation_source_auth_label"],
            "sequence_index_1based": position,
            "wt_residue": candidate["wt_residue"],
            "mutant_residue": candidate["mutant_residue"],
            "region": candidate["region"],
            "prepared_contact_sensitive": prepared_sensitive,
            "delta_dG_separated_median": dg_median,
            "delta_dG_separated_mad": float(summary["delta_dG_separated_mad"]),
            "delta_cross_interface_energy_median": cross_median,
            "delta_interface_fa_rep_median": fa_rep_median,
            "minimum_candidate_vs_paired_wt_vhh_contact_retention": paired_vhh,
            "minimum_candidate_vs_paired_wt_receptor_epitope_retention": paired_receptor,
            "minimum_experimental_vhh_contact_retention": float(
                summary["minimum_vhh_contact_retention"]
            ),
            "minimum_experimental_receptor_epitope_retention": float(
                summary["minimum_receptor_epitope_retention"]
            ),
            "maximum_interface_ca_rmsd": float(summary["maximum_interface_ca_rmsd"]),
            "negative_delta_dG_replicate_count": dg_negative_count,
            "negative_delta_cross_interface_replicate_count": cross_negative_count,
            "hard_validity_status": "pass" if not hard_blockers else "blocked",
            "hard_blockers": ";".join(hard_blockers),
            "tier": tier,
            "tier_reason": tier_reason,
            "risk_flags": ";".join(risk_flags),
            "strict_review_pool": tier in {"tier_1", "tier_2"},
            "pareto_front_within_tier": 0,
            "candidate_selection_performed": False,
        }
        output_rows.append(row)

    _assign_pareto_fronts(output_rows)
    output_rows.sort(
        key=lambda row: (
            TIERS.index(str(row["tier"])) if row["tier"] in TIERS else len(TIERS),
            int(row["pareto_front_within_tier"]),
            float(row["delta_dG_separated_median"]),
            str(row["candidate_id"]),
        )
    )
    tier_counts = Counter(str(row["tier"]) for row in output_rows)
    review_pool = [row.copy() for row in output_rows if row["strict_review_pool"]]
    return {
        "candidate_rows": output_rows,
        "review_pool_rows": review_pool,
        "tier_counts": dict(tier_counts),
        "tier_summary_rows": _tier_summary(output_rows),
        "position_summary_rows": _position_summary(output_rows),
        "region_summary_rows": _region_summary(output_rows),
        "thresholds": {
            "maximum_interface_ca_rmsd_angstrom": maximum_rmsd,
            "high_dg_mad_annotation_reu": high_mad_limit,
            "tier_1_requires_all_replicates_negative_for_both_energy_metrics": True,
            "tier_1_requires_paired_wt_contact_retention": 1.0,
            "tier_1_requires_nonpositive_median_delta_fa_rep": True,
        },
    }


def _assign_tier(
    *,
    hard_blockers: Sequence[str],
    dg_negative_count: int,
    cross_negative_count: int,
    dg_median: float,
    cross_median: float,
    paired_vhh: float,
    paired_receptor: float,
    fa_rep_median: float,
) -> tuple[str, str]:
    if hard_blockers:
        return "invalid_result", "one_or_more_hard_validity_checks_failed"
    all_both_negative = (
        dg_negative_count == EXPECTED_REPLICATES
        and cross_negative_count == EXPECTED_REPLICATES
    )
    if (
        all_both_negative
        and paired_receptor == 1.0
        and paired_vhh == 1.0
        and fa_rep_median <= 0.0
    ):
        return "tier_1", "replicate_consistent_dual_energy_support_low_observed_risk"
    if all_both_negative and paired_receptor == 1.0:
        return "tier_2", "replicate_consistent_dual_energy_support_needs_local_review"
    if dg_median < 0.0 and cross_median < 0.0:
        return "tier_3", "conditional_dual_median_support"
    if (dg_median < 0.0) != (cross_median < 0.0):
        return "tier_4", "median_energy_metrics_disagree"
    return "tier_5", "current_model_has_no_dual_favorable_median_support"


def _hard_blockers(
    replicates: Sequence[Mapping[str, object]],
    summary: Mapping[str, object],
    maximum_rmsd: float,
) -> list[str]:
    blockers: list[str] = []
    required_flags = (
        "mutant_mapping_pass",
        "mutant_breaks_pass",
        "mutant_disulfide_pass",
        "mutant_finite_metrics",
        "mutant_runtime_valid",
    )
    for field in required_flags:
        if not all(_as_bool(row[field]) for row in replicates):
            blockers.append(f"{field}_failed")
    if not all(str(row["status"]) == "pass" for row in replicates):
        blockers.append("replicate_status_failed")
    if str(summary["status"]) != "pass":
        blockers.append("summary_status_failed")
    numeric_fields = (
        "delta_dG_separated_median",
        "delta_dG_separated_mad",
        "delta_cross_interface_energy_median",
        "delta_interface_fa_rep_median",
        "minimum_candidate_vs_paired_wt_vhh_contact_retention",
        "minimum_candidate_vs_paired_wt_receptor_epitope_retention",
        "maximum_interface_ca_rmsd",
    )
    if not all(math.isfinite(float(summary[field])) for field in numeric_fields):
        blockers.append("nonfinite_summary_metric")
    elif float(summary["maximum_interface_ca_rmsd"]) > maximum_rmsd:
        blockers.append("interface_ca_rmsd_above_calibration_limit")
    return blockers


def _assign_pareto_fronts(rows: Sequence[dict[str, object]]) -> None:
    for tier in TIERS:
        remaining = [row for row in rows if row["tier"] == tier]
        front = 1
        while remaining:
            current = [
                row
                for row in remaining
                if not any(_dominates(other, row) for other in remaining if other is not row)
            ]
            if not current:
                raise AffinityPostScanError(f"Could not resolve Pareto front for {tier}")
            for row in current:
                row["pareto_front_within_tier"] = front
            current_ids = {id(row) for row in current}
            remaining = [row for row in remaining if id(row) not in current_ids]
            front += 1


def _dominates(first: Mapping[str, object], second: Mapping[str, object]) -> bool:
    fields = (
        "delta_dG_separated_median",
        "delta_cross_interface_energy_median",
        "delta_dG_separated_mad",
        "delta_interface_fa_rep_median",
        "maximum_interface_ca_rmsd",
    )
    first_values = [float(first[field]) for field in fields] + [
        -float(first["minimum_candidate_vs_paired_wt_vhh_contact_retention"]),
        -float(first["minimum_candidate_vs_paired_wt_receptor_epitope_retention"]),
    ]
    second_values = [float(second[field]) for field in fields] + [
        -float(second["minimum_candidate_vs_paired_wt_vhh_contact_retention"]),
        -float(second["minimum_candidate_vs_paired_wt_receptor_epitope_retention"]),
    ]
    return all(a <= b for a, b in zip(first_values, second_values)) and any(
        a < b for a, b in zip(first_values, second_values)
    )


def _tier_summary(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    counts = Counter(str(row["tier"]) for row in rows)
    return [
        {
            "tier": tier,
            "candidate_count": counts[tier],
            "strict_review_pool": tier in {"tier_1", "tier_2"},
        }
        for tier in TIERS
    ]


def _position_summary(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    grouped: dict[int, list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[int(row["sequence_index_1based"])].append(row)
    output = []
    for position in sorted(grouped):
        selected = grouped[position]
        counts = Counter(str(row["tier"]) for row in selected)
        output.append(
            {
                "sequence_index_1based": position,
                "wt_residue": selected[0]["wt_residue"],
                "region": selected[0]["region"],
                **{f"{tier}_count": counts[tier] for tier in TIERS},
            }
        )
    return output


def _region_summary(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["region"])].append(row)
    output = []
    for region in sorted(grouped):
        counts = Counter(str(row["tier"]) for row in grouped[region])
        output.append(
            {
                "region": region,
                "candidate_count": len(grouped[region]),
                **{f"{tier}_count": counts[tier] for tier in TIERS},
            }
        )
    return output


def _validate_releases(
    *,
    merge_gate: Mapping[str, object],
    scientific_review: Mapping[str, object],
    critical_residue_sets: Mapping[str, object],
    calibration_gate: Mapping[str, object],
) -> None:
    if (
        merge_gate.get("status") != "pass"
        or merge_gate.get("merged_candidate_set_complete") is not True
        or merge_gate.get("merged_replicate_keys_complete") is not True
        or merge_gate.get("candidate_filtering_applied") is not False
        or int(merge_gate.get("candidate_count", 0)) != EXPECTED_CANDIDATES
        or int(merge_gate.get("mutant_evaluation_count", 0))
        != EXPECTED_CANDIDATES * EXPECTED_REPLICATES
    ):
        raise AffinityPostScanError("Full-scan merge gate is not released")
    release = scientific_review.get("scientific_release", {})
    if (
        not isinstance(release, Mapping)
        or release.get("status") != "ready_for_post_scan_filter_implementation"
        or release.get("candidate_selection_performed") is not False
    ):
        raise AffinityPostScanError("Scientific review does not release post-scan filtering")
    if critical_residue_sets.get("status") != "pass":
        raise AffinityPostScanError("Critical residue sets are not released")
    if calibration_gate.get("pyrosetta_affinity_scoring_release") != "pass":
        raise AffinityPostScanError("PyRosetta scoring protocol is not released")


def _validate_identity(
    candidate: Mapping[str, object],
    summary: Mapping[str, object],
    replicates: Sequence[Mapping[str, object]],
) -> None:
    for field in (
        "candidate_id",
        "sequence_index_1based",
        "wt_residue",
        "mutant_residue",
        "mutation_reported_label",
        "mutation_numbering_label",
        "mutation_source_auth_label",
        "region",
    ):
        expected = str(candidate[field])
        if str(summary[field]) != expected or any(
            str(row[field]) != expected for row in replicates
        ):
            raise AffinityPostScanError(
                f"Candidate identity mismatch for {candidate['candidate_id']} field {field}"
            )


def _unique_rows(
    rows: Sequence[Mapping[str, object]], label: str
) -> dict[str, Mapping[str, object]]:
    output = {str(row["candidate_id"]): row for row in rows}
    if len(output) != len(rows):
        raise AffinityPostScanError(f"{label} contains duplicate candidate IDs")
    return output


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if str(value) in {"True", "true", "1"}:
        return True
    if str(value) in {"False", "false", "0"}:
        return False
    raise AffinityPostScanError(f"Invalid Boolean value {value!r}")
