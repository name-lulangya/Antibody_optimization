"""Contracts for the Nb252 property-candidate affinity non-inferiority review.

The module selects the fixed 30-member review pool from already computed
multi-tool evidence and validates paired PyRosetta runs.  It does not predict
expression, measured stability, or affinity, and it never filters candidates
from PyRosetta scores during the scan.
"""

from __future__ import annotations

from collections import Counter
import statistics
from typing import Mapping, Sequence


POOL_SIZE = 30
POSITION_COUNT = 10
REPLICATES = 3
MUTATION_NEIGHBORHOOD_ANGSTROM = 8.0
PILOT_CANDIDATES = (
    "Nb252_uni_seq001_Q1D",
    "Nb252_uni_seq005_Q5P",
    "Nb252_uni_seq023_A23R",
    "Nb252_uni_seq030_F30K",
    "Nb252_uni_seq055_S55G",
    "Nb252_uni_seq075_K75E",
)


class PropertyAffinityReviewError(ValueError):
    """Raised when the property-affinity review contract is violated."""


PROPERTY_FIELDS = (
    "property_pareto_layer",
    "netsolp_delta_usability_vs_wt",
    "netsolp_usability_magnitude",
    "netsolp_delta_solubility_vs_wt",
    "netsolp_solubility_magnitude",
    "nanomelt_delta_predicted_apparent_tm_c_vs_wt",
    "nanomelt_tm_magnitude",
    "experimental_complex_context_delta_log_probability",
    "material_favorable_count",
    "material_adverse_count",
    "chemical_risk_count",
    "new_liability_flags",
    "tnp_psh_delta_vs_wt",
    "tnp_flag_psh_vs_wt",
    "tnp_flag_regression_count",
    "tnp_flag_improvement_count",
    "tnp_developability_review",
)


def build_review_pool(
    evidence_rows: Sequence[Mapping[str, str]],
    candidate_rows: Sequence[Mapping[str, str]],
    mapping_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, object]]:
    """Build the exact 30-member review pool from released real-data records.

    Selection is restricted to property Pareto layer 1 with at least one
    material U/S/Tm improvement, no material U/S/Tm adverse change, and no new
    sequence-liability flag. AntiFold and TNP remain evidence columns rather
    than selection gates.
    """

    candidates = _unique(candidate_rows, "candidate_id")
    experimental_mapping = {
        int(row["sequence_index_1based"]): row
        for row in mapping_rows
        if row.get("source_model_role") == "experimental_nk2r_nb252"
    }
    selected = [
        row
        for row in evidence_rows
        if row.get("candidate_source") == "property_pareto_front_1"
        and int(row["material_favorable_count"]) > 0
        and int(row["material_adverse_count"]) == 0
        and int(row["chemical_risk_count"]) == 0
    ]
    if len(selected) != POOL_SIZE:
        raise PropertyAffinityReviewError(
            f"Expected {POOL_SIZE} property review candidates, found {len(selected)}"
        )

    output: list[dict[str, object]] = []
    for evidence in selected:
        candidate_id = evidence["candidate_id"]
        candidate = candidates.get(candidate_id)
        if candidate is None:
            raise PropertyAffinityReviewError(f"Missing unified candidate {candidate_id}")
        position = int(candidate["sequence_index_1based"])
        mapping = experimental_mapping.get(position)
        if mapping is None or mapping.get("coordinate_evaluable") != "True":
            raise PropertyAffinityReviewError(
                f"{candidate_id} lacks evaluable experimental coordinates"
            )
        if candidate.get("design_status") != "eligible_current_round":
            raise PropertyAffinityReviewError(f"{candidate_id} is not released")
        if candidate.get("sequence", "")[-4:] != "SSGS":
            raise PropertyAffinityReviewError(f"{candidate_id} does not preserve SSGS")
        if candidate["wt_residue"] != mapping["residue_aa"]:
            raise PropertyAffinityReviewError(f"WT mismatch for {candidate_id}")

        mutation_label = candidate["mutation_reported_label"]
        output.append(
            {
                "candidate_id": candidate_id,
                "sequence_index_1based": position,
                "wt_residue": candidate["wt_residue"],
                "mutant_residue": candidate["mutant_residue"],
                "mutation_reported_label": mutation_label,
                "mutation_numbering_label": (
                    f"Nb252 IMGT {mapping['numbering_position_label']} "
                    f"{candidate['wt_residue']}>{candidate['mutant_residue']}"
                ),
                "mutation_source_auth_label": (
                    f"Nb252 chain {mapping['auth_asym_id']} auth "
                    f"{candidate['wt_residue']}{mapping['auth_seq_id']}"
                    f"{mapping['insertion_code']}{candidate['mutant_residue']}"
                ),
                "region": candidate["region"],
                "experimental_model_name": mapping["source_model_name"],
                "experimental_auth_asym_id": mapping["auth_asym_id"],
                "experimental_auth_seq_id": int(mapping["auth_seq_id"]),
                "experimental_insertion_code": mapping["insertion_code"],
                "experimental_coordinate_status": mapping["coordinate_status"],
                "prepared_contact_sensitive": candidate["experimental_interface"],
                "sequence": candidate["sequence"],
                "pilot_selected": candidate_id in PILOT_CANDIDATES,
                "pilot_role": _pilot_role(candidate_id),
                **{field: evidence[field] for field in PROPERTY_FIELDS},
                "candidate_selection_performed": True,
                "selection_rule": (
                    "property_pareto_1_and_material_favorable_ge_1_and_"
                    "material_adverse_0_and_chemical_risk_0"
                ),
            }
        )
    output.sort(key=lambda row: (int(row["sequence_index_1based"]), str(row["mutant_residue"])))
    if len({row["candidate_id"] for row in output}) != POOL_SIZE:
        raise PropertyAffinityReviewError("Review pool candidate IDs are not unique")
    if len({int(row["sequence_index_1based"]) for row in output}) != POSITION_COUNT:
        raise PropertyAffinityReviewError("Review pool does not cover exactly 10 positions")
    if {row["candidate_id"] for row in output if row["pilot_selected"]} != set(PILOT_CANDIDATES):
        raise PropertyAffinityReviewError("Pilot candidate identity mismatch")
    return output


def combine_movable_indices(
    interface_indices: Sequence[int], mutation_neighborhood: Sequence[int]
) -> tuple[int, ...]:
    """Return the deterministic union used for both WT and mutant preparation."""

    combined = sorted({int(value) for value in interface_indices} | {int(value) for value in mutation_neighborhood})
    if not combined or any(value <= 0 for value in combined):
        raise PropertyAffinityReviewError("Movable pose indices must be positive")
    return tuple(combined)


def build_run_gate(
    *,
    run_kind: str,
    declared_candidate_ids: Sequence[str],
    declared_positions: Sequence[int],
    wt_controls: Sequence[Mapping[str, object]],
    paired_rows: Sequence[Mapping[str, object]],
    summaries: Sequence[Mapping[str, object]],
    expected_replicates: int,
) -> dict[str, object]:
    """Validate an unfiltered pilot or complete 30-candidate PyRosetta run."""

    if run_kind not in {"pilot", "full_scan"}:
        raise PropertyAffinityReviewError(f"Unsupported run kind: {run_kind}")
    expected_candidates = 6 if run_kind == "pilot" else POOL_SIZE
    expected_positions = 6 if run_kind == "pilot" else POSITION_COUNT
    blockers: list[str] = []
    if len(declared_candidate_ids) != expected_candidates or len(set(declared_candidate_ids)) != expected_candidates:
        blockers.append("candidate_identity_or_count_mismatch")
    if len(set(int(value) for value in declared_positions)) != expected_positions:
        blockers.append("position_count_mismatch")
    if len(wt_controls) != expected_positions * expected_replicates:
        blockers.append("position_specific_wt_count_mismatch")
    if len(paired_rows) != expected_candidates * expected_replicates:
        blockers.append("paired_row_count_mismatch")
    if len(summaries) != expected_candidates:
        blockers.append("summary_count_mismatch")
    if any(not bool(row.get("mutant_runtime_valid")) for row in paired_rows):
        blockers.append("mutant_runtime_failure")
    if any(
        not all(bool(row.get(field)) for field in ("mapping_pass", "breaks_pass", "disulfide_pass", "finite_metrics"))
        or float(row["dG_separated"]) >= 0
        or float(row["cross_interface_energy"]) >= 0
        for row in wt_controls
    ):
        blockers.append("position_specific_wt_failure")
    return {
        "schema_version": 1,
        "gate_name": "nb252_property_candidate_pyrosetta_affinity_noninferiority",
        "run_kind": run_kind,
        "status": "pass" if not blockers else "blocked",
        "release": (
            "ready_for_full_property_affinity_scan"
            if not blockers and run_kind == "pilot"
            else "ready_for_property_affinity_noninferiority_review"
            if not blockers
            else "blocked"
        ),
        "blockers": blockers,
        "candidate_count": len(summaries),
        "position_count": len(set(int(value) for value in declared_positions)),
        "replicate_count_per_candidate": expected_replicates,
        "wt_control_count": len(wt_controls),
        "paired_row_count": len(paired_rows),
        "candidate_status_counts": dict(Counter(str(row["status"]) for row in summaries)),
        "candidate_filtering_applied_during_scoring": False,
        "score_semantics": "mutant_minus_position_specific_paired_WT_Rosetta_ranking_signal",
    }


def review_completed_scan(
    *,
    summary_rows: Sequence[Mapping[str, str]],
    paired_rows: Sequence[Mapping[str, str]],
    wt_rows: Sequence[Mapping[str, str]],
    pilot_summary_rows: Sequence[Mapping[str, str]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Recompute and classify the completed 30-candidate paired scan.

    Classification uses energy direction only: ``directionally_favorable``
    requires negative medians for both interface-energy signals and at least
    two of three negative replicates for each signal. The symmetric rule marks
    ``directionally_adverse``; all other patterns are ``mixed``. These labels
    are computational evidence categories, not affinity measurements or final
    candidate selection.
    """

    if len(summary_rows) != POOL_SIZE or len(paired_rows) != POOL_SIZE * REPLICATES:
        raise PropertyAffinityReviewError("Completed scan must contain 30 summaries and 90 paired rows")
    if len(wt_rows) != POSITION_COUNT * REPLICATES:
        raise PropertyAffinityReviewError("Completed scan must contain 30 position-specific WT rows")
    summaries = _unique(summary_rows, "candidate_id")
    wt_by_id = _unique(wt_rows, "wt_control_id")
    grouped: dict[str, list[Mapping[str, str]]] = {}
    for row in paired_rows:
        grouped.setdefault(row["candidate_id"], []).append(row)
    if set(grouped) != set(summaries):
        raise PropertyAffinityReviewError("Paired and summary candidate identities differ")

    reviewed: list[dict[str, object]] = []
    for candidate_id, summary in summaries.items():
        replicates = sorted(grouped[candidate_id], key=lambda row: int(row["replicate"]))
        if len(replicates) != REPLICATES or {int(row["replicate"]) for row in replicates} != {1, 2, 3}:
            raise PropertyAffinityReviewError(f"{candidate_id} does not have three unique replicates")
        delta_dg: list[float] = []
        delta_cross: list[float] = []
        lost_receptor: set[int] = set()
        for row in replicates:
            wt = wt_by_id.get(row["wt_control_id"])
            if wt is None:
                raise PropertyAffinityReviewError(f"Missing paired WT for {candidate_id}")
            expected_dg = float(row["mutant_dG_separated"]) - float(wt["dG_separated"])
            expected_cross = float(row["mutant_cross_interface_energy"]) - float(wt["cross_interface_energy"])
            if abs(expected_dg - float(row["delta_dG_separated"])) > 1e-9 or abs(expected_cross - float(row["delta_cross_interface_energy"])) > 1e-9:
                raise PropertyAffinityReviewError(f"Stored paired delta mismatch for {candidate_id}")
            delta_dg.append(expected_dg)
            delta_cross.append(expected_cross)
            wt_contacts = _int_set(wt["receptor_contact_auth_positions"])
            mutant_contacts = _int_set(row["mutant_receptor_contact_auth_positions"])
            lost_receptor.update(wt_contacts - mutant_contacts)
        median_dg = statistics.median(delta_dg)
        median_cross = statistics.median(delta_cross)
        if abs(median_dg - float(summary["delta_dG_separated_median"])) > 1e-9 or abs(median_cross - float(summary["delta_cross_interface_energy_median"])) > 1e-9:
            raise PropertyAffinityReviewError(f"Stored summary median mismatch for {candidate_id}")
        negative_dg = sum(value < 0 for value in delta_dg)
        negative_cross = sum(value < 0 for value in delta_cross)
        both_negative = sum(dg < 0 and cross < 0 for dg, cross in zip(delta_dg, delta_cross, strict=True))
        if median_dg < 0 and median_cross < 0 and negative_dg >= 2 and negative_cross >= 2:
            direction = "directionally_favorable"
        elif median_dg > 0 and median_cross > 0 and negative_dg <= 1 and negative_cross <= 1:
            direction = "directionally_adverse"
        else:
            direction = "mixed"
        antifold = float(summary["experimental_complex_context_delta_log_probability"])
        intersection = (
            "rossetta_favorable_antifold_positive"
            if direction == "directionally_favorable" and antifold > 0
            else "rossetta_favorable_antifold_nonpositive"
            if direction == "directionally_favorable"
            else "not_rossetta_directionally_favorable"
        )
        reviewed.append(
            {
                **dict(summary),
                "short_mutation": f"{summary['wt_residue']}{summary['sequence_index_1based']}{summary['mutant_residue']}",
                "delta_dg_negative_replicate_count": negative_dg,
                "delta_cross_negative_replicate_count": negative_cross,
                "both_energy_negative_same_replicate_count": both_negative,
                "affinity_direction_class": direction,
                "all_three_replicates_both_energy_negative": both_negative == REPLICATES,
                "paired_contact_status": "preserved_all" if not lost_receptor else "one_or_more_receptor_contacts_not_retained",
                "lost_receptor_auth_positions": ";".join(map(str, sorted(lost_receptor))),
                "antifold_complex_direction": "positive" if antifold > 0 else "negative" if antifold < 0 else "zero",
                "multitool_intersection_class": intersection,
                "scientific_selection_performed": False,
            }
        )
    reviewed.sort(key=lambda row: (int(row["sequence_index_1based"]), str(row["mutant_residue"])))

    pilot = _unique(pilot_summary_rows, "candidate_id")
    comparison_fields = (
        "delta_dG_separated_median",
        "delta_cross_interface_energy_median",
        "delta_interface_fa_rep_median",
        "minimum_candidate_vs_paired_wt_vhh_contact_retention",
        "minimum_candidate_vs_paired_wt_receptor_epitope_retention",
        "maximum_interface_ca_rmsd",
    )
    if set(pilot) != set(PILOT_CANDIDATES):
        raise PropertyAffinityReviewError("Pilot summary identities differ from the fixed pilot")
    pilot_full_max_difference = max(
        abs(float(pilot[candidate_id][field]) - float(summaries[candidate_id][field]))
        for candidate_id in pilot
        for field in comparison_fields
    )
    counts = Counter(str(row["affinity_direction_class"]) for row in reviewed)
    intersection_counts = Counter(str(row["multitool_intersection_class"]) for row in reviewed)
    facts = {
        "candidate_count": len(reviewed),
        "paired_row_count": len(paired_rows),
        "wt_control_count": len(wt_rows),
        "direction_class_counts": dict(counts),
        "multitool_intersection_counts": dict(intersection_counts),
        "all_three_both_energy_negative_count": sum(bool(row["all_three_replicates_both_energy_negative"]) for row in reviewed),
        "paired_contact_preserved_all_count": sum(row["paired_contact_status"] == "preserved_all" for row in reviewed),
        "minimum_paired_vhh_contact_retention": min(float(row["minimum_candidate_vs_paired_wt_vhh_contact_retention"]) for row in reviewed),
        "minimum_paired_receptor_contact_retention": min(float(row["minimum_candidate_vs_paired_wt_receptor_epitope_retention"]) for row in reviewed),
        "maximum_interface_ca_rmsd_angstrom": max(float(row["maximum_interface_ca_rmsd"]) for row in reviewed),
        "pilot_full_max_absolute_difference": pilot_full_max_difference,
        "scientific_selection_performed": False,
    }
    return reviewed, facts


def _unique(rows: Sequence[Mapping[str, str]], key: str) -> dict[str, Mapping[str, str]]:
    result: dict[str, Mapping[str, str]] = {}
    for row in rows:
        value = row[key]
        if value in result:
            raise PropertyAffinityReviewError(f"Duplicate {key}: {value}")
        result[value] = row
    return result


def _int_set(value: object) -> set[int]:
    return {int(item) for item in str(value).split(";") if item}


def _pilot_role(candidate_id: str) -> str:
    roles = {
        PILOT_CANDIDATES[0]: "N_terminal_charge_change",
        PILOT_CANDIDATES[1]: "framework_proline_stress_case",
        PILOT_CANDIDATES[2]: "framework_charge_introduction",
        PILOT_CANDIDATES[3]: "CDR1_charge_introduction",
        PILOT_CANDIDATES[4]: "CDR2_flexibility_change",
        PILOT_CANDIDATES[5]: "framework_charge_reversal",
    }
    return roles.get(candidate_id, "not_in_pilot")
