"""Contracts for the Nb252 property-candidate affinity non-inferiority review.

The module selects the fixed 30-member review pool from already computed
multi-tool evidence and validates paired PyRosetta runs.  It does not predict
expression, measured stability, or affinity, and it never filters candidates
from PyRosetta scores during the scan.
"""

from __future__ import annotations

from collections import Counter
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


def _unique(rows: Sequence[Mapping[str, str]], key: str) -> dict[str, Mapping[str, str]]:
    result: dict[str, Mapping[str, str]] = {}
    for row in rows:
        value = row[key]
        if value in result:
            raise PropertyAffinityReviewError(f"Duplicate {key}: {value}")
        result[value] = row
    return result


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
