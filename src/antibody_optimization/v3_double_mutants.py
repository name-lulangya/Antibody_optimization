"""Build and analyze the complete V3 Nb252 double-mutant space.

The active V3 route combines only the 15 released parent single mutants.  All
unordered pairs are retained except alternatives at the same reported-sequence
position.  Stable-word occurrences and sequence liabilities are recomputed on
each complete 128-aa double sequence; they are never inferred by adding the two
single-mutant annotations.

AntiFold remains a negative single-position compatibility check.  This module
preserves the two constituent records but deliberately does not add their
``delta_logp`` values or invent a double-mutant AntiFold score.  NetSolP and
NanoMelt non-additivity residuals are calculated only after complete double
sequences have been scored and describe model output, not physical epistasis.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from math import comb, isfinite
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from Bio.PDB import MMCIFParser

from .expression_panel_selection import MAGNITUDE_THRESHOLDS, classify_change
from .stable_words import compare_stable_word_occurrences
from .unified_single_mutants import sequence_liability_deltas


EXPECTED_PARENT_COUNT = 15
EXPECTED_PARENT_POSITION_COUNT = 12
EXPECTED_THEORETICAL_PAIR_COUNT = 105
EXPECTED_INVALID_PAIR_COUNT = 3
EXPECTED_VALID_DOUBLE_COUNT = 102
EXPECTED_SCORE_SAMPLE_COUNT = 103
WT_SCORE_ID = "Nb252_v3_WT"
HYDROPHOBIC = frozenset("AVILMFWY")
SIDECHAIN_CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}


class V3DoubleMutantError(ValueError):
    """Raised when parent, pair, structure, or score identities disagree."""


def build_v3_double_mutant_space(
    selected_rows: Sequence[Mapping[str, object]],
    decision_rows: Sequence[Mapping[str, object]],
    stable_words: Sequence[str],
    mapping_rows: Sequence[Mapping[str, object]],
    experimental_cif: Path,
    af3_cif: Path,
) -> dict[str, object]:
    """Return the complete 102-double plan and its three invalid pairs.

    ``selected_rows`` is the concise released parent panel and ``decision_rows``
    is its 31-row audit, used to recover WT identities and structure evidence.
    Pair geometry is measured on WT coordinates.  Experimental coordinates are
    used when both sites are observed; otherwise both sites are measured in the
    independently predicted AF3 VHH so the distance remains internally
    consistent and explicitly prediction-only.
    """

    parents, parent_sequence = _normalize_parents(selected_rows, decision_rows)
    geometry = _pair_geometry_context(
        parents, mapping_rows, experimental_cif, af3_cif
    )
    candidates: list[dict[str, object]] = []
    invalid_pairs: list[dict[str, object]] = []
    stable_word_changes: list[dict[str, object]] = []
    theoretical_order = 0
    valid_order = 0
    for first, second in combinations(parents, 2):
        theoretical_order += 1
        mutation_a = _mutation_code(first)
        mutation_b = _mutation_code(second)
        position_a = int(first["reported_sequence_index_1based"])
        position_b = int(second["reported_sequence_index_1based"])
        if position_a == position_b:
            invalid_pairs.append(
                {
                    "theoretical_pair_order": theoretical_order,
                    "parent_a_candidate_id": first["candidate_id"],
                    "parent_b_candidate_id": second["candidate_id"],
                    "mutation_a": mutation_a,
                    "mutation_b": mutation_b,
                    "reported_sequence_index_1based": position_a,
                    "pair_status": "invalid_same_position",
                    "exclusion_reason": (
                        "same_position_alternative_substitutions_are_mutually_exclusive"
                    ),
                }
            )
            continue
        valid_order += 1
        sequence_list = list(parent_sequence)
        sequence_list[position_a - 1] = str(first["mutant_residue"])
        sequence_list[position_b - 1] = str(second["mutant_residue"])
        sequence = "".join(sequence_list)
        expected_positions = sorted((position_a, position_b))
        observed_positions = [
            index
            for index, (wt, mutant) in enumerate(
                zip(parent_sequence, sequence, strict=True), 1
            )
            if wt != mutant
        ]
        if observed_positions != expected_positions:
            raise V3DoubleMutantError(
                f"Combined sequence identity failed for {mutation_a}+{mutation_b}"
            )
        candidate_id = f"Nb252_v3_double_{mutation_a}__{mutation_b}"
        word_result = compare_stable_word_occurrences(
            parent_sequence, sequence, stable_words
        )
        created = word_result.pop("created_occurrences")
        lost = word_result.pop("lost_occurrences")
        for change_type, occurrences in (("created", created), ("lost", lost)):
            for occurrence in occurrences:
                stable_word_changes.append(
                    {
                        "double_candidate_id": candidate_id,
                        "mutation_set": f"{mutation_a};{mutation_b}",
                        "change_type": change_type,
                        **occurrence,
                    }
                )
        risks = sequence_risk_features(
            parent_sequence, sequence, expected_positions
        )
        pair_geometry = geometry[(min(position_a, position_b), max(position_a, position_b))]
        contains_t99f = mutation_a == "T99F" or mutation_b == "T99F"
        detailed_triggers: list[str] = []
        if contains_t99f:
            detailed_triggers.append("contains_user_directed_T99F_exploration")
        if pair_geometry["pair_experimental_coordinate_status"] != "both_observed":
            detailed_triggers.append("includes_experimental_missing_coordinate_site")
        if float(pair_geometry["pair_minimum_heavy_atom_distance_a"]) < 4.5:
            detailed_triggers.append("sites_share_direct_local_4p5A_neighborhood")
        if any(
            str(parent["near_interface_shell_status"])
            in {
                "within_4A_of_hard_interface_residue",
                "borderline_4_to_4p5A_from_hard_interface_residue",
            }
            for parent in (first, second)
        ):
            detailed_triggers.append("includes_hard_interface_shell_site")
        if int(risks["hard_sequence_risk_count"]) > 0:
            detailed_triggers.append("combined_sequence_hard_risk_flag")
        candidates.append(
            {
                "v3_double_plan_order_not_efficacy_rank": valid_order,
                "theoretical_pair_order": theoretical_order,
                "double_candidate_id": candidate_id,
                "mutation_set": f"{mutation_a};{mutation_b}",
                "mutation_a": mutation_a,
                "mutation_b": mutation_b,
                "parent_a_candidate_id": first["candidate_id"],
                "parent_b_candidate_id": second["candidate_id"],
                "parent_a_panel_order_not_efficacy_rank": first[
                    "v3_parent_panel_order_not_efficacy_rank"
                ],
                "parent_b_panel_order_not_efficacy_rank": second[
                    "v3_parent_panel_order_not_efficacy_rank"
                ],
                "position_a_reported_1based": position_a,
                "position_b_reported_1based": position_b,
                "region_a": first["region"],
                "region_b": second["region"],
                "sequence": sequence,
                "sequence_length_aa": len(sequence),
                "contains_t99f_stable_word_exploration_parent": contains_t99f,
                "parent_a_netsolp_delta_u": _finite(first, "netsolp_delta_u"),
                "parent_b_netsolp_delta_u": _finite(second, "netsolp_delta_u"),
                "parent_a_netsolp_delta_s": _finite(first, "netsolp_delta_s"),
                "parent_b_netsolp_delta_s": _finite(second, "netsolp_delta_s"),
                "parent_a_nanomelt_delta_tm_c": _finite(
                    first, "nanomelt_delta_tm_c"
                ),
                "parent_b_nanomelt_delta_tm_c": _finite(
                    second, "nanomelt_delta_tm_c"
                ),
                "antifold_component_a_source": first["antifold_selection_source"],
                "antifold_component_b_source": second["antifold_selection_source"],
                "antifold_component_a_delta_logp": _finite(
                    first, "antifold_delta_logp"
                ),
                "antifold_component_b_delta_logp": _finite(
                    second, "antifold_delta_logp"
                ),
                "antifold_component_a_rank_worst_first": first[
                    "antifold_mutant_rank_worst_first"
                ],
                "antifold_component_b_rank_worst_first": second[
                    "antifold_mutant_rank_worst_first"
                ],
                "antifold_component_a_veto_status": first["antifold_veto_status"],
                "antifold_component_b_veto_status": second["antifold_veto_status"],
                "antifold_constituent_gate": "pass",
                "antifold_double_mutant_scored": False,
                "antifold_component_values_combined": False,
                "antifold_double_mutant_score": "",
                **word_result,
                **risks,
                **pair_geometry,
                "machine_structure_triage_status": (
                    "detailed_review_triggered"
                    if detailed_triggers
                    else "routine_context_recorded"
                ),
                "machine_structure_triage_triggers": "|".join(detailed_triggers),
                "candidate_prefiltering_applied": False,
                "remote_property_scoring_status": "not_started",
                "final_double_selection_status": "not_performed",
            }
        )

    _validate_pair_space(parent_sequence, parents, candidates, invalid_pairs)
    facts = {
        "parent_single_count": len(parents),
        "parent_unique_position_count": len(
            {int(row["reported_sequence_index_1based"]) for row in parents}
        ),
        "theoretical_unordered_pair_count": comb(len(parents), 2),
        "invalid_same_position_pair_count": len(invalid_pairs),
        "valid_double_mutant_count": len(candidates),
        "score_sample_count_including_wt": len(candidates) + 1,
        "t99f_containing_valid_double_count": sum(
            bool(row["contains_t99f_stable_word_exploration_parent"])
            for row in candidates
        ),
        "stable_word_gain_candidate_count": sum(
            int(row["net_stable_word_occurrence_delta"]) > 0
            for row in candidates
        ),
        "stable_word_loss_candidate_count": sum(
            int(row["net_stable_word_occurrence_delta"]) < 0
            for row in candidates
        ),
        "hard_sequence_risk_candidate_count": sum(
            int(row["hard_sequence_risk_count"]) > 0 for row in candidates
        ),
        "detailed_structure_review_trigger_count": sum(
            row["machine_structure_triage_status"] == "detailed_review_triggered"
            for row in candidates
        ),
        "structure_distance_source_counts": dict(
            sorted(Counter(row["pair_structure_distance_source"] for row in candidates).items())
        ),
        "candidate_prefiltering_applied": False,
    }
    return {
        "parent_sequence": parent_sequence,
        "parents": parents,
        "candidates": candidates,
        "invalid_pairs": invalid_pairs,
        "stable_word_changes": stable_word_changes,
        "facts": facts,
    }


def build_v3_score_samples(
    parent_sequence: str, candidates: Sequence[Mapping[str, object]]
) -> list[dict[str, object]]:
    """Return the exact WT-plus-102 ID/sequence table shared by both tools."""

    rows: list[dict[str, object]] = [
        {
            "sample_uid": WT_SCORE_ID,
            "double_candidate_id": "WT",
            "sequence_raw": parent_sequence,
            "is_wt_control": True,
        }
    ]
    rows.extend(
        {
            "sample_uid": row["double_candidate_id"],
            "double_candidate_id": row["double_candidate_id"],
            "sequence_raw": row["sequence"],
            "is_wt_control": False,
        }
        for row in candidates
    )
    if len(rows) != EXPECTED_SCORE_SAMPLE_COUNT or len(
        {str(row["sample_uid"]) for row in rows}
    ) != EXPECTED_SCORE_SAMPLE_COUNT:
        raise V3DoubleMutantError("Expected WT plus 102 unique scoring samples")
    return rows


def merge_v3_property_scores(
    candidates: Sequence[Mapping[str, object]],
    netsolp_rows: Sequence[Mapping[str, object]],
    nanomelt_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Join complete scores and compute WT deltas and model non-additivity."""

    if len(candidates) != EXPECTED_VALID_DOUBLE_COUNT:
        raise V3DoubleMutantError("Expected exactly 102 double candidates")
    net = _score_index(netsolp_rows, "NetSolP")
    melt = _score_index(nanomelt_rows, "NanoMelt")
    wt_net = net[WT_SCORE_ID]
    wt_melt = melt[WT_SCORE_ID]
    output: list[dict[str, object]] = []
    for source in candidates:
        row = dict(source)
        identifier = str(row["double_candidate_id"])
        sequence = str(row["sequence"])
        if identifier not in net or identifier not in melt:
            raise V3DoubleMutantError(f"Missing complete-sequence score: {identifier}")
        net_row, melt_row = net[identifier], melt[identifier]
        if (
            str(net_row["sequence_raw"]) != sequence
            or str(melt_row["sequence_raw"]) != sequence
        ):
            raise V3DoubleMutantError(f"Score sequence mismatch: {identifier}")
        values = {
            "netsolp_u": _finite(net_row, "predicted_usability")
            - _finite(wt_net, "predicted_usability"),
            "netsolp_s": _finite(net_row, "predicted_solubility")
            - _finite(wt_net, "predicted_solubility"),
            "nanomelt_tm_c": _finite(
                melt_row, "nanomelt_predicted_apparent_tm_c"
            )
            - _finite(wt_melt, "nanomelt_predicted_apparent_tm_c"),
        }
        row.update(
            {
                "netsolp_predicted_usability": _finite(
                    net_row, "predicted_usability"
                ),
                "netsolp_predicted_solubility": _finite(
                    net_row, "predicted_solubility"
                ),
                "nanomelt_predicted_apparent_tm_c": _finite(
                    melt_row, "nanomelt_predicted_apparent_tm_c"
                ),
                "nanomelt_scoring_status": melt_row["scoring_status"],
                "nanomelt_scored_length_aa": melt_row.get("scored_length_aa", ""),
                "nanomelt_trimmed_c_terminal": melt_row.get(
                    "trimmed_c_terminal", ""
                ),
            }
        )
        parent_fields = {
            "netsolp_u": ("parent_a_netsolp_delta_u", "parent_b_netsolp_delta_u"),
            "netsolp_s": ("parent_a_netsolp_delta_s", "parent_b_netsolp_delta_s"),
            "nanomelt_tm_c": (
                "parent_a_nanomelt_delta_tm_c",
                "parent_b_nanomelt_delta_tm_c",
            ),
        }
        for metric, delta in values.items():
            row[f"{metric}_delta_vs_wt"] = delta
            row[f"{metric}_magnitude_band"] = classify_change(
                delta, MAGNITUDE_THRESHOLDS[metric]
            )
            first_field, second_field = parent_fields[metric]
            row[f"{metric}_model_nonadditivity_residual"] = (
                delta - _finite(row, first_field) - _finite(row, second_field)
            )
        row["model_nonadditivity_interpretation"] = (
            "predictor_output_residual_not_physical_epistasis"
        )
        row["remote_property_scoring_status"] = "complete"
        row["final_double_selection_status"] = "not_performed"
        output.append(row)
    return output


def sequence_risk_features(
    parent_sequence: str, mutant_sequence: str, positions: Sequence[int]
) -> dict[str, object]:
    """Recompute combined sequence liabilities and seven-residue windows."""

    liabilities = sequence_liability_deltas(parent_sequence, mutant_sequence)
    hard_flags: list[str] = []
    for position in positions:
        if (
            mutant_sequence[position - 1] == "P"
            and parent_sequence[position - 1] != "P"
        ):
            hard_flags.append("new_proline_backbone_constraint")
    parent_hydrophobic = _max_window_count(parent_sequence, HYDROPHOBIC)
    mutant_hydrophobic = _max_window_count(mutant_sequence, HYDROPHOBIC)
    parent_charge = _max_local_charge(parent_sequence)
    mutant_charge = _max_local_charge(mutant_sequence)
    if mutant_hydrophobic >= 6 and mutant_hydrophobic > parent_hydrophobic:
        hard_flags.append("new_dense_local_hydrophobic_window")
    if mutant_charge >= 5 and mutant_charge > parent_charge:
        hard_flags.append("new_extreme_local_charge_cluster")
    soft_flags = [
        value for value in str(liabilities["new_liability_flags"]).split("|") if value
    ]
    hard_unique = sorted(set(hard_flags))
    return {
        "hard_sequence_risk_flags": "|".join(hard_unique),
        "hard_sequence_risk_count": len(hard_unique),
        "soft_sequence_risk_flags": "|".join(soft_flags),
        "soft_sequence_risk_count": len(soft_flags),
        "max_local_hydrophobic_count_parent": parent_hydrophobic,
        "max_local_hydrophobic_count_mutant": mutant_hydrophobic,
        "max_abs_local_charge_parent": parent_charge,
        "max_abs_local_charge_mutant": mutant_charge,
        **liabilities,
    }


def _normalize_parents(
    selected_rows: Sequence[Mapping[str, object]],
    decision_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], str]:
    if len(selected_rows) != EXPECTED_PARENT_COUNT:
        raise V3DoubleMutantError("Expected exactly 15 released parent singles")
    decision_lookup = {
        str(row["candidate_id"]): dict(row)
        for row in decision_rows
        if str(row.get("v3_parent_selection_status", "")) == "selected"
    }
    if len(decision_lookup) != EXPECTED_PARENT_COUNT:
        raise V3DoubleMutantError("Decision audit must contain exactly 15 selected rows")
    parents: list[dict[str, object]] = []
    reconstructed: set[str] = set()
    seen_ids: set[str] = set()
    for concise in selected_rows:
        identifier = str(concise["candidate_id"])
        if identifier in seen_ids or identifier not in decision_lookup:
            raise V3DoubleMutantError(f"Parent identity mismatch: {identifier}")
        seen_ids.add(identifier)
        source = decision_lookup[identifier]
        if str(concise["sequence"]) != str(source["sequence"]):
            raise V3DoubleMutantError(f"Parent sequence mismatch: {identifier}")
        for field in (
            "netsolp_delta_u",
            "netsolp_delta_s",
            "nanomelt_delta_tm_c",
            "antifold_delta_logp",
        ):
            if _finite(concise, field) != _finite(source, field):
                raise V3DoubleMutantError(
                    f"Parent evidence mismatch for {identifier}: {field}"
                )
        if str(source["antifold_veto_status"]) != "pass":
            raise V3DoubleMutantError(
                f"Parent failed the frozen AntiFold veto: {identifier}"
            )
        position = int(source["reported_sequence_index_1based"])
        sequence = str(source["sequence"])
        wt = str(source["wt_residue"])
        mutant = str(source["mutant_residue"])
        if len(sequence) != 128 or sequence[position - 1] != mutant or wt == mutant:
            raise V3DoubleMutantError(f"Invalid released parent sequence: {identifier}")
        reconstructed.add(sequence[: position - 1] + wt + sequence[position:])
        parents.append(source)
    if len(reconstructed) != 1:
        raise V3DoubleMutantError("Released singles do not reconstruct one parent")
    parent_sequence = reconstructed.pop()
    if (
        not parent_sequence.endswith("SSGS")
        or parent_sequence[21] != "C"
        or parent_sequence[94] != "C"
        or parent_sequence.count("C") != 2
    ):
        raise V3DoubleMutantError("Authoritative Cys/SSGS parent contract failed")
    parents.sort(key=lambda row: int(row["v3_parent_panel_order_not_efficacy_rank"]))
    if [int(row["v3_parent_panel_order_not_efficacy_rank"]) for row in parents] != list(
        range(1, EXPECTED_PARENT_COUNT + 1)
    ):
        raise V3DoubleMutantError("Parent display orders must be exactly 1..15")
    if len({int(row["reported_sequence_index_1based"]) for row in parents}) != 12:
        raise V3DoubleMutantError("Expected the released parents to cover 12 positions")
    return parents, parent_sequence


def _pair_geometry_context(
    parents: Sequence[Mapping[str, object]],
    mapping_rows: Sequence[Mapping[str, object]],
    experimental_cif: Path,
    af3_cif: Path,
) -> dict[tuple[int, int], dict[str, object]]:
    positions = sorted(
        {int(row["reported_sequence_index_1based"]) for row in parents}
    )
    map_lookup: dict[tuple[str, int], Mapping[str, object]] = {}
    for row in mapping_rows:
        model = str(row.get("source_model_name", ""))
        if model not in {"NK2R-252.pdb", "fold_2r_252_nomg_model_0.cif"}:
            continue
        position = int(row["sequence_index_1based"])
        if position in positions:
            map_lookup[(model, position)] = row
    parser = MMCIFParser(QUIET=True)
    experimental = parser.get_structure("experimental", str(experimental_cif))[0]
    af3 = parser.get_structure("af3", str(af3_cif))[0]
    coordinates: dict[tuple[str, int], object] = {}
    for model_name, structure in (
        ("NK2R-252.pdb", experimental),
        ("fold_2r_252_nomg_model_0.cif", af3),
    ):
        for position in positions:
            mapping = map_lookup.get((model_name, position))
            if mapping is None:
                raise V3DoubleMutantError(
                    f"Structure mapping missing {model_name} reported position {position}"
                )
            if str(mapping["coordinate_status"]) != "observed":
                continue
            chain_id = str(mapping["auth_asym_id"])
            auth_seq_id = int(mapping["auth_seq_id"])
            insertion = str(mapping.get("insertion_code", "")).strip()
            residue = _find_residue(structure[chain_id], auth_seq_id, insertion)
            if str(mapping["structure_residue_aa"]) != _residue_aa(residue):
                raise V3DoubleMutantError(
                    f"Structure residue identity mismatch at {model_name}:{position}"
                )
            coordinates[(model_name, position)] = residue
    output: dict[tuple[int, int], dict[str, object]] = {}
    for first, second in combinations(positions, 2):
        both_experimental = all(
            ("NK2R-252.pdb", position) in coordinates
            for position in (first, second)
        )
        af3_a = coordinates[("fold_2r_252_nomg_model_0.cif", first)]
        af3_b = coordinates[("fold_2r_252_nomg_model_0.cif", second)]
        af3_ca_distance = float(np.linalg.norm(af3_a["CA"].coord - af3_b["CA"].coord))
        af3_minimum = _minimum_heavy_atom_distance(af3_a, af3_b)
        experimental_ca_distance: float | str = ""
        experimental_minimum: float | str = ""
        if both_experimental:
            experimental_a = coordinates[("NK2R-252.pdb", first)]
            experimental_b = coordinates[("NK2R-252.pdb", second)]
            experimental_ca_distance = float(
                np.linalg.norm(experimental_a["CA"].coord - experimental_b["CA"].coord)
            )
            experimental_minimum = _minimum_heavy_atom_distance(
                experimental_a, experimental_b
            )
            residue_a, residue_b = experimental_a, experimental_b
            ca_distance = experimental_ca_distance
            minimum = experimental_minimum
        else:
            residue_a, residue_b = af3_a, af3_b
            ca_distance = af3_ca_distance
            minimum = af3_minimum
        source = (
            "experimental_complex_vhh"
            if both_experimental
            else "af3_vhh_only_due_missing_experimental_coordinate"
        )
        neighbors_a = _neighbor_keys(residue_a, residue_a.get_parent())
        neighbors_b = _neighbor_keys(residue_b, residue_b.get_parent())
        shared = sorted(neighbors_a & neighbors_b)
        if minimum < 4.5:
            spatial_class = "direct_local_neighborhood_under_4p5A"
        elif ca_distance < 10.0:
            spatial_class = "nearby_ca_under_10A"
        else:
            spatial_class = "spatially_separated_ca_at_least_10A"
        output[(first, second)] = {
            "pair_experimental_coordinate_status": (
                "both_observed" if both_experimental else "includes_missing_coordinates"
            ),
            "experimental_pair_ca_distance_a": (
                round(float(experimental_ca_distance), 5)
                if experimental_ca_distance != ""
                else ""
            ),
            "experimental_pair_minimum_heavy_atom_distance_a": (
                round(float(experimental_minimum), 5)
                if experimental_minimum != ""
                else ""
            ),
            "af3_pair_ca_distance_a": round(af3_ca_distance, 5),
            "af3_pair_minimum_heavy_atom_distance_a": round(af3_minimum, 5),
            "pair_structure_distance_source": source,
            "pair_ca_distance_a": round(ca_distance, 5),
            "pair_minimum_heavy_atom_distance_a": round(minimum, 5),
            "pair_spatial_class": spatial_class,
            "pair_shared_local_neighbor_count": len(shared),
            "pair_shared_local_neighbors": ";".join(shared),
            "pair_geometry_role": "wt_structure_triage_not_mutant_effect_prediction",
        }
    return output


def _validate_pair_space(parent, parents, candidates, invalid_pairs) -> None:
    if (
        len(parents) != EXPECTED_PARENT_COUNT
        or comb(len(parents), 2) != EXPECTED_THEORETICAL_PAIR_COUNT
        or len(candidates) != EXPECTED_VALID_DOUBLE_COUNT
        or len(invalid_pairs) != EXPECTED_INVALID_PAIR_COUNT
    ):
        raise V3DoubleMutantError("V3 pair counts do not match 15/105/3/102")
    if len({str(row["double_candidate_id"]) for row in candidates}) != len(candidates):
        raise V3DoubleMutantError("Double candidate IDs are not unique")
    if len({str(row["sequence"]) for row in candidates}) != len(candidates):
        raise V3DoubleMutantError("Double candidate sequences are not unique")
    if sum(bool(row["contains_t99f_stable_word_exploration_parent"]) for row in candidates) != 14:
        raise V3DoubleMutantError("T99F must occur in exactly 14 valid pairs")
    for row in candidates:
        sequence = str(row["sequence"])
        if len(sequence) != 128 or not sequence.endswith("SSGS") or sequence.count("C") != 2:
            raise V3DoubleMutantError(
                f"Double sequence violates construct constraints: {row['double_candidate_id']}"
            )
        differences = [
            index
            for index, (wt, mutant) in enumerate(zip(parent, sequence, strict=True), 1)
            if wt != mutant
        ]
        if differences != sorted(
            (
                int(row["position_a_reported_1based"]),
                int(row["position_b_reported_1based"]),
            )
        ):
            raise V3DoubleMutantError(
                f"Candidate is not its declared double: {row['double_candidate_id']}"
            )
        if row["antifold_constituent_gate"] != "pass" or row[
            "antifold_component_values_combined"
        ] is not False:
            raise V3DoubleMutantError("AntiFold constituent-only contract failed")


def _score_index(rows, tool):
    if len(rows) != EXPECTED_SCORE_SAMPLE_COUNT:
        raise V3DoubleMutantError(f"{tool} must contain 103 rows")
    indexed: dict[str, Mapping[str, object]] = {}
    for row in rows:
        identifier = str(row["sample_uid"])
        if identifier in indexed or str(row.get("scoring_status")) != "pass":
            raise V3DoubleMutantError(f"{tool} score coverage is incomplete")
        indexed[identifier] = row
    if len(indexed) != EXPECTED_SCORE_SAMPLE_COUNT or WT_SCORE_ID not in indexed:
        raise V3DoubleMutantError(f"{tool} score identity set is incomplete")
    return indexed


def _mutation_code(row: Mapping[str, object]) -> str:
    return (
        f"{row['wt_residue']}{int(row['reported_sequence_index_1based'])}"
        f"{row['mutant_residue']}"
    )


def _finite(row: Mapping[str, object], field: str) -> float:
    value = float(row[field])
    if not isfinite(value):
        raise V3DoubleMutantError(f"Non-finite {field}")
    return value


def _max_window_count(sequence: str, residues: frozenset[str]) -> int:
    return max(
        sum(aa in residues for aa in sequence[start : start + 7])
        for start in range(len(sequence) - 6)
    )


def _max_local_charge(sequence: str) -> int:
    return max(
        abs(
            sum(
                SIDECHAIN_CHARGE.get(aa, 0)
                for aa in sequence[start : start + 7]
            )
        )
        for start in range(len(sequence) - 6)
    )


def _find_residue(chain, auth_seq_id: int, insertion_code: str):
    wanted_insertion = insertion_code or " "
    for residue in chain.get_residues():
        if residue.id[1] == auth_seq_id and residue.id[2].strip() == wanted_insertion.strip():
            return residue
    raise V3DoubleMutantError(
        f"Cannot find structure residue {chain.id}:{auth_seq_id}{insertion_code}"
    )


def _residue_aa(residue) -> str:
    mapping = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    try:
        return mapping[residue.resname]
    except KeyError as error:
        raise V3DoubleMutantError(f"Unsupported residue {residue.resname}") from error


def _heavy_atoms(residue):
    return [atom for atom in residue.get_atoms() if atom.element not in {"H", "D"}]


def _minimum_heavy_atom_distance(first, second) -> float:
    return min(
        float(np.linalg.norm(atom_a.coord - atom_b.coord))
        for atom_a in _heavy_atoms(first)
        for atom_b in _heavy_atoms(second)
    )


def _neighbor_keys(residue, chain) -> set[str]:
    output: set[str] = set()
    for other in chain.get_residues():
        if other is residue or other.id[0] != " ":
            continue
        if _minimum_heavy_atom_distance(residue, other) < 4.5:
            output.add(f"{_residue_aa(other)}{other.id[1]}{other.id[2].strip()}")
    return output


__all__ = [
    "EXPECTED_INVALID_PAIR_COUNT",
    "EXPECTED_PARENT_COUNT",
    "EXPECTED_SCORE_SAMPLE_COUNT",
    "EXPECTED_VALID_DOUBLE_COUNT",
    "V3DoubleMutantError",
    "WT_SCORE_ID",
    "build_v3_double_mutant_space",
    "build_v3_score_samples",
    "merge_v3_property_scores",
    "sequence_risk_features",
]
