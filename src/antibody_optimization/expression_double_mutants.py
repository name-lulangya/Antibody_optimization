"""Enumerate and summarize the active Nb252 19-to-162 double-mutant route.

The module combines only the 19 explicitly released single-mutant parents.
Same-position alternatives are mutually exclusive.  It recomputes sequence
liabilities and degenerate stable-word occurrences on each complete double
sequence.  AntiFold evidence remains position-wise fixed-backbone evidence;
values from experimental and AF3 views are never added across views.
"""

from __future__ import annotations

from itertools import combinations
from math import comb, isfinite
from typing import Mapping, Sequence

from .expression_panel_selection import MAGNITUDE_THRESHOLDS, classify_change
from .stable_words import compare_stable_word_occurrences
from .unified_single_mutants import sequence_liability_deltas


EXPECTED_PARENTS = 19
EXPECTED_VALID_DOUBLES = 162
EXPECTED_INVALID_PAIRS = 9
WT_SCORE_ID = "LTT__Nb252__WT"
HYDROPHOBIC = frozenset("AVILMFWY")
SIDECHAIN_CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
ADVERSE_ORDER = {
    "strong_adverse": 4,
    "moderate_adverse": 3,
    "weak_adverse": 2,
    "negligible": 1,
    "weak_favorable": 0,
    "moderate_favorable": -1,
    "strong_favorable": -2,
}


class ExpressionDoubleMutantError(ValueError):
    """Raised when active double-mutant inputs or score joins disagree."""


def build_double_mutant_space(
    parent_rows: Sequence[Mapping[str, object]],
    parent_gate: Mapping[str, object],
    stable_words: Sequence[str],
) -> dict[str, object]:
    """Return all 162 distinct-position double mutants and 9 invalid pairs."""

    if (
        parent_gate.get("status") != "pass"
        or parent_gate.get("release")
        != "parent_19_single_mutants_ready_for_double_enumeration"
    ):
        raise ExpressionDoubleMutantError("The 19-parent release gate is not open")
    if len(parent_rows) != EXPECTED_PARENTS:
        raise ExpressionDoubleMutantError("Expected exactly 19 parent single mutants")

    normalized: list[dict[str, object]] = []
    reconstructed: set[str] = set()
    for source in parent_rows:
        row = dict(source)
        if str(row.get("approved_as_double_mutant_parent")).lower() != "true":
            raise ExpressionDoubleMutantError("Every input row must be an approved parent")
        position = int(row["reported_sequence_index_1based"])
        wt, mutant, sequence = (
            str(row["wt_residue"]),
            str(row["mutant_residue"]),
            str(row["sequence"]),
        )
        if len(sequence) != 128 or sequence[position - 1] != mutant or wt == mutant:
            raise ExpressionDoubleMutantError(f"Invalid parent sequence: {row['candidate_id']}")
        reconstructed.add(sequence[: position - 1] + wt + sequence[position:])
        normalized.append(row)
    if len(reconstructed) != 1:
        raise ExpressionDoubleMutantError("The 19 singles do not reconstruct one parent")
    parent = reconstructed.pop()
    if not parent.endswith("SSGS") or parent[21] != "C" or parent[94] != "C":
        raise ExpressionDoubleMutantError("Authoritative Cys/terminal constraints failed")

    normalized.sort(key=lambda row: int(row["parent19_selection_order"]))
    candidates: list[dict[str, object]] = []
    invalid: list[dict[str, object]] = []
    word_changes: list[dict[str, object]] = []
    for first, second in combinations(normalized, 2):
        pos_a = int(first["reported_sequence_index_1based"])
        pos_b = int(second["reported_sequence_index_1based"])
        mutation_a = _mutation_code(first)
        mutation_b = _mutation_code(second)
        if pos_a == pos_b:
            invalid.append(
                {
                    "mutation_a": mutation_a,
                    "mutation_b": mutation_b,
                    "reported_sequence_index_1based": pos_a,
                    "exclusion_reason": "same_position_alternative_substitutions_are_mutually_exclusive",
                }
            )
            continue
        sequence_list = list(parent)
        sequence_list[pos_a - 1] = str(first["mutant_residue"])
        sequence_list[pos_b - 1] = str(second["mutant_residue"])
        sequence = "".join(sequence_list)
        differences = [
            index
            for index, (before, after) in enumerate(zip(parent, sequence, strict=True), 1)
            if before != after
        ]
        if differences != sorted((pos_a, pos_b)):
            raise ExpressionDoubleMutantError("A combined sequence is not the declared double")
        candidate_id = f"Nb252_expr_double_{mutation_a}__{mutation_b}"
        word = compare_stable_word_occurrences(parent, sequence, stable_words)
        created = word.pop("created_occurrences")
        lost = word.pop("lost_occurrences")
        for change_type, occurrences in (("created", created), ("lost", lost)):
            for occurrence in occurrences:
                word_changes.append(
                    {"candidate_id": candidate_id, "change_type": change_type, **occurrence}
                )
        risks = sequence_risks(parent, sequence, (pos_a, pos_b))
        source_a = str(first["antifold_selection_source"])
        source_b = str(second["antifold_selection_source"])
        delta_a = _finite(first, "antifold_selection_delta_log_probability")
        delta_b = _finite(second, "antifold_selection_delta_log_probability")
        same_view = source_a == source_b
        band_a = str(first["antifold_magnitude_band"])
        band_b = str(second["antifold_magnitude_band"])
        candidates.append(
            {
                "candidate_id": candidate_id,
                "mutation_set": f"{mutation_a};{mutation_b}",
                "mutation_a": mutation_a,
                "mutation_b": mutation_b,
                "position_a_reported_1based": pos_a,
                "position_b_reported_1based": pos_b,
                "region_a": str(first["region"]),
                "region_b": str(second["region"]),
                "source_single_candidate_a": str(first["candidate_id"]),
                "source_single_candidate_b": str(second["candidate_id"]),
                "sequence": sequence,
                "sequence_length_aa": len(sequence),
                "antifold_component_a_source": source_a,
                "antifold_component_b_source": source_b,
                "antifold_component_a_delta_log_probability": delta_a,
                "antifold_component_b_delta_log_probability": delta_b,
                "antifold_component_a_magnitude_band": band_a,
                "antifold_component_b_magnitude_band": band_b,
                "antifold_worst_component_band": max(
                    (band_a, band_b), key=lambda label: ADVERSE_ORDER[label]
                ),
                "antifold_same_view_additive_evaluable": same_view,
                "antifold_same_view_additive_delta_log_probability": delta_a + delta_b if same_view else "",
                "antifold_double_structure_rerun": False,
                "antifold_interaction_or_epistasis_evaluated": False,
                **risks,
                **word,
                "candidate_filtering_applied": False,
                "selection_status": "unfiltered_complete_double_space",
            }
        )

    if len(candidates) != EXPECTED_VALID_DOUBLES or len(invalid) != EXPECTED_INVALID_PAIRS:
        raise ExpressionDoubleMutantError(
            f"Unexpected pair counts: {len(candidates)} valid, {len(invalid)} invalid"
        )
    if len({row["candidate_id"] for row in candidates}) != EXPECTED_VALID_DOUBLES:
        raise ExpressionDoubleMutantError("Double candidate IDs are not unique")
    if len({row["sequence"] for row in candidates}) != EXPECTED_VALID_DOUBLES:
        raise ExpressionDoubleMutantError("Double candidate sequences are not unique")
    facts = {
        "parent_single_mutant_count": EXPECTED_PARENTS,
        "parent_position_count": len({int(row["reported_sequence_index_1based"]) for row in normalized}),
        "theoretical_all_pair_count": comb(EXPECTED_PARENTS, 2),
        "invalid_same_position_pair_count": len(invalid),
        "valid_double_mutant_count": len(candidates),
        "score_sample_count_including_wt": len(candidates) + 1,
        "hard_sequence_risk_count": sum(int(row["hard_sequence_risk_count"] > 0) for row in candidates),
        "stable_word_gain_candidate_count": sum(
            int(int(row["net_stable_word_occurrence_delta"]) > 0) for row in candidates
        ),
        "antifold_same_view_additive_evaluable_count": sum(
            int(bool(row["antifold_same_view_additive_evaluable"])) for row in candidates
        ),
        "candidate_filtering_applied": False,
    }
    return {
        "parent_sequence": parent,
        "candidates": candidates,
        "invalid_pairs": invalid,
        "stable_word_changes": word_changes,
        "facts": facts,
    }


def build_score_samples(
    parent_sequence: str, candidates: Sequence[Mapping[str, object]]
) -> list[dict[str, object]]:
    """Build the exact WT-plus-162 sequence order used by both predictors."""

    rows = [
        {
            "sample_uid": WT_SCORE_ID,
            "candidate_id": "WT",
            "sequence_raw": parent_sequence,
            "is_wt_control": True,
        }
    ]
    rows.extend(
        {
            "sample_uid": str(row["candidate_id"]),
            "candidate_id": str(row["candidate_id"]),
            "sequence_raw": str(row["sequence"]),
            "is_wt_control": False,
        }
        for row in candidates
    )
    if len(rows) != 163 or len({row["sample_uid"] for row in rows}) != 163:
        raise ExpressionDoubleMutantError("Expected WT plus 162 unique score samples")
    return rows


def merge_property_scores(
    candidates: Sequence[Mapping[str, object]],
    parent_rows: Sequence[Mapping[str, object]],
    netsolp_rows: Sequence[Mapping[str, object]],
    nanomelt_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Join complete remote scores and calculate WT deltas and interaction residuals."""

    if len(candidates) != 162 or len(netsolp_rows) != 163 or len(nanomelt_rows) != 163:
        raise ExpressionDoubleMutantError("Expected 162 candidates and 163 rows per tool")
    net = _score_index(netsolp_rows, "NetSolP")
    melt = _score_index(nanomelt_rows, "NanoMelt")
    wt_net, wt_melt = net[WT_SCORE_ID], melt[WT_SCORE_ID]
    singles = {str(row["candidate_id"]): row for row in parent_rows}
    output: list[dict[str, object]] = []
    for source in candidates:
        row = dict(source)
        candidate_id, sequence = str(row["candidate_id"]), str(row["sequence"])
        if candidate_id not in net or candidate_id not in melt:
            raise ExpressionDoubleMutantError(f"Missing property score: {candidate_id}")
        nrow, mrow = net[candidate_id], melt[candidate_id]
        if str(nrow["sequence_raw"]) != sequence or str(mrow["sequence_raw"]) != sequence:
            raise ExpressionDoubleMutantError(f"Score sequence mismatch: {candidate_id}")
        single_a = singles[str(row["source_single_candidate_a"])]
        single_b = singles[str(row["source_single_candidate_b"])]
        values = {
            "netsolp_u": _finite(nrow, "predicted_usability") - _finite(wt_net, "predicted_usability"),
            "netsolp_s": _finite(nrow, "predicted_solubility") - _finite(wt_net, "predicted_solubility"),
            "nanomelt_tm_c": _finite(mrow, "nanomelt_predicted_apparent_tm_c")
            - _finite(wt_melt, "nanomelt_predicted_apparent_tm_c"),
        }
        single_fields = {
            "netsolp_u": "netsolp_delta_usability_vs_current_wt",
            "netsolp_s": "netsolp_delta_solubility_vs_current_wt",
            "nanomelt_tm_c": "nanomelt_delta_predicted_apparent_tm_c_vs_current_wt",
        }
        row.update(
            {
                "netsolp_predicted_usability": _finite(nrow, "predicted_usability"),
                "netsolp_predicted_solubility": _finite(nrow, "predicted_solubility"),
                "nanomelt_predicted_apparent_tm_c": _finite(
                    mrow, "nanomelt_predicted_apparent_tm_c"
                ),
                "nanomelt_scoring_status": str(mrow["scoring_status"]),
                "nanomelt_scored_length_aa": mrow.get("scored_length_aa", ""),
                "nanomelt_trimmed_c_terminal": mrow.get("trimmed_c_terminal", ""),
            }
        )
        for metric, value in values.items():
            row[f"{metric}_delta_vs_wt"] = value
            row[f"{metric}_magnitude_band"] = classify_change(
                value, MAGNITUDE_THRESHOLDS[metric]
            )
            field = single_fields[metric]
            row[f"{metric}_interaction_residual"] = value - _finite(single_a, field) - _finite(single_b, field)
        row["candidate_selection_performed"] = False
        output.append(row)
    return output


def sequence_risks(
    parent: str, sequence: str, positions: Sequence[int]
) -> dict[str, object]:
    """Recompute combined sequence risks instead of adding single annotations."""

    liabilities = sequence_liability_deltas(parent, sequence)
    hard: list[str] = []
    for position in positions:
        if sequence[position - 1] == "P" and parent[position - 1] != "P":
            hard.append("new_proline_backbone_constraint")
    parent_hydro = _max_window_count(parent, HYDROPHOBIC)
    mutant_hydro = _max_window_count(sequence, HYDROPHOBIC)
    parent_charge = _max_local_charge(parent)
    mutant_charge = _max_local_charge(sequence)
    if mutant_hydro >= 6 and mutant_hydro > parent_hydro:
        hard.append("new_dense_local_hydrophobic_window")
    if mutant_charge >= 5 and mutant_charge > parent_charge:
        hard.append("new_extreme_local_charge_cluster")
    soft = [token for token in str(liabilities["new_liability_flags"]).split("|") if token]
    return {
        "hard_sequence_risk_flags": "|".join(sorted(set(hard))),
        "hard_sequence_risk_count": len(set(hard)),
        "soft_sequence_risk_flags": "|".join(soft),
        "soft_sequence_risk_count": len(soft),
        "max_local_hydrophobic_count_parent": parent_hydro,
        "max_local_hydrophobic_count_mutant": mutant_hydro,
        "max_abs_local_charge_parent": parent_charge,
        "max_abs_local_charge_mutant": mutant_charge,
        **liabilities,
    }


def _mutation_code(row: Mapping[str, object]) -> str:
    return f"{row['wt_residue']}{int(row['reported_sequence_index_1based'])}{row['mutant_residue']}"


def _max_window_count(sequence: str, residues: frozenset[str]) -> int:
    return max(sum(aa in residues for aa in sequence[start : start + 7]) for start in range(len(sequence) - 6))


def _max_local_charge(sequence: str) -> int:
    return max(
        abs(sum(SIDECHAIN_CHARGE.get(aa, 0) for aa in sequence[start : start + 7]))
        for start in range(len(sequence) - 6)
    )


def _finite(row: Mapping[str, object], field: str) -> float:
    value = float(row[field])
    if not isfinite(value):
        raise ExpressionDoubleMutantError(f"Non-finite {field}")
    return value


def _score_index(
    rows: Sequence[Mapping[str, object]], tool: str
) -> dict[str, Mapping[str, object]]:
    indexed = {str(row["sample_uid"]): row for row in rows}
    if len(indexed) != 163 or any(str(row.get("scoring_status")) != "pass" for row in rows):
        raise ExpressionDoubleMutantError(f"{tool} score coverage is incomplete")
    return indexed
