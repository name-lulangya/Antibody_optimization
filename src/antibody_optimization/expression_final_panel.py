"""Select the active 11 double mutants and assemble the 19+11 panel.

The selection consumes the frozen 162-row double-mutant property matrix.  It
uses only predeclared categorical magnitude bands: NetSolP U/S form one
predictor family, NanoMelt forms a second, and the worse of the two constituent
AntiFold bands forms a third.  Continuous values remain in the audit but never
break ties within a magnitude band.

Selection is a lexicographic binary optimization.  It first preserves the
strongest multi-family evidence and then applies risk, residual, stable-word,
and diversity criteria in the documented order.  This module does not run a
predictor and does not infer double-mutant AntiFold epistasis.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
from typing import Mapping, Sequence

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from .expression_panel_selection import GRADE, MAGNITUDE_THRESHOLDS, classify_change


EXPECTED_DOUBLE_COUNT = 162
EXPECTED_ELIGIBLE_COUNT = 84
EXPECTED_SELECTED_DOUBLE_COUNT = 11
EXPECTED_SINGLE_COUNT = 19
EXPECTED_FINAL_COUNT = 30


class ExpressionFinalPanelError(ValueError):
    """Raised when a frozen input or final-panel invariant is violated."""


def select_final_expression_panel(
    matrix_rows: Sequence[Mapping[str, object]],
    matrix_gate: Mapping[str, object],
    parent_rows: Sequence[Mapping[str, object]],
    parent_contract: Mapping[str, object],
    constraints: Mapping[str, object],
) -> dict[str, object]:
    """Return the complete audit, selected doubles, reserves, and final panel.

    The 11-double selection contract is:

    * reject hard sequence risk;
    * reject any moderate/strong adverse U, S, Tm, or AntiFold component band;
    * require at least one moderate/strong favorable predictor family;
    * select exactly 11 with at most one candidate per reported-position pair,
      at most two uses per exact component mutation, and at most three uses per
      reported position;
    * lexicographically maximize evidence class A, then B, strong-family count,
      absence of soft risk, absence of moderate/strong adverse interaction
      residuals, stable-word gain, exact-component coverage, and position
      coverage.  Candidate ID is the final deterministic tie-break only.
    """

    _validate_releases(matrix_rows, matrix_gate, parent_rows, parent_contract, constraints)
    parent_sequence = _parent_sequence(parent_rows)
    annotated = [_annotate(row) for row in matrix_rows]
    eligible = sorted(
        (row for row in annotated if row["double_selection_eligibility"] == "eligible"),
        key=lambda row: str(row["candidate_id"]),
    )
    if len(eligible) != EXPECTED_ELIGIBLE_COUNT:
        raise ExpressionFinalPanelError(
            f"Expected {EXPECTED_ELIGIBLE_COUNT} eligible doubles, observed {len(eligible)}"
        )

    selected_ids, optimizer = _lexicographic_select(eligible)
    selected = [row for row in eligible if str(row["candidate_id"]) in selected_ids]
    selected.sort(key=_display_order_key)
    selected = [
        {
            **row,
            "double_selection_status": "selected_final11",
            "double_selection_order": index,
            "double_selection_reason": _selection_reason(row),
        }
        for index, row in enumerate(selected, 1)
    ]
    _validate_diversity(selected)

    selected_by_id = {str(row["candidate_id"]): row for row in selected}
    selected_set = list(selected_by_id.values())
    reserves = []
    reserve_pool = []
    for row in eligible:
        if str(row["candidate_id"]) in selected_by_id:
            continue
        swaps = [
            str(chosen["candidate_id"])
            for chosen in selected_set
            if _valid_panel([item for item in selected_set if item is not chosen] + [row])
        ]
        reserve_pool.append((row, swaps))
    reserve_pool.sort(key=lambda item: (not bool(item[1]), _display_order_key(item[0])))
    for index, (row, swaps) in enumerate(reserve_pool, 1):
        reserves.append(
            {
                **row,
                "double_selection_status": (
                    "primary_reserve" if index <= 12 else "eligible_not_selected"
                ),
                "reserve_order": index,
                "direct_one_for_one_swap_feasible": bool(swaps),
                "feasible_replaced_selected_ids": "|".join(swaps),
                "double_selection_reason": (
                    "eligible_but_not_in_lexicographic_optimum; reserve requires the "
                    "listed one-for-one replacement to preserve all diversity caps"
                ),
            }
        )

    audit = []
    for row in annotated:
        candidate_id = str(row["candidate_id"])
        if candidate_id in selected_by_id:
            audit.append(dict(selected_by_id[candidate_id]))
        else:
            reserve = next(
                (item for item in reserves if str(item["candidate_id"]) == candidate_id),
                None,
            )
            if reserve is not None:
                audit.append(dict(reserve))
            else:
                audit.append(
                    {
                        **row,
                        "double_selection_status": "ineligible",
                        "double_selection_order": "",
                        "double_selection_reason": row["double_selection_exclusion_reasons"],
                    }
                )

    final_rows = _build_final_rows(parent_rows, selected)
    _validate_final_rows(final_rows, parent_sequence, constraints)
    facts = {
        "double_candidate_count": len(annotated),
        "eligible_double_count": len(eligible),
        "selected_double_count": len(selected),
        "reserve_count": len(reserves),
        "primary_reserve_count": min(12, len(reserves)),
        "selected_evidence_layer_counts": dict(
            sorted(Counter(str(row["evidence_layer"]) for row in selected).items())
        ),
        "selected_unique_component_mutation_count": len(
            {str(row[key]) for row in selected for key in ("mutation_a", "mutation_b")}
        ),
        "selected_unique_reported_position_count": len(
            {
                int(row[key])
                for row in selected
                for key in ("position_a_reported_1based", "position_b_reported_1based")
            }
        ),
        "selected_stable_word_gain_count": sum(
            int(int(row["net_stable_word_occurrence_delta"]) > 0) for row in selected
        ),
        "final_single_count": EXPECTED_SINGLE_COUNT,
        "final_double_count": EXPECTED_SELECTED_DOUBLE_COUNT,
        "final_candidate_count": len(final_rows),
    }
    return {
        "audit_rows": audit,
        "selected_double_rows": selected,
        "reserve_rows": reserves,
        "final_rows": final_rows,
        "optimizer": optimizer,
        "facts": facts,
        "parent_sequence": parent_sequence,
    }


def _annotate(source: Mapping[str, object]) -> dict[str, object]:
    row = dict(source)
    bands = {
        "netsolp_u": str(row["netsolp_u_magnitude_band"]),
        "netsolp_s": str(row["netsolp_s_magnitude_band"]),
        "nanomelt": str(row["nanomelt_tm_c_magnitude_band"]),
        "antifold": str(row["antifold_worst_component_band"]),
    }
    try:
        grades = {key: GRADE[value] for key, value in bands.items()}
    except KeyError as exc:
        raise ExpressionFinalPanelError(f"Unknown magnitude band: {exc}") from exc
    family_grades = {
        "netsolp": max(grades["netsolp_u"], grades["netsolp_s"]),
        "nanomelt": grades["nanomelt"],
        "antifold": grades["antifold"],
    }
    favorable = sum(value >= 1 for value in family_grades.values())
    strong = sum(value == 2 for value in family_grades.values())
    adverse_metrics = [key for key, value in grades.items() if value < 0]
    hard_risk = int(row["hard_sequence_risk_count"])
    reasons = []
    if hard_risk:
        reasons.append("hard_sequence_risk")
    if adverse_metrics:
        reasons.append("moderate_or_strong_adverse:" + "|".join(adverse_metrics))
    if favorable == 0:
        reasons.append("no_moderate_or_strong_favorable_family")
    residual_bands = {}
    for metric, field in (
        ("netsolp_u", "netsolp_u_interaction_residual"),
        ("netsolp_s", "netsolp_s_interaction_residual"),
        ("nanomelt_tm_c", "nanomelt_tm_c_interaction_residual"),
    ):
        residual_bands[metric] = classify_change(
            float(row[field]), MAGNITUDE_THRESHOLDS[metric]
        )
    residual_adverse_count = sum(GRADE[band] < 0 for band in residual_bands.values())
    layer = {3: "A_three_families", 2: "B_two_families", 1: "C_one_family"}.get(
        favorable, "ineligible"
    )
    return {
        **row,
        "netsolp_family_ordinal_grade": family_grades["netsolp"],
        "nanomelt_family_ordinal_grade": family_grades["nanomelt"],
        "antifold_family_ordinal_grade": family_grades["antifold"],
        "favorable_family_count": favorable,
        "strong_favorable_family_count": strong,
        "evidence_layer": layer,
        "netsolp_u_interaction_residual_band": residual_bands["netsolp_u"],
        "netsolp_s_interaction_residual_band": residual_bands["netsolp_s"],
        "nanomelt_tm_interaction_residual_band": residual_bands["nanomelt_tm_c"],
        "moderate_or_strong_adverse_interaction_residual_count": residual_adverse_count,
        "stable_word_gain_soft_tiebreak": int(row["net_stable_word_occurrence_delta"]) > 0,
        "double_selection_eligibility": "ineligible" if reasons else "eligible",
        "double_selection_exclusion_reasons": ";".join(reasons),
    }


def _lexicographic_select(
    eligible: Sequence[Mapping[str, object]],
) -> tuple[set[str], dict[str, object]]:
    candidates = list(eligible)
    components = sorted(
        {str(row[key]) for row in candidates for key in ("mutation_a", "mutation_b")}
    )
    positions = sorted(
        {
            int(row[key])
            for row in candidates
            for key in ("position_a_reported_1based", "position_b_reported_1based")
        }
    )
    n = len(candidates)
    component_index = {value: n + index for index, value in enumerate(components)}
    position_index = {
        value: n + len(components) + index for index, value in enumerate(positions)
    }
    variable_count = n + len(components) + len(positions)
    constraints: list[LinearConstraint] = []

    def add(coefficients: Mapping[int, float], lower: float, upper: float) -> None:
        vector = np.zeros(variable_count)
        for index, value in coefficients.items():
            vector[index] = value
        constraints.append(LinearConstraint(vector, lower, upper))

    add({index: 1 for index in range(n)}, EXPECTED_SELECTED_DOUBLE_COUNT, EXPECTED_SELECTED_DOUBLE_COUNT)
    by_pair: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, row in enumerate(candidates):
        by_pair[_position_pair(row)].append(index)
    for indices in by_pair.values():
        add({index: 1 for index in indices}, -np.inf, 1)
    for component in components:
        indices = [
            index
            for index, row in enumerate(candidates)
            if component in (str(row["mutation_a"]), str(row["mutation_b"]))
        ]
        add({**{index: 1 for index in indices}, component_index[component]: -2}, -np.inf, 0)
        add({component_index[component]: 1, **{index: -1 for index in indices}}, -np.inf, 0)
    for position in positions:
        indices = [
            index
            for index, row in enumerate(candidates)
            if position
            in (
                int(row["position_a_reported_1based"]),
                int(row["position_b_reported_1based"]),
            )
        ]
        add({**{index: 1 for index in indices}, position_index[position]: -3}, -np.inf, 0)
        add({position_index[position]: 1, **{index: -1 for index in indices}}, -np.inf, 0)

    zeros = [0] * (variable_count - n)
    objectives = [
        (
            "maximize_A_three_family_count",
            np.asarray([-int(row["favorable_family_count"] == 3) for row in candidates] + zeros),
            "maximize",
        ),
        (
            "maximize_B_two_family_count",
            np.asarray([-int(row["favorable_family_count"] == 2) for row in candidates] + zeros),
            "maximize",
        ),
        (
            "maximize_strong_favorable_family_count",
            np.asarray([-int(row["strong_favorable_family_count"]) for row in candidates] + zeros),
            "maximize",
        ),
        (
            "minimize_soft_sequence_risk_count",
            np.asarray([int(row["soft_sequence_risk_count"]) for row in candidates] + zeros),
            "minimize",
        ),
        (
            "minimize_adverse_interaction_residual_count",
            np.asarray(
                [
                    int(row["moderate_or_strong_adverse_interaction_residual_count"])
                    for row in candidates
                ]
                + zeros
            ),
            "minimize",
        ),
        (
            "maximize_stable_word_gain_count",
            np.asarray([-int(row["stable_word_gain_soft_tiebreak"]) for row in candidates] + zeros),
            "maximize",
        ),
        (
            "maximize_unique_component_mutations",
            np.asarray([0] * n + [-1] * len(components) + [0] * len(positions)),
            "maximize",
        ),
        (
            "maximize_unique_reported_positions",
            np.asarray([0] * (n + len(components)) + [-1] * len(positions)),
            "maximize",
        ),
    ]
    frozen: list[LinearConstraint] = []
    objective_results = []
    solution = None
    for name, vector, direction in objectives:
        solution = milp(
            vector,
            integrality=np.ones(variable_count),
            bounds=Bounds(np.zeros(variable_count), np.ones(variable_count)),
            constraints=[*constraints, *frozen],
        )
        if not solution.success:
            raise ExpressionFinalPanelError(f"MILP failed at {name}: {solution.message}")
        raw_value = int(round(float(vector @ np.rint(solution.x))))
        reported_value = -raw_value if direction == "maximize" else raw_value
        objective_results.append(
            {"objective": name, "direction": direction, "optimum": reported_value}
        )
        frozen.append(LinearConstraint(vector, raw_value, raw_value))

    deterministic = np.asarray(list(range(1, n + 1)) + zeros, dtype=float)
    solution = milp(
        deterministic,
        integrality=np.ones(variable_count),
        bounds=Bounds(np.zeros(variable_count), np.ones(variable_count)),
        constraints=[*constraints, *frozen],
    )
    if not solution.success:
        raise ExpressionFinalPanelError(f"MILP deterministic tie-break failed: {solution.message}")
    rounded = np.rint(solution.x)
    if not np.allclose(solution.x, rounded, atol=1e-7):
        raise ExpressionFinalPanelError("MILP did not return an integral solution")
    selected = {
        str(row["candidate_id"])
        for index, row in enumerate(candidates)
        if rounded[index] == 1
    }
    if len(selected) != EXPECTED_SELECTED_DOUBLE_COUNT:
        raise ExpressionFinalPanelError("MILP did not select exactly 11 doubles")
    return selected, {
        "solver": "scipy.optimize.milp",
        "selection_semantics": "lexicographic_categorical_objectives",
        "within_band_continuous_values_used_for_selection": False,
        "objectives_in_order": objective_results,
        "deterministic_final_tiebreak": "alphabetical_candidate_id_rank_sum",
        "deterministic_final_tiebreak_value": int(round(float(deterministic @ rounded))),
    }


def _build_final_rows(
    parent_rows: Sequence[Mapping[str, object]],
    selected: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    output = []
    for order, source in enumerate(
        sorted(parent_rows, key=lambda row: int(row["parent19_selection_order"])), 1
    ):
        output.append(
            {
                "final_panel_order": order,
                "candidate_id": source["candidate_id"],
                "candidate_kind": "single_mutant",
                "mutation_count": 1,
                "mutation_set": _single_mutation(source),
                "reported_positions_1based": source["reported_sequence_index_1based"],
                "regions": source["region"],
                "sequence": source["sequence"],
                "evidence_layer": "approved_parent19_single",
                "favorable_family_count": source["favorable_family_count"],
                "strong_favorable_family_count": source["strong_favorable_family_count"],
                "netsolp_family_ordinal_grade": max(
                    int(source["netsolp_u_ordinal_grade"]),
                    int(source["netsolp_s_ordinal_grade"]),
                ),
                "nanomelt_family_ordinal_grade": source["nanomelt_tm_ordinal_grade"],
                "antifold_family_ordinal_grade": source["antifold_ordinal_grade"],
                "hard_sequence_risk_flags": source["hard_sequence_risk_flags"],
                "soft_sequence_risk_flags": source["soft_sequence_risk_flags"],
                "stable_word_effect": source["stable_word_effect"],
                "selection_reason": source["parent19_selection_reason"],
                "source_artifact_role": "frozen_parent19_single_mutant",
            }
        )
    for offset, source in enumerate(selected, 1):
        output.append(
            {
                "final_panel_order": EXPECTED_SINGLE_COUNT + offset,
                "candidate_id": source["candidate_id"],
                "candidate_kind": "double_mutant",
                "mutation_count": 2,
                "mutation_set": source["mutation_set"],
                "reported_positions_1based": (
                    f"{source['position_a_reported_1based']};{source['position_b_reported_1based']}"
                ),
                "regions": f"{source['region_a']};{source['region_b']}",
                "sequence": source["sequence"],
                "evidence_layer": source["evidence_layer"],
                "favorable_family_count": source["favorable_family_count"],
                "strong_favorable_family_count": source["strong_favorable_family_count"],
                "netsolp_family_ordinal_grade": source["netsolp_family_ordinal_grade"],
                "nanomelt_family_ordinal_grade": source["nanomelt_family_ordinal_grade"],
                "antifold_family_ordinal_grade": source["antifold_family_ordinal_grade"],
                "hard_sequence_risk_flags": source["hard_sequence_risk_flags"],
                "soft_sequence_risk_flags": source["soft_sequence_risk_flags"],
                "stable_word_effect": source["stable_word_effect"],
                "selection_reason": source["double_selection_reason"],
                "source_artifact_role": "selected_from_complete_162_double_matrix",
            }
        )
    return output


def _validate_releases(
    matrix_rows: Sequence[Mapping[str, object]],
    matrix_gate: Mapping[str, object],
    parent_rows: Sequence[Mapping[str, object]],
    parent_contract: Mapping[str, object],
    constraints: Mapping[str, object],
) -> None:
    if len(matrix_rows) != EXPECTED_DOUBLE_COUNT:
        raise ExpressionFinalPanelError("Expected the complete 162-row double matrix")
    if matrix_gate.get("status") != "pass" or matrix_gate.get("release") != "ready_for_magnitude_aware_double_mutant_selection":
        raise ExpressionFinalPanelError("Double-mutant property-matrix gate is not open")
    if len(parent_rows) != EXPECTED_SINGLE_COUNT:
        raise ExpressionFinalPanelError("Expected exactly 19 approved parent singles")
    if list(parent_contract.get("selected_mutations_in_order", [])) != [
        _single_mutation(row)
        for row in sorted(parent_rows, key=lambda row: int(row["parent19_selection_order"]))
    ]:
        raise ExpressionFinalPanelError("Parent19 rows disagree with their frozen contract")
    if constraints.get("status") != "pass" or constraints.get("schema_version") != 2:
        raise ExpressionFinalPanelError("Expression design constraints are not released")
    if len({str(row["candidate_id"]) for row in matrix_rows}) != EXPECTED_DOUBLE_COUNT:
        raise ExpressionFinalPanelError("Double candidate IDs are not unique")
    if len({str(row["sequence"]) for row in matrix_rows}) != EXPECTED_DOUBLE_COUNT:
        raise ExpressionFinalPanelError("Double candidate sequences are not unique")


def _validate_final_rows(
    rows: Sequence[Mapping[str, object]],
    parent: str,
    constraints: Mapping[str, object],
) -> None:
    if len(rows) != EXPECTED_FINAL_COUNT:
        raise ExpressionFinalPanelError("Final panel must contain exactly 30 sequences")
    if Counter(str(row["candidate_kind"]) for row in rows) != {
        "single_mutant": EXPECTED_SINGLE_COUNT,
        "double_mutant": EXPECTED_SELECTED_DOUBLE_COUNT,
    }:
        raise ExpressionFinalPanelError("Final panel must contain 19 singles and 11 doubles")
    if len({str(row["candidate_id"]) for row in rows}) != EXPECTED_FINAL_COUNT:
        raise ExpressionFinalPanelError("Final candidate IDs are not unique")
    if len({str(row["sequence"]) for row in rows}) != EXPECTED_FINAL_COUNT:
        raise ExpressionFinalPanelError("Final sequences are not unique")
    frozen = {int(value) for value in constraints["hard_frozen_reported_indices_1based"]}
    consensus = {
        int(item["reported_sequence_index_1based"]): str(item["allowed_mutant_residue"])
        for item in constraints["consensus_reversion_only"]
    }
    expected_hash = str(constraints["authoritative_parent"]["sequence_sha256"])
    if sha256(parent.encode("ascii")).hexdigest() != expected_hash:
        raise ExpressionFinalPanelError("Reconstructed parent sequence hash disagrees with constraints")
    for row in rows:
        sequence = str(row["sequence"])
        if len(sequence) != len(parent):
            raise ExpressionFinalPanelError(f"Sequence length mismatch: {row['candidate_id']}")
        differences = [index for index, pair in enumerate(zip(parent, sequence, strict=True), 1) if pair[0] != pair[1]]
        if len(sequence) != 128 or len(differences) != int(row["mutation_count"]):
            raise ExpressionFinalPanelError(f"Mutation count mismatch: {row['candidate_id']}")
        if any(index in frozen for index in differences):
            raise ExpressionFinalPanelError(f"Frozen position changed: {row['candidate_id']}")
        if any(sequence[index - 1] == "C" for index in differences):
            raise ExpressionFinalPanelError(f"New cysteine introduced: {row['candidate_id']}")
        if sequence[21] != "C" or sequence[94] != "C" or not sequence.endswith("SSGS"):
            raise ExpressionFinalPanelError(f"Cys/terminal constraint failed: {row['candidate_id']}")
        for index in differences:
            if index in consensus and sequence[index - 1] != consensus[index]:
                raise ExpressionFinalPanelError(f"Consensus-only position violated: {row['candidate_id']}")


def _validate_diversity(rows: Sequence[Mapping[str, object]]) -> None:
    if not _valid_panel(rows):
        raise ExpressionFinalPanelError("Selected double panel violates a diversity cap")


def _valid_panel(rows: Sequence[Mapping[str, object]]) -> bool:
    if len({_position_pair(row) for row in rows}) != len(rows):
        return False
    component_counts = Counter(
        str(row[key]) for row in rows for key in ("mutation_a", "mutation_b")
    )
    position_counts = Counter(
        int(row[key])
        for row in rows
        for key in ("position_a_reported_1based", "position_b_reported_1based")
    )
    return max(component_counts.values(), default=0) <= 2 and max(position_counts.values(), default=0) <= 3


def _position_pair(row: Mapping[str, object]) -> tuple[int, int]:
    return tuple(
        sorted(
            (
                int(row["position_a_reported_1based"]),
                int(row["position_b_reported_1based"]),
            )
        )
    )


def _parent_sequence(parent_rows: Sequence[Mapping[str, object]]) -> str:
    reconstructed = set()
    for row in parent_rows:
        sequence = str(row["sequence"])
        position = int(row["reported_sequence_index_1based"])
        reconstructed.add(sequence[: position - 1] + str(row["wt_residue"]) + sequence[position:])
    if len(reconstructed) != 1:
        raise ExpressionFinalPanelError("Parent singles do not reconstruct one sequence")
    return reconstructed.pop()


def _single_mutation(row: Mapping[str, object]) -> str:
    return f"{row['wt_residue']}{int(row['reported_sequence_index_1based'])}{row['mutant_residue']}"


def _display_order_key(row: Mapping[str, object]) -> tuple[object, ...]:
    layer = {"A_three_families": 0, "B_two_families": 1, "C_one_family": 2}
    return (
        layer[str(row["evidence_layer"])],
        -int(row["strong_favorable_family_count"]),
        int(row["soft_sequence_risk_count"]),
        int(row["moderate_or_strong_adverse_interaction_residual_count"]),
        -int(bool(row["stable_word_gain_soft_tiebreak"])),
        str(row["candidate_id"]),
    )


def _selection_reason(row: Mapping[str, object]) -> str:
    return (
        f"{row['evidence_layer']}; {row['favorable_family_count']} predictor families "
        "have moderate/strong favorable bands; no moderate/strong adverse property "
        "band; retained by the frozen diversity-constrained lexicographic contract"
    )
