"""Select the V3 Nb252 expression-only 30-single-mutant panel.

NetSolP U, NetSolP S, and NanoMelt Tm are three separate positive evidence
metrics.  AntiFold never contributes positive selection credit: it only vetoes
a mutation when its log-probability change is at most -3 and the mutant amino
acid is among the four lowest-scoring of all 20 amino-acid states at that
position.  Raw values within a declared magnitude band never break ties.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Mapping, Sequence

from .expression_landscape_plot import build_expression_landscape_rows
from .expression_panel_selection import (
    EXPECTED_CANDIDATES,
    MAGNITUDE_THRESHOLDS,
    PANEL_SIZE,
    ExpressionPanelSelectionError,
    _finite,
    _sequence_risks,
    _validate_rows,
    _validate_upstream_gate,
    classify_change,
)


AA20 = frozenset("ACDEFGHIKLMNPQRSTVWY")
ANTIFOLD_DELTA_VETO = -3.0
ANTIFOLD_BOTTOM_COUNT = 4
TIER_ORDER = {
    "A_multi_metric": 0,
    "B_single_metric_strong": 1,
    "C_single_metric_moderate": 2,
    "D_controlled_tradeoff": 3,
}
GRADE = {
    "strong_adverse": -2,
    "moderate_adverse": -1,
    "weak_adverse": 0,
    "negligible": 0,
    "weak_favorable": 0,
    "moderate_favorable": 1,
    "strong_favorable": 2,
}


def build_expression_single_mutant_panel_v3(
    rows: Sequence[Mapping[str, object]],
    upstream_gate: Mapping[str, object],
    full_antifold_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Return the V3 audit, qualified pool, selected 30, reserves, and facts.

    ``rows`` is the released 847-row property/stable-word matrix.
    ``full_antifold_rows`` supplies every amino-acid substitution, including
    Cys substitutions that are forbidden as candidates but are required to
    define a complete 20-state within-position rank.  Selection is categorical
    and position-diverse; it does not compare raw values within a band.
    """

    _validate_upstream_gate(upstream_gate)
    parent = _validate_rows(rows)
    landscape_rows, _ = build_expression_landscape_rows(rows)
    anti_by_id = {str(row["candidate_id"]): row for row in landscape_rows}
    rank_evidence = _build_antifold_rank_evidence(rows, full_antifold_rows, anti_by_id)

    audit: list[dict[str, object]] = []
    for source in rows:
        identifier = str(source["candidate_id"])
        anti = anti_by_id[identifier]
        rank = rank_evidence[identifier]
        position = int(source["reported_sequence_index_1based"])
        wt = str(source["wt_residue"])
        mutant = str(source["mutant_residue"])
        sequence = str(source["sequence"])
        values = {
            "netsolp_u": _finite(source, "netsolp_delta_usability_vs_current_wt"),
            "netsolp_s": _finite(source, "netsolp_delta_solubility_vs_current_wt"),
            "nanomelt_tm_c": _finite(
                source, "nanomelt_delta_predicted_apparent_tm_c_vs_current_wt"
            ),
        }
        bands = {
            name: classify_change(value, MAGNITUDE_THRESHOLDS[name])
            for name, value in values.items()
        }
        grades = {name: GRADE[band] for name, band in bands.items()}
        favorable_count = sum(value >= 1 for value in grades.values())
        strong_favorable_count = sum(value == 2 for value in grades.values())
        moderate_adverse_count = sum(value == -1 for value in grades.values())
        strong_adverse_count = sum(value == -2 for value in grades.values())
        anti_delta = float(anti["antifold_landscape_delta_log_probability"])
        anti_veto = (
            anti_delta <= ANTIFOLD_DELTA_VETO
            and int(rank["antifold_mutant_rank_worst_first"]) <= ANTIFOLD_BOTTOM_COUNT
        )
        hard_flags, soft_flags, risk_details = _sequence_risks(
            parent, sequence, position, wt, mutant
        )

        if hard_flags:
            eligibility = "blocked_sequence_risk"
            tier = "not_eligible"
        elif anti_veto:
            eligibility = "blocked_antifold_lowest_20_percent"
            tier = "not_eligible"
        elif strong_adverse_count:
            eligibility = "blocked_property_strong_adverse"
            tier = "not_eligible"
        elif moderate_adverse_count == 0 and favorable_count >= 2:
            eligibility = "qualified"
            tier = "A_multi_metric"
        elif (
            moderate_adverse_count == 0
            and favorable_count == 1
            and strong_favorable_count == 1
        ):
            eligibility = "qualified"
            tier = "B_single_metric_strong"
        elif moderate_adverse_count == 0 and favorable_count == 1:
            eligibility = "qualified"
            tier = "C_single_metric_moderate"
        elif moderate_adverse_count == 1 and favorable_count >= 2:
            eligibility = "qualified"
            tier = "D_controlled_tradeoff"
        elif favorable_count == 0:
            eligibility = "no_moderate_or_strong_positive_metric"
            tier = "not_eligible"
        else:
            eligibility = "blocked_property_tradeoff"
            tier = "not_eligible"

        stable_word_gain = str(source["stable_word_effect"]) in {"gain_only", "net_gain"}
        audit.append(
            {
                **dict(source),
                "antifold_selection_source": anti["antifold_landscape_source"],
                "antifold_selection_delta_log_probability": anti_delta,
                **rank,
                "antifold_veto_delta_threshold": ANTIFOLD_DELTA_VETO,
                "antifold_veto_bottom_rank_maximum": ANTIFOLD_BOTTOM_COUNT,
                "antifold_veto_status": "veto" if anti_veto else "pass",
                "netsolp_u_magnitude_band_v3": bands["netsolp_u"],
                "netsolp_s_magnitude_band_v3": bands["netsolp_s"],
                "nanomelt_tm_magnitude_band_v3": bands["nanomelt_tm_c"],
                "netsolp_u_ordinal_grade_v3": grades["netsolp_u"],
                "netsolp_s_ordinal_grade_v3": grades["netsolp_s"],
                "nanomelt_tm_ordinal_grade_v3": grades["nanomelt_tm_c"],
                "positive_metric_count_v3": favorable_count,
                "strong_positive_metric_count_v3": strong_favorable_count,
                "moderate_adverse_property_count_v3": moderate_adverse_count,
                "strong_adverse_property_count_v3": strong_adverse_count,
                "hard_sequence_risk_flags_v3": ";".join(hard_flags),
                "soft_sequence_risk_flags_v3": ";".join(soft_flags),
                "hard_sequence_risk_count_v3": len(hard_flags),
                "soft_sequence_risk_count_v3": len(soft_flags),
                **{f"v3_{key}": value for key, value in risk_details.items()},
                "stable_word_gain_tiebreak_v3": stable_word_gain,
                "selection_eligibility_v3": eligibility,
                "selection_tier_v3": tier,
                "position_round_v3": "",
                "selection_status_v3": "not_selected",
                "selection_order_v3": "",
                "selection_reason_v3": "",
                "final_experimental_panel_released": False,
            }
        )

    qualified = [row for row in audit if row["selection_eligibility_v3"] == "qualified"]
    if len(qualified) < PANEL_SIZE:
        raise ExpressionPanelSelectionError(
            f"V3 qualified pool has only {len(qualified)} candidates; cannot select {PANEL_SIZE}"
        )
    selected, reserves = _select_position_diverse_panel(qualified)
    selected_ids = {str(row["candidate_id"]) for row in selected}
    reserve_ids = {str(row["candidate_id"]) for row in reserves}
    for order, row in enumerate(selected, 1):
        row["selection_status_v3"] = "selected_final30_single"
        row["selection_order_v3"] = order
        row["selection_reason_v3"] = (
            "selected_by_position_round_then_tier_and_categorical_evidence"
        )
    for order, row in enumerate(reserves, 1):
        row["selection_status_v3"] = "qualified_reserve"
        row["selection_order_v3"] = order
        row["selection_reason_v3"] = "qualified_but_outside_30_member_diversity_limit"
    for row in audit:
        if str(row["candidate_id"]) in selected_ids | reserve_ids:
            continue
        row["selection_reason_v3"] = str(row["selection_eligibility_v3"])

    facts = {
        "candidate_count": len(audit),
        "antifold_veto_count": sum(row["antifold_veto_status"] == "veto" for row in audit),
        "antifold_veto_experimental_complex_count": sum(
            row["antifold_veto_status"] == "veto"
            and row["antifold_selection_source"] == "experimental_complex_context"
            for row in audit
        ),
        "antifold_veto_af3_fallback_count": sum(
            row["antifold_veto_status"] == "veto"
            and str(row["antifold_selection_source"]).startswith("af3_")
            for row in audit
        ),
        "qualified_count": len(qualified),
        "qualified_tier_counts": dict(Counter(str(row["selection_tier_v3"]) for row in qualified)),
        "selected_count": len(selected),
        "reserve_count": len(reserves),
        "selected_tier_counts": dict(Counter(str(row["selection_tier_v3"]) for row in selected)),
        "selected_position_counts": dict(
            sorted(Counter(int(row["reported_sequence_index_1based"]) for row in selected).items())
        ),
        "selected_unique_position_count": len(
            {int(row["reported_sequence_index_1based"]) for row in selected}
        ),
        "selected_maximum_per_position": max(
            Counter(int(row["reported_sequence_index_1based"]) for row in selected).values()
        ),
        "selected_stable_word_gain_count": sum(
            bool(row["stable_word_gain_tiebreak_v3"]) for row in selected
        ),
        "raw_within_band_values_used_for_ranking": False,
        "antifold_positive_credit_used": False,
        "final_experimental_panel_released": True,
    }
    _validate_v3_result(audit, qualified, selected, reserves, facts)
    return {
        "audit_rows": audit,
        "qualified_rows": sorted(qualified, key=_quality_key),
        "panel_rows": selected,
        "reserve_rows": reserves,
        "facts": facts,
        "parent_sequence": parent,
    }


def _build_antifold_rank_evidence(
    rows: Sequence[Mapping[str, object]],
    full_antifold_rows: Sequence[Mapping[str, object]],
    anti_by_id: Mapping[str, Mapping[str, object]],
) -> dict[str, dict[str, object]]:
    full_by_position: dict[int, list[Mapping[str, object]]] = defaultdict(list)
    for row in full_antifold_rows:
        full_by_position[int(row["sequence_index_1based"])].append(row)
    current_by_position: dict[int, list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        current_by_position[int(row["reported_sequence_index_1based"])].append(row)

    evidence: dict[str, dict[str, object]] = {}
    for position, candidates in current_by_position.items():
        sources = {str(anti_by_id[str(row["candidate_id"])]["antifold_landscape_source"]) for row in candidates}
        if len(sources) != 1:
            raise ExpressionPanelSelectionError(f"Mixed AntiFold sources at position {position}")
        source = sources.pop()
        if source == "experimental_complex_context":
            mutant_key = "experimental_complex_context_mutant_log_probability"
            wt_key = "experimental_complex_context_wt_log_probability"
        elif source == "af3_vhh_only_fallback_for_missing_experimental_coordinates":
            mutant_key = "af3_vhh_only_mutant_log_probability"
            wt_key = "af3_vhh_only_wt_log_probability"
        else:
            raise ExpressionPanelSelectionError(f"Unexpected AntiFold source: {source}")
        wt_residues = {str(row["wt_residue"]) for row in candidates}
        if len(wt_residues) != 1:
            raise ExpressionPanelSelectionError(f"Mixed WT residues at position {position}")
        wt = wt_residues.pop()
        wt_values = {_finite(row, wt_key) for row in candidates}
        if len(wt_values) != 1:
            raise ExpressionPanelSelectionError(f"Mixed WT AntiFold scores at position {position}")
        scores = {wt: wt_values.pop()}
        full = full_by_position.get(position, [])
        if len(full) != 19:
            raise ExpressionPanelSelectionError(
                f"Expected 19 AntiFold substitutions at position {position}; observed {len(full)}"
            )
        for row in full:
            mutant = str(row["mutant_residue"])
            if mutant in scores:
                raise ExpressionPanelSelectionError(f"Duplicate AntiFold state at position {position}")
            scores[mutant] = _finite(row, mutant_key)
        if set(scores) != AA20:
            raise ExpressionPanelSelectionError(f"Incomplete 20-state AntiFold scores at position {position}")
        ordered = sorted(scores.items(), key=lambda item: (item[1], item[0]))
        cutoff = ordered[ANTIFOLD_BOTTOM_COUNT - 1][1]
        if ordered[ANTIFOLD_BOTTOM_COUNT][1] == cutoff:
            raise ExpressionPanelSelectionError(
                f"AntiFold bottom-20% boundary is tied at position {position}"
            )
        rank_by_aa = {amino_acid: rank for rank, (amino_acid, _) in enumerate(ordered, 1)}
        for row in candidates:
            identifier = str(row["candidate_id"])
            mutant = str(row["mutant_residue"])
            observed = float(anti_by_id[identifier]["antifold_landscape_delta_log_probability"])
            expected = scores[mutant] - scores[wt]
            if abs(observed - expected) > 1e-8:
                raise ExpressionPanelSelectionError(
                    f"AntiFold score mismatch for {identifier}: {observed} != {expected}"
                )
            evidence[identifier] = {
                "antifold_mutant_rank_worst_first": rank_by_aa[mutant],
                "antifold_position_state_count": len(scores),
                "antifold_bottom_20_percent_log_probability_cutoff": cutoff,
                "antifold_mutant_log_probability_for_rank": scores[mutant],
            }
    if len(evidence) != EXPECTED_CANDIDATES:
        raise ExpressionPanelSelectionError("AntiFold rank evidence does not cover 847 candidates")
    return evidence


def _select_position_diverse_panel(
    qualified: Sequence[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grouped: dict[int, list[dict[str, object]]] = defaultdict(list)
    for row in qualified:
        grouped[int(row["reported_sequence_index_1based"])].append(row)
    for group in grouped.values():
        group.sort(key=_quality_key)
        for rank, row in enumerate(group, 1):
            row["position_round_v3"] = rank

    selected: list[dict[str, object]] = []
    for position_round in range(1, 4):
        round_rows = [
            group[position_round - 1]
            for group in grouped.values()
            if len(group) >= position_round
        ]
        round_rows.sort(key=_quality_key)
        for row in round_rows:
            if len(selected) == PANEL_SIZE:
                break
            selected.append(row)
        if len(selected) == PANEL_SIZE:
            break
    if len(selected) != PANEL_SIZE:
        raise ExpressionPanelSelectionError("Three-per-position cap cannot fill V3 panel")
    selected_ids = {str(row["candidate_id"]) for row in selected}
    selected.sort(key=lambda row: (int(row["position_round_v3"]), *_quality_key(row)))
    reserves = sorted(
        [row for row in qualified if str(row["candidate_id"]) not in selected_ids],
        key=lambda row: (int(row["position_round_v3"]), *_quality_key(row)),
    )
    return selected, reserves


def _quality_key(row: Mapping[str, object]) -> tuple[object, ...]:
    return (
        TIER_ORDER[str(row["selection_tier_v3"])],
        -int(row["strong_positive_metric_count_v3"]),
        -int(row["positive_metric_count_v3"]),
        int(row["moderate_adverse_property_count_v3"]),
        int(row["soft_sequence_risk_count_v3"]),
        0 if bool(row["stable_word_gain_tiebreak_v3"]) else 1,
        str(row["candidate_id"]),
    )


def _validate_v3_result(
    audit: Sequence[Mapping[str, object]],
    qualified: Sequence[Mapping[str, object]],
    panel: Sequence[Mapping[str, object]],
    reserves: Sequence[Mapping[str, object]],
    facts: Mapping[str, object],
) -> None:
    if len(audit) != EXPECTED_CANDIDATES or len(panel) != PANEL_SIZE:
        raise ExpressionPanelSelectionError("Unexpected V3 audit or panel count")
    if len(qualified) != len(panel) + len(reserves):
        raise ExpressionPanelSelectionError("V3 qualified accounting mismatch")
    if len({str(row["sequence"]) for row in panel}) != PANEL_SIZE:
        raise ExpressionPanelSelectionError("V3 panel sequences are not unique")
    if any(row["antifold_veto_status"] != "pass" for row in panel):
        raise ExpressionPanelSelectionError("V3 panel contains an AntiFold-vetoed mutation")
    if any(int(row["hard_sequence_risk_count_v3"]) for row in panel):
        raise ExpressionPanelSelectionError("V3 panel contains a hard sequence risk")
    if any(int(row["strong_adverse_property_count_v3"]) for row in panel):
        raise ExpressionPanelSelectionError("V3 panel contains a strongly adverse property")
    if any(int(row["positive_metric_count_v3"]) < 1 for row in panel):
        raise ExpressionPanelSelectionError("V3 panel contains a candidate without positive evidence")
    position_counts = Counter(int(row["reported_sequence_index_1based"]) for row in panel)
    if max(position_counts.values()) > 3:
        raise ExpressionPanelSelectionError("V3 panel exceeds three candidates per position")
    if int(facts["antifold_veto_count"]) != 151:
        raise ExpressionPanelSelectionError("Unexpected AntiFold combined-veto count")
