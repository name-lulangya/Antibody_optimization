"""Magnitude-aware trial selection of 30 Nb252 expression single mutants.

The module consumes the released 847-row property and stable-word matrix.  It
converts continuous predictor changes into predeclared ordinal bands so that
small decimal differences cannot drive selection.  Hard sequence constraints,
strong adverse predictor bands, provenance, and position diversity remain
explicit.  The result is a trial panel for review, not experimental validation.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Mapping, Sequence

from .expression_landscape_plot import build_expression_landscape_rows
from .unified_single_mutants import sequence_liability_deltas


class ExpressionPanelSelectionError(ValueError):
    """Raised when the released matrix violates the trial-selection contract."""


EXPECTED_CANDIDATES = 847
PANEL_SIZE = 30
STABLE_WORD_EXPLORATORY_CANDIDATE_ID = "Nb252_expr_seq099_T99F"
STABLE_WORD_REPLACED_CANDIDATE_ID = "Nb252_expr_seq099_T99N"
HYDROPHOBIC = set("AVILMFWY")
SIDECHAIN_CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}

# (negligible/weak boundary, weak/moderate boundary, moderate/strong boundary)
MAGNITUDE_THRESHOLDS = {
    "netsolp_u": (0.005, 0.010, 0.015),
    "netsolp_s": (0.010, 0.020, 0.030),
    "nanomelt_tm_c": (0.5, 1.0, 1.5),
    "antifold_log_probability": (0.5, 1.5, 3.0),
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


def classify_change(value: float, thresholds: tuple[float, float, float]) -> str:
    """Classify a higher-is-better change without ranking within a band."""

    negligible, moderate, strong = thresholds
    magnitude = abs(float(value))
    if magnitude < negligible:
        return "negligible"
    direction = "favorable" if value > 0 else "adverse"
    if magnitude < moderate:
        return f"weak_{direction}"
    if magnitude < strong:
        return f"moderate_{direction}"
    return f"strong_{direction}"


def build_expression_trial_panel(
    rows: Sequence[Mapping[str, object]],
    upstream_gate: Mapping[str, object],
) -> dict[str, object]:
    """Return an 847-row audit, magnitude shortlist, trial 30, and reserves.

    NetSolP U and S form one predictor family.  Only moderate or strong bands
    count as favorable evidence.  Strong adverse bands block advancement.
    Candidates with no moderate adverse bands form the strict core; candidates
    with exactly one moderate adverse band form a separately labelled
    controlled-tradeoff layer.  New proline, a newly dense seven-residue
    hydrophobic window, or a newly extreme local charge cluster blocks the
    sequence before panel selection.

    Selection is deterministic and categorical.  Every strict-core candidate
    is retained first.  The remaining panel slots are filled from the
    controlled-tradeoff layer, prioritizing positions absent from the strict
    core and then categorical evidence.  A user-reviewed T99F stable-word
    hypothesis candidate then replaces T99N under an explicit weak-only/no-
    moderate-adverse exception.  Raw within-band decimals are unused.
    """

    _validate_upstream_gate(upstream_gate)
    parent = _validate_rows(rows)
    landscape_rows, _ = build_expression_landscape_rows(rows)
    anti_by_id = {str(row["candidate_id"]): row for row in landscape_rows}

    audit: list[dict[str, object]] = []
    for source in rows:
        identifier = str(source["candidate_id"])
        anti = anti_by_id[identifier]
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
            "antifold_log_probability": float(
                anti["antifold_landscape_delta_log_probability"]
            ),
        }
        bands = {
            name: classify_change(value, MAGNITUDE_THRESHOLDS[name])
            for name, value in values.items()
        }
        grades = {name: GRADE[label] for name, label in bands.items()}
        netsolp_family_grade = max(grades["netsolp_u"], grades["netsolp_s"])
        family_grades = (
            netsolp_family_grade,
            grades["nanomelt_tm_c"],
            grades["antifold_log_probability"],
        )
        favorable_family_count = sum(value >= 1 for value in family_grades)
        strong_favorable_family_count = sum(value == 2 for value in family_grades)
        strong_adverse_count = sum(value == -2 for value in grades.values())
        moderate_adverse_count = sum(value == -1 for value in grades.values())

        hard_flags, soft_flags, risk_details = _sequence_risks(
            parent, sequence, position, wt, mutant
        )
        magnitude_shortlist = strong_adverse_count == 0 and favorable_family_count >= 1
        if not magnitude_shortlist:
            eligibility = "not_in_magnitude_shortlist"
        elif hard_flags:
            eligibility = "blocked_sequence_risk"
        elif moderate_adverse_count == 0:
            eligibility = "strict_core"
        elif moderate_adverse_count == 1:
            eligibility = "controlled_tradeoff"
        else:
            eligibility = "blocked_multiple_moderate_adverse"

        if eligibility == "strict_core":
            if favorable_family_count >= 2:
                tier = "A_multi_family"
            elif strong_favorable_family_count >= 1:
                tier = "B_single_family_strong"
            else:
                tier = "C_single_family_moderate"
        elif eligibility == "controlled_tradeoff":
            tier = "D_controlled_tradeoff"
        else:
            tier = "not_eligible"

        stable_word_gain = str(source["stable_word_effect"]) in {"gain_only", "net_gain"}
        audit.append(
            {
                **dict(source),
                "antifold_selection_source": anti["antifold_landscape_source"],
                "antifold_selection_delta_log_probability": values[
                    "antifold_log_probability"
                ],
                "netsolp_u_magnitude_band": bands["netsolp_u"],
                "netsolp_s_magnitude_band": bands["netsolp_s"],
                "nanomelt_tm_magnitude_band": bands["nanomelt_tm_c"],
                "antifold_magnitude_band": bands["antifold_log_probability"],
                "netsolp_u_ordinal_grade": grades["netsolp_u"],
                "netsolp_s_ordinal_grade": grades["netsolp_s"],
                "nanomelt_tm_ordinal_grade": grades["nanomelt_tm_c"],
                "antifold_ordinal_grade": grades["antifold_log_probability"],
                "favorable_family_count": favorable_family_count,
                "strong_favorable_family_count": strong_favorable_family_count,
                "moderate_adverse_metric_count": moderate_adverse_count,
                "strong_adverse_metric_count": strong_adverse_count,
                "hard_sequence_risk_flags": ";".join(hard_flags),
                "soft_sequence_risk_flags": ";".join(soft_flags),
                "hard_sequence_risk_count": len(hard_flags),
                "soft_sequence_risk_count": len(soft_flags),
                **risk_details,
                "stable_word_gain_tiebreak": stable_word_gain,
                "magnitude_shortlist_status": "pass" if magnitude_shortlist else "fail",
                "selection_eligibility_class": eligibility,
                "selection_tier": tier,
                "position_round": "",
                "trial_selection_status": "not_selected",
                "trial_selection_reason": "",
                "trial_selection_order": "",
                "trial_panel_selection_performed": True,
                "final_experimental_panel_released": False,
            }
        )

    strict_core = [
        row
        for row in audit
        if row["selection_eligibility_class"] == "strict_core"
    ]
    controlled_tradeoffs = [
        row
        for row in audit
        if row["selection_eligibility_class"] == "controlled_tradeoff"
    ]
    eligible = [*strict_core, *controlled_tradeoffs]
    if len(strict_core) > PANEL_SIZE or len(eligible) < PANEL_SIZE:
        raise ExpressionPanelSelectionError(
            "Declared strict-core/control-tradeoff layers cannot form a 30-member panel"
        )
    _assign_position_rounds(eligible)
    strict_core.sort(key=_quality_key)
    strict_position_counts = Counter(
        int(row["reported_sequence_index_1based"]) for row in strict_core
    )
    controlled_tradeoffs.sort(
        key=lambda row: (
            0
            if int(row["reported_sequence_index_1based"]) not in strict_position_counts
            else 1,
            strict_position_counts[int(row["reported_sequence_index_1based"])],
            *_quality_key(row),
        )
    )
    remaining_slots = PANEL_SIZE - len(strict_core)
    selected = [*strict_core, *controlled_tradeoffs[:remaining_slots]]
    reserves = controlled_tradeoffs[remaining_slots:]
    selected_by_id = {str(row["candidate_id"]): row for row in selected}
    audit_by_id = {str(row["candidate_id"]): row for row in audit}
    if STABLE_WORD_REPLACED_CANDIDATE_ID not in selected_by_id:
        raise ExpressionPanelSelectionError("Expected T99N diversity candidate was not selected")
    exploratory = audit_by_id[STABLE_WORD_EXPLORATORY_CANDIDATE_ID]
    replaced = selected_by_id[STABLE_WORD_REPLACED_CANDIDATE_ID]
    if (
        exploratory["selection_eligibility_class"] != "not_in_magnitude_shortlist"
        or int(exploratory["strong_adverse_metric_count"]) != 0
        or int(exploratory["moderate_adverse_metric_count"]) != 0
        or not bool(exploratory["stable_word_gain_tiebreak"])
    ):
        raise ExpressionPanelSelectionError("T99F no longer satisfies its exploratory contract")
    selected = [row for row in selected if row is not replaced]
    selected.append(exploratory)
    reserves.append(replaced)
    exploratory["selection_tier"] = "E_stable_word_exploratory"
    selected_ids = {str(row["candidate_id"]) for row in selected}
    reserve_ids = {str(row["candidate_id"]) for row in reserves}

    for order, row in enumerate(selected, start=1):
        row["trial_selection_status"] = "trial_final30"
        row["trial_selection_order"] = order
        row["trial_selection_reason"] = (
            "strict_core_retained_before_diversity_layer"
            if row["selection_eligibility_class"] == "strict_core"
            else (
                "stable_word_exploratory_user_selected_replacing_T99N"
                if row["selection_tier"] == "E_stable_word_exploratory"
                else "controlled_tradeoff_selected_for_position_diversity"
            )
        )
    for order, row in enumerate(reserves, start=1):
        row["trial_selection_status"] = "reserve"
        row["trial_selection_order"] = order
        row["trial_selection_reason"] = "controlled_tradeoff_reserve_after_diversity_fill"
    for row in audit:
        if str(row["candidate_id"]) in selected_ids | reserve_ids:
            continue
        row["trial_selection_reason"] = _nonselection_reason(row)

    selected = sorted(selected, key=lambda row: int(row["trial_selection_order"]))
    reserves = sorted(reserves, key=lambda row: int(row["trial_selection_order"]))
    shortlist = [row for row in audit if row["magnitude_shortlist_status"] == "pass"]
    facts = {
        "candidate_count": len(audit),
        "magnitude_shortlist_count": len(shortlist),
        "magnitude_shortlist_stable_word_gain_count": sum(
            bool(row["stable_word_gain_tiebreak"]) for row in shortlist
        ),
        "strict_core_count": sum(
            row["selection_eligibility_class"] == "strict_core" for row in audit
        ),
        "controlled_tradeoff_count": sum(
            row["selection_eligibility_class"] == "controlled_tradeoff" for row in audit
        ),
        "blocked_sequence_risk_count": sum(
            row["selection_eligibility_class"] == "blocked_sequence_risk" for row in audit
        ),
        "blocked_multiple_moderate_adverse_count": sum(
            row["selection_eligibility_class"] == "blocked_multiple_moderate_adverse"
            for row in audit
        ),
        "trial_panel_count": len(selected),
        "reserve_count": len(reserves),
        "trial_panel_tier_counts": dict(
            Counter(str(row["selection_tier"]) for row in selected)
        ),
        "trial_panel_region_counts": dict(
            Counter(str(row["region"]) for row in selected)
        ),
        "trial_panel_antifold_source_counts": dict(
            Counter(str(row["antifold_selection_source"]) for row in selected)
        ),
        "trial_panel_position_counts": dict(
            sorted(
                Counter(
                    int(row["reported_sequence_index_1based"]) for row in selected
                ).items()
            )
        ),
        "trial_panel_unique_position_count": len(
            {int(row["reported_sequence_index_1based"]) for row in selected}
        ),
        "trial_panel_stable_word_gain_count": sum(
            bool(row["stable_word_gain_tiebreak"]) for row in selected
        ),
        "trial_panel_stable_word_exploratory_count": sum(
            row["selection_tier"] == "E_stable_word_exploratory" for row in selected
        ),
        "trial_panel_selection_performed": True,
        "final_experimental_panel_released": False,
    }
    _validate_result(audit, shortlist, selected, reserves, facts)
    return {
        "audit_rows": audit,
        "shortlist_rows": sorted(
            shortlist,
            key=lambda row: (
                int(row["reported_sequence_index_1based"]),
                str(row["mutant_residue"]),
            ),
        ),
        "panel_rows": selected,
        "reserve_rows": reserves,
        "facts": facts,
        "parent_sequence": parent,
    }


def _validate_upstream_gate(gate: Mapping[str, object]) -> None:
    if gate.get("status") != "pass" or int(gate.get("candidate_count", -1)) != EXPECTED_CANDIDATES:
        raise ExpressionPanelSelectionError("Stable-word upstream gate is not released")
    effect_counts = gate.get("stable_word_effect_counts")
    if not isinstance(effect_counts, Mapping) or int(effect_counts.get("gain_only", -1)) != 22:
        raise ExpressionPanelSelectionError("Unexpected stable-word gain count")


def _validate_rows(rows: Sequence[Mapping[str, object]]) -> str:
    if len(rows) != EXPECTED_CANDIDATES:
        raise ExpressionPanelSelectionError(f"Expected {EXPECTED_CANDIDATES} candidates")
    if len({str(row["candidate_id"]) for row in rows}) != EXPECTED_CANDIDATES:
        raise ExpressionPanelSelectionError("Candidate identifiers are not unique")
    if len({str(row["sequence"]) for row in rows}) != EXPECTED_CANDIDATES:
        raise ExpressionPanelSelectionError("Candidate sequences are not unique")
    parents: set[str] = set()
    for row in rows:
        sequence = str(row["sequence"])
        position = int(row["reported_sequence_index_1based"])
        wt = str(row["wt_residue"])
        mutant = str(row["mutant_residue"])
        if len(sequence) != 128 or not 1 <= position <= 128:
            raise ExpressionPanelSelectionError("Invalid candidate sequence or position")
        if sequence[position - 1] != mutant or wt == mutant:
            raise ExpressionPanelSelectionError("Candidate mutation identity mismatch")
        parent = sequence[: position - 1] + wt + sequence[position:]
        differences = [
            index
            for index, (left, right) in enumerate(zip(parent, sequence, strict=True), 1)
            if left != right
        ]
        if differences != [position] or sequence.count("C") != parent.count("C"):
            raise ExpressionPanelSelectionError("Candidate is not a valid non-Cys single mutant")
        parents.add(parent)
    if len(parents) != 1:
        raise ExpressionPanelSelectionError("Candidates do not share one parent")
    parent = parents.pop()
    if not parent.endswith("SSGS"):
        raise ExpressionPanelSelectionError("Authoritative terminal SSGS was not retained")
    return parent


def _sequence_risks(
    parent: str, sequence: str, position: int, wt: str, mutant: str
) -> tuple[list[str], list[str], dict[str, object]]:
    start = max(0, position - 4)
    stop = min(len(parent), position + 3)
    parent_window = parent[start:stop]
    mutant_window = sequence[start:stop]
    parent_hydrophobic = sum(residue in HYDROPHOBIC for residue in parent_window)
    mutant_hydrophobic = sum(residue in HYDROPHOBIC for residue in mutant_window)
    parent_charge = _max_local_charge(parent)
    mutant_charge = _max_local_charge(sequence)

    hard: list[str] = []
    if mutant == "P" and wt != "P":
        hard.append("new_proline_backbone_constraint")
    if mutant_hydrophobic >= 6 and mutant_hydrophobic > parent_hydrophobic:
        hard.append("new_dense_local_hydrophobic_window")
    if mutant_charge >= 5 and mutant_charge > parent_charge:
        hard.append("new_extreme_local_charge_cluster")

    liability = sequence_liability_deltas(parent, sequence)
    soft = [token for token in str(liability["new_liability_flags"]).split("|") if token]
    return hard, soft, {
        "local_sequence_parent": parent_window,
        "local_sequence_mutant": mutant_window,
        "local_hydrophobic_count_parent": parent_hydrophobic,
        "local_hydrophobic_count_mutant": mutant_hydrophobic,
        "max_abs_local_charge_parent": parent_charge,
        "max_abs_local_charge_mutant": mutant_charge,
        "new_n_linked_glycosylation_motif_delta": liability[
            "n_linked_glycosylation_motif_delta"
        ],
        "new_deamidation_motif_delta": liability["deamidation_motif_delta"],
        "new_isomerization_motif_delta": liability["isomerization_motif_delta"],
        "oxidation_susceptible_residue_delta": liability[
            "oxidation_susceptible_residue_delta"
        ],
    }


def _max_local_charge(sequence: str) -> int:
    return max(
        abs(sum(SIDECHAIN_CHARGE.get(residue, 0) for residue in sequence[start : start + 7]))
        for start in range(len(sequence) - 6)
    )


def _assign_position_rounds(rows: Sequence[dict[str, object]]) -> None:
    grouped: dict[int, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[int(row["reported_sequence_index_1based"])].append(row)
    for candidates in grouped.values():
        candidates.sort(key=_quality_key)
        for position_round, row in enumerate(candidates, start=1):
            row["position_round"] = position_round


def _quality_key(row: Mapping[str, object]) -> tuple[object, ...]:
    tier_order = {
        "A_multi_family": 0,
        "B_single_family_strong": 1,
        "C_single_family_moderate": 2,
        "D_controlled_tradeoff": 3,
        "E_stable_word_exploratory": 4,
    }
    source_order = {
        "experimental_complex_context": 0,
        "af3_vhh_only_fallback_for_missing_experimental_coordinates": 1,
    }
    return (
        tier_order[str(row["selection_tier"])],
        -int(row["strong_favorable_family_count"]),
        -int(row["favorable_family_count"]),
        int(row["moderate_adverse_metric_count"]),
        source_order[str(row["antifold_selection_source"])],
        int(row["soft_sequence_risk_count"]),
        0 if bool(row["stable_word_gain_tiebreak"]) else 1,
        str(row["candidate_id"]),
    )


def _nonselection_reason(row: Mapping[str, object]) -> str:
    eligibility = str(row["selection_eligibility_class"])
    if eligibility == "not_in_magnitude_shortlist":
        if int(row["strong_adverse_metric_count"]) > 0:
            return "blocked_by_at_least_one_strong_adverse_metric_band"
        return "no_moderate_or_strong_favorable_predictor_family"
    if eligibility == "blocked_sequence_risk":
        return "blocked_by_declared_sequence_risk"
    if eligibility == "blocked_multiple_moderate_adverse":
        return "blocked_by_multiple_moderate_adverse_metric_bands"
    return "not_selected"


def _validate_result(
    audit: Sequence[Mapping[str, object]],
    shortlist: Sequence[Mapping[str, object]],
    panel: Sequence[Mapping[str, object]],
    reserves: Sequence[Mapping[str, object]],
    facts: Mapping[str, object],
) -> None:
    if len(audit) != EXPECTED_CANDIDATES or len(panel) != PANEL_SIZE:
        raise ExpressionPanelSelectionError("Unexpected audit or trial-panel size")
    if len({str(row["sequence"]) for row in panel}) != PANEL_SIZE:
        raise ExpressionPanelSelectionError("Trial-panel sequences are not unique")
    if any(int(row["strong_adverse_metric_count"]) for row in panel):
        raise ExpressionPanelSelectionError("Trial panel contains a strong adverse metric")
    weak_only = [row for row in panel if int(row["favorable_family_count"]) < 1]
    if [str(row["candidate_id"]) for row in weak_only] != [
        STABLE_WORD_EXPLORATORY_CANDIDATE_ID
    ]:
        raise ExpressionPanelSelectionError("Unexpected weak-only trial candidate")
    exploratory = weak_only[0]
    if (
        not bool(exploratory["stable_word_gain_tiebreak"])
        or int(exploratory["moderate_adverse_metric_count"]) != 0
        or int(exploratory["strong_adverse_metric_count"]) != 0
    ):
        raise ExpressionPanelSelectionError("Exploratory stable-word candidate violates release rules")
    if any(int(row["hard_sequence_risk_count"]) for row in panel):
        raise ExpressionPanelSelectionError("Trial panel contains a hard sequence risk")
    if len(panel) + len(reserves) != int(facts["strict_core_count"]) + int(
        facts["controlled_tradeoff_count"]
    ) + 1:
        raise ExpressionPanelSelectionError("Eligible candidate accounting mismatch")
    if len(shortlist) != int(facts["magnitude_shortlist_count"]):
        raise ExpressionPanelSelectionError("Magnitude-shortlist accounting mismatch")


def _finite(row: Mapping[str, object], key: str) -> float:
    value = float(row[key])
    if value != value or value in (float("inf"), float("-inf")):
        raise ExpressionPanelSelectionError(f"Non-finite value for {key}")
    return value
