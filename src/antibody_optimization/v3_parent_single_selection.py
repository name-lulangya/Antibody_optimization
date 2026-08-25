"""Freeze the user-approved V3 Nb252 parent-single selection.

The input is the released 31-candidate expert-review pool plus the complete
V3 single-mutant audit.  The module never recalculates predictors and never
uses within-band decimal differences as a global score.  It records a curated
decision for every reviewed candidate while preserving the immutable upstream
selection status, including the deliberately exceptional T99F stable-word
hypothesis.

This module selects the 15 parent single mutants only.  It reports the size of
the future valid unordered pair space but does not enumerate, score, or select
double mutants.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from math import comb


V3_REVIEW_CANDIDATE_COUNT = 31
V3_SELECTED_PARENT_COUNT = 15
V3_T99F_ID = "Nb252_expr_seq099_T99F"


class V3ParentSingleSelectionError(ValueError):
    """Raised when the frozen V3 selection and its upstream evidence disagree."""


def _candidate_id(position: int, mutation: str) -> str:
    return f"Nb252_expr_seq{position:03d}_{mutation}"


# Display order follows the released review-pool order and is explicitly not an
# efficacy rank.  The set contains all nine non-hard-risk candidates with at
# least one strong-favorable U/S/Tm band, five complementary moderate-evidence
# candidates, and the user-directed T99F stable-word exploration exception.
SELECTED_PARENT_IDS = (
    _candidate_id(11, "L11Y"),
    _candidate_id(30, "F30S"),
    _candidate_id(86, "K86S"),
    _candidate_id(23, "A23R"),
    _candidate_id(5, "Q5V"),
    _candidate_id(55, "S55G"),
    _candidate_id(75, "K75A"),
    _candidate_id(29, "F29Q"),
    _candidate_id(43, "K43A"),
    _candidate_id(76, "N76G"),
    _candidate_id(30, "F30N"),
    _candidate_id(75, "K75E"),
    _candidate_id(11, "L11M"),
    _candidate_id(1, "Q1D"),
    V3_T99F_ID,
)

HARD_EXPERT_EXCLUSION_IDS = frozenset(
    {
        _candidate_id(49, "A49F"),
        _candidate_id(50, "S50F"),
        _candidate_id(71, "R71G"),
        _candidate_id(96, "A96R"),
        _candidate_id(49, "A49M"),
    }
)


_SELECTED_DECISIONS = {
    _candidate_id(11, "L11Y"): (
        "selected_strong_property_evidence",
        "NetSolP S为强改善且U为中等改善；实验坐标缺失带来的低置信度仅作标注，未达到专家硬排条件。",
    ),
    _candidate_id(30, "F30S"): (
        "selected_strong_property_evidence",
        "NetSolP S为强改善且U为中等改善；NanoMelt预测Tm仅弱下降，CDR1重排风险为中等置信度。",
    ),
    _candidate_id(86, "K86S"): (
        "selected_strong_property_evidence",
        "NanoMelt预测Tm强改善、NetSolP S中等改善且U弱改善；表面Ser保留极性，结构判断合理。",
    ),
    _candidate_id(23, "A23R"): (
        "selected_moderate_property_evidence",
        "NetSolP U和NanoMelt预测Tm均为中等改善；Cys22邻域的大正电侧链风险为中等置信度，按现行规则不作硬排。",
    ),
    _candidate_id(5, "Q5V"): (
        "selected_strong_property_evidence",
        "NanoMelt预测Tm强改善，且Q5V是天然VHH保守性合同唯一允许的亲本到共识回变；暴露疏水风险保留为谨慎项。",
    ),
    _candidate_id(55, "S55G"): (
        "selected_strong_property_evidence",
        "NanoMelt预测Tm强改善；CDR2引入Gly可能增加柔性的专家担忧为中等置信度，因此保留为模型与机理分歧的实验候选。",
    ),
    _candidate_id(75, "K75A"): (
        "selected_strong_property_evidence",
        "NanoMelt预测Tm强改善；NetSolP U弱下降被完整保留，实验结构显示该表面替换可容纳且无高置信度结构风险。",
    ),
    _candidate_id(29, "F29Q"): (
        "selected_moderate_property_evidence",
        "NetSolP S为中等改善且专家水化方向有利；虽仅有AF3局部证据，但没有性质弱下降档，保留该独立CDR1位置假设。",
    ),
    _candidate_id(43, "K43A"): (
        "selected_moderate_property_evidence",
        "NetSolP S为中等改善且NanoMelt预测Tm弱改善；表面电荷平衡风险为中等置信度，未达到硬排标准。",
    ),
    _candidate_id(76, "N76G"): (
        "selected_moderate_property_evidence",
        "NanoMelt预测Tm中等改善且NetSolP U弱改善；正phi主链构象为Gly提供明确几何依据。",
    ),
    _candidate_id(30, "F30N"): (
        "selected_strong_property_evidence",
        "NetSolP U为强改善且S为中等改善；新增Asn-Gly脱酰胺基序属于软化学风险，记录但不直接排除。",
    ),
    _candidate_id(75, "K75E"): (
        "selected_strong_property_evidence",
        "NanoMelt预测Tm强改善且U/S无明确下降档；保留为与K75A不同的表面电荷反转假设。",
    ),
    _candidate_id(11, "L11M"): (
        "selected_strong_property_evidence",
        "NetSolP U为强改善且S为弱改善；实验坐标缺失和Met氧化均为软风险，不构成高置信度硬排。",
    ),
    _candidate_id(1, "Q1D"): (
        "selected_moderate_property_evidence",
        "NetSolP S为中等改善且U弱改善，专家水化方向有利；真实N端加工语境保留为不确定性。",
    ),
    V3_T99F_ID: (
        "selected_user_directed_stable_word_exploration",
        "用户指定的稳定词gain_only探索例外；U无改善、S和Tm均弱下降且存在中等置信度结构风险，不得描述为预测器合格候选。",
    ),
}


_REJECTED_DECISIONS = {
    _candidate_id(49, "A49F"): (
        "rejected_high_confidence_expert_risk",
        "实验结构显示该Ala深埋；替换为大芳香Phe的核心过度包装风险明确，且专家结构、溶解度和热稳定判断均不利，置信度高。",
    ),
    _candidate_id(50, "S50F"): (
        "rejected_high_confidence_expert_risk",
        "该位点部分埋藏并靠近受体，Ser到大芳香Phe同时带来过度包装和疏水风险；专家不利判断置信度高。",
    ),
    _candidate_id(71, "R71G"): (
        "rejected_high_confidence_expert_risk",
        "实验结构支持内部极性网络作用；R到G同时丢失该网络并增加非转角主链柔性，两个风险机制同向且置信度高。",
    ),
    _candidate_id(96, "A96R"): (
        "rejected_high_confidence_expert_risk",
        "深埋Ala替换为大带电Arg并紧邻Cys95，存在明确的埋藏电荷、过度包装和二硫键邻域风险，置信度高。",
    ),
    _candidate_id(49, "A49M"): (
        "rejected_high_confidence_expert_risk",
        "实验结构显示该Ala深埋；替换为大Met的核心过度包装风险明确，并附加Met氧化风险，专家判断置信度高。",
    ),
    _candidate_id(1, "Q1A"): (
        "rejected_stronger_same_position_selected",
        "Q1A主要只有Tm中等改善且U弱下降；同位点Q1D具有S中等和U弱改善并有更清晰的水化解释。",
    ),
    _candidate_id(28, "I28Y"): (
        "rejected_limited_slots_weaker_evidence",
        "仅NetSolP S达到中等改善，同时U弱下降；实验坐标缺失且AF3芳香过度包装方向无法确定，优先级低于入选项。",
    ),
    _candidate_id(32, "Y32L"): (
        "rejected_limited_slots_weaker_evidence",
        "只有NanoMelt预测Tm中等改善，U/S均无改善；芳香和氢键网络丢失为中等置信度谨慎项。",
    ),
    _candidate_id(40, "A40G"): (
        "rejected_limited_slots_weaker_evidence",
        "只有NanoMelt预测Tm中等改善；部分埋藏转角处的柔性与小空腔权衡使其综合证据弱于入选候选。",
    ),
    _candidate_id(60, "A60D"): (
        "rejected_limited_slots_weaker_evidence",
        "只有NetSolP S中等改善；部分埋藏负电簇为中等置信度风险，综合证据不足以占用有限名额。",
    ),
    _candidate_id(69, "I69V"): (
        "rejected_limited_slots_weaker_evidence",
        "只有NanoMelt预测Tm中等改善且U弱下降；保守核心缩小仍可能形成小空腔，综合证据弱于入选项。",
    ),
    _candidate_id(79, "Y79T"): (
        "rejected_limited_slots_weaker_evidence",
        "NetSolP S中等改善伴NanoMelt预测Tm弱下降；部分埋藏芳香空腔风险使其优先级降低。",
    ),
    _candidate_id(83, "N83A"): (
        "rejected_limited_slots_weaker_evidence",
        "NanoMelt预测Tm中等改善伴NetSolP S弱下降；还可能丢失beta片层边缘极性网络。",
    ),
    _candidate_id(99, "T99N"): (
        "rejected_stronger_same_position_selected",
        "T99N仅S中等改善且Tm弱下降；同位点名额按用户决定用于T99F稳定词探索，不再额外保留T99N。",
    ),
    _candidate_id(86, "K86A"): (
        "rejected_stronger_same_position_selected",
        "同位点K86S具有Tm强改善、S中等和U弱改善，并保留表面极性；K86A的S和Tm仅为中等改善。",
    ),
    _candidate_id(23, "A23Q"): (
        "rejected_stronger_same_position_selected",
        "A23Q只有U中等改善；同位点A23R同时具有U和Tm中等改善，且其中等置信度风险按现行规则不作硬排。",
    ),
}


def build_v3_parent_single_selection(
    review_rows: Sequence[Mapping[str, object]],
    complete_audit_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Return the 31-row decision audit, selected 15 rows, and stage facts."""

    review = [dict(row) for row in review_rows]
    if len(review) != V3_REVIEW_CANDIDATE_COUNT:
        raise V3ParentSingleSelectionError("Expected exactly 31 expert-review rows")
    review_ids = [str(row.get("candidate_id", "")).strip() for row in review]
    if any(not value for value in review_ids) or len(set(review_ids)) != len(review_ids):
        raise V3ParentSingleSelectionError("Expert-review candidate IDs are not unique")

    decisions = {**_SELECTED_DECISIONS, **_REJECTED_DECISIONS}
    if set(review_ids) != set(decisions):
        raise V3ParentSingleSelectionError(
            "Frozen selection decisions do not exactly cover the 31-candidate review pool"
        )
    if len(SELECTED_PARENT_IDS) != V3_SELECTED_PARENT_COUNT or set(
        SELECTED_PARENT_IDS
    ) != set(_SELECTED_DECISIONS):
        raise V3ParentSingleSelectionError("Selected-parent identity contract is invalid")

    complete_lookup: dict[str, dict[str, object]] = {}
    for source in complete_audit_rows:
        identifier = str(source.get("candidate_id", "")).strip()
        if identifier in review_ids:
            if identifier in complete_lookup:
                raise V3ParentSingleSelectionError(
                    f"Complete V3 audit contains duplicate candidate {identifier}"
                )
            complete_lookup[identifier] = dict(source)
    if set(complete_lookup) != set(review_ids):
        raise V3ParentSingleSelectionError(
            "Complete V3 audit does not cover every expert-review candidate"
        )

    selected_order = {
        identifier: index for index, identifier in enumerate(SELECTED_PARENT_IDS, start=1)
    }
    output: list[dict[str, object]] = []
    for source in review:
        identifier = str(source["candidate_id"])
        audit_source = complete_lookup[identifier]
        if str(audit_source.get("sequence", "")) != str(source.get("sequence", "")):
            raise V3ParentSingleSelectionError(
                f"Expert review and complete audit sequence mismatch for {identifier}"
            )
        selected = identifier in selected_order
        decision_class, reason = decisions[identifier]
        hard_exclusion = identifier in HARD_EXPERT_EXCLUSION_IDS
        if selected and hard_exclusion:
            raise V3ParentSingleSelectionError(
                f"High-confidence expert-risk candidate cannot be selected: {identifier}"
            )
        row = dict(source)
        row.update(
            {
                "v3_parent_selection_status": "selected" if selected else "not_selected",
                "v3_parent_panel_order_not_efficacy_rank": selected_order.get(identifier, ""),
                "v3_parent_decision_class": decision_class,
                "v3_parent_high_confidence_expert_risk_exclusion": hard_exclusion,
                "v3_parent_decision_reason_cn": reason,
                "upstream_selection_eligibility_v3": audit_source.get(
                    "selection_eligibility_v3", ""
                ),
                "upstream_selection_tier_v3": audit_source.get("selection_tier_v3", ""),
                "upstream_selection_status_v3": audit_source.get(
                    "selection_status_v3", ""
                ),
                "upstream_selection_reason_v3": audit_source.get(
                    "selection_reason_v3", ""
                ),
                "upstream_strong_positive_metric_count_v3": audit_source.get(
                    "strong_positive_metric_count_v3", ""
                ),
            }
        )
        output.append(row)

    selected_rows = sorted(
        (row for row in output if row["v3_parent_selection_status"] == "selected"),
        key=lambda row: int(row["v3_parent_panel_order_not_efficacy_rank"]),
    )
    _validate_selection(output, selected_rows)
    position_counts = Counter(
        int(row["reported_sequence_index_1based"]) for row in selected_rows
    )
    invalid_same_position_pairs = sum(comb(count, 2) for count in position_counts.values())
    all_pairs = comb(len(selected_rows), 2)
    facts = {
        "review_candidate_count": len(output),
        "selected_parent_single_count": len(selected_rows),
        "not_selected_candidate_count": len(output) - len(selected_rows),
        "high_confidence_expert_risk_exclusion_count": len(HARD_EXPERT_EXCLUSION_IDS),
        "competitive_not_selected_count": len(output)
        - len(selected_rows)
        - len(HARD_EXPERT_EXCLUSION_IDS),
        "selected_with_strong_favorable_metric_count": sum(
            any(
                row[field] == "strong_favorable"
                for field in (
                    "netsolp_u_band_v3",
                    "netsolp_s_band_v3",
                    "nanomelt_tm_band_v3",
                )
            )
            for row in selected_rows
        ),
        "selected_user_directed_exploration_count": sum(
            row["v3_parent_decision_class"]
            == "selected_user_directed_stable_word_exploration"
            for row in selected_rows
        ),
        "selected_unique_position_count": len(position_counts),
        "selected_same_position_groups": {
            str(position): count
            for position, count in sorted(position_counts.items())
            if count > 1
        },
        "theoretical_unordered_pair_count": all_pairs,
        "invalid_same_position_pair_count": invalid_same_position_pairs,
        "valid_unordered_double_mutant_count": all_pairs - invalid_same_position_pairs,
        "double_mutant_enumeration_performed": False,
    }
    return {"audit_rows": output, "selected_rows": selected_rows, "facts": facts}


def selected_parent_export_rows(
    selected_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Return a report-ready projection of the approved 15 parent singles."""

    fields = (
        "v3_parent_panel_order_not_efficacy_rank",
        "candidate_id",
        "mutation_reported_label",
        "reported_sequence_index_1based",
        "imgt_position_label",
        "region",
        "sequence",
        "netsolp_delta_u",
        "netsolp_u_band_v3",
        "netsolp_delta_s",
        "netsolp_s_band_v3",
        "nanomelt_delta_tm_c",
        "nanomelt_tm_band_v3",
        "antifold_selection_source",
        "antifold_delta_logp",
        "antifold_mutant_rank_worst_first",
        "antifold_veto_status",
        "stable_word_effect",
        "expert_structural_assessment",
        "expert_solubility_expectation",
        "expert_thermal_stability_expectation",
        "expert_confidence",
        "expert_primary_concern",
        "expert_rationale_cn",
        "v3_parent_decision_class",
        "v3_parent_decision_reason_cn",
        "upstream_selection_eligibility_v3",
        "upstream_selection_status_v3",
    )
    return [{field: row[field] for field in fields} for row in selected_rows]


def _validate_selection(
    audit_rows: Sequence[Mapping[str, object]],
    selected_rows: Sequence[Mapping[str, object]],
) -> None:
    if len(selected_rows) != V3_SELECTED_PARENT_COUNT:
        raise V3ParentSingleSelectionError("Expected exactly 15 selected parents")
    selected_ids = {str(row["candidate_id"]) for row in selected_rows}
    if selected_ids != set(SELECTED_PARENT_IDS):
        raise V3ParentSingleSelectionError("Selected-parent identities changed")
    hard_rows = {
        str(row["candidate_id"])
        for row in audit_rows
        if row["v3_parent_high_confidence_expert_risk_exclusion"] is True
    }
    if hard_rows != HARD_EXPERT_EXCLUSION_IDS:
        raise V3ParentSingleSelectionError("High-confidence risk exclusion set changed")
    for row in audit_rows:
        hard = str(row["candidate_id"]) in HARD_EXPERT_EXCLUSION_IDS
        if hard and not (
            row["expert_structural_assessment"] == "structurally_concerning"
            and row["expert_confidence"] == "high"
        ):
            raise V3ParentSingleSelectionError(
                f"Hard expert exclusion lacks high-confidence structural concern: {row['candidate_id']}"
            )
    if any(row["antifold_veto_status"] != "pass" for row in audit_rows):
        raise V3ParentSingleSelectionError("The 31-candidate pool must pass AntiFold veto")
    if len({str(row["sequence"]) for row in selected_rows}) != len(selected_rows) or any(
        len(str(row["sequence"])) != 128 for row in selected_rows
    ):
        raise V3ParentSingleSelectionError("Selected parent sequences are not unique 128-aa constructs")
    reconstructed_parents: set[str] = set()
    for row in selected_rows:
        sequence = str(row["sequence"])
        position = int(row["reported_sequence_index_1based"])
        wt = str(row["wt_residue"])
        mutant = str(row["mutant_residue"])
        if sequence[position - 1] != mutant or mutant == wt:
            raise V3ParentSingleSelectionError(
                f"Selected sequence does not encode its declared substitution: {row['candidate_id']}"
            )
        parent = sequence[: position - 1] + wt + sequence[position:]
        reconstructed_parents.add(parent)
        if not sequence.endswith("SSGS") or sequence.count("C") != 2:
            raise V3ParentSingleSelectionError(
                f"Selected sequence violates SSGS/Cys construct constraints: {row['candidate_id']}"
            )
    if len(reconstructed_parents) != 1:
        raise V3ParentSingleSelectionError(
            "Selected substitutions do not reconstruct one common Nb252 parent"
        )
    regular_selected = [
        row for row in selected_rows if row["candidate_id"] != V3_T99F_ID
    ]
    if any(
        not any(
            row[field] in {"moderate_favorable", "strong_favorable"}
            for field in (
                "netsolp_u_band_v3",
                "netsolp_s_band_v3",
                "nanomelt_tm_band_v3",
            )
        )
        for row in regular_selected
    ):
        raise V3ParentSingleSelectionError(
            "Every regular selected parent must retain moderate/strong property evidence"
        )
    strong_pool = {
        str(row["candidate_id"])
        for row in audit_rows
        if any(
            row[field] == "strong_favorable"
            for field in (
                "netsolp_u_band_v3",
                "netsolp_s_band_v3",
                "nanomelt_tm_band_v3",
            )
        )
    }
    if len(strong_pool) != 11 or strong_pool - selected_ids != {
        _candidate_id(49, "A49F"),
        _candidate_id(49, "A49M"),
    }:
        raise V3ParentSingleSelectionError(
            "Strong-property candidate retention/exclusion contract changed"
        )
    t99f = next(row for row in audit_rows if row["candidate_id"] == V3_T99F_ID)
    if not (
        t99f["v3_parent_selection_status"] == "selected"
        and t99f["v3_parent_decision_class"]
        == "selected_user_directed_stable_word_exploration"
        and t99f["stable_word_effect"] == "gain_only"
        and t99f["upstream_selection_tier_v3"] == "not_eligible"
        and t99f["upstream_selection_status_v3"] == "not_selected"
        and all(
            t99f[field] not in {"moderate_favorable", "strong_favorable"}
            for field in (
                "netsolp_u_band_v3",
                "netsolp_s_band_v3",
                "nanomelt_tm_band_v3",
            )
        )
    ):
        raise V3ParentSingleSelectionError("T99F exploratory exception contract changed")


__all__ = [
    "HARD_EXPERT_EXCLUSION_IDS",
    "SELECTED_PARENT_IDS",
    "V3ParentSingleSelectionError",
    "V3_REVIEW_CANDIDATE_COUNT",
    "V3_SELECTED_PARENT_COUNT",
    "V3_T99F_ID",
    "build_v3_parent_single_selection",
    "selected_parent_export_rows",
]
