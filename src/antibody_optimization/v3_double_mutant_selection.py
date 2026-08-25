"""Build the V3 double-mutant review, final 15 doubles, and 15+15 panel.

The module consumes the released 102-row V3 property matrix, the released
15-parent panel, the 31-row parent decision audit, and the post-sync annotation
erratum.  It does not run a predictor, model a double-mutant side chain, or use
historical V1/V2 selection semantics.

All 102 doubles receive the same expert-review fields and the same eligibility
rules.  ``enhanced`` versus ``standard`` controls review detail only.  T99F and
stable-word gains do not trigger a quota, bonus, penalty, veto, or forced slot.
The selected set is an explicit human-reviewed decision over predeclared
magnitude bands; raw decimals are retained as evidence but are not combined
into a continuous ranking score.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from itertools import combinations


SEPARATED_PAIR_CLASS = "spatially_separated_ca_at_least_10A"
T99F_MUTATION = "T99F"

SELECTED_DOUBLE_MUTATION_SETS: tuple[str, ...] = (
    "F30S;Q5V",
    "S55G;K43A",
    "K86S;Q5V",
    "L11Y;K86S",
    "S55G;F30N",
    "N76G;L11M",
    "F30S;K75E",
    "K86S;K43A",
    "A23R;S55G",
    "K43A;N76G",
    "K75E;Q1D",
    "L11Y;K75A",
    "Q5V;N76G",
    "L11Y;Q1D",
    "F30S;Q1D",
)

SELECTED_MECHANISM_RATIONALE_CN: dict[str, str] = {
    "F30S;Q5V": "CDR1表面去疏水与FR1天然共识回变属于不同机制，两个位点空间分离且无责任基序风险。",
    "S55G;K43A": "CDR2局部柔性假设与FR2表面电荷调整空间独立；保留为预测Tm与溶解度信号一致的组合。",
    "K86S;Q5V": "FR3表面Lys到Ser与FR1共识回变空间分离，兼顾表面极性和预测热稳定性。",
    "L11Y;K86S": "两个外露位点的水化与侧链熵假设互补；L11仅有AF3坐标，因此保留低置信度结构限制。",
    "S55G;F30N": "三项性质幅度档均为强改善，虽有F30N新增NG脱酰胺基序，仍作为高收益、明确软风险假设保留。",
    "N76G;L11M": "转角Gly适配假设与FR1保守疏水替换空间分离；L11坐标缺失和Met氧化敏感性均保留为限制。",
    "F30S;K75E": "CDR1去疏水与FR3表面电荷反转空间分离，S和预测Tm提供互补支持。",
    "K86S;K43A": "两个外露框架电荷调整位点相距较远，未发现责任基序或局部耦联触发。",
    "A23R;S55G": "补充reported 23机制并保留强预测Tm信号；A23R邻近二硫键的担忧为中等置信度，未作为硬排除。",
    "K43A;N76G": "两个FR位点空间分离，表面电荷调整与转角构象假设互补且无责任基序风险。",
    "K75E;Q1D": "FR3电荷反转与N端水化假设空间分离；保留电荷环境依赖这一专家不确定性。",
    "L11Y;K75A": "FR1外露水化和FR3表面去电荷假设互补；L11实验坐标缺失，结构结论限于AF3。",
    "Q5V;N76G": "天然共识回变与正phi转角Gly假设空间分离，并提供明显预测Tm支持。",
    "L11Y;Q1D": "两个水化相关假设位于远距离位置，U和S提供互补支持；L11仍为AF3-only证据。",
    "F30S;Q1D": "CDR1去疏水与N端负电荷假设空间分离，U和S均达到中等改善档。",
}

BAND_CN = {
    "strong_favorable": "强改善",
    "moderate_favorable": "中等改善",
    "weak_favorable": "轻微改善",
    "negligible": "近似中性",
    "weak_adverse": "轻微不利",
    "moderate_adverse": "中等不利",
    "strong_adverse": "明显不利",
}

PARENT_EXPERT_FIELDS = (
    "expert_structural_assessment",
    "expert_solubility_expectation",
    "expert_thermal_stability_expectation",
    "expert_confidence",
    "expert_primary_concern",
    "expert_rationale_cn",
    "expert_uncertainty_cn",
    "primary_structure_source",
    "experimental_coordinate_status",
    "near_interface_shell_status",
    "manual_visual_review_status",
)


def build_v3_double_mutant_selection(
    matrix_rows: Sequence[Mapping[str, object]],
    parent_selected_rows: Sequence[Mapping[str, object]],
    parent_audit_rows: Sequence[Mapping[str, object]],
    post_sync_review: Mapping[str, object],
) -> dict[str, object]:
    """Return the complete V3 review audit, selected doubles, and final panel.

    Args:
        matrix_rows: Released 102-row complete V3 double-mutant property matrix.
        parent_selected_rows: Released ordered 15-parent export.
        parent_audit_rows: Released 31-row parent decision and expert audit.
        post_sync_review: Released post-sync JSON containing annotation errata.

    Returns:
        A mapping with ``audit_rows``, ``selected_double_rows``,
        ``final_panel_rows``, and compact ``facts``.

    This function deliberately does not infer physical epistasis from model
    non-additivity and does not claim that double-mutant side chains were
    structurally modeled.
    """

    matrix = [dict(row) for row in matrix_rows]
    selected_parents = [dict(row) for row in parent_selected_rows]
    parent_audit = [dict(row) for row in parent_audit_rows]
    _validate_inputs(matrix, selected_parents, parent_audit, post_sync_review)

    audit_by_id = {str(row["candidate_id"]): row for row in parent_audit}
    errata = _parse_errata(post_sync_review)
    selected_order = {
        mutation_set: index
        for index, mutation_set in enumerate(SELECTED_DOUBLE_MUTATION_SETS, start=1)
    }
    available_sets = {str(row["mutation_set"]) for row in matrix}
    missing_selected = set(selected_order) - available_sets
    if missing_selected:
        raise ValueError(f"Selected V3 doubles missing from matrix: {sorted(missing_selected)}")

    audit_rows: list[dict[str, object]] = []
    for source in matrix:
        row = dict(source)
        candidate_id = str(row["double_candidate_id"])
        parent_a = audit_by_id[str(row["parent_a_candidate_id"])]
        parent_b = audit_by_id[str(row["parent_b_candidate_id"])]
        effective_flags, effective_count, erratum_applied = _effective_soft_risk(
            row, errata.get(candidate_id)
        )
        review_depth, review_triggers = _review_depth(row)
        property_layer = _property_evidence_layer(row)
        expert = _combined_expert_assessment(
            row,
            parent_a,
            parent_b,
            effective_flags=effective_flags,
            review_depth=review_depth,
            review_triggers=review_triggers,
        )
        mutation_set = str(row["mutation_set"])
        is_selected = mutation_set in selected_order
        if is_selected:
            decision_class = "selected_final15_double"
            decision_reason = _selected_reason(row, expert, effective_flags)
            selection_order: object = selected_order[mutation_set]
        else:
            decision_class, decision_reason = _not_selected_reason(
                row, expert, effective_flags
            )
            selection_order = ""

        row.update(
            {
                "source_soft_sequence_risk_flags": row["soft_sequence_risk_flags"],
                "source_soft_sequence_risk_count": row["soft_sequence_risk_count"],
                "soft_sequence_risk_flags": effective_flags,
                "soft_sequence_risk_count": effective_count,
                "effective_soft_sequence_risk_flags": effective_flags,
                "effective_soft_sequence_risk_count": effective_count,
                "post_sync_annotation_erratum_applied": erratum_applied,
                "expert_review_depth": review_depth,
                "expert_review_triggers": "|".join(review_triggers),
                "expert_review_completed": True,
                "property_evidence_layer": property_layer,
                "double_expert_assessment": expert["structural"],
                "double_expert_solubility_interpretation": expert["solubility"],
                "double_expert_thermal_interpretation": expert["thermal"],
                "double_expert_confidence": expert["confidence"],
                "double_expert_primary_concern": expert["concern"],
                "double_expert_rationale_cn": expert["rationale"],
                "double_expert_uncertainty_cn": expert["uncertainty"],
                "final_double_selection_status": (
                    "selected" if is_selected else "not_selected"
                ),
                "final_double_panel_order_not_efficacy_rank": selection_order,
                "final_double_decision_class": decision_class,
                "final_double_decision_reason_cn": decision_reason,
                "t99f_specific_selection_rule_applied": False,
                "stable_word_specific_selection_rule_applied": False,
            }
        )
        for prefix, parent in (("parent_a", parent_a), ("parent_b", parent_b)):
            for field in PARENT_EXPERT_FIELDS:
                row[f"{prefix}_{field}"] = parent[field]
        audit_rows.append(row)

    selected_rows = sorted(
        (row for row in audit_rows if row["final_double_selection_status"] == "selected"),
        key=lambda row: int(row["final_double_panel_order_not_efficacy_rank"]),
    )
    final_panel = _build_final_panel(selected_parents, selected_rows)
    facts = _validate_outputs(audit_rows, selected_rows, final_panel)
    return {
        "audit_rows": audit_rows,
        "selected_double_rows": selected_rows,
        "final_panel_rows": final_panel,
        "selection_policy": selection_policy(),
        "facts": facts,
    }


def selection_policy() -> dict[str, object]:
    """Return the frozen V3 double-selection semantics."""

    return {
        "positive_metrics": ["NetSolP U", "NetSolP S", "NanoMelt predicted Tm"],
        "magnitude_bands_drive_selection": True,
        "within_band_raw_decimals_used_as_rank": False,
        "mutation_specific_quota_or_exception": False,
        "stable_word_role": "uniform_soft_evidence_only_not_a_selection_requirement",
        "review_depth_role": "documentation_depth_only_not_eligibility_or_rank",
        "expert_hard_exclusion_role": (
            "only_high_confidence_concrete_physical_risk; none newly asserted from "
            "unmodeled double side chains"
        ),
        "antifold_role": (
            "constituent_negative_veto_only_no_double_score_no_positive_rank"
        ),
        "predictors_rerun": False,
        "double_sidechain_modeling_performed": False,
        "portfolio_constraints": {
            "maximum_exact_parent_usage": 3,
            "maximum_reported_position_usage": 4,
            "maximum_unordered_position_pair_usage": 1,
            "minimum_parent_component_coverage": 10,
            "minimum_reported_position_coverage": 9,
            "maximum_af3_only_pairs": 4,
            "maximum_soft_risk_pairs": 2,
            "maximum_local_pairs": 1,
        },
    }


def selected_double_export_rows(
    selected_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Return a compact, report-ready projection of the selected 15 doubles."""

    fields = (
        "final_double_panel_order_not_efficacy_rank",
        "double_candidate_id",
        "mutation_set",
        "mutation_a",
        "mutation_b",
        "parent_a_candidate_id",
        "parent_b_candidate_id",
        "position_a_reported_1based",
        "position_b_reported_1based",
        "region_a",
        "region_b",
        "sequence",
        "netsolp_u_delta_vs_wt",
        "netsolp_u_magnitude_band",
        "netsolp_s_delta_vs_wt",
        "netsolp_s_magnitude_band",
        "nanomelt_tm_c_delta_vs_wt",
        "nanomelt_tm_c_magnitude_band",
        "moderate_or_strong_favorable_metric_count",
        "antifold_constituent_gate",
        "effective_soft_sequence_risk_flags",
        "stable_word_effect",
        "pair_structure_distance_source",
        "pair_ca_distance_a",
        "pair_spatial_class",
        "expert_review_depth",
        "double_expert_assessment",
        "double_expert_solubility_interpretation",
        "double_expert_thermal_interpretation",
        "double_expert_confidence",
        "double_expert_rationale_cn",
        "double_expert_uncertainty_cn",
        "final_double_decision_reason_cn",
    )
    return [{field: row[field] for field in fields} for row in selected_rows]


def _validate_inputs(
    matrix: Sequence[Mapping[str, object]],
    selected_parents: Sequence[Mapping[str, object]],
    parent_audit: Sequence[Mapping[str, object]],
    post_sync_review: Mapping[str, object],
) -> None:
    if len(matrix) != 102:
        raise ValueError("Expected exactly 102 V3 double-mutant matrix rows")
    ids = [str(row["double_candidate_id"]) for row in matrix]
    sequences = [str(row["sequence"]) for row in matrix]
    if len(set(ids)) != 102 or len(set(sequences)) != 102:
        raise ValueError("V3 double-mutant IDs and sequences must each be unique")
    if len(selected_parents) != 15 or len(parent_audit) != 31:
        raise ValueError("Expected 15 selected parents and the complete 31-row audit")
    selected_ids = {str(row["candidate_id"]) for row in selected_parents}
    if len(selected_ids) != 15:
        raise ValueError("Selected-parent IDs must be unique")
    audit_by_id = {str(row["candidate_id"]): row for row in parent_audit}
    if len(audit_by_id) != 31 or not selected_ids.issubset(audit_by_id):
        raise ValueError("Parent audit does not contain the exact selected-parent set")
    for row in selected_parents:
        audit = audit_by_id[str(row["candidate_id"])]
        if str(row["sequence"]) != str(audit["sequence"]):
            raise ValueError(f"Parent sequence identity mismatch: {row['candidate_id']}")
    for row in matrix:
        if str(row["parent_a_candidate_id"]) not in selected_ids or str(
            row["parent_b_candidate_id"]
        ) not in selected_ids:
            raise ValueError("Double includes a component outside the selected parent set")
        if int(row["position_a_reported_1based"]) == int(
            row["position_b_reported_1based"]
        ):
            raise ValueError("Same-position alternatives cannot form a V3 double")
        if str(row["antifold_constituent_gate"]) != "pass":
            raise ValueError("Every double must retain two passing constituent AntiFold gates")
        if str(row["antifold_double_mutant_scored"]).lower() != "false":
            raise ValueError("Double-mutant AntiFold scores are outside the V3 contract")
        if str(row["antifold_component_values_combined"]).lower() != "false":
            raise ValueError("Constituent AntiFold values must not be combined")
        if str(row["antifold_double_mutant_score"]) != "":
            raise ValueError("Unexpected double-mutant AntiFold score")
    if str(post_sync_review.get("status", "")) != "pass_with_one_annotation_erratum":
        raise ValueError("Expected the released post-sync annotation review")


def _parse_errata(review: Mapping[str, object]) -> dict[str, Mapping[str, object]]:
    entries = review.get("annotation_errata")
    if not isinstance(entries, list) or len(entries) != 1:
        raise ValueError("Expected exactly one released annotation erratum")
    output: dict[str, Mapping[str, object]] = {}
    for entry in entries:
        if not isinstance(entry, Mapping):
            raise ValueError("Invalid annotation erratum")
        candidate_id = str(entry["double_candidate_id"])
        output[candidate_id] = entry
    return output


def _effective_soft_risk(
    row: Mapping[str, object], erratum: Mapping[str, object] | None
) -> tuple[str, int, bool]:
    source_flags = str(row["soft_sequence_risk_flags"])
    source_count = int(row["soft_sequence_risk_count"])
    if erratum is None:
        return source_flags, source_count, False
    if source_flags != str(erratum["source_value"]) or source_count != int(
        erratum["source_soft_sequence_risk_count"]
    ):
        raise ValueError("Annotation erratum does not match its source matrix row")
    return (
        str(erratum["corrected_value"]),
        int(erratum["corrected_soft_sequence_risk_count"]),
        True,
    )


def _review_depth(row: Mapping[str, object]) -> tuple[str, tuple[str, ...]]:
    favorable = int(row["moderate_or_strong_favorable_metric_count"])
    adverse = int(row["moderate_adverse_metric_count"]) + int(
        row["strong_adverse_metric_count"]
    )
    triggers: list[str] = []
    if favorable >= 2 and adverse == 0:
        triggers.append("multi_metric_property_evidence")
    if adverse > 0:
        triggers.append("moderate_or_strong_property_tradeoff")
    if str(row["pair_spatial_class"]) != SEPARATED_PAIR_CLASS:
        triggers.append("local_pair_geometry")
    return (
        ("enhanced", tuple(triggers))
        if triggers
        else ("standard", ("standard_common_review_no_enhanced_trigger",))
    )


def _property_evidence_layer(row: Mapping[str, object]) -> str:
    favorable = int(row["moderate_or_strong_favorable_metric_count"])
    moderate_adverse = int(row["moderate_adverse_metric_count"])
    strong_adverse = int(row["strong_adverse_metric_count"])
    if strong_adverse:
        return "strong_adverse_tradeoff"
    if moderate_adverse:
        return "moderate_adverse_tradeoff"
    return {
        3: "three_moderate_or_strong_favorable",
        2: "two_moderate_or_strong_favorable",
        1: "one_moderate_or_strong_favorable",
        0: "no_moderate_or_strong_favorable",
    }[favorable]


def _combined_expert_assessment(
    row: Mapping[str, object],
    parent_a: Mapping[str, object],
    parent_b: Mapping[str, object],
    *,
    effective_flags: str,
    review_depth: str,
    review_triggers: Sequence[str],
) -> dict[str, str]:
    parent_structural = {
        str(parent_a["expert_structural_assessment"]),
        str(parent_b["expert_structural_assessment"]),
    }
    missing = str(row["pair_experimental_coordinate_status"]) != "both_observed"
    local = str(row["pair_spatial_class"]) != SEPARATED_PAIR_CLASS
    if local or "structurally_concerning" in parent_structural:
        structural = "reasonable_with_caution"
    elif parent_structural == {"reasonable"}:
        structural = "reasonable"
    else:
        structural = "reasonable_with_caution"
    confidence = "low" if missing or "low" in {
        str(parent_a["expert_confidence"]),
        str(parent_b["expert_confidence"]),
    } else "medium"

    u_band = str(row["netsolp_u_magnitude_band"])
    s_band = str(row["netsolp_s_magnitude_band"])
    tm_band = str(row["nanomelt_tm_c_magnitude_band"])
    sol_adverse = u_band in {"moderate_adverse", "strong_adverse"} or s_band in {
        "moderate_adverse",
        "strong_adverse",
    }
    sol_positive = u_band in {"moderate_favorable", "strong_favorable"} or s_band in {
        "moderate_favorable",
        "strong_favorable",
    }
    parent_sol_unfavorable = "unfavorable" in {
        str(parent_a["expert_solubility_expectation"]),
        str(parent_b["expert_solubility_expectation"]),
    }
    if sol_adverse:
        solubility = "adverse_or_tradeoff_predicted"
    elif sol_positive and parent_sol_unfavorable:
        solubility = "predictor_supported_with_expert_caveat"
    elif sol_positive:
        solubility = "predictor_supported"
    else:
        solubility = "neutral_or_uncertain"

    parent_tm_unfavorable = "unfavorable" in {
        str(parent_a["expert_thermal_stability_expectation"]),
        str(parent_b["expert_thermal_stability_expectation"]),
    }
    if tm_band in {"moderate_adverse", "strong_adverse"}:
        thermal = "adverse_or_tradeoff_predicted"
    elif tm_band in {"moderate_favorable", "strong_favorable"} and parent_tm_unfavorable:
        thermal = "predictor_supported_with_expert_caveat"
    elif tm_band in {"moderate_favorable", "strong_favorable"}:
        thermal = "predictor_supported"
    else:
        thermal = "neutral_or_uncertain"

    concerns: list[str] = []
    if missing:
        concerns.append("experimental_coordinates_missing_for_at_least_one_site")
    if local:
        concerns.append("local_pair_coupling_not_double_mutant_modeled")
    if effective_flags:
        concerns.append(effective_flags.replace("|", "+"))
    for parent in (parent_a, parent_b):
        concern = str(parent["expert_primary_concern"])
        if concern and concern not in concerns:
            concerns.append(concern)
    if not concerns:
        concerns.append("no_specific_pair_level_risk_beyond_model_uncertainty")

    band_text = (
        f"U={BAND_CN[u_band]}、S={BAND_CN[s_band]}、预测Tm={BAND_CN[tm_band]}"
    )
    structure_text = (
        f"WT位点几何采用{row['pair_structure_distance_source']}，空间类别为"
        f"{row['pair_spatial_class']}；该信息不是突变后双侧链结构。"
    )
    parent_text = (
        f"组成单突专家判断分别为{parent_a['expert_structural_assessment']}和"
        f"{parent_b['expert_structural_assessment']}。"
    )
    risk_text = (
        f"有效软责任基序：{effective_flags}。" if effective_flags else "未检出声明的软责任基序。"
    )
    detail_text = (
        f"增强审查触发：{'、'.join(review_triggers)}。" if review_depth == "enhanced"
        else "按统一字段完成标准审查；标准深度不表示自动淘汰。"
    )
    return {
        "structural": structural,
        "solubility": solubility,
        "thermal": thermal,
        "confidence": confidence,
        "concern": "|".join(concerns),
        "rationale": (
            f"{row['mutation_set']}：{band_text}；{parent_text}{structure_text}"
            f"{risk_text}{detail_text}两个组成单突AntiFold否决门均通过，但没有双突AntiFold分数。"
        ),
        "uncertainty": (
            "未进行双突侧链重建、ChimeraX双突rotamer检查或实验验证；结构推断仅来自WT位点几何、"
            "父单突专家审查和完整双突序列责任基序。预测非加和残差只描述模型输出，不能称为物理上位性。"
        ),
    }


def _selected_reason(
    row: Mapping[str, object], expert: Mapping[str, str], effective_flags: str
) -> str:
    mutation_set = str(row["mutation_set"])
    mechanism = SELECTED_MECHANISM_RATIONALE_CN[mutation_set]
    favorable = int(row["moderate_or_strong_favorable_metric_count"])
    risk = f"有效软风险为{effective_flags}；" if effective_flags else "无有效软责任基序；"
    return (
        f"入选最终15条双突：{favorable}个独立性质指标达到中等或强改善，且无中等或明显不利；"
        f"{risk}{mechanism}综合结构判断为{expert['structural']}、置信度{expert['confidence']}。"
    )


def _not_selected_reason(
    row: Mapping[str, object], expert: Mapping[str, str], effective_flags: str
) -> tuple[str, str]:
    favorable = int(row["moderate_or_strong_favorable_metric_count"])
    adverse = int(row["moderate_adverse_metric_count"]) + int(
        row["strong_adverse_metric_count"]
    )
    local = str(row["pair_spatial_class"]) != SEPARATED_PAIR_CLASS
    missing = str(row["pair_experimental_coordinate_status"]) != "both_observed"
    if adverse:
        decision = "not_selected_property_tradeoff"
        reason = "存在中等或明显不利性质档；在有足够无该权衡的候选时未进入最终15条。"
    elif favorable == 0:
        decision = "not_selected_no_moderate_or_strong_positive"
        reason = "三个性质指标均未达到中等或强改善，综合正向证据不足。"
    elif favorable == 1:
        decision = "not_selected_less_complete_property_support"
        reason = "仅一个性质指标达到中等或强改善，弱于入选候选的双指标或三指标支持。"
    elif effective_flags:
        decision = "not_selected_soft_liability_or_portfolio_balance"
        reason = "虽有多指标正向证据，但含软责任基序；同机制已有证据相当且风险更低的组合。"
    elif local:
        decision = "not_selected_local_pair_coupling_uncertainty"
        reason = "虽有多指标正向证据，但两位点局部邻近且未完成双突侧链建模；优先选择空间分离组合。"
    elif missing:
        decision = "not_selected_structural_uncertainty_or_portfolio_balance"
        reason = "虽有多指标正向证据，但至少一个位点缺少实验坐标；AF3-only名额用于证据更完整的互补组合。"
    else:
        decision = "not_selected_portfolio_diversity_or_redundancy"
        reason = "具备多指标正向证据，但在父单突使用上限、位置覆盖和位置对去冗余后未进入15个名额。"
    return (
        decision,
        f"未入选最终15条双突：{reason}专家综合判断为{expert['structural']}、置信度{expert['confidence']}。",
    )


def _build_final_panel(
    selected_parents: Sequence[Mapping[str, object]],
    selected_doubles: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    parents = sorted(
        selected_parents,
        key=lambda row: int(row["v3_parent_panel_order_not_efficacy_rank"]),
    )
    output: list[dict[str, object]] = []
    for panel_order, row in enumerate(parents, start=1):
        output.append(
            {
                "final_panel_order_not_efficacy_rank": panel_order,
                "candidate_kind": "single_mutant",
                "candidate_id": row["candidate_id"],
                "mutation_set": str(row["mutation_reported_label"]).replace(
                    "Nb252 reported_seq ", ""
                ),
                "reported_positions_1based": row["reported_sequence_index_1based"],
                "regions": row["region"],
                "component_candidate_ids": row["candidate_id"],
                "sequence": row["sequence"],
                "selection_reason_cn": row["v3_parent_decision_reason_cn"],
            }
        )
    for index, row in enumerate(selected_doubles, start=1):
        output.append(
            {
                "final_panel_order_not_efficacy_rank": 15 + index,
                "candidate_kind": "double_mutant",
                "candidate_id": row["double_candidate_id"],
                "mutation_set": row["mutation_set"],
                "reported_positions_1based": (
                    f"{row['position_a_reported_1based']};{row['position_b_reported_1based']}"
                ),
                "regions": f"{row['region_a']};{row['region_b']}",
                "component_candidate_ids": (
                    f"{row['parent_a_candidate_id']};{row['parent_b_candidate_id']}"
                ),
                "sequence": row["sequence"],
                "selection_reason_cn": row["final_double_decision_reason_cn"],
            }
        )
    return output


def _validate_outputs(
    audit: Sequence[Mapping[str, object]],
    selected: Sequence[Mapping[str, object]],
    final_panel: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    depth_counts = Counter(str(row["expert_review_depth"]) for row in audit)
    trigger_counts = Counter(
        trigger
        for row in audit
        for trigger in str(row["expert_review_triggers"]).split("|")
        if trigger
    )
    if depth_counts != {"enhanced": 58, "standard": 44}:
        raise ValueError(f"Unexpected V3 review-depth counts: {depth_counts}")
    expected_triggers = {
        "multi_metric_property_evidence": 42,
        "moderate_or_strong_property_tradeoff": 7,
        "local_pair_geometry": 12,
        "standard_common_review_no_enhanced_trigger": 44,
    }
    if dict(trigger_counts) != expected_triggers:
        raise ValueError(f"Unexpected V3 review triggers: {trigger_counts}")
    t99f_rows = [
        row
        for row in audit
        if T99F_MUTATION in {str(row["mutation_a"]), str(row["mutation_b"])}
    ]
    t99f_depth = Counter(str(row["expert_review_depth"]) for row in t99f_rows)
    if len(t99f_rows) != 14 or t99f_depth != {"enhanced": 2, "standard": 12}:
        raise ValueError("T99F doubles must follow the same generic 2/12 review-depth split")
    if len(selected) != 15:
        raise ValueError("Expected exactly 15 selected V3 doubles")
    if tuple(str(row["mutation_set"]) for row in selected) != SELECTED_DOUBLE_MUTATION_SETS:
        raise ValueError("Selected V3 double order does not match the approved explicit decision")
    if any(int(row["moderate_or_strong_favorable_metric_count"]) < 2 for row in selected):
        raise ValueError("Every selected double must have at least two positive property bands")
    if any(
        int(row["moderate_adverse_metric_count"]) + int(row["strong_adverse_metric_count"])
        for row in selected
    ):
        raise ValueError("Selected doubles cannot carry moderate or strong property adversity")
    parent_usage = Counter(
        mutation
        for row in selected
        for mutation in (str(row["mutation_a"]), str(row["mutation_b"]))
    )
    position_usage = Counter(
        int(position)
        for row in selected
        for position in (
            row["position_a_reported_1based"],
            row["position_b_reported_1based"],
        )
    )
    position_pairs = [
        tuple(
            sorted(
                (
                    int(row["position_a_reported_1based"]),
                    int(row["position_b_reported_1based"]),
                )
            )
        )
        for row in selected
    ]
    if max(parent_usage.values()) > 3 or max(position_usage.values()) > 4:
        raise ValueError("Selected-double diversity limits were exceeded")
    if len(set(position_pairs)) != 15:
        raise ValueError("Selected doubles must use unique unordered position pairs")
    if len(parent_usage) < 10 or len(position_usage) < 9:
        raise ValueError("Selected doubles do not meet parent or position coverage")
    af3_only = sum(
        str(row["pair_experimental_coordinate_status"]) != "both_observed"
        for row in selected
    )
    soft_risk = sum(bool(str(row["effective_soft_sequence_risk_flags"])) for row in selected)
    local = sum(
        str(row["pair_spatial_class"]) != SEPARATED_PAIR_CLASS for row in selected
    )
    if af3_only > 4 or soft_risk > 2 or local > 1:
        raise ValueError("Selected doubles exceed structural or soft-risk portfolio limits")
    if len(final_panel) != 30 or Counter(
        str(row["candidate_kind"]) for row in final_panel
    ) != {"single_mutant": 15, "double_mutant": 15}:
        raise ValueError("Expected a final 15-single plus 15-double panel")
    sequences = [str(row["sequence"]) for row in final_panel]
    if len(set(sequences)) != 30 or any(
        len(sequence) != 128 or not sequence.endswith("SSGS") or sequence.count("C") != 2
        for sequence in sequences
    ):
        raise ValueError("Final V3 panel sequences failed identity or construct checks")
    reconstructed_parents: set[str] = set()
    for row in audit:
        sequence = list(str(row["sequence"]))
        for mutation_key, position_key in (
            ("mutation_a", "position_a_reported_1based"),
            ("mutation_b", "position_b_reported_1based"),
        ):
            mutation = str(row[mutation_key])
            position = int(row[position_key])
            if sequence[position - 1] != mutation[-1]:
                raise ValueError(f"Declared mutation does not match {row['mutation_set']}")
            sequence[position - 1] = mutation[0]
        reconstructed_parents.add("".join(sequence))
    if len(reconstructed_parents) != 1:
        raise ValueError("V3 doubles do not reconstruct one common parent sequence")
    parent_sequence = reconstructed_parents.pop()
    for row in final_panel:
        difference_count = sum(
            wt != mutant
            for wt, mutant in zip(parent_sequence, str(row["sequence"]), strict=True)
        )
        expected = 1 if str(row["candidate_kind"]) == "single_mutant" else 2
        if difference_count != expected:
            raise ValueError("Final-panel substitution count is inconsistent")
    selected_t99f = sum(
        T99F_MUTATION in {str(row["mutation_a"]), str(row["mutation_b"])}
        for row in selected
    )
    selected_three = sum(
        int(row["moderate_or_strong_favorable_metric_count"]) == 3 for row in selected
    )
    return {
        "source_double_candidate_count": 102,
        "enhanced_expert_review_count": depth_counts["enhanced"],
        "standard_expert_review_count": depth_counts["standard"],
        "multi_metric_review_trigger_count": trigger_counts[
            "multi_metric_property_evidence"
        ],
        "adverse_review_trigger_count": trigger_counts[
            "moderate_or_strong_property_tradeoff"
        ],
        "local_geometry_review_trigger_count": trigger_counts["local_pair_geometry"],
        "standard_common_review_trigger_count": trigger_counts[
            "standard_common_review_no_enhanced_trigger"
        ],
        "t99f_double_count": len(t99f_rows),
        "t99f_enhanced_review_count": t99f_depth["enhanced"],
        "t99f_standard_review_count": t99f_depth["standard"],
        "selected_double_mutant_count": 15,
        "selected_three_metric_positive_count": selected_three,
        "selected_two_metric_positive_count": 15 - selected_three,
        "selected_parent_component_count": len(parent_usage),
        "selected_reported_position_count": len(position_usage),
        "maximum_exact_parent_usage": max(parent_usage.values()),
        "maximum_reported_position_usage": max(position_usage.values()),
        "selected_unique_position_pair_count": len(set(position_pairs)),
        "selected_af3_only_pair_count": af3_only,
        "selected_soft_risk_pair_count": soft_risk,
        "selected_local_pair_count": local,
        "selected_t99f_pair_count": selected_t99f,
        "t99f_specific_selection_rule_applied": False,
        "stable_word_specific_selection_rule_applied": False,
        "final_single_count": 15,
        "final_double_count": 15,
        "final_panel_candidate_count": 30,
    }


__all__ = [
    "SELECTED_DOUBLE_MUTATION_SETS",
    "build_v3_double_mutant_selection",
    "selection_policy",
    "selected_double_export_rows",
]
