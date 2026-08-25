"""Curated, non-ranking VHH expert judgements for the active V3 shortlist.

The entries in :data:`V3_EXPERT_ASSESSMENTS` interpret verified WT coordinate
contexts and amino-acid chemistry.  They are hypotheses for review, not
experimental observations, predictor outputs, stability calculations, or
candidate-selection decisions.  Coordinate-derived numeric facts remain the
responsibility of :mod:`antibody_optimization.vhh_expert_review`.
"""

from __future__ import annotations

from collections.abc import Iterable

from .v3_expert_review_pool import V3_REVIEW_POOL_COUNT


REQUIRED_EXPERT_FIELDS = frozenset(
    {
        "structural_facts_cn",
        "chimerax_single_rotamer_observation_cn",
        "vhh_expert_inference_cn",
        "expert_structural_assessment",
        "expert_solubility_expectation",
        "expert_thermal_stability_expectation",
        "expert_confidence",
        "expert_primary_concern",
        "expert_rationale_cn",
        "expert_uncertainty_cn",
        "expert_rule_flags",
    }
)

STRUCTURAL_CATEGORY = {
    "structurally_plausible": "reasonable",
    "structurally_plausible_core_conservative": "reasonable",
    "structurally_plausible_positive_phi_glycine": "reasonable",
    "structurally_plausible_terminal_context": "reasonable",
    "structurally_plausible_with_hallmark_neighborhood_caution": "reasonable_with_caution",
    "plausible_with_model_only_caution": "reasonable_with_caution",
    "plausible_with_loop_caution": "reasonable_with_caution",
    "plausible_with_tradeoff": "reasonable_with_caution",
    "plausible_with_glycine_caution": "reasonable_with_caution",
    "plausible_with_network_caution": "reasonable_with_caution",
    "plausible_with_cdr3_model_sensitivity": "reasonable_with_caution",
    "plausible_with_loop_and_liability_caution": "reasonable_with_caution",
    "plausible_with_electrostatic_caution": "reasonable_with_caution",
    "plausible_with_disulfide_neighborhood_caution": "reasonable_with_caution",
    "structure_sensitive": "structurally_concerning",
    "structure_sensitive_high_concern": "structurally_concerning",
    "high_structural_concern": "structurally_concerning",
    "very_high_structural_concern": "structurally_concerning",
    "model_only_structure_sensitive": "indeterminate",
}
SOLUBILITY_CATEGORY = {
    "likely_favorable": "favorable",
    "possibly_favorable": "favorable",
    "neutral_to_possibly_favorable": "neutral_or_uncertain",
    "possibly_favorable_but_weakly_supported": "neutral_or_uncertain",
    "neutral_or_uncertain": "neutral_or_uncertain",
    "neutral_to_possibly_adverse": "neutral_or_uncertain",
    "neutral_to_possibly_adverse_via_folding_risk": "unfavorable",
    "possibly_adverse": "unfavorable",
    "likely_adverse": "unfavorable",
    "likely_adverse_via_folding_risk": "unfavorable",
}
THERMAL_CATEGORY = {
    "possibly_favorable_or_neutral": "neutral_or_uncertain",
    "neutral_to_possibly_favorable": "neutral_or_uncertain",
    "neutral_or_uncertain": "neutral_or_uncertain",
    "neutral_to_possibly_adverse": "neutral_or_uncertain",
    "possibly_adverse": "unfavorable",
    "likely_adverse": "unfavorable",
    "likely_strongly_adverse": "unfavorable",
}


def _assessment(
    *,
    mutation: str,
    structural_facts_cn: str,
    vhh_expert_inference_cn: str,
    structural: str,
    solubility: str,
    thermal: str,
    confidence: str,
    concern: str,
    rationale_cn: str,
    uncertainty_cn: str,
    flags: tuple[str, ...],
) -> dict[str, object]:
    return {
        "mutation": mutation,
        "structural_facts_cn": structural_facts_cn,
        "vhh_expert_inference_cn": vhh_expert_inference_cn,
        "expert_structural_assessment": STRUCTURAL_CATEGORY[structural],
        "expert_solubility_expectation": SOLUBILITY_CATEGORY[solubility],
        "expert_thermal_stability_expectation": THERMAL_CATEGORY[thermal],
        "expert_confidence": confidence,
        "expert_primary_concern": concern,
        "expert_rationale_cn": rationale_cn,
        "expert_uncertainty_cn": uncertainty_cn,
        "expert_rule_flags": flags,
    }


V3_EXPERT_ASSESSMENTS: dict[str, dict[str, object]] = {
    "Nb252_expr_seq011_L11Y": _assessment(
        mutation="L11Y",
        structural_facts_cn="实验结构中reported 11缺失坐标；AF3中该位点位于FR1的外露β样主链位置。",
        vhh_expert_inference_cn="Tyr保留较大疏水芳香表面但增加羟基，可能改善相对Leu的表面水化；较大芳香侧链仍需检查局部rotamer空间。",
        structural="plausible_with_model_only_caution",
        solubility="possibly_favorable",
        thermal="neutral_or_uncertain",
        confidence="low",
        concern="experimental_coordinates_missing",
        rationale_cn="表面疏水残基转为带羟基芳香残基具有可解释的水化收益，但没有实验坐标支持其侧链环境。",
        uncertainty_cn="结论仅来自AF3单体；IMGT12具有VHH框架背景，但Y并不是可直接假定有利的典型共识回变。",
        flags=("af3_only", "beta_like", "exposed", "aromatic_introduced", "imgt12_vhh_context"),
    ),
    "Nb252_expr_seq030_F30S": _assessment(
        mutation="F30S",
        structural_facts_cn="实验结构中reported 30有坐标，位于外露CDR1 loop并靠近实验缺失片段边界；该位点不属于释放的严格界面集合。",
        vhh_expert_inference_cn="移除外露Phe的芳香疏水面通常有利于水化，但Ser缩小侧链并提高loop局部柔性，可能牺牲原有构象预组织。",
        structural="plausible_with_loop_caution",
        solubility="likely_favorable",
        thermal="neutral_to_possibly_adverse",
        confidence="medium",
        concern="cdr1_loop_repacking",
        rationale_cn="表面去疏水的溶解度机制清楚，热稳定性方向则受到loop包装和折叠熵的竞争影响。",
        uncertainty_cn="相邻未解析残基限制了对完整CDR1构象和受体邻近关系的判断；WT结构不能直接证明突变rotamer无碰撞。",
        flags=("experimental_coordinates", "cdr1", "loop", "surface_dehydrophobization", "gap_boundary"),
    ),
    "Nb252_expr_seq086_K86S": _assessment(
        mutation="K86S",
        structural_facts_cn="实验结构中reported 86有坐标，位于FR3外露β样位置并紧邻P87；WT Lys未显示为埋藏核心残基。",
        vhh_expert_inference_cn="Lys变Ser可降低长侧链构象熵并保留极性，但移除正电荷对溶解度的作用取决于整体pI和局部电荷斑。",
        structural="structurally_plausible",
        solubility="neutral_or_uncertain",
        thermal="neutral_to_possibly_favorable",
        confidence="medium",
        concern="surface_charge_removal",
        rationale_cn="外露β样位点可容纳较小极性替换，且没有证据表明必须保留Lys作为核心包装残基。",
        uncertainty_cn="仅凭静态WT距离不能排除弱盐桥或pH依赖电荷作用，也不能将Ser-Pro上下文解释为已验证的稳定化。",
        flags=("experimental_coordinates", "fr3", "beta_like", "exposed", "surface_charge_removed", "pre_proline"),
    ),
    "Nb252_expr_seq023_A23R": _assessment(
        mutation="A23R",
        structural_facts_cn="实验结构中reported 23有坐标且较外露，其主链紧邻形成天然二硫键的Cys22，并位于实验缺失片段的N端边界；因后续残基缺失，实验构象不能完整赋值。",
        vhh_expert_inference_cn="Arg可能提高表面水化，却将小Ala替换为长而带正电的侧链，需优先排查二硫键邻域碰撞和未补偿电荷。",
        structural="structure_sensitive",
        solubility="possibly_favorable",
        thermal="neutral_to_possibly_adverse",
        confidence="medium",
        concern="disulfide_adjacent_large_charge",
        rationale_cn="表面加电荷和二硫键邻域过度包装形成明确权衡，不能只按溶解度方向判断。",
        uncertainty_cn="风险取决于Arg可采用的朝溶剂rotamer；静态WT结构没有突变后侧链和局部松弛信息。",
        flags=("experimental_coordinates", "adjacent_cys22", "gap_boundary", "large_volume_increase", "positive_charge_introduced", "model_sensitive"),
    ),
    "Nb252_expr_seq005_Q5V": _assessment(
        mutation="Q5V",
        structural_facts_cn="实验结构中reported 5有坐标，位于FR1外露的β样位置并靠近Cys22空间邻域；该突变是冻结合同允许的天然共识回变。",
        vhh_expert_inference_cn="Val可能改善局部β片层疏水包装，但会移除Gln极性作用并增加暴露疏水性，因此预期是溶解度与局部稳定性的权衡。",
        structural="plausible_with_tradeoff",
        solubility="possibly_adverse",
        thermal="possibly_favorable_or_neutral",
        confidence="medium",
        concern="exposed_hydrophobic_consensus_reversion",
        rationale_cn="天然共识支持结构可容纳性，却不能证明在Nb252特定表面环境中提升溶解度或产量。",
        uncertainty_cn="实验与AF3暴露程度存在差异；天然频率、结构相容性和本构建体表达收益是不同证据。",
        flags=("experimental_coordinates", "fr1", "exposed", "hydrophobicity_increase", "consensus_reversion", "near_disulfide_space"),
    ),
    "Nb252_expr_seq049_A49F": _assessment(
        mutation="A49F",
        structural_facts_cn="实验结构中reported 49有坐标，位于FR2 β链的深埋位置，周围为密集VHH核心包装；它不是界面位点，但距冻结界面残基N58约3.6埃。",
        vhh_expert_inference_cn="将完全埋藏的小Ala替换为大芳香Phe极易产生过度包装和主链/侧链碰撞；除非存在未见空腔，否则不符合保守核心设计经验。",
        structural="high_structural_concern",
        solubility="likely_adverse_via_folding_risk",
        thermal="likely_adverse",
        confidence="high",
        concern="buried_large_aromatic_substitution",
        rationale_cn="深埋、体积大幅增加和邻居密集三项事实同时指向核心包装风险。",
        uncertainty_cn="尚未显式生成并局部松弛全部Phe rotamer；因此不能声称必然碰撞，但风险机制明确。",
        flags=("experimental_coordinates", "fr2", "beta_like", "buried", "large_volume_increase", "aromatic_introduced", "core_overpacking_risk", "near_interface_shell"),
    ),
    "Nb252_expr_seq055_S55G": _assessment(
        mutation="S55G",
        structural_facts_cn="实验结构中reported 55有坐标，位于外露CDR2 loop；WT主链不是需要Gly才能容纳的明确正phi构象。",
        vhh_expert_inference_cn="Gly去除Ser羟基并提高局部及未折叠态柔性，可能改变CDR2 loop预组织；没有强表面疏水变化可支持明确溶解度收益。",
        structural="plausible_with_glycine_caution",
        solubility="neutral_or_uncertain",
        thermal="possibly_adverse",
        confidence="medium",
        concern="cdr2_glycine_flexibility",
        rationale_cn="Gly适合某些紧转角，但此处缺乏正phi这一强几何依据，柔性代价应显式保留。",
        uncertainty_cn="loop动态不能由单一静态构象量化，Ser侧链水桥也未通过显式溶剂模拟评估。",
        flags=("experimental_coordinates", "cdr2", "loop", "glycine_introduced", "sidechain_hydroxyl_removed"),
    ),
    "Nb252_expr_seq075_K75A": _assessment(
        mutation="K75A",
        structural_facts_cn="实验结构中reported 75有坐标，位于FR3外露turn/loop；WT Lys侧链不属于深埋核心。",
        vhh_expert_inference_cn="Ala可降低长Lys侧链熵，但同时移除正电荷并增加小片表面非极性区域，溶解度方向不能仅凭去电荷判断。",
        structural="structurally_plausible",
        solubility="neutral_to_possibly_adverse",
        thermal="neutral_to_possibly_favorable",
        confidence="medium",
        concern="surface_charge_removal",
        rationale_cn="外露loop的Ala替换通常可容纳，但电荷和构象熵效应方向相反。",
        uncertainty_cn="局部电荷斑、构建体pI及缓冲液pH将影响实际溶解度，WT静态结构不足以确定。",
        flags=("experimental_coordinates", "fr3", "loop", "exposed", "surface_charge_removed", "alanine_substitution"),
    ),
    "Nb252_expr_seq001_Q1A": _assessment(
        mutation="Q1A",
        structural_facts_cn="实验结构中reported 1有坐标并高度暴露，位于Nb252报告序列N端。",
        vhh_expert_inference_cn="Ala替换结构上容易容纳，但移除Gln极性侧链；其表达影响可能更多来自真实成熟N端加工而非VHH核心稳定性。",
        structural="structurally_plausible_terminal_context",
        solubility="neutral_to_possibly_adverse",
        thermal="neutral_or_uncertain",
        confidence="medium",
        concern="construct_n_terminus_context",
        rationale_cn="高度外露的末端不构成核心过度包装风险，主要问题是表达构建体中的实际N端身份。",
        uncertainty_cn="若存在信号肽、标签或切割残基，reported Q1可能并非细胞内真实N端；不能据此推断N端降解或焦谷氨酸效应。",
        flags=("experimental_coordinates", "n_terminus", "exposed", "construct_context_uncertain", "polar_sidechain_removed"),
    ),
    "Nb252_expr_seq028_I28Y": _assessment(
        mutation="I28Y",
        structural_facts_cn="实验结构中reported 28缺失坐标；AF3中该位点位于CDR1部分埋藏loop并邻近F29、F30和Y32芳香环境。",
        vhh_expert_inference_cn="Tyr体积大于Ile且扩展局部芳香簇，可能改善特定包装，也可能造成过度包装或形成暴露疏水斑。",
        structural="model_only_structure_sensitive",
        solubility="neutral_to_possibly_adverse",
        thermal="neutral_or_uncertain",
        confidence="low",
        concern="af3_only_aromatic_cluster_overpacking",
        rationale_cn="AF3给出部分埋藏和邻近芳香簇的风险线索，但缺乏实验loop与受体环境。",
        uncertainty_cn="实验接触状态为不可评价；AF3单体CDR1构象和Tyr rotamer均不能作为实验事实。",
        flags=("af3_only", "cdr1", "partially_buried", "large_volume_increase", "aromatic_cluster", "experimental_contact_not_evaluable"),
    ),
    "Nb252_expr_seq029_F29Q": _assessment(
        mutation="F29Q",
        structural_facts_cn="实验结构中reported 29缺失坐标；AF3中该位点位于外露CDR1 loop并邻近局部芳香残基。",
        vhh_expert_inference_cn="Gln可移除外露Phe疏水芳香面并提高水化，但也可能破坏CDR1芳香预组织或改变loop构象。",
        structural="plausible_with_model_only_caution",
        solubility="likely_favorable",
        thermal="neutral_to_possibly_adverse",
        confidence="low",
        concern="experimental_coordinates_missing",
        rationale_cn="表面去芳香疏水具有直接溶解度机制，但结构稳定性取决于实验未解析的CDR1网络。",
        uncertainty_cn="实验受体接触和完整CDR1构象不可评价；结论仅来自AF3单体。",
        flags=("af3_only", "cdr1", "exposed", "surface_dehydrophobization", "aromatic_removed", "experimental_contact_not_evaluable"),
    ),
    "Nb252_expr_seq032_Y32L": _assessment(
        mutation="Y32L",
        structural_facts_cn="实验结构中reported 32有坐标，位于CDR1外露β样位置并处在跨区域紧密邻居环境；它不是界面位点，但距冻结界面残基D98约3.4埃，且AF3给出更低的暴露度。",
        vhh_expert_inference_cn="Leu移除Tyr羟基和芳香几何，同时仍保持较强疏水性，既不构成明确表面水化改良，也可能削弱局部包装。",
        structural="structure_sensitive",
        solubility="possibly_adverse",
        thermal="possibly_adverse",
        confidence="medium",
        concern="aromatic_and_hydrogen_bond_network_loss",
        rationale_cn="从带羟基芳香残基变为脂肪族疏水残基，缺乏清晰的溶解度收益并存在局部结构网络损失。",
        uncertainty_cn="需要突变rotamer和局部松弛才能区分芳香包装损失与可能的去张力收益。",
        flags=("experimental_coordinates", "cdr1", "exposed", "aromatic_removed", "hydroxyl_removed", "local_network_risk", "near_interface_shell", "model_sensitive"),
    ),
    "Nb252_expr_seq040_A40G": _assessment(
        mutation="A40G",
        structural_facts_cn="实验结构中reported 40有坐标，位于FR2部分埋藏的loop/β边界并紧邻P41。",
        vhh_expert_inference_cn="Gly-Pro上下文可能有利于局部转角，但Gly也去除Ala甲基包装并提高未折叠态熵，稳定性方向具有竞争性。",
        structural="plausible_with_glycine_caution",
        solubility="neutral_or_uncertain",
        thermal="neutral_to_possibly_adverse",
        confidence="medium",
        concern="glycine_at_partly_buried_turn_boundary",
        rationale_cn="该位置有转角合理性，却不是可仅凭Gly偏好就宣称稳定化的正phi位点。",
        uncertainty_cn="静态主链不能量化Gly导致的构象集合扩展；邻近Pro的顺反异构影响未评估。",
        flags=("experimental_coordinates", "fr2", "partially_buried", "turn_boundary", "glycine_introduced", "pre_proline"),
    ),
    "Nb252_expr_seq043_K43A": _assessment(
        mutation="K43A",
        structural_facts_cn="实验结构中reported 43有坐标，位于外露FR2 loop及经典VHH hallmark邻域，但并非IMGT 42、49、50或52本身。",
        vhh_expert_inference_cn="移除Lys可能减弱正电荷斑或局部电荷拥挤，也可能降低总体水化并改变邻近hallmark表面的电荷平衡。",
        structural="structurally_plausible_with_hallmark_neighborhood_caution",
        solubility="neutral_or_uncertain",
        thermal="neutral_or_uncertain",
        confidence="medium",
        concern="hallmark_neighborhood_charge_balance",
        rationale_cn="外露β样位置可容纳Ala，但VHH former-VL表面的电荷效应不能用单一去电荷规则预测。",
        uncertainty_cn="未确认Lys侧链的pH依赖水桥和局部静电能；hallmark邻域不等同于hallmark位点。",
        flags=("experimental_coordinates", "fr2", "exposed", "surface_charge_removed", "hallmark_neighborhood", "not_direct_hallmark"),
    ),
    "Nb252_expr_seq050_S50F": _assessment(
        mutation="S50F",
        structural_facts_cn="实验结构中reported 50有坐标，位于FR2 β链的部分埋藏位置、邻居密集，且处于严格4埃界面外的近受体壳层。",
        vhh_expert_inference_cn="小极性Ser变为大芳香Phe会增加局部体积和疏水表面，可能同时造成核心过度包装、表面聚集倾向和近界面构象扰动。",
        structural="high_structural_concern",
        solubility="likely_adverse",
        thermal="likely_adverse",
        confidence="high",
        concern="partly_buried_aromatic_overpacking_near_receptor",
        rationale_cn="部分埋藏、显著体积增加、芳香疏水引入和近受体位置构成相互一致的风险证据。",
        uncertainty_cn="该位点不是严格界面残基，不能据此宣称一定改变结合；显式rotamer检查仍用于确认碰撞程度。",
        flags=("experimental_coordinates", "fr2", "beta_like", "partially_buried", "large_volume_increase", "aromatic_introduced", "near_interface_shell"),
    ),
    "Nb252_expr_seq060_A60D": _assessment(
        mutation="A60D",
        structural_facts_cn="实验结构中reported 60有坐标，位于FR3部分埋藏β样位置并邻接D61；突变后形成局部DD序列，且该位点距冻结界面残基F47约3.8埃。",
        vhh_expert_inference_cn="Asp可能提高暴露部分的水化，却可能在部分埋藏环境形成未补偿负电荷并与D61产生局部排斥。",
        structural="structure_sensitive",
        solubility="neutral_or_uncertain",
        thermal="possibly_adverse",
        confidence="medium",
        concern="partly_buried_negative_charge_cluster",
        rationale_cn="表面加负电与部分埋藏电荷代价相互竞争，局部DD使稳定性风险更加具体。",
        uncertainty_cn="Asp质子化状态、溶剂可达性和突变后水网络需要显式结构或实验验证。",
        flags=("experimental_coordinates", "fr3", "partially_buried", "negative_charge_introduced", "adjacent_acidic_residue", "buried_charge_risk", "near_interface_shell"),
    ),
    "Nb252_expr_seq069_I69V": _assessment(
        mutation="I69V",
        structural_facts_cn="实验结构中reported 69有坐标，位于FR3 β链深埋核心；Ile到Val为同类疏水残基的保守缩小。",
        vhh_expert_inference_cn="替换通常具有较高结构可容纳性，但较小Val可能留下轻微空腔；没有直接表面机制支持明显溶解度改变。",
        structural="structurally_plausible_core_conservative",
        solubility="neutral_or_uncertain",
        thermal="neutral_to_possibly_adverse",
        confidence="high",
        concern="small_core_cavity",
        rationale_cn="化学类别保守且体积变化有限，主要剩余风险是深埋核心的欠包装。",
        uncertainty_cn="是否存在恰好适配Val的局部压缩需rotamer/局部松弛或实验Tm确认。",
        flags=("experimental_coordinates", "fr3", "beta_like", "buried", "conservative_hydrophobic", "core_cavity_risk"),
    ),
    "Nb252_expr_seq071_R71G": _assessment(
        mutation="R71G",
        structural_facts_cn="实验结构中reported 71有坐标，主链呈扩展β样构象且位于FR3部分埋藏环境；WT Arg邻近多处VHH内部极性接触，并距冻结界面残基D33约4.0002埃。",
        vhh_expert_inference_cn="Gly会同时移除长Arg侧链的极性网络并增加扩展主链处的构象自由度，可能削弱β夹层预组织。",
        structural="high_structural_concern",
        solubility="neutral_to_possibly_adverse_via_folding_risk",
        thermal="likely_adverse",
        confidence="high",
        concern="internal_polar_network_and_backbone_flexibility_loss",
        rationale_cn="极性网络丢失和非转角Gly引入是两个独立且同向的稳定性风险。",
        uncertainty_cn="部分接触涉及主链和静态几何，实际能量贡献尚未由突变结构或实验验证。",
        flags=("experimental_coordinates", "fr3", "extended_backbone", "partially_buried", "glycine_introduced", "internal_polar_network_loss", "borderline_interface_shell"),
    ),
    "Nb252_expr_seq076_N76G": _assessment(
        mutation="N76G",
        structural_facts_cn="实验结构中reported 76有坐标，位于FR3紧转角并采用正phi主链构象；实验与AF3对暴露程度有差异。",
        vhh_expert_inference_cn="Gly对正phi构象具有明确几何相容性，可能缓解局部主链应变；代价是失去Asn侧链极性作用并增加未折叠态熵。",
        structural="structurally_plausible_positive_phi_glycine",
        solubility="neutral_or_uncertain",
        thermal="possibly_favorable_or_neutral",
        confidence="medium",
        concern="loss_of_asparagine_sidechain_network",
        rationale_cn="正phi是支持Gly替换的具体结构证据，但不能抵消全部侧链和熵代价。",
        uncertainty_cn="实验与AF3暴露分类不一致，且静态几何不能量化净自由能变化。",
        flags=("experimental_coordinates", "fr3", "positive_phi_turn", "glycine_introduced", "model_sensitive", "asparagine_removed"),
    ),
    "Nb252_expr_seq079_Y79T": _assessment(
        mutation="Y79T",
        structural_facts_cn="实验结构中reported 79有坐标并部分埋藏；实验与AF3均支持β样主链，WT Tyr参与密集VHH内部包装并邻近K75。",
        vhh_expert_inference_cn="Thr可降低芳香疏水面并保留羟基，但显著缩小侧链可能在部分埋藏环境形成空腔并丢失芳香包装。",
        structural="structure_sensitive_high_concern",
        solubility="possibly_favorable_but_weakly_supported",
        thermal="likely_adverse",
        confidence="medium",
        concern="partly_buried_aromatic_cavity",
        rationale_cn="表面去疏水收益受到部分埋藏事实限制，而芳香核心欠包装是更直接的结构风险。",
        uncertainty_cn="二级结构模型敏感，需突变rotamer和局部松弛区分空腔风险与可能的去张力收益。",
        flags=("experimental_coordinates", "fr3", "partially_buried", "aromatic_removed", "core_cavity_risk"),
    ),
    "Nb252_expr_seq083_N83A": _assessment(
        mutation="N83A",
        structural_facts_cn="实验结构中reported 83有坐标，位于FR3 β边缘的外露位置；实验与AF3暴露度和β样主链描述一致，WT Asn邻近跨片层极性残基。",
        vhh_expert_inference_cn="Ala可降低侧链熵，但会移除Asn的氢键能力并略增表面非极性，净稳定性取决于该极性网络是否真实占优。",
        structural="plausible_with_network_caution",
        solubility="neutral_to_possibly_adverse",
        thermal="neutral_to_possibly_adverse",
        confidence="medium",
        concern="polar_network_loss_at_sheet_edge",
        rationale_cn="侧链缩小容易容纳，但跨片层极性作用和表面水化可能同时受损。",
        uncertainty_cn="实验与AF3局部二级结构不一致，Asn水桥的占有率未知。",
        flags=("experimental_coordinates", "fr3", "sheet_edge", "exposed", "polar_sidechain_removed"),
    ),
    "Nb252_expr_seq096_A96R": _assessment(
        mutation="A96R",
        structural_facts_cn="实验结构中reported 96有坐标，位于CDR3起始β链的深埋核心并紧邻形成天然二硫键的Cys95，周围包装密集；它不是界面位点，但距冻结界面残基Y115约3.2埃。",
        vhh_expert_inference_cn="在此处将小Ala替换为大而带正电的Arg，极易产生过度包装和未补偿埋藏电荷，并可能扰动二硫键邻域。",
        structural="very_high_structural_concern",
        solubility="likely_adverse_via_folding_risk",
        thermal="likely_strongly_adverse",
        confidence="high",
        concern="buried_charged_substitution_adjacent_disulfide",
        rationale_cn="深埋、巨大体积增加、引入电荷和紧邻Cys95四项结构事实给出一致的高风险机制。",
        uncertainty_cn="尚无突变后局部松弛结构，因此不声称必然失折叠；但静态结构风险置信度高。",
        flags=("experimental_coordinates", "cdr3", "beta_like", "buried", "adjacent_cys95", "large_volume_increase", "buried_charge_risk", "core_overpacking_risk", "near_interface_shell"),
    ),
    "Nb252_expr_seq099_T99N": _assessment(
        mutation="T99N",
        structural_facts_cn="实验结构中reported 99有坐标并位于CDR3 loop/turn；AF3将该位置置于β样构象，WT Thr参与局部极性接触，且距冻结界面残基I112约2.4埃。",
        vhh_expert_inference_cn="Thr到Asn保持极性但改变侧链几何和氢键供受体，通常比跨化学类别替换温和，却可能重排CDR3内部网络。",
        structural="plausible_with_cdr3_model_sensitivity",
        solubility="neutral_or_uncertain",
        thermal="neutral_or_uncertain",
        confidence="medium",
        concern="cdr3_hydrogen_bond_repatterning",
        rationale_cn="化学类别相近支持可容纳性，实验与AF3构象差异则限制稳定性方向判断。",
        uncertainty_cn="CDR3是已知模型差异最大的区域；单一结构不能给出突变后loop构象分布。",
        flags=("experimental_coordinates", "cdr3", "polar_to_polar", "hydrogen_bond_repatterning", "model_sensitive", "near_interface_shell"),
    ),
    "Nb252_expr_seq030_F30N": _assessment(
        mutation="F30N",
        structural_facts_cn="实验结构中reported 30有坐标，位于外露CDR1 loop及实验缺失边界；完整突变序列新增NG脱酰胺敏感上下文。",
        vhh_expert_inference_cn="Asn可移除Phe芳香疏水面并提供氢键，可能改善水化；同时侧链缩小和NG化学不稳定性限制其总体可开发性。",
        structural="plausible_with_loop_and_liability_caution",
        solubility="likely_favorable",
        thermal="neutral_or_uncertain",
        confidence="medium",
        concern="new_asparagine_glycine_deamidation_motif",
        rationale_cn="表面去疏水机制明确，但不能忽略新NG序列及CDR1预组织变化。",
        uncertainty_cn="脱酰胺影响取决于制备、pH和储存条件；相邻实验缺失区域限制完整loop判断。",
        flags=("experimental_coordinates", "cdr1", "exposed", "surface_dehydrophobization", "asparagine_introduced", "new_deamidation_motif", "gap_boundary"),
    ),
    "Nb252_expr_seq086_K86A": _assessment(
        mutation="K86A",
        structural_facts_cn="实验结构中reported 86有坐标，位于FR3外露β样位置并紧邻P87；WT Lys不属于深埋核心。",
        vhh_expert_inference_cn="Ala可降低长Lys侧链熵，但相比Ser不保留极性且移除正电荷，可能略降低表面水化。",
        structural="structurally_plausible",
        solubility="neutral_to_possibly_adverse",
        thermal="neutral_to_possibly_favorable",
        confidence="medium",
        concern="surface_charge_and_polarity_removal",
        rationale_cn="结构容纳性较好，但溶解度和热稳定性机制方向相反。",
        uncertainty_cn="整体pI、局部电荷斑及Ala对turn熵的净效应未由实验确定。",
        flags=("experimental_coordinates", "fr3", "beta_like", "exposed", "surface_charge_removed", "polar_surface_reduced", "pre_proline"),
    ),
    "Nb252_expr_seq075_K75E": _assessment(
        mutation="K75E",
        structural_facts_cn="实验结构中reported 75有坐标，位于FR3外露turn/loop；替换导致正电到负电的电荷反转，局部邻近D72等酸性环境。",
        vhh_expert_inference_cn="Glu可降低正电荷斑，却可能形成新的负电荷簇或排斥；其溶解度与稳定性均强烈依赖pH和局部溶剂屏蔽。",
        structural="plausible_with_electrostatic_caution",
        solubility="neutral_or_uncertain",
        thermal="neutral_to_possibly_adverse",
        confidence="medium",
        concern="surface_charge_reversal",
        rationale_cn="外露位置可容纳Glu，但两单位形式电荷变化不能作为普适溶解度改良规则。",
        uncertainty_cn="未计算pKa、整体pI或盐浓度依赖静电；静态距离不能量化电荷反转自由能。",
        flags=("experimental_coordinates", "fr3", "loop", "exposed", "charge_reversal", "local_acidic_environment"),
    ),
    "Nb252_expr_seq011_L11M": _assessment(
        mutation="L11M",
        structural_facts_cn="实验结构中reported 11缺失坐标；AF3中该位点位于FR1外露β样位置。",
        vhh_expert_inference_cn="Met与Leu体积和疏水类别接近，结构替换较保守，但表面Met不提供明确水化收益并新增氧化敏感性。",
        structural="plausible_with_model_only_caution",
        solubility="neutral_to_possibly_adverse",
        thermal="neutral_or_uncertain",
        confidence="low",
        concern="experimental_coordinates_missing_and_methionine_oxidation",
        rationale_cn="保守疏水替换降低了几何风险，但AF3-only证据和新增Met限制了可开发性解释。",
        uncertainty_cn="没有实验结构侧链环境；Met氧化程度取决于表达、纯化和储存条件。",
        flags=("af3_only", "fr1", "beta_like", "exposed", "conservative_hydrophobic", "methionine_introduced", "oxidation_liability", "imgt12_vhh_context"),
    ),
    "Nb252_expr_seq049_A49M": _assessment(
        mutation="A49M",
        structural_facts_cn="实验结构中reported 49有坐标，位于FR2 β链深埋且邻居密集的VHH核心；它不是界面位点，但距冻结界面残基N58约3.6埃。",
        vhh_expert_inference_cn="Ala变Met显著增大柔性疏水侧链，可能造成核心过度包装；Met还新增氧化敏感性。",
        structural="high_structural_concern",
        solubility="likely_adverse_via_folding_risk",
        thermal="likely_adverse",
        confidence="high",
        concern="buried_large_methionine_substitution",
        rationale_cn="深埋位置的大体积增加是主要风险，氧化是额外而非替代性的可开发性问题。",
        uncertainty_cn="若WT存在真实空腔，Met可能填充该空腔；目前静态结构和未松弛rotamer没有支持这一例外。",
        flags=("experimental_coordinates", "fr2", "beta_like", "buried", "large_volume_increase", "core_overpacking_risk", "methionine_introduced", "oxidation_liability", "near_interface_shell"),
    ),
    "Nb252_expr_seq001_Q1D": _assessment(
        mutation="Q1D",
        structural_facts_cn="实验结构中reported 1有坐标并高度暴露，位于Nb252报告序列N端；替换新增负电荷。",
        vhh_expert_inference_cn="暴露Asp结构上易容纳且可能提高水化或调整pI，但其表达影响取决于真实成熟N端和构建体加工。",
        structural="structurally_plausible_terminal_context",
        solubility="possibly_favorable",
        thermal="neutral_or_uncertain",
        confidence="medium",
        concern="construct_n_terminus_and_charge_context",
        rationale_cn="不存在核心包装障碍，主要潜在收益来自表面负电和水化。",
        uncertainty_cn="若存在信号肽、标签或额外残基，reported D1并非实际N端；不能在未知构建体背景下推断N端降解规则。",
        flags=("experimental_coordinates", "n_terminus", "exposed", "negative_charge_introduced", "construct_context_uncertain"),
    ),
    "Nb252_expr_seq023_A23Q": _assessment(
        mutation="A23Q",
        structural_facts_cn="实验结构中reported 23有坐标且较外露，其主链紧邻Cys22二硫键，并位于实验缺失片段的N端边界；因后续残基缺失，实验构象不能完整赋值。",
        vhh_expert_inference_cn="Gln增加极性和侧链体积，可能改善水化；相比Arg不引入形式电荷，但仍需检查二硫键邻域rotamer和局部包装。",
        structural="plausible_with_disulfide_neighborhood_caution",
        solubility="neutral_to_possibly_favorable",
        thermal="neutral_or_uncertain",
        confidence="medium",
        concern="disulfide_adjacent_volume_increase",
        rationale_cn="外露极性替换具备可解释的水化机制，二硫键相邻位置要求比普通表面Ala更谨慎。",
        uncertainty_cn="静态WT结构不包含Gln rotamer和局部松弛；局部二级结构存在模型敏感性。",
        flags=("experimental_coordinates", "adjacent_cys22", "gap_boundary", "polar_sidechain_introduced", "volume_increase", "model_sensitive"),
    ),
    "Nb252_expr_seq099_T99F": _assessment(
        mutation="T99F",
        structural_facts_cn="实验结构和AF3均在reported 99具有坐标；该CDR3位点部分埋藏且不是直接界面残基，但距冻结界面残基I112约2.4埃。实验与AF3的局部主链及Thr极性接触伙伴不同。",
        vhh_expert_inference_cn="Thr到Phe使侧链体积增加约42立方埃、疏水性显著升高并移除羟基。芳环或可参与局部芳香包装，但在邻居密集、构象敏感的CDR3中更直接的风险是过度包装、极性网络丢失和外露疏水斑增加。",
        structural="structure_sensitive_high_concern",
        solubility="likely_adverse",
        thermal="possibly_adverse",
        confidence="medium",
        concern="partly_buried_aromatic_overpacking_and_polar_network_loss_in_model_sensitive_cdr3",
        rationale_cn="实验复合物与AF3均显示该位点只有有限空间，且WT Thr羟基具有距离相容的局部极性接触；稳定词新增仅是软偏好，不能抵消体积、疏水性和局部网络三项同向风险。",
        uncertainty_cn="只检查了未松弛固定骨架上的单个Dunbrack rotamer；实验与AF3的CDR3构象不一致，因此不声称该突变必然失稳，也不据此执行选择或淘汰。",
        flags=("experimental_coordinates", "cdr3", "partially_buried", "large_volume_increase", "aromatic_introduced", "hydrophobicity_increase", "polar_network_loss_risk", "near_interface_shell", "model_sensitive", "stable_word_gain"),
    ),
}


# These observations were made from the candidate-specific ChimeraX 1.12
# ``swapaa`` views released with this review.  They describe only the selected
# Dunbrack rotamer on the unchanged source backbone; they are not clash scores,
# relaxed mutant structures, or selection decisions.
V3_CHIMERAX_VISUAL_OBSERVATIONS: dict[str, str] = {
    "Nb252_expr_seq011_L11Y": "AF3-only视图中Tyr侧链朝向溶剂，所选rotamer未见肉眼可辨的严重原子穿插；实验局部环境仍不可评价。",
    "Nb252_expr_seq030_F30S": "实验视图中Ser为朝外的小侧链，所选rotamer无明显过度包装；该静态视图不能评估移除Phe后CDR1预组织的变化。",
    "Nb252_expr_seq086_K86S": "实验视图中Ser位于表面开放空间，所选rotamer未见明显空间冲突。",
    "Nb252_expr_seq023_A23R": "实验视图中Arg可向外伸展，但起点紧邻Cys22且位于缺失片段边界；AF3敏感性视图显示局部骨架背景不同。",
    "Nb252_expr_seq005_Q5V": "实验视图中Val体积紧凑，所选rotamer未见明显穿插；其较疏水侧链仍处于部分可接触溶剂的环境。",
    "Nb252_expr_seq049_A49F": "实验视图中Phe芳环被置入邻居密集的核心，周围可用空间有限，呈现明显的过度包装警示。",
    "Nb252_expr_seq055_S55G": "实验视图中突变后仅保留Gly主链，不产生侧链冲突；柔性与折叠态熵代价无法从固定主链图片直接观察。",
    "Nb252_expr_seq075_K75A": "实验视图中Ala位于表面开放空间，未见明显包装问题；主要改变是移除Lys长侧链和正电荷。",
    "Nb252_expr_seq001_Q1A": "实验视图中Ala位于开放的报告序列N端，未见局部空间障碍。",
    "Nb252_expr_seq028_I28Y": "AF3-only视图中Tyr芳环进入已有芳香/疏水邻居附近，局部空间偏紧；没有实验坐标可确认该构象。",
    "Nb252_expr_seq029_F29Q": "AF3-only视图中Gln侧链可朝外伸展，所选rotamer未见明显穿插；缺失实验坐标限制判断。",
    "Nb252_expr_seq032_Y32L": "实验视图中Leu可被局部空间容纳，但会移除原Tyr芳环和羟基；该位置靠近冻结界面残基网络。",
    "Nb252_expr_seq040_A40G": "实验视图中Gly不产生侧链冲突，但Ala甲基被移除后留下的局部空腔及主链柔性无法由固定骨架修复。",
    "Nb252_expr_seq043_K43A": "实验视图中Ala位于表面开放空间，未见明显包装冲突；主要代价是移除表面正电与极性。",
    "Nb252_expr_seq050_S50F": "实验视图中Phe芳环进入部分埋藏且邻近受体的拥挤区域，呈现过度包装与表面疏水增加警示。",
    "Nb252_expr_seq060_A60D": "实验视图中Asp可采用一个不明显穿插的rotamer，但该位置部分埋藏，固定骨架视图不能消除埋藏电荷风险。",
    "Nb252_expr_seq069_I69V": "实验视图中Val可直接容纳且无穿插，但相对Ile缩小后的核心空腔不会在固定骨架视图中闭合。",
    "Nb252_expr_seq071_R71G": "实验视图中Gly本身无侧链冲突，但原Arg占据的局部接触空间被完全清空，支持局部网络丢失警示。",
    "Nb252_expr_seq076_N76G": "实验视图中Gly与该转角主链几何相容且不产生侧链冲突；原Asn侧链接触被移除。",
    "Nb252_expr_seq079_Y79T": "实验视图中Thr可被容纳，但相对Tyr显著缩小并移除芳环，留下局部包装与接触损失警示。",
    "Nb252_expr_seq083_N83A": "实验视图中Ala体积紧凑且无明显穿插；原Asn的极性接触能力不再存在。",
    "Nb252_expr_seq096_A96R": "实验视图中Arg被置入深埋、邻居密集且靠近Cys95的区域，长带电侧链显示强烈过度包装和埋藏电荷警示；AF3视图不消除该风险。",
    "Nb252_expr_seq099_T99N": "实验视图中Asn可置入局部网络而无肉眼可辨的严重穿插，但实验与AF3敏感性视图的CDR3骨架背景明显不同。",
    "Nb252_expr_seq030_F30N": "实验视图中Asn朝向外侧且无明显过度包装；AF3敏感性视图提示缺失片段边界的局部构象仍不确定。",
    "Nb252_expr_seq086_K86A": "实验视图中Ala位于表面开放空间且无明显穿插，主要改变为移除Lys电荷和极性。",
    "Nb252_expr_seq075_K75E": "实验视图中Glu可向溶剂伸展且无明显空间冲突；图片不能判定电荷反转在实际pH和盐条件下的净效应。",
    "Nb252_expr_seq011_L11M": "AF3-only视图中Met朝向溶剂且几何上可容纳；实验环境缺失且新增Met氧化风险不由图片体现。",
    "Nb252_expr_seq049_A49M": "实验视图中Met被置入深埋核心并占据有限空间，呈现过度包装警示；固定骨架视图未显示支持空腔填充的明确余量。",
    "Nb252_expr_seq001_Q1D": "实验视图中Asp位于开放N端且无明显局部空间障碍；真实构建体N端加工仍不由结构图决定。",
    "Nb252_expr_seq023_A23Q": "实验视图中Gln可向外伸展且未见明显穿插，但仍邻近Cys22和缺失片段边界；AF3敏感性视图显示局部模型依赖。",
    "Nb252_expr_seq099_T99F": "实验复合物视图中Phe芳环进入邻居密集的CDR3局部空间，呈现过度包装警示；AF3敏感性视图中的骨架与朝向不同，但同样未提供开放空间或保留WT极性网络的证据。",
}

if set(V3_CHIMERAX_VISUAL_OBSERVATIONS) != set(V3_EXPERT_ASSESSMENTS):
    raise ValueError("ChimeraX visual-observation identities must match V3 assessments")
for _candidate_id, _observation in V3_CHIMERAX_VISUAL_OBSERVATIONS.items():
    V3_EXPERT_ASSESSMENTS[_candidate_id][
        "chimerax_single_rotamer_observation_cn"
    ] = _observation


def validate_v3_expert_assessments(expected_candidate_ids: Iterable[str] | None = None) -> None:
    """Validate identity and field completeness without performing selection."""

    if len(V3_EXPERT_ASSESSMENTS) != V3_REVIEW_POOL_COUNT:
        raise ValueError(
            f"Expected {V3_REVIEW_POOL_COUNT} expert assessments, "
            f"found {len(V3_EXPERT_ASSESSMENTS)}"
        )
    if expected_candidate_ids is not None:
        expected = {str(value) for value in expected_candidate_ids}
        observed = set(V3_EXPERT_ASSESSMENTS)
        if expected != observed:
            raise ValueError(
                "Expert-assessment candidate mismatch: "
                f"missing={sorted(expected - observed)}, extra={sorted(observed - expected)}"
            )
    for candidate_id, row in V3_EXPERT_ASSESSMENTS.items():
        missing = REQUIRED_EXPERT_FIELDS - set(row)
        if missing:
            raise ValueError(f"{candidate_id} lacks expert fields: {sorted(missing)}")
        if any(not str(row[field]).strip() for field in REQUIRED_EXPERT_FIELDS):
            raise ValueError(f"{candidate_id} contains an empty expert field")
        flags = row["expert_rule_flags"]
        if not isinstance(flags, tuple) or not flags:
            raise ValueError(f"{candidate_id} expert_rule_flags must be a non-empty tuple")


def get_v3_expert_assessment(candidate_id: str) -> dict[str, object]:
    """Return a shallow copy of one curated assessment."""

    try:
        return dict(V3_EXPERT_ASSESSMENTS[candidate_id])
    except KeyError as exc:
        raise KeyError(f"Unknown V3 candidate_id: {candidate_id}") from exc


def get_all_v3_expert_assessments() -> dict[str, dict[str, object]]:
    """Return copies of all assessments; insertion order is not a rank."""

    return {candidate_id: dict(row) for candidate_id, row in V3_EXPERT_ASSESSMENTS.items()}


validate_v3_expert_assessments()


__all__ = [
    "REQUIRED_EXPERT_FIELDS",
    "V3_EXPERT_ASSESSMENTS",
    "get_all_v3_expert_assessments",
    "get_v3_expert_assessment",
    "validate_v3_expert_assessments",
]
