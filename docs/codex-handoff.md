# Codex 项目交接

Last updated: 2026-08-18 23:30:00

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是当前状态快照，不是追加式历史。现行路线以 Nb252 在 BL21 体系中的表达量优化为唯一设计目标。

## 当前目标与硬约束

- 最终向合作者提供 30 条 Nb252 单突序列；单突实验完成前不设计或交付组合突变。
- 冻结实验复合物中已复现的 24 个 VHH 界面位点，不再显式优化亲和力，也不在当前路线使用 Rosetta 排序。
- 冻结天然 VHH 邻域中的高置信保守位点、Cys22/Cys95 和末端 SSGS（reported 125–128）。
- 候选必须保持完整 128-aa 父序列长度、末端 SSGS、两枚原有 Cys，且不得引入新 Cys。
- 现行核心预测工具仅为 NetSolP、NanoMelt 和实验复合物视图 AntiFold。新工具必须先在 47 条可比产量数据上验证，证明有独立且可重复的样本外信息后才能纳入筛选。
- RP3Net 0.0.2 已完成环境和权重身份固定，47条验证计划及连续/离散统计实现已就绪；远程正式评分尚未运行，因此当前不能用于候选排序。

## 权威输入基线

- Nb252 权威设计父本为 collaborator-confirmed 128-aa reported sequence；末端 SSGS 是构建体组成而非 linker。
- 实验结构 VHH 为链 C，NK2R 为链 R；24 个界面残基采用项目已复现的严格 `<4.0 Å` 聚合物重原子距离定义。
- AF3 VHH 是预测结构，只补充结构视图；不得把预测 CDR3 当作实验界面证据。
- 47 条产量记录来自 BL21 表达体系，LTT/WCC 个体近似值可作为数值，LLJ 仍保持分组/删失语义。

## Natural-VHH conservation contract

- 数据源固定为 TNP 论文仓库提交 `a9ba3edc3d967ecf8a2b9b5c2c29bf7495bbc9a0` 的最终 VHH-OAS 描述符表，共 4,059 条非冗余天然 VHH；输入文件 SHA-256 为 `D87A2E66CE0E46D34547D25DF10BF07ABC10B06CCF0D9D0C4A304A36A9D0EBE5`。
- 采用项目已固定的 ANARCII 2.0.8 / IMGT 编号。4,057 条通过 H 链与 framework coverage 审核，2 条因重复编号失败而排除。
- 先按完整 IMGT 域 90% identity 单连接聚类，每簇总权重为 1；再以 framework-only IMGT identity `>=0.80` 且 coverage `>=0.80` 定义 Nb252 邻域。
- 得到 3,784 个去冗余簇；Nb252 邻域含 1,564 条序列、1,532 个有效簇。
- reported 128 个位置分为 55 个 `hard_conserved`、33 个 `cautious`、33 个 `variable` 和 7 个 `insufficient`。硬保守要求邻域 dominant frequency `>=0.90`、coverage `>=0.80`、有效簇数 `>=50`，且全局/邻域优势残基一致并且全局频率 `>=0.80`。
- 保守位点与界面、Cys22/Cys95、末端 SSGS 合并后冻结 81 个 reported positions；剩余 47 个位置允许生成 846 条非 Cys 单突。这是待预测的完整约束空间，不是最终 30 条。
- 权威机器可读合同与结果位于 `docs/result_artifacts/input_baseline/vhh_conservation_20260818/`；下游必须读取合同，不得从本文复制残基列表。
- 已生成带IMGT FR/CDR标注的全局天然VHH、Nb252邻域及项目表达序列Logo。项目Logo以47条源序列为审核范围，仅纳入45条编号成功的H链序列；编号失败和非H链各1条保持显式排除，且产量不作为频率权重。

## 当前工具证据

- NetSolP：保留 S（solubility）与 U（usability）原始连续值；在 47 条产量数据中仅显示有限关联，不能单独决定候选。
- NanoMelt：预测 apparent melting temperature；与产量关系有限，只作为稳定性约束信号。
- AntiFold：使用实验复合物视图评价结构条件下序列相容性；缺失实验坐标位置不得伪装成可评价位置。
- TNP 与 nanoBERT 已完成探索性验证，但不在现行精简筛选工具集中。
- 所有工具验证须同时报告连续关联和预注册方向下的离散分类性能；阈值必须在训练折内选择，外层交叉验证报告 ROC-AUC、PR-AUC、MCC、balanced accuracy、sensitivity、specificity 及阈值稳定性。

## 下一执行路线

1. 以 846 条允许单突为统一输入，分别运行 NetSolP、NanoMelt 和实验复合物视图 AntiFold；未解析坐标导致 AntiFold 不可评价的候选保留明确缺失状态。
2. 远程运行已冻结的RP3Net 47条验证并按连续关联、逐样本留一和序列簇留一分类证据决定是否纳入；随后用同一离散合同补齐现行三工具验证。
3. 对 846 条单突进行硬风险审核，包括新增糖基化基序、未配对 Cys、强疏水/电荷斑块、Pro/Gly 结构风险和其他明显表达风险。
4. 在硬约束通过者中，以经过验证的表达预测证据、天然保守性等级、工具一致性和位置/机制多样性形成 30 条单突面板；WT 作为独立实验对照，不占 30 条名额。
5. 产量实验完成后才讨论组合；组合资格由真实单突效应决定，而不是当前预测分数。

## 阶段门

- `structure_and_interface_identity=pass`
- `natural_vhh_conservation_contract=pass`
- `expression_single_mutant_constraint_space=pass`：846 条仅表示允许进入预测的单突空间。
- `predictor_continuous_and_classification_validation=blocked`：RP3Net计划已就绪但正式评分未运行，现行三工具分类合同仍待补齐。
- `new_30_single_mutant_panel_release=blocked`：需要完成全空间预测、风险审核与分层选择。
- `combination_design_release=blocked`：等待单突实验结果。

## 本轮验证状态

- 本地 CPU 真实数据运行约 226 秒，无需 Slurm、checkpoint 或 resume。
- 全套测试`242 passed, 1 skipped, 4 subtests passed`；`pip check`、Python编译检查、Bash语法检查与`git diff --check`通过。
- 4,059 条输入、4,057 条合格序列、1,564 条邻域序列、128 个位置和 846 条单突均已回读核对；结果图已人工检查。
- RP3Net 47条计划已生成，固定checkpoint SHA-256为`443743bd031689aaf17dc6f7c22c5da3d23cf87b38e10341f114b27d651e6d2b`；远程正式分数和结论尚不存在。
- 本轮尚未提交或推送。
