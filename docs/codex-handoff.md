# Codex 项目交接

Last updated: 2026-08-19 22:45:00

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是当前状态快照，不是追加式历史。现行路线以 Nb252 在 BL21 体系中的表达量优化为唯一设计目标。

## 当前目标与硬约束

- 最终向合作者提供 30 条 Nb252 单突序列；单突实验完成前不设计或交付组合突变。
- 冻结实验复合物中已复现的 24 个 VHH 界面位点，不再显式优化亲和力，也不在当前路线使用 Rosetta 排序。
- 冻结天然 VHH 邻域中“Nb252亲本残基等于全局/邻域共同优势残基”的高置信保守位点、Cys22/Cys95 和末端 SSGS（reported 125–128）；对高保守但亲本偏离共识的位置只开放共识回变。
- 候选必须保持完整 128-aa 父序列长度、末端 SSGS、两枚原有 Cys，且不得引入新 Cys。
- 现行核心预测工具仅为 NetSolP、NanoMelt 和实验复合物视图 AntiFold。新工具必须先在 47 条可比产量数据上验证，证明有独立且可重复的样本外信息后才能纳入筛选。
- RP3Net 0.0.2 已完成47条正式验证，最终证据等级为`no_supported_use`，不得加入847条候选的生成、筛选或排序。
- PLM_Sol V1.0已完成47条正式验证，最终证据等级为`no_supported_use`，不得加入847条候选的生成、筛选或排序；Nb252单序列smoke分数只证明调用链可运行，不是实测溶解度或产量。

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
- v2将reported 128个位置分为54个`hard_conserved`、1个`conserved_nonconsensus`、33个`cautious`、33个`variable`和7个`insufficient_evidence`。硬保守除要求邻域dominant frequency `>=0.90`、coverage `>=0.80`、有效簇数`>=50`、全局/邻域优势残基一致且全局频率`>=0.80`外，还要求Nb252亲本残基等于该共同优势残基。
- 唯一`conserved_nonconsensus`为reported Q5：邻域/全局优势残基均为V，频率分别为0.9980/0.9976，因此不把Q5亲本状态称为天然硬保守，也不开放任意扫描，只允许`Q5V`共识回变。
- 54个亲本匹配的硬保守位点与界面、Cys22/Cys95、末端SSGS合并后冻结80个reported positions；47个常规可扫描位置产生846条非Cys单突，另加Q5V，共847条。这是待预测的完整约束空间，不是最终30条。
- 权威机器可读合同与结果位于 `docs/result_artifacts/input_baseline/vhh_conservation_consensus_v2_20260819/`；下游必须读取合同，不得从本文复制残基列表。旧`vhh_conservation_20260818`仅保留为历史v1 provenance。
- 下游阶段边界preflight已核对847条候选、48个可变位置、80个硬冻结位置、24个界面位置及CSV/FASTA/亲本一致性；后续三个评分入口复用该通过结果，不再重复同一检查。
- 已完成旧结果的精确回映计划：当前847条中，721条具有可复用NetSolP/NanoMelt及三个AntiFold视图，126条仅缺NetSolP/NanoMelt且AntiFold为AF3-only。缺口严格位于reported 11、14、24、26–29，每个位点18条；Q5V已有三工具旧结果。
- 已生成带IMGT FR/CDR标注的全局天然VHH、Nb252邻域及项目表达序列Logo。项目Logo以47条源序列为审核范围，仅纳入45条编号成功的H链序列；编号失败和非H链各1条保持显式排除，且产量不作为频率权重。

## 当前工具证据

- NetSolP：保留 S（solubility）与 U（usability）原始连续值；在 47 条产量数据中仅显示有限关联，不能单独决定候选。
- NanoMelt：43/47条可评分，其中27条具有LTT/WCC个体数值产量。嵌套逐样本/序列簇留一结果相同：ROC-AUC 0.571、PR-AUC 0.554、MCC 0.408、balanced accuracy 0.701、sensitivity 0.786、specificity 0.615；连续与分类联合证据仍为`no_supported_use`，不得作为BL21产量排序器，只保留为预测稳定性约束。
- AntiFold：47条产量序列中只有Nb252具有匹配实验NK2R复合物，另外46条没有可比实验复合物，而且输出是逐位点结构条件概率而不是统一表达分数；正式分类状态为`not_applicable`，不得报告伪AUC/MCC。其唯一现行用途是实验复合物视图中坐标可评价位点的突变相容性约束。
- TNP 与 nanoBERT 已完成探索性验证，但不在现行精简筛选工具集中。
- RP3Net：31条数值记录的直接合并Spearman为0.476，但来源内分层Spearman仅0.198且95% bootstrap区间跨0；LLJ有序Kendall为-0.451，与预声明方向相反。嵌套分类ROC-AUC为0.621、PR-AUC为0.620、MCC为0.313，未达到预声明综合门，因此不支持候选使用。
- PLM_Sol：47/47条评分成功。31条数值记录的来源内分层Spearman为0.473，95% bootstrap区间为0.096–0.749，但WCC内部Spearman为-0.096；嵌套分类ROC-AUC为0.638、PR-AUC为0.662、MCC为0.313。PLM_Sol与NetSolP U/S高度重叠，在NetSolP S基础上的序列簇外增量为-0.140，因此gate为`plm_sol_not_supported_for_candidate_use`，不纳入候选使用。固定5 mg结果仍仅供展示。
- 固定5 mg探索图：31条数值记录显示为高产14条、低产17条；RP3Net、NetSolP U和NetSolP S均展示训练折最大MCC阈值及留出指标。该图仅用于直观展示，不作为工具准入、候选筛选或5 mg阈值有效性的正式证据。
- 所有工具验证须同时报告连续关联和预注册方向下的离散分类性能；阈值必须在训练折内选择，外层交叉验证报告 ROC-AUC、PR-AUC、MCC、balanced accuracy、sensitivity、specificity 及阈值稳定性。

## 下一执行路线

1. 远程运行`bash scripts/candidate_design/submit_expression_property_completion_v2.sh`。作业先以WT、Q5V、性质极值及FR/CDR代表组成的12条复核面板验证NetSolP/NanoMelt旧值，并对三个固定WT结构视图各重跑一次AntiFold；只有逐值精度门通过才继续。首次复核中NetSolP和NanoMelt全部通过，AntiFold仅出现最大`6.04e-6`的GPU浮点复现差异；其绝对容差已按实测精度固定为`1e-5`，等待重新运行确认。
2. 仅补算reported 11、14、24、26–29的126条NetSolP/NanoMelt结果；复用全部847条AntiFold结果。实验缺失位置的实验视图继续为`not_evaluable`，AF3-only值作为预测结构补充证据独立保留。
3. 合并并审核847行完整原始评分矩阵；只记录绝对值、相对当前WT变化、评价状态和来源，不筛选、不计算Tier或综合分。
4. 完成用户计划中的后续独立步骤后，再制定风险审核和30条单突选择合同；当前不得提前淘汰候选。

## 阶段门

- `structure_and_interface_identity=pass`
- `natural_vhh_conservation_contract=pass`
- `expression_single_mutant_constraint_space=pass`：847 条仅表示允许进入预测的单突空间。
- `predictor_continuous_and_classification_validation=pass_with_tool_specific_roles`：NanoMelt正式判为仅稳定性约束；AntiFold因结构覆盖不足正式判为分类不适用、仅实验复合物相容性约束；不得把该门解释为产量预测器已验证。
- `expression_property_completion_plan_v2=pass`：721条精确复用、126条补算和12条抽样复核计划已冻结。
- `expression_property_complete_matrix_v2=blocked`：等待远程重复一致性门和126条补算完成。
- `new_30_single_mutant_panel_release=blocked`：需要完成全空间预测、风险审核与分层选择。
- `combination_design_release=blocked`：等待单突实验结果。

## 本轮验证状态

- 本地 CPU 真实数据运行约 226 秒，无需 Slurm、checkpoint 或 resume。
- RP3Net正式运行覆盖47/47条序列；连续、分类、逐折预测、结果图和gate均已生成并完成schema及计数核验。
- 4,059条输入、4,057条合格序列、1,564条邻域序列、128个位置和847条单突均已回读核对；v2相对v1只新增Q5V，未删除或意外开放其他候选。
- v2下游preflight为pass：847条候选中846条来自常规非Cys扫描、1条为Q5V共识回变，无多突、新Cys、冻结或界面突变。
- v2性质补全计划为pass：847/847条AntiFold旧证据精确回映；NetSolP/NanoMelt为721条可复用、126条待补算；连接同时核对位置、WT、mutant和完整128-aa序列，不依赖旧candidate ID。
- NanoMelt分类制品覆盖27条数值样本并保留4条未评分状态；正式图已人工检查。AntiFold分类不适用合同确认1条匹配结构、46条结构缺失且不生成任何分类指标。
- RP3Net gate固定为`rp3net_not_supported_for_candidate_use`；其模型分数不能解释为mg/L或通用表达阈值。
- PLM_Sol正式运行覆盖47/47条序列；31条数值记录和16条LLJ有序/删失记录语义保持不变，结果图、gate和run summary已回读核对。其固定5 mg展示的序列簇留一ROC-AUC/PR-AUC/MCC为0.761/0.685/0.411，但该展示不参与工具准入。
