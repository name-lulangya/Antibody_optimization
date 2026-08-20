# Codex 项目交接

Last updated: 2026-08-20 10:40:00

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

1. 三工具全空间补全已经完成，不在输入和工具合同不变时重复运行。重复一致性门共183项、0失败；NetSolP和NanoMelt逐值一致，AntiFold最大绝对差为`5.72e-6`并通过固定`1e-5`容差。
2. 847行完整原始评分矩阵已经释放：721条NetSolP/NanoMelt沿用经重复验证的旧值，126条缺口采用本轮新算值；AntiFold为721条三视图和126条AF3-only证据。实验缺失位置的实验视图保持`not_evaluable`。
3. 1,336条简并稳定词已按固定12符号、大小写敏感、允许重叠的精确连续子串合同加入847行矩阵。Nb252 WT没有稳定词命中；22条单突新增24个命中，825条不变，无减少项。该步骤没有筛选、排名或生成Tier。
4. 47条BL21产量验证不支持“稳定词密度越高、产量越高”：主指标来源分层Spearman为-0.185，95% bootstrap区间为-0.572–0.216，序列簇留一ROC-AUC为0.431、MCC为-0.277，证据等级为`no_supported_use`。稳定词只作为用户指定的可解释软偏好，不能覆盖硬约束、明确性质恶化或AntiFold/NanoMelt的既定用途。
5. 847条单突的四指标全景图已经生成：四张位置×替换残基热图分别展示NetSolP ΔU、NetSolP ΔS、NanoMelt预测ΔTm和AntiFold ΔlogP，并另有ΔU–ΔS及AntiFold–ΔTm两张散点图。AntiFold优先采用721条实验复合物视图结果，126条实验坐标不可评价候选以独立AF3 VHH-only结果补充；来源在逐候选表、图例和散点点形中显式区分。22条稳定词新增候选均用星号标记。
6. 幅度分档试选已经完成：847条中40条满足“无明显恶化且至少一个工具家族中等改善”；其中2条新Pro风险和2条含两个中等恶化的候选不放行，得到26条严格核心和10条单中等恶化的受控权衡候选。用户复核后，以新增稳定词的`T99F`替换原受控权衡候选`T99N`：现试选30条包含26条严格核心、3条受控权衡和1条稳定词假设探索候选，另有7条受控权衡替补。`T99F`不在40条幅度短名单中，其U/S/Tm/AntiFold均无中等改善，也无中等或明显恶化；纳入只用于实验检验稳定词假设，不得写成多工具支持的优化候选。
7. 该30条仍是待用户复核的计算试选，不是已冻结实验面板。当前仍需重点审查F30共9条、Q1共5条、T27共4条的位点集中是否可接受；不得为了凑多样性重新让微弱变化参与排名。`T99F`在图表中以稳定词探索类别单独标记。

## 阶段门

- `structure_and_interface_identity=pass`
- `natural_vhh_conservation_contract=pass`
- `expression_single_mutant_constraint_space=pass`：847 条仅表示允许进入预测的单突空间。
- `predictor_continuous_and_classification_validation=pass_with_tool_specific_roles`：NanoMelt正式判为仅稳定性约束；AntiFold因结构覆盖不足正式判为分类不适用、仅实验复合物相容性约束；不得把该门解释为产量预测器已验证。
- `expression_property_completion_plan_v2=pass`：721条精确复用、126条补算和12条抽样复核计划已冻结。
- `expression_property_complete_matrix_v2=pass`：847条NetSolP、NanoMelt和AntiFold证据已按视图可评价范围完成并通过完整性审核；该门不代表候选已筛选。
- `stable_word_feature_evaluation_v1=pass`：847条单突已完成稳定词新增/减少审计，47条产量验证也已完成；`pass`只表示制品完整，不表示稳定词已被验证为产量预测器。
- `expression_single_mutant_four_metric_landscape_v1=pass`：847条、48个位置、721条实验复合物AntiFold值和126条AF3补充值均已准确绘图并保留来源；该门不执行候选选择。
- `expression_single_mutant_trial_selection_v2=pass`：40条幅度短名单、26条严格核心、10条受控权衡、30条试选和7条替补已生成；试选含1条用户指定的`T99F`稳定词探索例外，原始同档小数没有参与排序。
- `new_30_single_mutant_panel_release=blocked`：试选30条尚待用户复核位点集中、3条受控权衡和1条稳定词探索候选，不能直接视为最终交付。
- `combination_design_release=blocked`：等待单突实验结果。

## 本轮验证状态

- 本地 CPU 真实数据运行约 226 秒，无需 Slurm、checkpoint 或 resume。
- RP3Net正式运行覆盖47/47条序列；连续、分类、逐折预测、结果图和gate均已生成并完成schema及计数核验。
- 4,059条输入、4,057条合格序列、1,564条邻域序列、128个位置和847条单突均已回读核对；v2相对v1只新增Q5V，未删除或意外开放其他候选。
- v2下游preflight为pass：847条候选中846条来自常规非Cys扫描、1条为Q5V共识回变，无多突、新Cys、冻结或界面突变。
- v2性质完整矩阵为pass：847条ID和序列均唯一且均为合法128-aa单突；721条复用值与126条新算值来源明确，全部保留末端SSGS，无冻结位点突变或新增Cys。NetSolP、NanoMelt和AF3 AntiFold均覆盖847条；实验单体/复合物AntiFold仅覆盖坐标可评价的721条，其余126条明确为`not_evaluable`。
- 稳定词输入共1,336条且均符合固定12符号合同；847条单突与完整性质矩阵逐ID联接，24条长表变化均覆盖实际突变位点。47条验证保留31条个体数值和16条LLJ有序/删失语义，未训练高容量模型，也未执行候选选择。
- 四指标景观图的数据表和gate已回读，847条候选、48个位置、22条稳定词新增、721条实验复合物AntiFold值与126条AF3补充值计数一致；热图和散点PNG均已人工检查，位置标签、FR/CDR标记、独立色标、星号及AntiFold来源点形均可辨识。
- 试选流程逐条回读847条审计、40条短名单、30条试选和7条替补；30条均为唯一128-aa单突且无明显不利档位和声明的硬序列风险。29条至少有一个中等/明显有利工具家族；唯一例外`T99F`新增一个长度5的稳定词且四指标仅为中性/微弱负向。试选包含24条实验复合物AntiFold和6条AF3补充来源，覆盖13个位置；PNG已人工检查，`T99F`以金色星号和边框单独标识。
- NanoMelt分类制品覆盖27条数值样本并保留4条未评分状态；正式图已人工检查。AntiFold分类不适用合同确认1条匹配结构、46条结构缺失且不生成任何分类指标。
- RP3Net gate固定为`rp3net_not_supported_for_candidate_use`；其模型分数不能解释为mg/L或通用表达阈值。
- PLM_Sol正式运行覆盖47/47条序列；31条数值记录和16条LLJ有序/删失记录语义保持不变，结果图、gate和run summary已回读核对。其固定5 mg展示的序列簇留一ROC-AUC/PR-AUC/MCC为0.761/0.685/0.411，但该展示不参与工具准入。
