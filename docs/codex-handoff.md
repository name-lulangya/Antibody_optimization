# Codex 项目交接

Last updated: 2026-08-18 17:10:00

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是 replace-in-place 的当前状态快照，不是追加式历史。2026-08-18导师讨论后的路线替代此前亲和力、多工具和组合突变路线。

## Project Goal

当前唯一显式优化目标是提高NK2R纳米抗体Nb252在BL21体系中的表达产量。项目不再显式优化亲和力；通过冻结实验界面、保守位点和既有不可变位点，降低表达优化破坏结合功能与VHH基本结构的风险。

## Current Non-Negotiable Design Decisions

- 实验复现的24个Nb252界面残基全部冻结，不再生成任何界面突变。权威集合仍来自`docs/result_artifacts/input_baseline/reviews/nb252_critical_residue_sets.json`，reported indices为`33,37,45,46,47,58,98,100-116`。
- 通过批准的多序列比对确定保守位点；保守位点冻结。比对输入范围、序列清理、编号映射、gap处理、去冗余方法、保守性指标和阈值必须先形成机器可读合同，不能从图或聊天重建。
- Cys22、Cys95和reported positions 125–128的末端`SSGS`继续不可突变；不得新增未配对Cys或删除母本残基。
- 新的实验交付面板必须恰好包含30条单突序列。当前阶段不得生成或交付双突/多突；只有单突表达实验完成并复核后，才讨论有实验依据的组合。
- PyRosetta、Flex ddG及既有Rosetta能量不再用于现行候选生成、排序或淘汰。既有结果只保留为历史provenance，不运行新的Rosetta任务。
- 现行评分工具仅为NetSolP、NanoMelt和AntiFold。可以依次测试一至两个新的表达相关工具，但只有在47条产量数据上完成覆盖、连续关联和分类验证后，才决定是否替代或补充NetSolP；不保留同类型重复工具作简单投票。
- 所有预测分数都是模型信号，不是实测产量、溶解度、Tm或结合功能。最终结论必须来自BL21同体系实验。

## Authoritative Inputs and Baseline

- 母本为合作者确认的完整128-aa reported Nb252；末端`SSGS`属于真实构建体而非linker。
- 47条reported序列及产量来自BL21体系，合作者确认不同来源数值可直接比较。LTT/WCC的31条individual approximate记录可作数值分析；LLJ的16条仍保持group anchor/lower-bound或有序/删失语义，不能伪造个体点估计。
- 实验结构、链身份、128-aa序列映射和严格`<4.0 Å`界面已释放并通过；实验缺失坐标为reported indices`9-15,24-29`，仍保持`not_evaluable`结构语义。
- 24个界面位点以前属于“谨慎可突变”，本次路线变更后已提升为硬冻结集合。候选生成代码和新合同必须显式更新这一语义，不能只修改文档。

## Active Tool Contracts

- NetSolP：远程环境`/data/software/env/luly25/netsolp`，固定从`/homes/Tianlab/luly25/software/netsolp`调用Distilled模型，输出Solubility（S）和Usability（U）。既有连续关联较弱，只能在新分类验证前作为待验证信号。
- NanoMelt：远程环境`/data/software/env/luly25/nanomelt`，输出预测表观Tm。既有结果不支持把Tm作为连续yield代理；现路线仅把它作为潜在表达/折叠相容性辅助指标并重新进行分类验证。
- AntiFold：远程环境`/data/software/env/luly25/antifold`，AntiFold 0.3.1；固定模型位于`/homes/Tianlab/luly25/software/AntiFold/models/model.pt`，CUDA评分使用`num_threads=0`。现路线只评价允许位点单突对既定Nb252结构的条件相容性。
- PyRosetta、Flex ddG、TNP、nanoBERT和批量AF3均退出当前执行路线。环境和历史制品保留，但不得作为新候选依据。

## Predictor–Yield Validation Contract

- 每个工具先保留原始连续指标，报告来源分层或来源内的Spearman/Kendall关联、置信区间和置换检验；连续相关性不再是唯一评价。
- 同时建立预先声明的离散评价：将reported yield划分为高/低（如数据支持可增加中间不确定区），将每个工具指标按训练数据中的候选阈值离散化，计算ROC-AUC、PR-AUC、MCC、balanced accuracy、sensitivity、specificity及混淆矩阵。
- 产量标签阈值、预测方向、指标阈值优化目标和tie-break规则必须在查看最终分类表现前冻结。阈值应在训练折内选择，外层使用LOOCV、重复分层交叉验证或按序列簇留出进行评估；不得在同一47条上选择阈值后再把表观AUC/MCC称为泛化性能。
- 类别不平衡时PR-AUC和MCC与ROC-AUC并列报告。报告阈值在重采样中的分布/稳定性；如果阈值不稳定或样本外性能不足，则该工具不进入候选排序。
- LLJ记录只在其分组/删失语义与所选产量分类阈值能够确定类别时纳入；跨越阈值或无法唯一分类的记录标为不确定并从对应二分类评价中排除，不作点估计填补。
- 新工具与NetSolP在相同样本、标签和交叉验证切分下比较。只有覆盖充分、样本外分类性能有增益且输出与项目目标相符时，才替换或补充NetSolP。

## Current Expression-Only Workflow

1. **更新设计合同**：把24个界面位点从谨慎集合提升为硬冻结集合，并使所有候选入口读取同一机器可读约束；旧的界面亲和力候选全部失效。
2. **建立MSA与保守位点合同**：确定用于保守性判断的序列范围和清理规则，对齐到Nb252 reported 1-based索引，输出逐位点频率、gap/覆盖、保守性得分和最终冻结集合。对低覆盖、截短、ANARCII失败或非同类链记录必须显式处理。
3. **重新验证工具**：在47条BL21产量数据上，对NetSolP、NanoMelt、AntiFold同时完成连续关联和预声明二分类评价；随后至多逐一测试两个新工具，采用完全相同的评价切分和指标。
4. **生成完整允许单突空间**：仅在非界面、非保守、非Cys22/Cys95、非末端`SSGS`且序列身份可确认的位置枚举单替换；先应用化学liability和VHH结构硬风险，再运行三个现行工具。
5. **形成30条单突面板**：依据经验证的表达相关指标、AntiFold结构相容性、风险控制、位点和替换类型多样性形成30条；保留原始各指标和入选理由，不生成组合突变。
6. **BL21实验验证**：WT与30条单突在同一构建体和实验条件下平行表达、纯化和定量，设置重复及预定义成功判据。
7. **实验后组合**：只有单突结果返回后，才根据真实表达改善、机制互补和风险决定是否组合；组合序列必须重新评价并单独实验验证。

## Superseded Historical Results

- 既有456个界面单突、Flex ddG复核、Rosetta亲和力核心、86条双突及原`final_candidate_panel_20260817`均属于已完成的旧路线结果，不删除、不改写，但不再作为当前订购或实验面板。
- 旧的合作者交付包包含亲和力和双突候选，已被本次路线变更撤销；在新的30条表达单突面板生成前，不得继续发送或订购其中序列。
- 旧周报可用于说明已完成的探索工作，但面向导师或合作者的后续材料必须明确其候选结论已被新目标替代。

## Current Gates

- `structure_and_interface_identity=pass`：结构身份、链角色、映射和24位界面集合可复用。
- `expression_single_mutant_design_release=blocked`：等待机器可读的界面硬冻结更新、MSA保守位点合同和新的工具分类验证合同。
- `new_30_single_mutant_panel_release=blocked`：旧30条已撤销；新面板尚未生成。
- `combination_design_release=blocked`：必须等待30条单突的真实BL21表达实验结果。

## Required Next Steps

1. 更新机器可读关键残基/候选合同及相关测试，使任何界面突变和多突候选都被主动阻断。
2. 在实施MSA前确认保守性分析的输入范围和判定口径，特别是47条序列中的截短、失败编号、轻链样记录和高度相似序列如何处理。
3. 实现统一的连续+离散工具验证入口，先重评NetSolP、NanoMelt和AntiFold，再决定是否测试新工具。
4. 三项前置合同通过后重新生成允许单突空间和新的30条表达单突面板；不得复用旧最终30条作为订单源。

## Verification State

- 本轮只替换项目规则和交接计划，尚未修改候选生成代码、机器可读关键残基集合或任何生成结果。
- 既有测试和历史制品仍保持原状态；新路线在完成合同迁移和新增测试前保持blocked。
