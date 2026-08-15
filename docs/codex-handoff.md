# Codex 项目交接

Last updated: 2026-08-15 19:15:00

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是 replace-in-place 的当前状态快照，不是追加式历史。

## Project Goal

优化 NK2R 纳米抗体 Nb252 的亲和力、稳定性和表达相关性质，同时保留结构、序列、实验和计算证据的可追溯性。

## Environment and Git

- 本地环境：`D:\miniconda\envs\ab_optim`，Python 3.11.15、ANARCII 2.0.8、Gemmi 0.7.5。
- 远程项目环境：`/data/software/env/luly25/ab_optim`；计划检出父目录 `/homes/Tianlab/luly25/`，登录别名尚未建立。
- PyRosetta：`/data/software/env/luly25/multi_ligand`，Python 3.10.20，PyRosetta 2026.03，Rosetta commit `5e498f1409c68ade56c8ce5842bf79e1b02e8db4`。
- nanoBERT：`/data/software/env/luly25/vhh-lm`，`NaturalAntibody/nanoBERT` revision `edc8182ad89a827f8737fa572c6b5fac6197e6b0`，使用已记录离线缓存。
- NetSolP：`/data/software/env/luly25/netsolp`；官方5.63 GB发行包解压在`/homes/Tianlab/luly25/software/netsolp`。固定从该顶层工作目录调用`predict.py`，使用`MODEL_TYPE=Distilled`、`PREDICTION_TYPE=SU`，以Usability为主指标、Solubility为辅指标；环境快照保存在远程软件目录的`environment.freeze.txt`。官方3序列测试同时完成S/U预测，实测57.812787秒。
- TNP：`/data/software/env/luly25/tnp`，TNP 0.0.1 commit `29dcac72f1380e8538e8870f45a699d3c6156162`、ImmuneBuilder 1.2、ANARCI 2024.05.21、Biopython 1.77、OpenMM 8.5.2、DSSP 4.6.1、Torch 2.7.1+cu126。固定入口为该环境的`bin/TNP`，并设置`PYTHONPATH=/homes/Tianlab/luly25/software/TNP`。ImmuneBuilder安装版`refine.py`的OpenMM `Threads` set-literal错误已按上游正确dict语义修复并由`LTT__Nb294`失败前/成功后smoke验证；V2评分入口会在运行前核对该补丁。Nb252完整128-aa输入建模126 aa并仅裁掉末端`GS`，不改变authoritative母本或`SSGS`冻结规则。
- AntiFold：独立环境`/data/software/env/luly25/antifold`，AntiFold 0.3.1、Torch 2.2.0+cu121、PyG 2.4.0；模型位于`/homes/Tianlab/luly25/software/AntiFold/models/model.pt`。A100模型加载、官方纳米抗体评分及CDR1采样smoke均通过。0.3.1加载器忽略普通`checkpoint_path`，因此环境内固定模型路径使用符号链接，CUDA评分必须`num_threads=0`；完整版本和模型身份只记录在`antifold_validation_plan_20260815/antifold_environment_contract.json`。
- Git：`main`，远程 `git@github.com:name-lulangya/Antibody_optimization.git`；TNP V2真实结果提交`1a6e908`已同步，当前`HEAD=origin/main`。
- 阶段2的本地阶段0、远程WT安全导入、WT评分校准、456候选全量扫描、post-scan分层、50候选×20样本Flex ddG复核及三类yield关联验证均已完成。AntiFold环境与官方smoke已通过，项目最小验证计划/派生结构已生成，但Nb252三个真实视图尚未提交远程评分；亲和力组合与AF3终选复核也尚未运行。

## Frozen Inputs

- `Nb252-optimization.cxs`：537670 bytes，SHA-256 `1bc636c28f66ae60edc658d2e1c4aad0b07f4141ca5411c78662aa19da793c4d`，原文件未被覆盖。
- `nb序列及产量（1L）.docx`：14172 bytes，SHA-256 `a6e4022f0978fbd70a0e04dc78f479140ab6f55caaa90b467fb77a62eb5db5d1`。
- 47 条 reported 序列的 provisional IMGT 审核：46 pass、1 failed。`WCC__4-28` 为 ANARCII `Score less than cut off.`；Nb252的前126 aa被编号，末端`GS`未编号。
- 合作者已确认完整128-aa reported Nb252就是authoritative design parent，末端`SSGS`不是linker；用户进一步明确“保留”是候选设计硬约束：reported-sequence positions 125–128 的`SSGS`不得替换或删除。未编号的末端`GS`仍属于构建体，该约束必须用reported-sequence index表达，不能虚构IMGT编号。机器可读约束位于`input_baseline/reviews/nb252_design_constraints.json`。

## Structure and Interface Baseline

- `data/structures/cxs_exports/` 是 ChimeraX 1.12 的真实导出：3 个目标原子模型、5 个 mmCIF、颜色清单、manifest 和 run summary 均已验证；5 个 mmCIF 可由 Gemmi 0.7.5 读回。
- CXS 中的精确模型名为 `NK2R-252.pdb`、`NK2R-NKA` 和 `fold_2r_252_nomg_model_0.cif`。导出入口已修复 ChimeraX 非 `__main__` sandbox 调用问题。
- 用户已在 ChimeraX 1.12 中确认链角色：实验 `C=Nb252`、`R=NK2R`；NKA 模型 `L=NKA`、`R=NK2R`；AF3 `A=Nb252`。
- 用户确认实验 chain C 的 24 个精确 `[255,165,0,255]`、`atom+ribbon` 位点就是合作者橙色区域。机器可读记录：`docs/result_artifacts/input_baseline/structure_review_20260810/baseline_review.json`。
- `structure_released_20260810/` 已在authoritative 128-aa确认后重建可逆reported sequence/IMGT/实验/AF3映射和FR-only Cα对齐。实验结构115个有坐标VHH残基通过source auth exact-WT映射；AF3 126个有坐标残基通过source label ID和536-aa polymer中唯一的128-aa Nb252连续片段映射。
- 82 个共同 framework Cα、Kabsch、无 outlier rejection：RMSD 0.631994 Å。拟合后 FR aggregate RMSD 0.631994 Å；CDR3 RMSD 6.490853 Å、最大位移 10.733119 Å。实验结构仍是结合构象证据，AF3 始终是预测。
- `interface_released_20260810/`：严格 polymer heavy-atom center `<4.0 Å`，排除 H/D、非正 occupancy、水/配体/糖/离子、晶体/NCS images，并遵守 altloc 兼容性。结果为246个原子接触对、24个VHH界面残基。
- 严格 `<4.0 Å` 残基集合与确认橙色 24 位点完全相同；保护并集为 reported sequence index `33,37,45,46,47,58,98,100-116`。该集合仅用于保守保护，不是能量热点或突变效应结论。
- 关键残基集合已集中绑定到机器可读文件`docs/result_artifacts/input_baseline/reviews/nb252_critical_residue_sets.json`：实验缺失坐标为reported indices `9-15,24-29`，在实验界面语义中必须保持`not_evaluable`；24个复现界面位点为谨慎突变而非禁区；positions `125-128`的`SSGS`不可突变。所有候选生成和评分入口必须重新读取并核对该文件及其上游哈希，不得从聊天、handoff或记忆重建残基集合。

## Expression Audit

- 合作者已确认LTT、WCC、LLJ的reported yield可直接相互比较；机器可读证据位于`input_baseline/reviews/expression_cross_provider_confirmation.json`。
- 合作者已确认这47条产量来自BL21表达体系；该事实使E. coli可溶表达/纯化可用性预测与本项目标签具有直接体系相关性，但预测分数仍不能等同于实测产量。
- Expression audit 1.1.0允许47条记录进行保留观测语义的跨来源联合探索：31条LTT/WCC individual approximate按数值处理，LLJ的9条group lower bound和7条group approximate仍按删失/分组信息处理。
- 跨来源pooling为`pass`，但普通连续回归、把LLJ分组插值为个体点估计，以及未经验证地向新Nb252突变转移仍被禁止。

## Current Gates

- `input_freeze_manifest.status=pass`。
- `local_baseline_build=pass`：结构导出、清单、链身份、可逆映射和临时界面安全均已完成。
- `candidate_design_release=pass`：authoritative 128-aa Nb252、结构身份、映射和界面安全均已确认；设计必须维持实验表位和结合构象。
- `pooled_expression_model_release=pass`：允许保留删失/分组语义的联合建模；独立的`nb252_expression_transfer`仍为`blocked`，直到模型经适当验证。
- `stage0_local_contract=pass`、`candidate_manifest_release=pass`：关键来源哈希、母本、映射、链身份、缺失坐标、界面、SSGS和实验二硫键均通过本地重检。
- `pyrosetta_wt_import_release=pass`：远程WT导入门的396个polymer残基、4个断点、PDBInfo映射、Cys22–Cys95二硫键和有限raw score均通过；该阶段当时只释放协议校准，其`ready_for_scoring_protocol_calibration`状态现已被下述v2 gate取代。
- source结构有20个侧链不完整的标准残基（Nb252 3、NK2R 17），按标准重原子计数共缺90个；Nb252 auth 102 TYR属于已确认界面。PyRosetta补建这些原子产生的warning不单独阻断，但校准输出必须记录该事实，且界面邻域使用统一repack，避免WT/候选准备不对称。
- v2 `pyrosetta_affinity_scoring_release=pass`，唯一选择`interface_repack_constrained_min`；repack-only因8/8重复的`dG_separated`和跨界面能均为正而被阻断。选中协议8/8重复两项均为负：中位`dG_separated=-58.147156 REU`（MAD 2.357017）、中位跨界面能`-82.876253 REU`，界面`fa_rep`中位22.466594 REU，较raw 183.483686下降87.755536%；最大界面Cα RMSD 0.071804 Å。
- 代表prepared WT（seed 8112027）保持21/24个原实验VHH接触并新增2个，保持34/37个NK2R表位残基并新增4个。VHH source-auth E46、D101、I103跨过4 Å阈值成为`lost`，S50、T53成为`gained`；其中I103为4.000093 Å的阈值边界变化，D101变化更明显。候选仍按原实验24位界面和原表位约束生成，prepared接触集合只作评分/QC，不能反向改写实验界面；E46/D101/I103候选不得仅凭单一prepared构象淘汰。
- 科学review状态为`released_for_paired_relative_candidate_scoring`：候选必须从该prepared WT按同一受约束局部协议比较，并保留映射、断点、二硫键、接触/表位和构象门。Rosetta分数仍只是相对排序信号，不是实测亲和力或膜蛋白绝对稳定性；无需RosettaMP、缺失区补全或全局relax。
- pilot v1已作为历史诊断保留：其3 pass/9 blocked来自候选绝对0.90表位门低于配对WT自身33/37表现，不能作为候选分类。
- pilot V2已远程完成：1146.969247秒，3个共享WT、36个突变体评估、12个候选摘要全部通过运行安全，所有摘要均为`not_applied_scan_stage`，没有执行候选筛选。接触集合、计数和相对同seed配对WT保持率已逐行重算一致；NK2R表位配对WT保持率全部为1.0，VHH配对WT保持率范围0.958333–1.0。V1/V2共同数值字段逐字一致，证明V2修正的是门和记录语义而非评分轨迹。
- `affinity_pilot_v2_scientific_review.json`将路线发布为`released_for_full_456_scan_implementation`。全量合同固定为先计算全部456候选并核验合并完整性，再统一筛选；pilot能量只作未筛选诊断，不产生实验候选推荐。
- 456×3全量扫描已完成并通过merge gate：12片累计计算42645.718076秒，456个候选、1368个唯一重复键、3个去重WT均完整，0个runtime failure，全部摘要为`not_applied_scan_stage`。接触集合/保持率和`mutant-WT`能量差已逐行复算一致；`affinity_full_scan_scientific_review.json`发布`ready_for_post_scan_filter_implementation`，尚未选择或淘汰任何候选。
- 统一post-scan分层已在本地完成：456条全部通过硬有效性门并互斥分为Tier 1/2/3/4/5=`18/30/39/82/287`；Tier 1/2合计48条形成严格复核池。分层使用三重复双能量方向、配对WT接触保持和`Δfa_rep`，风险标签与层内Pareto保持独立；`candidate_selection_performed=false`，尚未形成最终实验面板。
- 亲和力ensemble核心筛选已完成：50条完整证据中8个单点通过“两项能量均至少18/20样本为负且中位数均为负”的唯一核心门，分别为`R45C/R45V/D101W/I103W/E105F/E105L/N107A/S114M`，覆盖6个位置；R45和E105的同位点替换互斥。风险、接触保持、prepared敏感性与Pareto层仍独立保留；未生成双突或最终实验序列。
- 稳定性/表达WT发现合同已完成并作为历史范围证据保留：其81个非界面framework位点不再整体冻结，也不能直接视为可突变清单。当前允许在FR和非关键CDR中提出有明确结构、序列和性质依据的稳定性/可开发性候选，但必须重新读取硬约束，并通过亲和力非劣、表位/接触保持和完整多目标复核。
- nanoBERT—reported yield验证gate为`pass`，但证据等级为`no_supported_use`、release为`nanobert_not_supported_for_candidate_use`。主指标完整reported序列mean PLL在31条LTT/WCC上的来源分层Spearman为-0.2295，长度调整partial Spearman为-0.2020，LTT/WCC内部为-0.2413/0.0602；95% bootstrap CI为[-0.5982, 0.2219]，来源内置换p=0.2799，加入nanoBERT后的LOOCV Spearman增益为-0.0507。16条LLJ有序/删失Kendall tau-b为-0.2308。不得用nanoBERT为Nb252稳定性、表达或候选序列打分/排序。
- NetSolP—BL21 yield验证gate为`pass`，科学证据等级为`compatibility_filter_only`、release为`netsolp_compatibility_filter_only`。主指标U在31条LTT/WCC上的来源分层Spearman为0.3943、长度调整partial Spearman为0.2858，但LTT/WCC内部为0.4154/-0.1205，95% bootstrap CI为[-0.0665, 0.7478]、来源内置换p=0.05039；普通和leave-cluster-out相对provider-only增益为0.3600/0.2661，但其绝对预测Spearman均仅约0.108。LLJ有序/删失tau-b为-0.1648。U/S只能作为完整候选的辅助相容性信号，不能单独排名、硬过滤或解释为产量；辅助S的跨来源/CV表现更一致，但未经过主指标同等级的预声明重采样门，保持探索性。
- TNP V1依次尝试47条、实测1416.280578秒，38 pass/9 failed；诊断确认5条触发ImmuneBuilder OpenMM set-literal错误，`LTT__Nb294`补丁后smoke已恢复。另4条WCC为TNP不适用：`WCC__4-1/4-28/4-11`被TNP ANARCI拒绝，`WCC__4-42`被NanoBodyBuilder2判定缺失过多残基；不得补写或改造其真实序列。
- TNP V2验证gate为`pass`，43/43适用序列通过、4条WCC保留`not_applicable`，科学证据等级为`compatibility_filter_only`。唯一预声明主指标PSH在27条数值记录中的来源分层Spearman为-0.0856、长度调整partial Spearman为0.0665，95% bootstrap CI为[-0.5377, 0.3820]、来源内置换p=0.6901；LTT/WCC内部为-0.0841/-0.4000，16条LLJ tau-b为-0.0769。TNP PSH单独的普通/簇外CV绝对Spearman为-0.5925/-0.3867，加入NetSolP U后也低于U单独模型；不得用TNP预测或排序yield。Nb252的正式单次结果为PSH 136.6287、amber，在43条中第三高；其余5项flag均为green。TNP仅保留为完整候选的developability风险/相容性证据，优先避免相对WT进一步恶化，不作单独硬淘汰。
- NanoMelt—BL21 yield验证gate为`pass`，科学证据等级为`no_supported_use`、release为`nanomelt_not_supported_for_yield_use`。43条真实评分中27条为LTT/WCC数值记录、16条为LLJ有序/删失记录；分层Spearman为0.0404，95% bootstrap CI为[-0.4329, 0.5511]，置换p=0.8519，LTT/WCC内部为0.0396/0.2108，LLJ tau-b为0.0110。加入Tm后普通/序列簇外CV Spearman为-0.7481/-0.7469，较provider-only进一步降低0.1070/0.2562；不得用NanoMelt预测或排序yield。NanoMelt仍作为独立VHH热稳定性预测器，只在完整候选上相对WT评价，不把预测Tm称为实验Tm或设为产量代理。Nb252预测表观Tm为65.18 °C，在43条中处于18.6百分位。
- `finalize_input_baseline.py` 已用真实 structure mapping/interface 重建 canonical 128 位点图和 `stage1_gate.json`。程序 `status=pass` 不代表科学 release gate 自动通过。

## Stage-2 Phase 0 Result

- 正式制品位于`docs/result_artifacts/candidate_design/stage0_contract_20260810/`，run summary位于`docs/run_summaries/candidate_design/stage0_contract_20260810.json`。
- 128个位点中：实验缺失坐标13位、实验界面24位、硬冻结6位、首轮亲和力放行24位。硬冻结为reported positions 125–128的`SSGS`和实验SG–SG几何支持的Cys 22/Cys 95；界面位点是“谨慎但可突变”。
- 首轮亲和力以未批量补全的实验complex为主；缺失位点不进入首轮亲和力扫描。只有终选涉及缺失CDR1、或原实验结构与完整VHH预测轨道冲突时，才触发定点补全敏感性分析。

## Current Optimization Plan

1. **远程WT结构门（已完成，实测14.123538秒）**：PyRosetta 2026.03直接导入未补全实验complex，396个polymer残基和4个预期断点均安全，映射与二硫键通过；未做relax、补全、突变或候选评分。权威结果位于`structure_preparation/pyrosetta_wt_import_20260810/pyrosetta_wt_import_gate.json`。
2. **WT评分协议校准（已完成，实测465.154733秒）**：schema v2 gate通过并唯一选择受约束局部最小化；代表WT、逐位置接触变化、重复表、逐残基能量和QC图均已同步。该步骤只释放成对相对候选评分，不生成突变或给出实验亲和力结论。
3. **总体多目标架构（当前有效）**：不再按CDR/FR严格划分亲和力与稳定性/表达相关优化区域，也不再默认冻结全部framework。候选分为亲和力驱动和稳定性/可开发性驱动两个来源，但每条完整序列从单突阶段起都同步评价亲和力、AntiFold结构相容性、NetSolP U/S、NanoMelt相对WT预测Tm、TNP developability和化学风险。亲和力驱动候选必须满足性质非严重恶化，性质驱动候选必须满足亲和力、表位和构象非劣；完整组合逐条重算，不机械相加分数，也不把预测值解释为实测产量、Tm或亲和力。
4. **亲和力单点空间（本地已完成，实测约2秒）**：`affinity_single_mutants_20260811/`已从阶段0合同、实验24位界面、可逆编号和v2评分gate生成每个位点19种非WT替换，共456个单点；每条均保留完整128-aa序列、reported/IMGT/source-auth编号、实验接触和prepared WT敏感标记。12个pilot候选覆盖FR/CDR、保守/跨类替换及E46/D101/I103，未按prepared接触或主观化学偏好预删候选。
5. **成对PyRosetta亲和力pilot V2（已完成，实测19.116154分钟）**：12候选×3重复和共享WT运行全部有效，配对WT相对接触指标复算一致，扫描阶段未筛选候选；路线已释放456全量扫描实现。
6. **456单点全量扫描（已完成）**：12片累计11.846小时计算量，完整得到456候选×3重复且没有扫描中筛选。未筛选景观中113个候选`ΔdG<0`、143个跨界面能变化<0、87个两者均<0，但这些只是描述性计数；82个候选两指标方向不一致，说明后续不能使用单一能量列筛选。
7. **统一筛选（已完成，本地实测约1.5秒）**：完整456候选按统一规则分为`18/30/39/82/287`，48条Tier 1/2进入严格复核池；未删除候选、未形成最终实验推荐。
8. **Flex ddG生产参数计时pilot（已完成）**：4代表候选×2独立样本共8任务全部pass，累计1.666723 job-hours；单任务中位717.980720秒、P90 874.753734秒，backrub占累计任务时间94.7212%，峰值内存975.65–1000.62 MiB。全部WT/突变映射、断点和二硫键检查通过。两重复只验证生产参数、耗时和协议可行性；代表候选存在样本间能量方向翻转，因此不得用pilot分数排名。
9. **亲和力严格复核与核心选择（已完成）**：50候选×20样本=1000任务全部pass；随后按唯一双指标18/20方向门从完整50条中选择8个核心单点、覆盖6个位置。该选择不使用加权综合分，仍保留全部50条证据和每个核心的`fa_rep`、接触及prepared风险；当前只释放核心模块使用，不代表8条均等安全，也未生成亲和力双突。
10. **候选、组合与最终三类面板（当前有效）**：单突是可组合的证据模块而非设计上限。组合池包含亲和力×亲和力、性质×性质和亲和力×性质三类双突，必要时才保留极少量三突；同位点替换互斥。最终30条分为亲和力优先、稳定性/可开发性优先和平衡组合三类，类别表示计算证据优先级而非实验成功声明；具体配额按合格候选数确定，不为凑数降低安全门。现有`R45C/R45V`仅作高风险单突对照：前者因额外未配对Cys退出组合，后者因VHH hallmark和接触重排风险不进入常规组合。
11. **工具分工（当前选择）**：PyRosetta/Flex ddG提供相对亲和力及接触/表位非劣证据；AntiFold评价候选对既定Nb252结构的条件相容性，并可为允许区域的性质候选提供序列先验。NetSolP U/S只作表达/溶解性辅助相容性信号，TNP只作developability风险工具，NanoMelt只提供相对WT的VHH热稳定性预测；三者已验证不能单独或联合承担yield排序。nanoBERT已退出，不训练现有47条数据的Nb252局部产量排序器，也不部署同类型重复工具投票。AF3只复核少量终选完整组合构象。
12. **framework与硬约束政策（当前有效）**：FR可以提出稳定性/可开发性候选，但不是无约束全扫描。reported positions 125–128的`SSGS`、Cys22–Cys95实验二硫键和实验缺失坐标的首轮结构不可评估状态继续保持；额外未配对Cys、VHH hallmark、CDR支撑网络、核心疏水包装、实验界面和表位变化均作为明确风险或阻断条件。现有81位合同只提供历史范围，47条序列只提供同框架频率和背景先验，不得据此把某个FR替换称为表达因果突变；所有FR候选必须在完整序列上重新评价亲和力和结构安全。
13. **nanoBERT—yield验证（已完成）**：47条真实序列均完成固定revision逐残基single-mask评分；主指标未通过方向、跨来源稳定性、不确定性或LOOCV门，证据等级为`no_supported_use`。全样本合并ρ=0.1887与来源分层ρ=-0.2295方向相反，显示来源混杂；nanoBERT从稳定性/表达候选评价路线中移除，不作为排序或过滤信号。理化特征只产生未校正的探索性关联，当前也不得作为表达优化方向或训练模型依据。

## Verification

- 阶段0专项：`3 passed`；覆盖真实128位合同、过期哈希拒绝、CSV BOM/LF、拒绝覆盖、固定时间戳双跑六制品逐字节一致。
- 阶段0图、亲和力候选空间图、post-scan四面板图、ensemble核心图和稳定性/表达合同图均已人工检查，坐标轴、图例和说明无遮挡。当前全套验收为`187 passed, 1 skipped, 4 subtests passed`；新增5项AntiFold计划、映射、ΔlogP与结果门测试，唯一skip仍是Windows真实symlink权限测试。
- `pip check`、`python -m compileall -q src scripts tests`、`git diff --check` 均通过。
- NetSolP真实结果含47条样本、15项指标和完整gate；本地从紧凑样本表重算主指标全部10个关键统计及证据等级，与gate逐项一致。600 dpi PNG已人工检查，四个面板、坐标轴、图例和说明无遮挡。
- TNP V2真实结果含47条身份、43条pass、4条`not_applicable`、6项指标和4种CV模型；本地从紧凑证据表复算全部关联、5000次bootstrap/permutation及CV，最大差异不超过浮点舍入误差。结果PNG已人工检查，四面板、坐标轴、图例和说明无遮挡。
- NanoMelt真实结果含47条身份、43条pass、4条`nanomelt_not_scored`；本地从样本证据重算主要相关、10000次bootstrap/permutation、CV和证据等级，与gate逐项一致。未评分为`WCC__4-1/4-28/4-11/4-40`，与TNP集合只重合前三条；`WCC__4-42`由NanoMelt评分且预测61.99 °C。结果PNG已人工检查，四面板、覆盖说明、坐标轴和图例无遮挡。
- v2结果一致性复核：schema 2 gate为pass、16行重复数据和67行接触状态与gate一致、代表PDB含396个polymer残基；评分校准专项测试`8 passed`。
- pilot V2一致性复核：3行WT、36行候选重复、12行候选汇总与machine gate一致；0个runtime failure，36/36重复和12/12摘要为pass，全部摘要未筛选。接触集合计数及配对WT保持率复算一致；V1/V2共同数值字段逐字一致。
- 全量扫描plan为pass：456候选、1368预期突变体评估、12片×38、每片两个完整位置、最大并发4且`candidate_filtering_applied=false`。实现测试覆盖真实分片、完整合并、缺片/提前筛选阻断、最终PNG/SVG和Slurm依赖合同。
- 全量结果复核：456摘要、1368重复、3个WT与manifest/merge gate一致，所有candidate×replicate×seed唯一，接触和能量差复算通过。最终PNG已人工检查；右侧重复Y轴标签已从绘图代码和已有plot data修正，不改评分数据。
- 修复覆盖：ChimeraX sandbox 入口、实际模型名、无 polymer sequence 时严格 source-auth exact-WT 映射、较长 polymer 中唯一 authoritative segment、`not_evaluable` summary 状态。
- Flex ddG修复后真实precheck和8任务均成功；结果commit为`0bfcc88`，gate为pass。8/8任务身份与manifest一致、结构安全全pass，时间投影独立复算一致，PNG已人工检查无遮挡；科学review位于`flex_ddg_pilot_scientific_review.json`。

## Required Next Steps

1. 远程pull后运行`bash scripts/candidate_design/submit_antifold_validation.sh`；该单GPU作业依次评分实验VHH-only、实验complex-context和AF3 VHH，并在项目环境生成8个核心单突的三视图`ΔlogP`、perplexity、方向一致性、gate和图。该阶段不采样新序列、不设AntiFold淘汰阈值。
2. 同步真实结果并复核后建立统一单突评价表：为已有亲和力候选补齐AntiFold、NetSolP、NanoMelt、TNP和化学风险；同时在通过硬约束的FR及非关键CDR位置提出一批有结构/序列依据的稳定性/可开发性单突，并用PyRosetta接触、表位和亲和力非劣门复核。
3. 从双向过门的单突模块生成亲和力×亲和力、性质×性质和亲和力×性质双突；先对完整组合运行快速统一评分，再只对入围组合运行20构象Flex ddG，必要时以AF3复核少量终选构象。
4. 按硬约束加Pareto关系形成最终30条：分别标为亲和力优先、稳定性/可开发性优先和平衡组合，保留必要WT/单突对照及全部原始分项，不使用yield综合分，也不把计算类别写成实验改善结论。
