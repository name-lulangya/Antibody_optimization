# Codex 项目交接

Last updated: 2026-08-11 15:40:32

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是 replace-in-place 的当前状态快照，不是追加式历史。

## Project Goal

优化 NK2R 纳米抗体 Nb252 的亲和力、稳定性和表达相关性质，同时保留结构、序列、实验和计算证据的可追溯性。

## Environment and Git

- 本地环境：`D:\miniconda\envs\ab_optim`，Python 3.11.15、ANARCII 2.0.8、Gemmi 0.7.5。
- 远程项目环境：`/data/software/env/luly25/ab_optim`；计划检出父目录 `/homes/Tianlab/luly25/`，登录别名尚未建立。
- PyRosetta：`/data/software/env/luly25/multi_ligand`，Python 3.10.20，PyRosetta 2026.03，Rosetta commit `5e498f1409c68ade56c8ce5842bf79e1b02e8db4`。
- nanoBERT：`/data/software/env/luly25/vhh-lm`，`NaturalAntibody/nanoBERT` revision `edc8182ad89a827f8737fa572c6b5fac6197e6b0`，使用已记录离线缓存。
- Git：`main`，远程 `git@github.com:name-lulangya/Antibody_optimization.git`；pilot v1结果及科学复核均已同步，pilot V2实现随当前`main`提交同步。
- 阶段2的本地阶段0、远程WT安全导入、WT评分校准v1/v2和12候选PyRosetta pilot均已运行。尚未运行456候选全量评分、AF3候选复核、nanoBERT/AntiFold或表达模型训练。

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
- 12候选×3重复pilot的计算执行完整：3个共享WT、36个突变体评估、12个候选汇总，耗时1157.355071秒；所有突变体映射、断点、二硫键和有限数值均通过，WT的`dG_separated`及跨界面能也全部为负。生成machine gate为pass，但科学release由`affinity_pilot_scientific_review.json`覆盖为`held_pending_candidate_structure_gate_revision`。
- hold原因：再次运行局部准备后的配对WT本身只保持33/37个实验NK2R表位残基，即0.8918918918918919，低于候选硬编码的0.90绝对阈值；9个候选仅与WT相同却被标为blocked。因此生成的“3 pass/9 blocked”和`full_affinity_scan_release=pass`不可用于启动456扫描。候选结构门必须先改为显式相对配对WT，同时继续核对原实验表位集合，再以新目录重跑pilot。
- pilot V2代码已实现但尚未远程运行：扫描阶段不按能量、接触保持率或RMSD筛选；逐重复保存WT和突变体精确接触集合，同时计算相对原实验集合及相对同seed配对WT的保持率。V2 gate只验证候选运行安全和WT对照有效性，并显式声明`score_all_declared_candidates_then_filter_once`。v1制品保持历史状态、不覆盖。
- `finalize_input_baseline.py` 已用真实 structure mapping/interface 重建 canonical 128 位点图和 `stage1_gate.json`。程序 `status=pass` 不代表科学 release gate 自动通过。

## Stage-2 Phase 0 Result

- 正式制品位于`docs/result_artifacts/candidate_design/stage0_contract_20260810/`，run summary位于`docs/run_summaries/candidate_design/stage0_contract_20260810.json`。
- 128个位点中：实验缺失坐标13位、实验界面24位、硬冻结6位、首轮亲和力放行24位。硬冻结为reported positions 125–128的`SSGS`和实验SG–SG几何支持的Cys 22/Cys 95；界面位点是“谨慎但可突变”。
- 首轮亲和力以未批量补全的实验complex为主；缺失位点不进入首轮亲和力扫描。只有终选涉及缺失CDR1、或原实验结构与完整VHH预测轨道冲突时，才触发定点补全敏感性分析。

## Current Optimization Plan

1. **远程WT结构门（已完成，实测14.123538秒）**：PyRosetta 2026.03直接导入未补全实验complex，396个polymer残基和4个预期断点均安全，映射与二硫键通过；未做relax、补全、突变或候选评分。权威结果位于`structure_preparation/pyrosetta_wt_import_20260810/pyrosetta_wt_import_gate.json`。
2. **WT评分协议校准（已完成，实测465.154733秒）**：schema v2 gate通过并唯一选择受约束局部最小化；代表WT、逐位置接触变化、重复表、逐残基能量和QC图均已同步。该步骤只释放成对相对候选评分，不生成突变或给出实验亲和力结论。
3. **总体多目标架构（保留，但调整执行顺序）**：最终仍采用亲和力与稳定性/表达两条主轨从同一authoritative WT独立发现、风险作为横向监测层、汇合后少量组合、实验后再串行迭代。当前先完整完成证据最确定的亲和力轨道；稳定性/表达候选、47条yield小模型、nanoBERT PLL相关性和风险修复候选均暂停，待亲和力单点结果形成后另行讨论，不让未验证模型阻断亲和力进度。
4. **亲和力单点空间（本地已完成，实测约2秒）**：`affinity_single_mutants_20260811/`已从阶段0合同、实验24位界面、可逆编号和v2评分gate生成每个位点19种非WT替换，共456个单点；每条均保留完整128-aa序列、reported/IMGT/source-auth编号、实验接触和prepared WT敏感标记。12个pilot候选覆盖FR/CDR、保守/跨类替换及E46/D101/I103，未按prepared接触或主观化学偏好预删候选。
5. **成对PyRosetta亲和力pilot V2（实现完成、待远程运行）**：沿用12候选×3重复和共享WT，但计算阶段不筛选候选；同时输出配对WT相对接触保持、原实验集合保持、能量差、RMSD和运行安全。预计约20分钟，写入全新V2目录。
6. **456单点全量扫描（V2通过后实现）**：所有456候选先按同一协议完成3重复评分，不在分片、任务或扫描中按分数/结构指标提前淘汰。结果合并且完整性核验后，再一次性采用统一规则筛选，避免不同批次或局部门改变候选空间。
7. **亲和力复核与面板（约1–3天，不含实验）**：统一筛选综合重复稳定性、相对配对WT及原实验表位/接触保持、界面Cα RMSD、排斥和局部几何，不以单一能量值决定。E46/D101/I103不得因单一prepared构象接触丢失而直接淘汰。对小规模入围者再运行更严格的flex-ddG或等价多构象复核，并用AF3检查最终候选是否维持总体VHH构象；输出WT、机制多样的亲和力单点实验面板。此阶段不组合突变。

## Verification

- 阶段0专项：`3 passed`；覆盖真实128位合同、过期哈希拒绝、CSV BOM/LF、拒绝覆盖、固定时间戳双跑六制品逐字节一致。
- 阶段0图及亲和力候选空间图均已人工检查，坐标轴、图例和说明无遮挡。当前全套验收为`126 passed, 1 skipped, 4 subtests passed`；唯一skip仍是Windows真实symlink权限测试。
- `pip check`、`python -m compileall -q src scripts tests`、`git diff --check` 均通过。
- v2结果一致性复核：schema 2 gate为pass、16行重复数据和67行接触状态与gate一致、代表PDB含396个polymer残基；评分校准专项测试`8 passed`。
- pilot结果一致性复核：3行WT、36行候选重复、12行候选汇总与machine gate计数一致；0个突变体runtime failure。WT表位保持0.8918918918918919低于候选0.90绝对门，生成候选状态不接受为科学分类。
- 修复覆盖：ChimeraX sandbox 入口、实际模型名、无 polymer sequence 时严格 source-auth exact-WT 映射、较长 polymer 中唯一 authoritative segment、`not_evaluable` summary 状态。

## Required Next Steps

1. 在远程pull当前pilot V2实现；通过`sbatch scripts/candidate_design/submit_affinity_scoring_pilot.slurm`写入全新`affinity_pyrosetta_pilot_v2_20260811/`目录，不覆盖v1。
2. 同步V2结果，核验3个WT、36个突变体评估、12个未筛选摘要、精确接触集合、两套接触保持率及machine gate；此时只判断计算路线是否可用于全量扫描，不进行候选优劣筛选。
3. V2运行门通过后实现456候选显式分片与Slurm array；全部候选完成并合并后，再实现一次统一筛选入口。每任务控制在5小时内，不增加复杂checkpoint。
4. 表达/稳定性建模、nanoBERT/AntiFold/AbMPNN和风险修复候选继续暂停，亲和力全量筛选结果完成后重新讨论。
