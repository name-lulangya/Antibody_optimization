# Codex 项目交接

Last updated: 2026-08-16 23:48:51

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
- Git：`main`，远程 `git@github.com:name-lulangya/Antibody_optimization.git`；V2.1真实结果提交`ba654be`已拉取，开始本轮复核时`HEAD=origin/main`。
- 统一单突AntiFold复用已完成且gate为`pass/ready_for_unified_property_scoring`：2318条均有记录，1962条本轮候选三视图全部可评价，247条实验缺失坐标候选仅AF3可评价；模型没有重新运行。1962条中66条三视图均为正、1733条均为负、163条方向不一致。界面432条中30条三视图均正，性质发现1530条中36条三视图均正；结合既有20构象严格亲和力门后仍只有R45V和E105L为兼容性方向一致的亲和力核心，没有新增核心。
- 统一性质评分已远程完成并复核：NetSolP与NanoMelt均为1963/1963 pass，1962条候选无丢失；432条界面轨道分为Pareto 1/2/background=`65/154/213`，1530条性质发现轨道为`49/111/1370`。全体候选在U、S、预测表观Tm和实验复合物AntiFold四个方向同时为正者仅17条；性质发现轨道仅5条，其中`F30D`新增异构化motif，余4条无当前化学风险标记。Pareto表示轨道内非支配关系，不等于所有性质改善或最终入选；无yield预测、加权总分、TNP、双突或最终候选选择。
- 统一TNP复核已远程完成并复核：WT+95候选为96/96 pass，49条性质来源和46条亲和力来源身份完整；74条官方flag不变、21条改善、0条恶化、0条新增red。21条改善全部是PSH由WT amber变为green，其他5类flag对全部候选均为green且不变。本次WT PSH为130.9926，既有同序列正式运行WT为136.6287，两次均为amber但相差5.6361；因此单次连续PSH及临界flag变化只作辅助风险证据，不单独排序或入选。release为`ready_for_multitool_shortlist_with_tnp_as_supporting_risk_only`。

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

- `unified_tnp_candidate_review=pass`：固定96条输入全部通过；没有候选因TNP新增硬风险，TNP只作为多工具短名单中的辅助developability风险证据。49条性质来源中9条同时满足“至少一个U/S/Tm非微小改善、无这些指标的明显恶化、实验complex AntiFold为正、无化学风险”，其中`Q1D/F30K`同时出现PSH flag改善；46条亲和力来源中该交集有5条，仅`V108Y`同时出现PSH flag改善。这些是证据交集而非最终候选选择。

- `unified_single_mutant_property_scoring=pass`：NetSolP/NanoMelt覆盖完整，科学复核release为`ready_for_preliminary_property_pool_definition`。既有8个亲和力核心中7个进入本轮评分，`R45C`因新增未配对Cys继续阻断；只有`R45V`在四个性质方向均为正，其他核心保留亲和力证据但存在至少一个性质或兼容性权衡。该gate只允许定义小型复核池，不允许把Pareto 1或任一预测器当成单独淘汰/入选门。

- `unified_single_mutant_antifold_landscape=pass`：2318条证据身份完整，1962条本轮候选全部具有实验VHH-only、实验complex-context和AF3三视图；release为`ready_for_unified_property_scoring`。AntiFold不得单独筛选，缺失坐标247条仍不释放。

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

- 性质候选专用PyRosetta的pilot与完整扫描均已完成：pilot 6条×3重复、实测1081.568448秒；完整扫描30条×3重复、10个位点各3个位置特异WT、实测3516.793389秒。30条/90个配对结果均完整且扫描中未筛选，共同6条候选在pilot/full的汇总数值最大绝对差为0。独立复算得到方向一致有利/混合/方向一致不利=`9/12/9`；9条为`Q1D/Q3Y/Q5A/F30A/F30P/F30Q/S62D/K86S/K86T`，其中只有`Q1D`在3/3重复中两项能量同时为负。该方向分类不设REU效应阈值，也不代表实验亲和力或最终入选。
- 30条中26条在三个配对重复中完整保持既定VHH/NK2R接触；`Q5V`丢失NK2R auth 106，`K86A/K86S/K86T`在部分或全部重复中丢失NK2R auth 37。VHH接触最低保持率为1.0、NK2R为0.972973，最大界面Cα RMSD为0.045763 Å。方向一致有利且实验complex AntiFold `ΔlogP>0`的交集只有`Q1D/F30A/F30P`；这是跨工具方向交集，不是最终实验序列选择。科学review release为`ready_for_property_module_shortlist_review`。
- 统一单突安全门已在50条Flex ddG亲和力候选和30条性质候选上完成，release为`ready_for_targeted_structure_review_not_combination_generation`。80条分为`combination_ready/single_mutant_test_only/targeted_alternative_review/blocked_pending_structure/not_prioritized/blocked=1/19/1/10/45/4`；唯一直接组合就绪为`Q1D`，较温和亲和力备选为`R45T`，现有8条亲和力能量核心中7条为单突测试、`R45C`因新增游离Cys硬阻断。A23/F30的10条候选因邻接实验缺失positions 24–29而暂缓，其中`F30P`另有Pro主链风险。该门使用链C单独SASA和局部结构/序列风险；所有阈值仅是项目保守分诊，不是实验表达、聚集或亲和力界限。
- 非冗余定点结构复核已远程完成：9条非Pro A23/F30候选×3重复共27次AF3完整VHH局部repack全部通过运行完整性门，实测20.075577秒且评分阶段未筛选。V2释放为`ready_for_combination_module_review`；80条总池更新为`combination_ready/single_mutant_test_only/targeted_alternative_review/not_prioritized/do_not_advance=5/24/1/43/7`。组合模块审阅池为`Q1D/A23S/F30A/F30S/F30T`；其中本轮新增4条均满足三重复有效，且AF3 VHH总分与突变位点8 Å局部加权`fa_rep`的中位变化均不大于0。`A23Q/A23R`保留强负向AntiFold/既有亲和力不利风险，`F30K/F30Q/F30R`因局部`fa_rep`中位上升，均仅为`single_mutant_test_only`。7条不可补偿Cys/Pro风险继续为`do_not_advance`；局部repack不清除序列或developability软风险，也不代表实测稳定性、表达或亲和力。
- 性质单突空间已按综合风险进一步缩小，gate为`pass/ready_for_small_combination_contract`。V2的30条活跃单突降为14条：亲和力8条全部保留，性质22条降为6条`Q1D/A23S/F30A/F30S/F30T/S55G`。16条性质候选因受体接触变化4条、AF3局部非劣门失败3条、亲和力方向不利2条、强负向AntiFold伴暴露疏水3条或单独强负向AntiFold 4条而降级；它们和既有43条未优先、7条硬排除仍保留在80行审计表，不从历史数据删除。该缩减只复用V2证据，没有运行新模型或生成组合。

## Stage-2 Phase 0 Result

- 正式制品位于`docs/result_artifacts/candidate_design/stage0_contract_20260810/`，run summary位于`docs/run_summaries/candidate_design/stage0_contract_20260810.json`。
- 128个位点中：实验缺失坐标13位、实验界面24位、硬冻结6位、首轮亲和力放行24位。硬冻结为reported positions 125–128的`SSGS`和实验SG–SG几何支持的Cys 22/Cys 95；界面位点是“谨慎但可突变”。
- 首轮亲和力以未批量补全的实验complex为主；缺失位点不进入首轮亲和力扫描。只有终选涉及缺失CDR1、或原实验结构与完整VHH预测轨道冲突时，才触发定点补全敏感性分析。

## Current Optimization Plan

1. **远程WT结构门（已完成，实测14.123538秒）**：PyRosetta 2026.03直接导入未补全实验complex，396个polymer残基和4个预期断点均安全，映射与二硫键通过；未做relax、补全、突变或候选评分。权威结果位于`structure_preparation/pyrosetta_wt_import_20260810/pyrosetta_wt_import_gate.json`。
2. **WT评分协议校准（已完成，实测465.154733秒）**：schema v2 gate通过并唯一选择受约束局部最小化；代表WT、逐位置接触变化、重复表、逐残基能量和QC图均已同步。该步骤只释放成对相对候选评分，不生成突变或给出实验亲和力结论。
3. **总体多目标架构（当前有效）**：不再按CDR/FR严格划分亲和力与稳定性/表达相关优化区域，也不再默认冻结全部framework。候选分为亲和力驱动和稳定性/可开发性驱动两个来源，但每条完整序列从单突阶段起都同步评价亲和力、AntiFold结构相容性、NetSolP U/S、NanoMelt相对WT预测Tm、TNP developability和化学风险。亲和力驱动候选必须满足性质非严重恶化，性质驱动候选必须满足亲和力、表位和构象非劣；完整组合逐条重算，不机械相加分数，也不把预测值解释为实测产量、Tm或亲和力。
4. **亲和力单点空间（本地已完成，实测约2秒）**：`affinity_single_mutants_20260811/`已从阶段0合同、实验24位界面、可逆编号和v2评分gate生成每个位点19种非WT替换，共456个单点；每条均保留完整128-aa序列、reported/IMGT/source-auth编号、实验接触和prepared WT敏感标记。12个pilot候选覆盖FR/CDR、保守/跨类替换及E46/D101/I103，未按prepared接触或主观化学偏好预删候选。
5. **456单点全量扫描（已完成）**：12片累计11.846小时计算量，完整得到456候选×3重复且没有扫描中筛选。未筛选景观中113个候选`ΔdG<0`、143个跨界面能变化<0、87个两者均<0，但这些只是描述性计数；82个候选两指标方向不一致，说明后续不能使用单一能量列筛选。
6. **统一筛选（已完成，本地实测约1.5秒）**：完整456候选按统一规则分为`18/30/39/82/287`，48条Tier 1/2进入严格复核池；未删除候选、未形成最终实验推荐。
7. **Flex ddG生产参数计时pilot（已完成）**：4代表候选×2独立样本共8任务全部pass，累计1.666723 job-hours；单任务中位717.980720秒、P90 874.753734秒，backrub占累计任务时间94.7212%，峰值内存975.65–1000.62 MiB。全部WT/突变映射、断点和二硫键检查通过。两重复只验证生产参数、耗时和协议可行性；代表候选存在样本间能量方向翻转，因此不得用pilot分数排名。
8. **亲和力严格复核与证据核心（已完成）**：50候选×20样本=1000任务全部pass；按双指标18/20方向门选择的8条现在明确称为“亲和力证据核心”，不是组合许可。统一安全门将其中7条保留为`single_mutant_test_only`，`R45C`因新增游离Cys硬阻断；没有亲和力证据核心直接成为`combination_ready`。
9. **完整双突合同与扫描（已完成）**：86条均完成四工具评分；V2.1真实分析用258条突变重复和135条WT逐行复算配对接触，3.698408秒完成并以`pass/ready_for_scientific_shortlist_definition`释放。86条全部通过配对WT VHH/NK2R保持率与界面RMSD主结构门；最低保持率分别为0.958333/0.972973，最大界面Cα RMSD为0.047716 Å。绝对实验4 Å参考敏感/不敏感为57/29，只作准备敏感性标记。证据类为平衡支持/亲和力支持且性质非不利/性质支持且亲和力非不利/权衡或不清晰=`28/12/2/44`，即42条进入科学短名单定义范围。40条存在任一配对接触变化、46条三重复均不变；变化侧为VHH 18条、NK2R 37条（15条两侧均变），但均未跌破主门。42条支持候选中31条接触不变、11条有局部变化，且无新增当前化学liability或TNP flag回退。AntiFold仍只是两个固定骨架单点分数的数学加和，未评价双突上位性。
10. **工具分工（当前选择）**：PyRosetta/Flex ddG提供相对亲和力及接触/表位非劣证据；AntiFold评价候选对既定Nb252结构的条件相容性，并可为允许区域的性质候选提供序列先验。NetSolP U/S只作表达/溶解性辅助相容性信号，TNP只作developability风险工具，NanoMelt只提供相对WT的VHH热稳定性预测；三者已验证不能单独或联合承担yield排序。nanoBERT已退出，不训练现有47条数据的Nb252局部产量排序器，也不部署同类型重复工具投票。AF3只复核少量终选完整组合构象。
11. **framework与硬约束政策（当前有效）**：FR可以提出稳定性/可开发性候选，但不是无约束全扫描。reported positions 125–128的`SSGS`、Cys22–Cys95实验二硫键和实验缺失坐标的首轮结构不可评估状态继续保持；额外未配对Cys、VHH hallmark、CDR支撑网络、核心疏水包装、实验界面和表位变化均作为明确风险或阻断条件。现有81位合同只提供历史范围，47条序列只提供同框架频率和背景先验，不得据此把某个FR替换称为表达因果突变；所有FR候选必须在完整序列上重新评价亲和力和结构安全。

## Verification

- 当前全套为`223 passed, 1 skipped, 4 subtests passed`，唯一skip为Windows真实symlink权限测试；V2.1专项7项通过。Python编译、Slurm Bash语法和`git diff --check`通过。真实V2.1已生成schema 3的86行联合表、258行接触审计、gate、PNG/SVG和run summary；图已人工检查，三面板、图例、黑色接触变化圈和图注均清晰无遮挡。

## Required Next Steps

1. 从V2.1释放的42条支持双突、保留的安全单突及WT对照中建立亲和力驱动、性质驱动和组合驱动三类非冗余候选池；将11条支持但存在配对接触变化的双突单独标记，不把局部4 Å集合变化自动当作淘汰。
2. 按证据强度、性质非劣、化学风险、接触变化、位点/机制多样性和三类配额形成30条初选；不得把42条支持类别或31条接触不变候选机械截取为最终名单。
3. 仅当近邻候选选择或机制解释确实需要时，对少量优先组合补做统一四态PyRosetta分解；不对86条普遍重跑。随后仅让综合非劣且具有明确收益的少量组合进入AF3完整构象复核，再冻结30条待实验序列。
