# Codex 项目交接

Last updated: 2026-08-13 16:30:00

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是 replace-in-place 的当前状态快照，不是追加式历史。

## Project Goal

优化 NK2R 纳米抗体 Nb252 的亲和力、稳定性和表达相关性质，同时保留结构、序列、实验和计算证据的可追溯性。

## Environment and Git

- 本地环境：`D:\miniconda\envs\ab_optim`，Python 3.11.15、ANARCII 2.0.8、Gemmi 0.7.5。
- 远程项目环境：`/data/software/env/luly25/ab_optim`；计划检出父目录 `/homes/Tianlab/luly25/`，登录别名尚未建立。
- PyRosetta：`/data/software/env/luly25/multi_ligand`，Python 3.10.20，PyRosetta 2026.03，Rosetta commit `5e498f1409c68ade56c8ce5842bf79e1b02e8db4`。
- nanoBERT：`/data/software/env/luly25/vhh-lm`，`NaturalAntibody/nanoBERT` revision `edc8182ad89a827f8737fa572c6b5fac6197e6b0`，使用已记录离线缓存。
- Git：`main`，远程 `git@github.com:name-lulangya/Antibody_optimization.git`；亲和力ensemble核心与稳定性/表达合同提交`a243399`已同步。
- 阶段2的本地阶段0、远程WT安全导入、WT评分校准、456候选全量扫描、post-scan分层、50候选×20样本Flex ddG复核和nanoBERT—yield验证均已完成。尚未运行亲和力组合、AF3候选复核、AntiFold或表达模型训练。

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
- pilot v1已作为历史诊断保留：其3 pass/9 blocked来自候选绝对0.90表位门低于配对WT自身33/37表现，不能作为候选分类。
- pilot V2已远程完成：1146.969247秒，3个共享WT、36个突变体评估、12个候选摘要全部通过运行安全，所有摘要均为`not_applied_scan_stage`，没有执行候选筛选。接触集合、计数和相对同seed配对WT保持率已逐行重算一致；NK2R表位配对WT保持率全部为1.0，VHH配对WT保持率范围0.958333–1.0。V1/V2共同数值字段逐字一致，证明V2修正的是门和记录语义而非评分轨迹。
- `affinity_pilot_v2_scientific_review.json`将路线发布为`released_for_full_456_scan_implementation`。全量合同固定为先计算全部456候选并核验合并完整性，再统一筛选；pilot能量只作未筛选诊断，不产生实验候选推荐。
- 456×3全量扫描已完成并通过merge gate：12片累计计算42645.718076秒，456个候选、1368个唯一重复键、3个去重WT均完整，0个runtime failure，全部摘要为`not_applied_scan_stage`。接触集合/保持率和`mutant-WT`能量差已逐行复算一致；`affinity_full_scan_scientific_review.json`发布`ready_for_post_scan_filter_implementation`，尚未选择或淘汰任何候选。
- 统一post-scan分层已在本地完成：456条全部通过硬有效性门并互斥分为Tier 1/2/3/4/5=`18/30/39/82/287`；Tier 1/2合计48条形成严格复核池。分层使用三重复双能量方向、配对WT接触保持和`Δfa_rep`，风险标签与层内Pareto保持独立；`candidate_selection_performed=false`，尚未形成最终实验面板。
- 亲和力ensemble核心筛选已完成：50条完整证据中8个单点通过“两项能量均至少18/20样本为负且中位数均为负”的唯一核心门，分别为`R45C/R45V/D101W/I103W/E105F/E105L/N107A/S114M`，覆盖6个位置；R45和E105的同位点替换互斥。风险、接触保持、prepared敏感性与Pareto层仍独立保留；未生成双突或最终实验序列。
- 稳定性/表达WT发现合同已完成：128位中72个有实验坐标的非界面framework正常放行，9个实验缺失framework仅在完整VHH AF3证据下谨慎放行；24个实验界面、17个CDR/terminal flank、Cys22/Cys95和末端SSGS共47位冻结。合同不生成突变，也不声称性质改善。
- nanoBERT—reported yield验证gate为`pass`，但证据等级为`no_supported_use`、release为`nanobert_not_supported_for_candidate_use`。主指标完整reported序列mean PLL在31条LTT/WCC上的来源分层Spearman为-0.2295，长度调整partial Spearman为-0.2020，LTT/WCC内部为-0.2413/0.0602；95% bootstrap CI为[-0.5982, 0.2219]，来源内置换p=0.2799，加入nanoBERT后的LOOCV Spearman增益为-0.0507。16条LLJ有序/删失Kendall tau-b为-0.2308。不得用nanoBERT为Nb252稳定性、表达或候选序列打分/排序。
- `finalize_input_baseline.py` 已用真实 structure mapping/interface 重建 canonical 128 位点图和 `stage1_gate.json`。程序 `status=pass` 不代表科学 release gate 自动通过。

## Stage-2 Phase 0 Result

- 正式制品位于`docs/result_artifacts/candidate_design/stage0_contract_20260810/`，run summary位于`docs/run_summaries/candidate_design/stage0_contract_20260810.json`。
- 128个位点中：实验缺失坐标13位、实验界面24位、硬冻结6位、首轮亲和力放行24位。硬冻结为reported positions 125–128的`SSGS`和实验SG–SG几何支持的Cys 22/Cys 95；界面位点是“谨慎但可突变”。
- 首轮亲和力以未批量补全的实验complex为主；缺失位点不进入首轮亲和力扫描。只有终选涉及缺失CDR1、或原实验结构与完整VHH预测轨道冲突时，才触发定点补全敏感性分析。

## Current Optimization Plan

1. **远程WT结构门（已完成，实测14.123538秒）**：PyRosetta 2026.03直接导入未补全实验complex，396个polymer残基和4个预期断点均安全，映射与二硫键通过；未做relax、补全、突变或候选评分。权威结果位于`structure_preparation/pyrosetta_wt_import_20260810/pyrosetta_wt_import_gate.json`。
2. **WT评分协议校准（已完成，实测465.154733秒）**：schema v2 gate通过并唯一选择受约束局部最小化；代表WT、逐位置接触变化、重复表、逐残基能量和QC图均已同步。该步骤只释放成对相对候选评分，不生成突变或给出实验亲和力结论。
3. **总体多目标架构（已确认）**：采用“并行发现、条件化组合”。亲和力与稳定性/表达模块先从同一authoritative WT独立发现，以保留因果归因和统一基线；随后不做机械拼接，而是选择多个亲和力核心作为新序列上下文，固定其亲和力突变、实验表位/界面关键约束、Cys22/Cys95和末端SSGS，再重新设计和筛选稳定性/表达模块。最终完整组合必须重新计算亲和力、结构相容性、序列相容性和开发性指标，不能相加单模块分数，也不锁定唯一亲和力母本。
4. **亲和力单点空间（本地已完成，实测约2秒）**：`affinity_single_mutants_20260811/`已从阶段0合同、实验24位界面、可逆编号和v2评分gate生成每个位点19种非WT替换，共456个单点；每条均保留完整128-aa序列、reported/IMGT/source-auth编号、实验接触和prepared WT敏感标记。12个pilot候选覆盖FR/CDR、保守/跨类替换及E46/D101/I103，未按prepared接触或主观化学偏好预删候选。
5. **成对PyRosetta亲和力pilot V2（已完成，实测19.116154分钟）**：12候选×3重复和共享WT运行全部有效，配对WT相对接触指标复算一致，扫描阶段未筛选候选；路线已释放456全量扫描实现。
6. **456单点全量扫描（已完成）**：12片累计11.846小时计算量，完整得到456候选×3重复且没有扫描中筛选。未筛选景观中113个候选`ΔdG<0`、143个跨界面能变化<0、87个两者均<0，但这些只是描述性计数；82个候选两指标方向不一致，说明后续不能使用单一能量列筛选。
7. **统一筛选（已完成，本地实测约1.5秒）**：完整456候选按统一规则分为`18/30/39/82/287`，48条Tier 1/2进入严格复核池；未删除候选、未形成最终实验推荐。
8. **Flex ddG生产参数计时pilot（已完成）**：4代表候选×2独立样本共8任务全部pass，累计1.666723 job-hours；单任务中位717.980720秒、P90 874.753734秒，backrub占累计任务时间94.7212%，峰值内存975.65–1000.62 MiB。全部WT/突变映射、断点和二硫键检查通过。两重复只验证生产参数、耗时和协议可行性；代表候选存在样本间能量方向翻转，因此不得用pilot分数排名。
9. **亲和力严格复核与核心选择（已完成）**：50候选×20样本=1000任务全部pass；随后按唯一双指标18/20方向门从完整50条中选择8个核心单点、覆盖6个位置。该选择不使用加权综合分，仍保留全部50条证据和每个核心的`fa_rep`、接触及prepared风险；当前只释放核心模块使用，不代表8条均等安全，也未生成亲和力双突。
10. **模块与组合层级（已确认）**：单突是证据模块而非最终设计上限。亲和力先从20样本结果形成约6–8个机制多样单突核心，再主要探索双亲和力组合；稳定性/表达模块允许1–3个协调framework突变，必要时保留极少数4突变探索项；跨轨道完整候选以总计2–4个突变为主，极少超过4个。最终30条以重新评价的多目标组合为主体，同时保留少量WT/单模块对照以解释实验结果。
11. **工具分工（当前选择）**：AntiFold作为结构条件稳定性/表达模块的主要生成器，暂不并行部署同类AbMPNN；nanoBERT已因47条yield验证不支持而退出本项目候选排序/过滤路线；PyRosetta用于完整亲和力组合和跨轨道组合的界面/结构复核；理化、聚集和化学风险作为所有完整序列的横向过滤；AF3仅用于少量终选完整组合的构象复核。
12. **稳定性/表达设计合同（已完成）**：WT发现空间固定为81个非界面framework位点，其中9个缺失坐标位点要求AF3证据；47位冻结。主要模块允许1–3突变、极少数探索模块最多4突变；此处只有范围和证据合同，AntiFold尚未运行。
13. **nanoBERT—yield验证（已完成）**：47条真实序列均完成固定revision逐残基single-mask评分；主指标未通过方向、跨来源稳定性、不确定性或LOOCV门，证据等级为`no_supported_use`。全样本合并ρ=0.1887与来源分层ρ=-0.2295方向相反，显示来源混杂；nanoBERT从稳定性/表达候选评价路线中移除，不作为排序或过滤信号。理化特征只产生未校正的探索性关联，当前也不得作为表达优化方向或训练模型依据。

## Verification

- 阶段0专项：`3 passed`；覆盖真实128位合同、过期哈希拒绝、CSV BOM/LF、拒绝覆盖、固定时间戳双跑六制品逐字节一致。
- 阶段0图、亲和力候选空间图、post-scan四面板图、ensemble核心图和稳定性/表达合同图均已人工检查，坐标轴、图例和说明无遮挡。当前全套验收为`166 passed, 1 skipped, 4 subtests passed`；唯一skip仍是Windows真实symlink权限测试。
- `pip check`、`python -m compileall -q src scripts tests`、`git diff --check` 均通过。
- v2结果一致性复核：schema 2 gate为pass、16行重复数据和67行接触状态与gate一致、代表PDB含396个polymer残基；评分校准专项测试`8 passed`。
- pilot V2一致性复核：3行WT、36行候选重复、12行候选汇总与machine gate一致；0个runtime failure，36/36重复和12/12摘要为pass，全部摘要未筛选。接触集合计数及配对WT保持率复算一致；V1/V2共同数值字段逐字一致。
- 全量扫描plan为pass：456候选、1368预期突变体评估、12片×38、每片两个完整位置、最大并发4且`candidate_filtering_applied=false`。实现测试覆盖真实分片、完整合并、缺片/提前筛选阻断、最终PNG/SVG和Slurm依赖合同。
- 全量结果复核：456摘要、1368重复、3个WT与manifest/merge gate一致，所有candidate×replicate×seed唯一，接触和能量差复算通过。最终PNG已人工检查；右侧重复Y轴标签已从绘图代码和已有plot data修正，不改评分数据。
- 修复覆盖：ChimeraX sandbox 入口、实际模型名、无 polymer sequence 时严格 source-auth exact-WT 映射、较长 polymer 中唯一 authoritative segment、`not_evaluable` summary 状态。
- Flex ddG修复后真实precheck和8任务均成功；结果commit为`0bfcc88`，gate为pass。8/8任务身份与manifest一致、结构安全全pass，时间投影独立复算一致，PNG已人工检查无遮挡；科学review位于`flex_ddg_pilot_scientific_review.json`。

## Required Next Steps

1. 建立AntiFold最小运行路线，先在WT完整VHH结构上验证81位合同、冻结位点和1–3突变模块输出；不同时部署同类型AbMPNN，也不再接入nanoBERT评分。
2. 依据8个核心模块的风险和同位点互斥关系，设计少量机制互补亲和力双突并用完整组合重新评分；不是把单点REU相加。
3. 对AntiFold模块使用结构相容性、保守位点、聚集/化学风险及基础理化边界筛选；现有47条理化关联仅作探索背景，不能作为表达量训练标签或单向优化目标。
4. 将稳定性/表达模块分别放入多个亲和力核心背景重新评价，生成受控的2–4突变完整组合；只对有限组合运行多突变PyRosetta/Flex ddG，并对少量终选使用AF3。
5. 在完整组合证据形成后确定最终30条及必要对照；不把REU、AntiFold输出或未验证理化趋势解释为实验性质。
