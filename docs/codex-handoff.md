# Codex 项目交接

Last updated: 2026-08-10 10:30:00

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是 replace-in-place 的当前状态快照，不是追加式历史。

## Project Goal

优化 NK2R 纳米抗体 Nb252 的亲和力、稳定性和表达相关性质，同时保留结构、序列、实验和计算证据的可追溯性。

## Environment and Git

- 本地环境：`D:\miniconda\envs\ab_optim`，Python 3.11.15、ANARCII 2.0.8、Gemmi 0.7.5。
- 远程项目环境：`/data/software/env/luly25/ab_optim`；计划检出父目录 `/homes/Tianlab/luly25/`，登录别名尚未建立。
- PyRosetta：`/data/software/env/luly25/multi_ligand`，Python 3.10.20，PyRosetta 2026.03，Rosetta commit `5e498f1409c68ade56c8ce5842bf79e1b02e8db4`。
- nanoBERT：`/data/software/env/luly25/vhh-lm`，`NaturalAntibody/nanoBERT` revision `edc8182ad89a827f8737fa572c6b5fac6197e6b0`，使用已记录离线缓存。
- Git：`main`，远程 `git@github.com:name-lulangya/Antibody_optimization.git`；已同步基线为 `facfe848`。当前有未提交的真实 CXS 导出、review、结构/界面结果、修复、测试和文档。
- 本阶段未运行 AF3、PyRosetta、nanoBERT、候选生成或表达模型训练。

## Frozen Inputs

- `Nb252-optimization.cxs`：537670 bytes，SHA-256 `1bc636c28f66ae60edc658d2e1c4aad0b07f4141ca5411c78662aa19da793c4d`，原文件未被覆盖。
- `nb序列及产量（1L）.docx`：14172 bytes，SHA-256 `a6e4022f0978fbd70a0e04dc78f479140ab6f55caaa90b467fb77a62eb5db5d1`。
- 47 条 reported 序列的 provisional IMGT 审核：46 pass、1 failed。`WCC__4-28` 为 ANARCII `Score less than cut off.`；Nb252 编号 126 aa，末端 `GS` 未编号，边界仍是 provisional。

## Structure and Interface Baseline

- `data/structures/cxs_exports/` 是 ChimeraX 1.12 的真实导出：3 个目标原子模型、5 个 mmCIF、颜色清单、manifest 和 run summary 均已验证；5 个 mmCIF 可由 Gemmi 0.7.5 读回。
- CXS 中的精确模型名为 `NK2R-252.pdb`、`NK2R-NKA` 和 `fold_2r_252_nomg_model_0.cif`。导出入口已修复 ChimeraX 非 `__main__` sandbox 调用问题。
- 用户已在 ChimeraX 1.12 中确认链角色：实验 `C=Nb252`、`R=NK2R`；NKA 模型 `L=NKA`、`R=NK2R`；AF3 `A=Nb252`。
- 用户确认实验 chain C 的 24 个精确 `[255,165,0,255]`、`atom+ribbon` 位点就是合作者橙色区域。机器可读记录：`docs/result_artifacts/input_baseline/structure_review_20260810/baseline_review.json`。
- `structure_reviewed_20260810/` 已生成可逆 reported sequence/IMGT/实验/AF3 映射和 FR-only Cα 对齐。实验结构 115 个有坐标 VHH 残基通过 source auth exact-WT 映射；AF3 126 个有坐标残基通过 source label ID 和 536-aa polymer 中唯一的 128-aa Nb252 连续片段映射。
- 82 个共同 framework Cα、Kabsch、无 outlier rejection：RMSD 0.631994 Å。拟合后 FR aggregate RMSD 0.631994 Å；CDR3 RMSD 6.490853 Å、最大位移 10.733119 Å。实验结构仍是结合构象证据，AF3 始终是预测。
- `interface_reviewed_20260810/`：严格 polymer heavy-atom center `<4.0 Å`，排除 H/D、非正 occupancy、水/配体/糖/离子、晶体/NCS images，并遵守 altloc 兼容性。结果为 246 个原子接触对、24 个 VHH 界面残基。
- 严格 `<4.0 Å` 残基集合与确认橙色 24 位点完全相同；保护并集为 reported sequence index `33,37,45,46,47,58,98,100-116`。该集合仅用于保守保护，不是能量热点或突变效应结论。

## Expression Audit

- LTT/WCC 31 条：仅允许来源内 exploratory numeric，状态 conditional。
- LLJ 16 条：仅允许来源内 ordinal/censored 探索，状态 conditional。
- 47 条跨来源 pooling 和向 Nb252 候选转移均 blocked；“同一体系”是 user-provided 概括证据，不能替代具体协议字段。

## Current Gates

- `input_freeze_manifest.status=pass`。
- `local_baseline_build=pass`：结构导出、清单、链身份、可逆映射和临时界面安全均已完成。
- `candidate_design_release=blocked`：唯一阻断为 authoritative Nb252 mature/full construct 确认，包括准确边界和末端 `GS` 解释。
- `pooled_expression_model_release=blocked`：缺少跨 assay 协议等价性证据。
- `finalize_input_baseline.py` 已用真实 structure mapping/interface 重建 canonical 128 位点图和 `stage1_gate.json`。程序 `status=pass` 不代表科学 release gate 自动通过。

## Verification

- 全套：`96 passed, 1 skipped, 4 subtests passed`；唯一 skip 为 Windows 真实 symlink 权限测试。
- `pip check`、`python -m compileall -q src scripts tests`、`git diff --check` 均通过。
- 修复覆盖：ChimeraX sandbox 入口、实际模型名、无 polymer sequence 时严格 source-auth exact-WT 映射、较长 polymer 中唯一 authoritative segment、`not_evaluable` summary 状态。

## Required Next Steps

1. 向合作者确认 Nb252 的成熟 VHH/完整表达构建体准确序列边界、末端 `GS` 的来源与是否属于构建体；确认后更新 `authoritative_construct_review`，重跑结构/finalizer，才可释放候选设计。
2. 向合作者补齐 LTT/WCC/LLJ 的宿主、载体、信号肽、标签、诱导、表达区室、纯化、定量、批次、重复、误差和跨来源协议等价性；证据不足时不得 pooled training 或向 Nb252 外推。
3. release gate 通过后再进入第二阶段候选设计：先冻结保护位点和允许突变空间，再以 nanoBERT、结构/developability 过滤、PyRosetta 与 AF3 分层评估，不把计算分数称为实验亲和力或表达量。
