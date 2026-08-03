# Codex 项目交接

Last updated: 2026-08-03 21:35:46

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是 replace-in-place 的当前状态快照；状态变化时重写过期内容，不追加日记。

## Project Goal

优化实验筛选获得的 NK2R 纳米抗体 Nb252 的亲和力、稳定性和 reported yield/表达相关性质，同时保留可追溯的结构、序列、实验依据和发布门。

## Current Status

- Git 分支为 `main`，远程为 `git@github.com:name-lulangya/Antibody_optimization.git`；本轮第一阶段实现纳入当前 Git HEAD，并以 `origin/main` 为同步目标。
- 本地项目环境：`D:\miniconda\envs\ab_optim`，Python 3.11.15、ANARCII 2.0.8、Gemmi 0.7.5。完整本地依赖见 `requirements-local*.txt`。
- 远程项目环境：`/data/software/env/luly25/ab_optim`；计划检出父目录为 `/homes/Tianlab/luly25/`，登录别名和实际仓库根仍待核验。
- PyRosetta 工具环境：`/data/software/env/luly25/multi_ligand`，Python 3.10.20、PyRosetta 2026.03、Rosetta commit `5e498f1409c68ade56c8ce5842bf79e1b02e8db4`。
- nanoBERT 工具环境：`/data/software/env/luly25/vhh-lm`，模型 `NaturalAntibody/nanoBERT` revision `edc8182ad89a827f8737fa572c6b5fac6197e6b0`，使用已记录的离线缓存。
- 本阶段没有运行 AF3、PyRosetta、nanoBERT、候选生成、突变设计或表达量模型训练。

## Frozen Inputs and Completed Local Artifacts

- `Nb252-optimization.cxs`：537670 bytes，SHA-256 `1BC636C28F66AE60EDC658D2E1C4AAD0B07F4141CA5411C78662AA19DA793C4D`，mtime 未改变。
- `nb序列及产量（1L）.docx`：14172 bytes，SHA-256 `a6e4022f0978fbd70a0e04dc78f479140ab6f55caaa90b467fb77a62eb5db5d1`，mtime 未改变。
- `docs/result_artifacts/input_baseline/sequence/`：47 条真实序列 provisional IMGT 审核；46 pass、1 failed。唯一失败为 `WCC__4-28`（`Score less than cut off.`）。成功结果为 H=45、L=1；chain type 只是工具输出。
- Nb252：ANARCII H，原索引 0–125，共 126 aa 进入编号，IMGT 1–128，末端 `GS` 未编号；仍是 provisional 边界。
- `WCC__4-11`：ANARCII L、止于 IMGT117，必须人工复核，不能据此宣称真实轻链。`WCC__4-1` 仅编号 87 aa、止于 IMGT102。
- `docs/result_artifacts/input_baseline/expression/`：78 行 assay 字段审核和 47 行样本审核。LTT/WCC 31 条只允许来源内 exploratory numeric；LLJ 16 条只允许来源内 ordinal/censored；47 条全部禁止跨来源 pooling 和向 Nb252 转移。
- `docs/result_artifacts/input_baseline/summary/`：输入冻结、128 位点绘图数据、状态计数、600 dpi PNG、SVG、summary manifest 和 `stage1_gate.json` 已生成。图中结构轨道当前全为 `not_available`，橙色/临时界面为空，不能解释为实测结构覆盖。
- 固定时间戳真实数据双跑时，序列、表达审核和 summary 的正式数据/图制品逐字节一致；run summary 的耗时字段属于运行遥测，不要求逐字节固定。
- 所有人工 CSV 为 UTF-8 with BOM、LF、固定列序；PowerShell `Import-Csv` 可读，全部新增制品继承项目 ACL 且不是只读。

## Current Gates

- `input_freeze_manifest.status=pass`。其中 `construct_boundary_confirmed=false` 表示冻结源输入本身不能证明边界，不是动态 review 状态。
- `local_baseline_build=blocked`：缺 structure export、inventory、chain identity、residue mapping 和 interface safety。
- `candidate_design_release=blocked`：上述五项之外还缺 authoritative Nb252 sequence/construct confirmation。
- `pooled_expression_model_release=blocked`：缺跨 assay 协议等价证据。
- Finalizer 的 `status=pass` 只表示汇总程序成功，不表示任何 blocked 科学发布门已通过。

## Structure and Interface State

- CXS 中三个目标模型名为 `NK2R-252.pdb`、`NK2R-NKA.pdb`、`fold_2r_252_nomg_model_0.cif`；当前 `data/structures/cxs_exports/` 尚不存在。
- 源会话由 ChimeraX 1.9/macOS 生成是 `user_provided` 信息；本轮尚未从会话自动验证。导出必须使用本机 ChimeraX 1.12。
- `structure_precheck/structure_baseline_manifest.json` 和 `interface_precheck/interface_manifest.json` 是诚实的 blocked precheck；尚未导出结构、确认链角色、建立真实结构映射或得到界面残基。
- 导出器会清单化全部 session model（包括 surface/group/child），但仅要求三个目标 `AtomicStructure` 恰好唯一匹配；保存 native/reference-frame mmCIF、变换、计数、颜色/显示/选择和 surface 回映信息。
- source-aware mapping 直接读取 mmCIF `_entity*`、`_struct_asym` 和 `_pdbx_poly_seq_scheme`，优先源 polymer/entity 与源 label/auth 编号；源证据存在但冲突时阻断，只在对应证据缺失时允许固定参数序列 fallback。
- 单一 `baseline_review.json` 分别审核链角色、精确 orange RGBA/渠道和 authoritative construct。chain/orange 可先确认而 construct 保持 pending；前者允许完成结构/界面基线，后者继续阻断候选设计。
- 临时保护集合定义为“confirmed orange ∪ strict non-H/D atom-center `<4.0 Å`”；`not_evaluable` 不等于 false，该集合不是能量热点。功能已实现但当前没有实际残基集合。

## Active Entry Points

- `build_sequence_review.py`：ANARCII/IMGT 审核。
- `build_expression_audit.py`：assay 与样本可比性审核。
- `export_cxs_session_chimerax.py`：ChimeraX 1.12 会话导出，无覆盖。
- `build_structure_baseline.py`：Gemmi 清单、source-aware 映射和 FR-only Cα 定量对齐；每轮用新输出目录。
- `calculate_temporary_interface.py`：review 后的严格 `<4.0 Å` 接触和橙色比较；每轮用新输出目录。
- `finalize_input_baseline.py`：冻结、stage gates、紧凑绘图数据和正式图。
- `plot_input_baseline.py`：只从紧凑 CSV 复现图。
- 详细接口、假设和不支持范围见 `src/README.md`；实际运行摘要在 `docs/run_summaries/input_baseline/`。

## Verification

- 全套测试：91 passed、1 skipped、4 subtests passed；跳过项仍是 Windows 真实 symlink 权限测试。
- `pip check`、`python -m compileall -q src scripts tests`、`git diff --check` 均通过。
- CXS 真实集成验收尚未执行，因为它要求用户保存其他 ChimeraX 工作、视觉核验会话并明确确认链/橙色；不得用自动猜测补过该门。

## Required Next Local Steps

1. 用户先保存任何现有 ChimeraX 工作；在 ChimeraX 1.12 中打开原始 CXS，视觉确认三个原子模型和橙色注释，不重新保存会话。
2. 在 ChimeraX 命令行运行：

```text
runscript "C:\Users\16217\Desktop\Codex Projects\Antibody_optimization\scripts\input_baseline\export_cxs_session_chimerax.py" --source-cxs "C:\Users\16217\Desktop\Codex Projects\Antibody_optimization\Nb252-optimization.cxs" --expected-source-sha256 1BC636C28F66AE60EDC658D2E1C4AAD0B07F4141CA5411C78662AA19DA793C4D --output-dir "C:\Users\16217\Desktop\Codex Projects\Antibody_optimization\data\structures\cxs_exports"
```

3. 在全新目录首跑 `build_structure_baseline.py`，取得 inventory 和 `baseline_review_template.json`；复制为 `baseline_review.json`，确认所有链角色及唯一 orange RGBA/渠道。没有合作者证据时，`authoritative_construct_review` 保持 pending。
4. 在另一全新目录以 `--confirmed-review` 重跑 structure baseline；其通过后在新目录运行 interface 入口。然后给 finalizer 增加真实 structure mapping、interface manifest 和 `orange_vs_4A.csv`，用 `--overwrite` 更新 canonical summary/gates。
5. 向合作者确认 Nb252 成熟 VHH/完整表达构建体边界、末端 `GS` 来源，以及 LTT/WCC/LLJ 的宿主、载体、信号肽、标签、诱导、表达区室、纯化、定量、批次、重复、误差和跨来源协议等价性。
