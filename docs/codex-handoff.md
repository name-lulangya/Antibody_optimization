# Codex 项目交接

Last updated: 2026-08-03 15:09:46

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是 replace-in-place 的当前状态快照；状态变化时应重写过期内容，而不是追加日记。

## Project Goal

优化实验筛选获得的 NK2R 纳米抗体 Nb252，重点改善亲和力、稳定性和 reported yield/表达相关性质，同时保留可追溯的结构、序列和实验依据。

## Current Project Status

- Git 分支为 `main`，同步目标为 `origin/main`；远程地址为 `git@github.com:name-lulangya/Antibody_optimization.git`。本轮表达量数据提取与权限修复也以该分支为同步目标。
- 计划在远程服务器父目录 `/homes/Tianlab/luly25/` 下检出项目；实际仓库根目录和登录别名仍待核验。
- 项目主 Conda 环境为 `/data/software/env/luly25/ab_optim`。
- PyRosetta 工具环境为 `/data/software/env/luly25/multi_ligand`：Python 3.10.20、PyRosetta 2026.03、Rosetta commit `5e498f1409c68ade56c8ce5842bf79e1b02e8db4`。
- nanoBERT 工具环境为 `/data/software/env/luly25/vhh-lm`：`NaturalAntibody/nanoBERT` revision `edc8182ad89a827f8737fa572c6b5fac6197e6b0`，使用已记录的本地缓存和离线模式。
- Slurm 默认使用 `batch`，至少 1 GPU、每 GPU 12 CPU，默认不显式指定内存；多 GPU 只用 `n1`/`n2`，单节点多 GPU 设置 `--exclude=n3`。
- 结构输入 `Nb252-optimization.cxs` 为 537670 bytes，SHA-256 `1BC636C28F66AE60EDC658D2E1C4AAD0B07F4141CA5411C78662AA19DA793C4D`。
- 新增表达量源文件 `nb序列及产量（1L）.docx`，14172 bytes，SHA-256 `a6e4022f0978fbd70a0e04dc78f479140ab6f55caaa90b467fb77a62eb5db5d1`；源文件未被改写。
- DOCX 已自动提取为 47 条记录：LTT 23、WCC 8、LLJ 16。全部序列逐段回读、长度和 SHA-256 校验通过，47 条均唯一且 mismatch 为 0。
- Git 跟踪制品位于 `docs/result_artifacts/nb_expression/`，包括分表 CSV、宽表、原文转录、FASTA、manifest、validation、数据字典和 QC SVG；run summaries 位于 `docs/run_summaries/nb_expression/`。
- 已修复 Windows 文件访问权限：事务安装先复制到最终父目录内的候选文件，使结果继承项目 ACL。当前全部结果制品和 run summary 均允许交互账号 `Tian_lab_luly25\\16217` 访问，且不是只读；修复前后文件内容哈希一致。
- `LTT__Nb252` 序列长 128 aa，SHA-256 `df5b83ddde8a3486383c12afe45e22af6a358f507eab5503d5dbd4430710288d`，reported yield 为个体近似 `~0.5 mg`。

## Active Workflows

- `scripts/data_preparation/prepare_nb_expression_dataset.py`：从冻结哈希的 DOCX 和显式文档标题体积元数据生成可审计数据集；默认拒绝覆盖，覆盖时使用备份/rollback 事务。
- `scripts/data_preparation/verify_nb_expression_outputs.py`：不导入生产解析器，以独立状态机复核全部序列与产量字段。
- `tests/`：28 个测试中 27 个通过；1 个真实符号链接测试因当前 Windows 权限跳过。已覆盖真实源文件、所有序列、CLI 路径碰撞、词法符号链接终点、事务 rollback、同目录安装候选、Windows ACL 继承、固定时间双跑字节一致和隐式时间戳重放。
- 尚未启动突变设计、结构打分、nanoBERT/PyRosetta 预测或实验候选排序。

## Data Semantics and Cautions

- 表型统一命名为 `reported_yield`，不是已校准的表达速率，也未自动换算为 `expression_mg_per_l`。
- LTT/WCC 的 31 条是个体近似值；LLJ 是共享分档：9 条 `>20 mg` 只存下界，6 条 `~10 mg` 和 1 条 `~2 mg` 只存组级锚点。LLJ 个体点估计全部留空。
- WCC 明确记录 1 L TB 和纯化后近似产量；LTT/LLJ 的 1 L 仅来自原文档标题，并作为显式 CLI 元数据传入，不从文件名推断。跨来源协议、构建体、批次、重复、回收率和误差尚未核实。
- 所有序列均按文档完整保留，`sequence_scope=unknown`；未裁剪 VHH。WCC 的来源特异末端可能造成建模泄漏。
- 原始段落与解析段落分列保存；只有标题段落 40、42 含尾随空白，序列段落均未做空白归一化。可能隐藏/拆分内容的 Word 结构会触发失败。
- 6 条序列短于 115 aa，`Nb257` 与 WCC `4-11` 各仅有一个字面 Cys，9 条不含 `WGQGT`；这些是人工复核标记，不是自动错误判定。
- CSV 使用安全 `sample_uid`；导入 Excel 时仍需把 `source_sample_id` 设为文本，避免 `1E2`/`4-7` 自动转型。
- 当前会话未提供合规的 `@oai/artifact-tool` 运行时，因此未生成 XLSX，也未使用未经授权的替代 Excel 库。

## Structural Context and Cautions

- `.cxs` 内模型名：`NK2R-252.pdb`（实验 NK2R–Nb252）、`NK2R-NKA.pdb`（NKA 结合背景）、`fold_2r_252_nomg_model_0.cif`（AF3 VHH）；它们尚不是仓库中的独立结构文件。
- 实验 VHH 有未搭建区域；用户/合作者目视观察称除 CDR3 外 AF3 与实验结构总体对齐良好，尚未定量核验，AF3 不能替代实验构象证据。
- `NK2R-252` 中橙色 VHH 残基是合作者留下的推定界面注释；“<4 Å”的原子/距离定义、确切残基和编号映射仍未知，突变前必须重新计算并谨慎处理。

## Recent Changes

- 完成 47 条合作者序列与 reported-yield 数据的无手工转录提取。
- 将样本、产量观察和实验上下文分表，并提供方便使用的宽表与 FASTA。
- 为每条序列加入长度、完整 SHA-256、源段落索引和非破坏性人工复核标记。
- 严格区分 LTT/WCC 个体近似值与 LLJ 组级下界/近似锚点。
- 将解析与制品写出拆成聚焦模块，并新增独立验证器、完整字段回读、28 个安全/回归测试、manifest、validation 和 run summaries。
- 加入输出路径碰撞检查、词法符号链接终点拒绝和多文件备份/rollback 事务；run summary 改为结构化 argv 与 PowerShell 重放命令，省略 `--generated-at` 时也会记录实际时间以便确定性重放。
- 修复私有 staging 文件经 Windows `os.replace` 后无法由 VSCode/Excel 读取的问题；最终文件现在通过同目录候选安装并继承目标目录 ACL，全部数据哈希保持不变。
- 新增数据字典与只展示提取计数/语义的可复现 QC SVG。

## Suggested Next Steps

1. 向合作者确认 LTT/WCC/LLJ 是否共享宿主、载体、标签、诱导、纯化、定量和批次条件，以及每条数据的重复与误差。
2. 回查短序列、单 Cys 和缺少 `WGQGT` 的记录，确认测序与构建体边界；在不改原序列的前提下生成经验证的统一 VHH 区域。
3. 建立抗体编号/CDR 注释，并按来源或批次分层设计表达相关建模与验证拆分，避免构建体/来源泄漏。
4. 导出并核验 `.cxs` 中的独立结构、链 ID、缺失残基、编号映射和橙色界面残基定义，再开展 Nb252 突变设计。
5. 在 `@oai/artifact-tool` 合规运行时可用后，补充 ID 列强制文本格式、分表和数据字典工作表的 XLSX。
