# 纳米抗体 reported yield 数据字典

## 数据来源与范围

- 源文件：`nb序列及产量（1L）.docx`
- 源文件 SHA-256：`a6e4022f0978fbd70a0e04dc78f479140ab6f55caaa90b467fb77a62eb5db5d1`
- 记录数：47（LTT 23、WCC 8、LLJ 16）
- 序列：47 条原文字面序列，未裁剪、未改大小写、未修复；每条均保存长度和 SHA-256。
- 表型名称使用 `reported_yield`。当前数据是文档所报产量/产量分档，仅 WCC 明确为纯化后产量；它不等同于已校准的表达速率，也不自动换算为 `mg/L` 表达量。

## 文件用途

| 文件 | 用途 |
|---|---|
| `nb_expression_records.csv` | 方便浏览和分析的宽表；序列、产量语义与上下文均在同一行。 |
| `samples.csv` | 样本、原始序列、序列哈希、来源段落和序列复核标记。 |
| `yield_observations.csv` | 产量观察长表；个体近似值、组级近似锚点和下界严格分列。 |
| `assay_context.csv` | LTT/WCC/LLJ 三个来源的已知实验上下文和缺失信息。 |
| `raw_transcription.csv` | DOCX 中 100 个非空正文段落的顺序转录；分别保存字面 `raw_text` 和仅供语法解析的 `parse_text`。 |
| `nb_expression_sequences.fasta` | 使用稳定 `sample_uid` 的原序列 FASTA；每个 header 带序列 SHA-256。 |
| `manifest.json` | 数据来源、字段、语义、限制和各输出哈希。 |
| `validation_report.json` | 主解析器的计数、字符集、哈希和回读校验结果。 |
| `qc_plot_data.csv` / `nb_expression_qc.svg` | QC 图的精确源数据和可复现矢量图，仅展示提取计数与产量语义。 |

主数据 CSV 为 UTF-8 with BOM；仅供绘图脚本读取的 `qc_plot_data.csv` 为无 BOM 的 UTF-8。用 Excel 导入时，仍应把 `source_sample_id` 指定为文本，避免 `1E2` 被解释为科学计数、`4-7` 被解释为日期；跨文件关联优先使用不会自动转型的 `sample_uid`，例如 `LLJ__1E2`。

## 关键标识与序列字段

| 字段 | 定义 |
|---|---|
| `sample_uid` | `provider_code__source_sample_id`，跨表稳定主键。 |
| `provider_code` / `source_section` | 来源段：`LTT`、`WCC` 或 `LLJ`。 |
| `source_sample_id` | 文档中的原始 clone ID；WCC 行首的 `>` 不属于 clone ID。 |
| `source_header_raw` | 含 clone ID 的原始标题；LLJ 因 clone ID 独立成段，此字段为原 clone ID 段。 |
| `sequence_raw` | 从对应 DOCX 正文段落直接提取的完整字面序列。 |
| `sequence_length_aa` | `sequence_raw` 的字符数。 |
| `sequence_sha256` | 对 `sequence_raw` 的 ASCII 字节计算的 SHA-256。 |
| `sequence_scope` | 当前固定为 `unknown`；尚未判定每条字符串是成熟 VHH、完整构建体或带末端延伸的构建体。 |
| `vhh_region_sequence` | 当前留空；完成构建体核验和抗体编号后才能生成。 |
| `source_*_paragraph_index` | 在 `raw_transcription.csv` 中使用的 1-based 非空正文段落索引。 |

`raw_text` 不删首尾空白；`parse_text` 只为识别 section、clone 和产量标题而去除段落首尾空白。序列段落若两者不完全相同，解析器会直接失败，不会静默清理序列。本源文件只有非序列标题段落 40、42 含尾随空白。

`sequence_review_flags` 只提示人工复核，不代表序列错误，也不会触发自动修改或删除：

- `short_literal_sequence_lt115`：字面长度小于 115 aa；
- `single_cysteine_literal_sequence`：字面序列仅含一个 Cys；
- `wgqgt_motif_absent`：字面序列不含字符串 `WGQGT`。

## 产量字段与 LLJ 分档

| 字段 | 定义 |
|---|---|
| `reported_text` | 原始产量标题，例如 `Nb252  ~0.5 mg`、`>20 mg`。 |
| `observation_semantics` | `individual_approximate`、`group_lower_bound` 或 `group_approximate`。 |
| `value_relation` | `approx` 或 `gt`。 |
| `point_estimate_mg` | 仅 LTT/WCC 的个体近似值有数值；LLJ 全部留空。 |
| `group_anchor_mg` | 仅 LLJ 的 `~10 mg`、`~2 mg` 分档分别填 10、2；不是个体点估计。 |
| `lower_bound_mg` | 仅 LLJ 的 `>20 mg` 分档填 20；`lower_bound_inclusive=False`。 |
| `upper_bound_mg` | 当前全部留空；文档没有提供分档上界。 |
| `assignment_level` | `individual` 或 `group`。 |
| `group_id` | LLJ 分档标识，例如 `LLJ_GT20`、`LLJ_APPROX10`。 |
| `individual_numeric_available` | LLJ 全部为 `False`。 |
| `censoring_type` | `>20 mg` 标为 `right_censored`；这只是数值约束描述，删失原因未知。 |
| `replicate_count` / uncertainty 字段 | 当前留空；不能把未报告的重复或误差填成 0 或 1。 |

不得把 LLJ `>20 mg` 当作精确 20 mg，也不得为 `~10 mg`、`~2 mg` 臆造分档区间或个体数值。普通连续回归应排除 LLJ 的个体点估计；若后续使用删失模型或来源内序数任务，必须显式保留上述语义。

## 实验上下文限制

- WCC 每条记录明确写有 1 L、TB 和“纯化得到约”；因此记录为 `volume_evidence=entry_text`、`medium=TB`、`yield_stage=post_purification`。
- LTT/LLJ 的 1 L 仅来自原文档标题，并通过 CLI 元数据参数显式传入，不从可变文件名推断；培养基、纯化阶段和方法均未提供。
- 宿主、诱导条件、构建体边界、标签、批次、纯化回收率、定量方法、重复数和不确定度尚未记录。
- WCC 序列存在来源特异的 N/C 端特征；在核验构建体并提取统一 VHH 区域前，模型可能学习到来源或构建体差异，而不是目标性质。
- 跨来源建模前应确认是否确为同一实验体系，并按来源/批次分层拆分数据，避免来源泄漏。
