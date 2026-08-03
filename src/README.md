# 可复用源码索引

经验证、具有稳定输入输出契约并可跨分析流程复用的 Python 逻辑放在：

`src/antibody_optimization/`

实验编排、一次性探索和特定批次参数不应仅为复用形式而机械迁入该包。新增、删除、重命名或实质修改共享工具时，应在同一任务中更新本索引，并记录其用途、输入、返回值、算法假设及明确不支持的范围。

## 当前工具

- `antibody_optimization.nb_expression`
  - 用途：从表达量 DOCX 的 `word/document.xml` 中按正文段落顺序解析 LTT、WCC 和 LLJ，执行逐序列 SHA-256、计数、来源原文和产量语义校验。
  - 主要输入：沿用当前三段格式的 `.docx`，以及显式传入的文档标题培养体积元数据；当前冻结输入 SHA-256 为 `a6e4022f0978fbd70a0e04dc78f479140ab6f55caaa90b467fb77a62eb5db5d1`。
  - 主要返回：`ExpressionRecord` 列表和 `DocxParagraph` 列表；后者同时保留 `raw_text` 与仅供语法解析的 `parse_text`。
  - 算法假设：序列位于对应标题或 clone ID 后一个非空正文段落；序列原文必须等于解析文本且只含 20 种标准大写氨基酸字母；LTT/WCC 是个体近似值，LLJ 是共享分档且不生成个体点估计。
  - 明确不支持：Word 表格、文本框、tab/换行节点、修订或隐藏文字、未知分段、自动修复序列、VHH 边界裁剪、抗体编号和构建体/标签识别。遇到可能隐藏或拆分序列的 DOCX 结构会直接失败，不静默跳过。
- `antibody_optimization.nb_expression_artifacts`
  - 用途：把解析记录写为样本、产量观察、实验上下文、宽表、原文转录和 FASTA，生成 manifest/validation/QC SVG，并逐字段回读验证全部 CSV 与序列输出。
  - 主要输入：已通过 `validate_records` 的记录和段落、输出路径、生成时间及源文件路径。
  - 主要返回：UTF-8 with BOM CSV、ASCII FASTA、UTF-8 JSON/SVG，以及结构化写出校验报告。
  - 算法假设：实验上下文从记录元数据归并，不从可变文件名推断；QC 图只展示来源和语义计数。
  - 明确不支持：XLSX 编写、表达性能推断、分档插值、VHH 裁剪或任何序列变换。
- `antibody_optimization.file_transaction`
  - 用途：在写出前拒绝 source/target、target/target 的精确或祖先路径碰撞，并用同目录安装候选、备份与 rollback 事务替换一组 staged 文件；同目录候选确保 Windows 最终文件继承目标目录 ACL，而不是保留私有临时目录权限。
  - 主要输入：项目根、受保护源路径、staged→final 路径对；测试可注入符合 `os.replace` 合同的替换函数。
  - 主要返回：通过后的规范化路径，或完成全部替换；不返回部分成功状态。
  - 算法假设：最终目标必须位于项目根内，已存在目标必须是普通非符号链接文件，父目录必须存在；提交前需为全部 staged 文件保留一份同目录候选空间。Windows 候选继承目标父目录 ACL，POSIX 候选保留 staged 文件 mode。
  - 明确不支持：跨项目根写出、目录替换、符号链接目标、时间戳/xattr/显式 ACL 复制和操作系统崩溃级持久性保证；运行期替换失败会尽力完整回滚并显式报告不完整 rollback。
