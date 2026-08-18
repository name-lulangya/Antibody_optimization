# TNP natural-VHH reference set

本目录冻结 TNP 论文分析实际使用的 4,059 条非冗余天然 VHH 序列，来源为 TNP 仓库提交 `a9ba3edc3d967ecf8a2b9b5c2c29bf7495bbc9a0` 中的最终 VHH-OAS 描述符表：

- 文件：`VHH_OAS_all_properties_FINAL.csv`
- SHA-256：`D87A2E66CE0E46D34547D25DF10BF07ABC10B06CCF0D9D0C4A304A36A9D0EBE5`
- 样本数：4,059；序列 ID 和氨基酸序列均唯一。

这些序列源自 Li 等 2016 年报道的三只双峰驼天然 VHH repertoire，并经 OAS/TNP 汇集。它们用于估计天然序列保守性，不含本项目可直接使用的 BL21 表达量标签。

TNP 仓库当前同时存在 4,383 条的原始 `vhh_oas.fasta`；它与论文最终分析表的 4,059 条范围不一致，因此本项目没有将其作为权威输入。完整来源、版本和选择理由记录在保守性结果目录的 `source_manifest.json` 与 `nb252_vhh_conservation_contract.json` 中。
