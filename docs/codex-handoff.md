# Codex 项目交接

Last updated: 2026-07-31 16:48:33

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是 replace-in-place 的当前状态快照；状态变化时应重写过期内容，而不是追加日记。

## Project Goal

优化实验筛选获得的 NK2R 纳米抗体 Nb252，重点改善亲和力、稳定性和表达量，同时保留可追溯的结构依据与设计记录。

## Current Project Status

- 当前可验证的科学/结构输入仅有根目录下的 `Nb252-optimization.cxs`。
- Git 仓库已在根目录初始化，当前分支为 `main`，同步目标为 `origin/main`。
- Git 远程 `origin` 已配置为 `git@github.com:name-lulangya/Antibody_optimization.git`。
- 用户计划将项目检出到远程服务器的 `/homes/Tianlab/luly25/` 目录下；实际仓库根目录尚待服务器端检出后核验。
- 项目自身的远程 Conda 环境为 `/data/software/env/luly25/ab_optim`。
- PyRosetta 工具环境为 `/data/software/env/luly25/multi_ligand`：Python 3.10.20、PyRosetta 2026.03、Rosetta commit `5e498f1409c68ade56c8ce5842bf79e1b02e8db4`。
- nanoBERT 工具环境为 `/data/software/env/luly25/vhh-lm`：模型 `NaturalAntibody/nanoBERT`，固定 revision `edc8182ad89a827f8737fa572c6b5fac6197e6b0`；使用本地 Hugging Face 缓存和离线模式。完整激活与检查命令见 `docs/history/2026-W31.md`。
- 远程计算使用与参考项目一致的 Slurm 规则：默认 `batch` 分区、至少 1 GPU、每 GPU 配 12 CPU，默认不显式指定内存。
- 多 GPU 作业仅允许使用 `n1` 或 `n2`；单节点多 GPU 作业须请求一个节点并设置 `--exclude=n3`。
- 文件大小：537670 bytes。
- SHA256：`1BC636C28F66AE60EDC658D2E1C4AAD0B07F4141CA5411C78662AA19DA793C4D`。
- 根据合作者说明，该 ChimeraX 会话包含三个 model：
  - `NK2R-252.pdb`：实验解析的 NK2R–Nb252 结合构象。
  - `NK2R-NKA.pdb`：NK2R–NKA 结合结构。
  - `fold_2r_252_nomg_model_0.cif`：AF3 预测的 Nb252 VHH 结构。
- 实验结构中的 VHH 有部分区域未搭建；补充的 AF3 VHH 与实验 VHH 对齐后，除 CDR3 外其余区域总体对齐良好。
- `NK2R-252` 中合作者将推定界面 VHH 残基标为橙色，据称筛选条件为距离小于 4 Å；距离定义、确切残基列表和编号映射尚未记录或复核。
- 尚未记录缺失残基范围、链 ID、编号映射或定量结构比较结果。

## Active Experiments

- 当前没有已启动或已记录结果的计算实验。
- 尚未建立突变候选、打分流程或实验验证批次。

## Recent Changes

- 建立 NK2R/Nb252 项目规则、handoff、周历史和可复用源码索引。
- 将 ChimeraX 会话标记为 Git 二进制，并在 `main` 分支初始化仓库。
- 将 GitHub SSH 地址配置为本地仓库的 `origin`，并记录计划使用的远程服务器父目录与环境待办。
- 增加与参考项目一致的 Slurm 分区、GPU/CPU、内存及多 GPU 节点约束。
- 区分项目主环境 `ab_optim` 与 PyRosetta、nanoBERT 两个独立工具环境，并记录工具固定版本/revision 及离线缓存要求。
- 记录 `NK2R-252` 橙色 VHH 界面残基注释及其 4 Å 定义未明确的风险。

## Decision-Relevant Cautions

- `NK2R-252.pdb`、`NK2R-NKA.pdb` 和 `fold_2r_252_nomg_model_0.cif` 是当前 `.cxs` 会话内显示的 model 名称，不是项目目录中的独立 PDB/CIF 文件。
- 未从会话显式导出并核验前，不应假定上述 model 可由普通结构解析工具直接读取。
- CDR3 是实验结构与 AF3 预测结构的主要不一致区域；不能在未核验实验密度、缺失区段和编号映射的情况下直接据此设计突变。
- 橙色 VHH 残基只是合作者在 ChimeraX 中留下的界面注释；“距离小于 4 Å”尚未说明是任意原子、重原子、主链或其他距离定义。重新计算并形成明确残基表前，不得将其当作已验证的精确界面集合，后续突变需谨慎。
- 当前结构描述来自用户与合作者说明，尚无脚本化检查或定量分析结果。
- 项目环境及两个工具环境的信息来自用户提供，尚未通过服务器会话独立检查；不得把工具环境当作项目环境，也不得假定工具依赖同时存在于 `ab_optim`。实际项目检出路径和登录别名仍待核验。

## Deferred Report Items

- 暂无已指定延期至后续周报的内容。

## Suggested Next Steps

1. 在 ChimeraX 中核验三个 model 的 model ID、链 ID、残基编号和缺失区段，并导出 `NK2R-252` 中橙色 VHH 残基列表及其选择定义。
2. 将需要程序化分析的结构显式导出为独立 PDB/mmCIF，并记录导出命令、文件哈希和来源 model。
3. 建立实验 VHH 与 AF3 VHH 的残基级编号映射，单独标注 CDR1–CDR3、框架区和未建模残基。
4. 确定亲和力、稳定性与表达量优化的优先级、允许突变范围及实验筛选约束。
5. 在真实导出结构上定义界面与基线指标后，再制定可复现的候选突变生成和排序流程。
