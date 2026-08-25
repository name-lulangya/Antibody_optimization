# Nb252 V3计算链与报告就绪性审计

审计状态：**Complete**  
审计日期范围：2026-08-19至2026-08-25  
审计提交范围：`862394d78d229618d13048229457f9be1ed2f759`至`e44c2ae87f215510ae15729c4dac6553656fc645`

## 1. 审计结论

审计对象为Nb252表达优化V3完整计算链及现行最终面板，审计提交为`e44c2ae87f215510ae15729c4dac6553656fc645`，审计范围从V3路线开始前的`862394d78d229618d13048229457f9be1ed2f759`延伸至该提交。现行唯一候选权威来源是：

`docs/result_artifacts/candidate_design/v3_final_15plus15_panel_20260825/v3_final_panel_manifest.json`

总体判定为：**通过，但存在必须披露的实质性科学与溯源限制**（`pass_with_material_caveats`）。

- 未发现V3候选计算链、序列身份、约束执行或报告准备层面的Critical问题。
- 未发现亲本或候选序列错配、突变位置重建错误、候选遗漏、冻结位点违规、直接界面位点违规、Cys22/Cys95改变、末端`SSGS`改变、CSV/FASTA错配或manifest哈希漂移。
- 847条单突至最终30条的身份链完整闭合；机器审计18/18项通过。
- 全套项目验证为320 passed、1 skipped、4 subtests passed；`compileall`、`py_compile`、`pip check`、`git diff --check`及证据JSON重复生成一致性均通过。
- 当前最终面板包含15条单突和15条双突。展示顺序不是效力排名，全部结果仍是计算优先的实验假设，不是已验证的表达改善。
- W34报告、PPT及其生成入口对应历史V2的19单突＋11双突，这是预期的历史版本状态，不构成V3审计缺陷。现行V3项目报告已从final manifest、final CSV/FASTA及本审计重新生成DOCX/PDF并完成manifest绑定；历史V2材料保持原样。

本审计的机器可读支持证据为`docs/result_artifacts/weekly_report_result/Nb252_V3_expression_report/Nb252_V3_audit_evidence.json`；该JSON由`scripts/reporting/audit_v3_release.py`从实体CSV、FASTA、JSON和结构映射重新联接生成。V3项目报告与审计报告统一存放在`docs/result_artifacts/weekly_report_result/Nb252_V3_expression_report/`。本轮范围仅包括项目报告DOCX/PDF及审计制品，不制作PPT或压缩交付包，其缺失不构成发布阻断。

## 2. 审计范围与方法

### 2.1 检查范围

本次审计覆盖：

1. 4,059条公开天然VHH来源、编号资格、去冗余簇和Nb252邻域保守性合同；
2. 128-aa权威亲本、reported-sequence位置、IMGT映射、实验结构和AF3结构映射；
3. 80个硬冻结位置、24个实验界面位置、二硫键Cys22/Cys95和末端125–128位`SSGS`；
4. 847条允许单突、完整NetSolP/NanoMelt/AntiFold证据和V3单突审计；
5. 30条上游单突短名单、31条专家审查池和15条父单突；
6. 105个理论父单突对、3个同位点互斥对、102条有效双突、共同评分矩阵和责任基序勘误；
7. 102条双突统一审查、15条双突终选和最终30条CSV/FASTA；
8. 各阶段manifest、输入/输出SHA-256、运行门、图表及V3报告DOCX/PDF绑定状态。

### 2.2 独立复核方法

审计没有只接受各阶段已有的`pass`字段，而是重新执行了以下独立检查：

- 从31条单突完整序列及其WT/Mut字段反向重建亲本，并与权威SHA-256核对；
- 按candidate ID和完整序列联接847条允许空间、性质矩阵及V3审计；
- 对每条单突和双突逐字符重算与亲本的差异位置；
- 从15条父单突独立枚举无序组合，验证`C(15,2)=105`、去除3个同位点对后恰为102条；
- 将最终15条双突回接102条性质矩阵，将最终30条回接父单突和双突序列；
- 对最终CSV与FASTA逐ID、逐序列比对；
- 回读reported-sequence、IMGT、实验结构和AF3映射，检查WT身份与坐标状态；
- 重算性质幅度档、AntiFold组成门、责任基序和`N76G+F30N`勘误overlay；
- 复核4个V3阶段manifest中的47个本地可用哈希绑定；
- 分别检查历史V2材料与现行V3报告生成器的来源合同，并核对V3 DOCX/PDF及manifest绑定；
- 执行全套项目测试、Python编译、依赖一致性与diff检查，并以相同输入重复生成证据JSON核对字节确定性。

### 2.3 未验证或未评价事项

- 双突远程运行的两个103行normalized score表和两个model-run JSON未随本地仓库归档；它们在manifest中的SHA-256可追溯，但本次不能从本地原始文件独立复算联接或核验实际远程环境记录。
- 当前路线没有进行双突侧链重建、双突AntiFold或双突能量最小化，因此突变后局部构象与远程结构耦联属于未评价事项，不得由WT位点距离代替。
- 最终30条尚无本项目BL21表达量、结合保持、聚集或热稳定性实验结果；所有改善与风险方向仍是计算或专家假设。

## 3. 完整性核验结果

| 检查层 | 结果 | 结论 |
|---|---:|---|
| 天然VHH来源 | 4,059条；4,057条编号合格；3,784个去冗余簇 | 与冻结来源合同一致 |
| Nb252邻域 | 1,564条序列；1,532个有效去冗余簇 | 与保守性合同一致 |
| 权威亲本 | 128 aa；SHA-256 `df5b83dd…710288d` | 唯一重建，无分歧 |
| 允许单突空间 | 847条 | ID、序列唯一；全部为单突 |
| 上游V3短名单 | 30条 | 是847条的精确子集 |
| 父单突专家池 | 31条 | 精确等于30条短名单＋T99F |
| 父单突终选 | 15条、12个位点 | CSV/FASTA与31条审计一致 |
| 双突空间 | 105个理论对－3个同位点对＝102条 | 无缺失、无额外、无重复 |
| 双突共同评分 | WT＋102条＝103条 | NetSolP/NanoMelt门均记录103/103通过 |
| 双突终选 | 15条 | 6条有3项中/强改善，9条有2项 |
| 最终面板 | 15单突＋15双突＝30条 | CSV/FASTA逐条一致 |
| 序列硬约束 | 30/30为128 aa、保留SSGS、恰有2个Cys | 冻结/界面重叠均为0 |
| manifest哈希 | 47个本地可用绑定全部一致 | 哈希漂移0 |
| 自动审计 | 18/18项通过 | 无机器完整性失败 |

### 3.1 结构映射

- 实验复合物与AF3各有128行reported-sequence映射；序列身份与亲本一致。
- 实验结构中113位有坐标、13位缺失坐标、末端2位属于未编号侧翼；AF3中126位有坐标、末端2位属于未编号侧翼。
- 最终15条父单突中，L11Y、F29Q和L11M仅能使用AF3结构上下文；最终15条双突中有4条因为包含L11替换而使用AF3-only位点距离。最终30条中合计7条构建体依赖实验缺失位点的AF3补充证据。

### 3.2 最终面板已登记的风险

- 4/30条带明确化学软责任：F30N和S55G＋F30N新增`NG`脱酰胺基序；L11M和N76G＋L11M增加Met氧化敏感残基。
- 6/30条包含reported position 30的F30S或F30N，属于同一CDR1缺口边界假设家族，不能视为六份独立结构证据。
- 2个父单突的专家结构判断为`structurally_concerning`且置信度为medium：A23R和T99F。
- 专家溶解度方向判断为不利的父单突有Q5V和T99F；专家热稳定性方向判断为不利的父单突有S55G和T99F。
- 最终15条双突全部使用两个在WT骨架上Cα距离至少10 Å的位点，没有局部/邻近位点对；这降低了直接局部耦联风险，但不能证明没有远程非加和效应。

## 4. 分级发现

### Critical

未发现Critical级问题。

### High

#### H-01：上游30单突制品仍自称“最终实验面板”

- **证据**：`expression_single_mutant_v3_gate.json`仍含`final_experimental_panel_released=true`，release字符串为`v3_final_30_single_mutants_ready_for_experimental_testing`；现行权威manifest则释放15单突＋15双突。
- **影响范围**：报告、PPT、交付文件夹和未来自动化读取；不影响已释放V3最终序列本身。
- **影响**：若报告脚本误读旧gate，会把历史上游30单突错交付为当前最终面板。
- **行动**：所有V3报告和交付只读取现行final manifest/CSV/FASTA；旧30单突只能称“不可改写的上游短名单”。建议后续新增轻量supersession索引，而不改写历史结果。
- **状态**：Open，但已由AGENTS、handoff和本审计控制。**在新V3报告完成前阻止对外发布，不阻止报告书写。**

#### H-02：T99F是已接受但风险明确的稳定词探索例外

- **证据**：T99F上游为`not_eligible/not_selected`；ΔU近似中性，ΔS为弱不利，预测ΔTm为−0.50°C；专家判断为`structurally_concerning`、溶解度/热稳定性方向均不利、置信度medium。该位点部分埋藏，距离最近硬界面残基约2.39 Å，并位于模型敏感的CDR3。
- **影响范围**：最终面板中的T99F单突；14条T99F双突均按通用规则审查且0条入选。
- **影响**：T99F不能与性质支持候选合并表述，也不能作为“预计提高表达”的证据。
- **行动**：V3报告将其单列为稳定词假设探索对照，完整保留弱不利预测与结构风险；实验中与WT直接比较，并设置结合/功能保持质控。
- **状态**：Accepted scientific exception；不属于遗漏，也不触发面板回滚。

### Medium

#### M-01：三类工具没有验证为Nb252突变产量预测器

- **证据**：NetSolP yield gate的证据等级为`compatibility_filter_only`且`nb252_expression_prediction_validated=false`；NanoMelt为`no_supported_use`且未验证Nb252表达预测；AntiFold产量分类为`not_applicable`，仅允许作结构条件序列相容性负向门。
- **影响**：最终候选不能描述为“预测提高BL21产量”或“已优化表达量”。
- **行动**：报告统一使用“表达相关性质的计算优先实验假设”；NetSolP U、S和NanoMelt预测Tm保留各自定义，不换算为mg/L。
- **状态**：Open methodological limitation；只能通过本项目实验验证收敛。

NetSolP U和S作为两个输出指标分别保留，但不是两个独立模型。847条单突中U/S Spearman为0.269，102条双突中为0.301；报告不得写成“两个统计独立证据”。

#### M-02：远程双突原始评分表未随本地仓库归档

- **证据**：双突property manifest绑定的两个103行normalized score表和两个model-run JSON在本地checkout中不存在；四个SHA-256仍保留，派生102行矩阵和全部本地输出哈希一致。
- **影响**：可验证派生矩阵的内部一致性和冻结身份，但本地不能从原始远程输出重新执行103行join或独立回读实际环境记录。
- **行动**：在V3报告发布前，优先将四个小型原始文件同步到受控归档，或在交付manifest中记录可访问的外部归档位置；若未完成，报告必须明确该溯源限制。
- **状态**：Open provenance gap；不导致现有序列或数值失效。

#### M-03：双突“专家审查”不是双突侧链结构验证

- **证据**：`double_sidechain_modeling_performed=false`；没有双突AntiFold分数、双rotamer建模、能量最小化或逐条ChimeraX双突构象检查。审查使用父单突人工判断、WT位点几何、完整序列责任基序和U/S/Tm幅度档。
- **影响**：不能声称102条或15条双突已经通过结构验证；远程耦联和突变后侧链相互作用仍未知。
- **行动**：报告准确称为“规则化专家复核与WT几何分层”，将15条双突的性质和风险定位为实验优先级证据。
- **状态**：Open scientific limitation；当前合同明确接受。

#### M-04：实验缺失坐标和F30缺口边界造成结构不确定性

- **证据**：最终30中7条涉及实验缺失坐标位置；L11相关双突的结构置信度为low。F30虽有实验坐标，但紧邻reported 24–29缺失片段；最终面板有6条F30家族构建体。
- **影响**：这些候选的局部包装和loop预组织判断弱于具有完整实验坐标的候选。
- **行动**：报告逐项标注AF3-only和F30 gap-boundary，不把AF3模型写成实验结构，也不把同一位点多个替换写成独立结构重复。
- **状态**：Open and disclosed。

#### M-05：软责任基序和专家风险不能被“hard risk=0”掩盖

- **证据**：4条最终构建体有明确化学软责任；A23R在Cys22邻位引入大体积正电侧链；S55G有loop柔性/热稳定不利假设；Q5V虽是唯一合法天然共识回变，但在暴露表面增加疏水性。
- **影响**：这些风险可能影响制备、储存或折叠，且部分不由NetSolP/NanoMelt直接覆盖。
- **行动**：最终报告建立候选级风险栏，分别列出脱酰胺、Met氧化、二硫键邻域、loop柔性和表面疏水假设；只能写“未触发既定硬规则”，不能写“无风险”。
- **状态**：Open experimental risks；均未达到现行“高置信明确物理风险”硬排门。

#### M-06：AntiFold门通过不等于没有AntiFold负向信号

- **证据**：AntiFold否决要求同时满足`ΔlogP≤−3`和该位置20种氨基酸中最差4名。15个父单突中有8个`ΔlogP≤−3`但因rank大于4而通过；它们至少出现在最终20/30条构建体中。L11Y为−11.545、rank 5，距离否决边界最近。
- **影响**：只报告`pass`会隐藏单条件显著负向；尤其L11Y还依赖AF3-only坐标。
- **行动**：报告解释双条件门，并在候选证据表保留ΔlogP、position rank和来源；AntiFold只作负向排除，不给正向信用。
- **状态**：Open interpretation risk；没有合同违规。

#### M-07：冻结直接界面不等于证明结合保持

- **证据**：最终30没有突变24个直接界面残基，但T99F位于界面壳层，F30位于CDR1且靠近缺失片段；当前V3没有亲和力/结合保持计算。
- **影响**：非界面位点仍可能通过loop构象或远程耦联影响结合。
- **行动**：报告只能写“直接界面位点未突变”，不得写“保证亲和力不变”；表达筛选后至少设置binding/function保持检测。
- **状态**：Open experimental requirement。

#### M-08：关键残基机器文件存在旧界面语义

- **证据**：较早的`nb252_critical_residue_sets.json`把24位界面写为`cautious_not_forbidden`；现行保守性约束已把同一24位完整列入`experimental_interface_frozen`。实际847条和最终30均无界面突变。
- **影响**：未来代码若只读取旧关键事实文件，可能使用错误的突变许可语义。
- **行动**：现行流程以`vhh_conservation_consensus_v2_20260819/nb252_expression_design_constraints.json`为约束权威；后续新增supersession/统一关键事实索引，避免改写历史字节。
- **状态**：Open metadata conflict；当前结果未受影响。

#### M-09：最终15条双突是显式人工冻结面板，不是唯一全局最优解

- **证据**：选择在冻结幅度档、风险、结构证据和组合多样性约束下人工明确确定；同档27条落选包含组合冗余和机制覆盖的人工裁量。原始连续小数不用于档内排名。
- **影响**：结果可重复释放，但不能声称数学上或生物学上“最优15条”。
- **行动**：报告使用“综合证据后选定/优先测试”，并保留102行逐条选择与淘汰理由。
- **状态**：Accepted decision method。

#### M-10：BL21真实表达的若干决定因素未进入当前计算

- **证据**：当前序列合同评价蛋白氨基酸序列，不包含DNA/codon、mRNA二级结构、载体、信号肽、标签、真实N端加工、培养/诱导参数、蛋白酶敏感性及实验聚集测定；也没有计算双突的solvent-exposed hydrophobic/charge patch。
- **影响**：即使U/S/Tm方向有利，也不能排除构建与工艺层面的低产原因。Q1D的真实N端语义尤其依赖表达构建体。
- **行动**：不再追加无验证的同类预测器；在实验设计和报告局限中记录这些未覆盖因素，WT作为30条候选之外的独立对照。
- **状态**：Open experimental-context limitation。

#### M-11：部分上游图的标题或分层语义已不适合V3报告直接复用

- **证据**：上游30单突图仍写`Final 30`和`Three independent ordinal selection metrics`，但该30条现在只是父单突上游短名单，U/S又是同一NetSolP模型的两个输出；双突计划图的53/49和性质矩阵图的pre-selection review strata都不是终选阶段58/44。最终父单突与双突总览图虽清晰，但纵向尺寸较高，整图缩入Word会使标签过小。
- **影响**：直接复制整图可能造成版本、指标独立性或复核深度误读，也可能降低对外材料可读性。
- **行动**：V3报告不直接复用上游单突“Final 30”整图；终选复核深度只用final manifest的58/44。报告已采用V3专用拆图，并在图注中说明U/S是分别评价的同模型输出及相应结构/责任标记。
- **状态**：Closed for the V3 report-only release；不影响图源数据。

### Low

#### L-01：NanoMelt评分域为126 aa，而交付构建为128 aa

- **证据**：WT和全部102条双突均以126-aa numbered domain评分，仅统一裁去末端`GS`；NetSolP使用完整128 aa。
- **影响**：同一工具内WT/突变比较一致，但不同工具的输入域不完全相同。
- **行动**：报告方法和表头注明“预测Tm基于前126 aa；末端GS未进入NanoMelt评分”，不得写成128-aa实测Tm。
- **状态**：Accepted tool contract。

#### L-02：FASTA中的裸突变标签可能脱离reported编号语义

- **证据**：最终CSV记录`reported_positions_1based`，FASTA header使用`L11Y`或`F30S;Q5V`等简写。
- **影响**：若FASTA脱离CSV传播，可能被误读为IMGT编号。
- **行动**：交付README和报告统一说明所有突变标签均指Nb252 128-aa reported-sequence位置，并同时提供CSV。
- **状态**：Open documentation task。

#### L-03：单突V3选择阶段的溯源元数据弱于后续阶段

- **证据**：该阶段有contract、gate、847行审计和30条输出，但没有与后续阶段同等完整的tracked run summary/输入哈希manifest；当前代码仍能从tracked输入精确复现98个共享字段和30条顺序。
- **影响**：不影响结果身份，但使单一阶段的运行环境追溯弱于父单突和双突阶段。
- **行动**：V3报告引用下游已绑定该847审计SHA的parent manifest；不回写历史结果。
- **状态**：Accepted minor provenance limitation。

### Observations

- 历史报告、PPT和交付包属于V2的19单突＋11双突版本，这是预期的版本沿革，不是V3计算链或报告准备缺陷。它们继续原样保留为历史provenance；本目录中的V3项目报告仅读取现行15＋15 final manifest、CSV/FASTA与本审计，并已完成DOCX/PDF绑定。
- `N76G+F30N`在源矩阵中因删除旧N76T并新增F30N的NG基序而净计数抵消。post-sync overlay已经正确保留“新增NG”责任信息；该组合未入选最终15双突。
- 15条入选双突中4条为AF3-only结构距离来源、2条带软责任、0条为局部邻近组合；6条有3个中/强正向性质档，9条有2个。
- 14条含T99F的双突自然分为2条enhanced和12条standard复核，未设T99F专属奖励、配额或否决；0条入选。
- 模型非加和残差仅描述预测器输出的非加和，未用于终选，也不能称为物理epistasis。

## 5. V3报告发布门

### 5.1 报告制品状态

- 现行V3项目报告DOCX与渲染PDF已在本目录生成，并由V3 report manifest绑定文件身份与SHA-256；
- 报告候选身份回读为15条单突＋15条双突，且只使用现行final manifest、final CSV/FASTA和V3审计证据；
- AntiFold在报告中仅作为负向风险排除：不提议、不奖励、不排序候选，双突不计算AntiFold联合分数；
- 历史V2模板字节身份保持不变，历史V2报告、PPT和delivery不参与V3生成；
- 本轮是report-only release，不制作PPT或压缩交付包；二者缺失不是finalization阻断条件。

### 5.2 保留的科学披露与实验要求

1. 当前面板为15单突＋15双突，旧30单突仅为不可改写的上游短名单。
2. 报告逐项披露T99F、A23R、F30缺口边界、AF3-only、F30N脱酰胺和L11M氧化风险。
3. NetSolP、NanoMelt和AntiFold均为预测/相容性证据，不是Nb252突变产量或实测Tm；AntiFold没有正向入选信用。
4. 冻结直接界面不能证明结合保持，实验中仍需binding/function质控。
5. 15条双突未做双侧链结构建模，最终选择不是唯一全局最优解。
6. 四个远程原始评分文件尚未本地归档；当前派生矩阵和本地输出哈希一致，报告保留该溯源限制。
7. WT作为独立实验对照，不计入30条候选。

## 6. 最终审计判定

| 门 | 判定 |
|---|---|
| V3序列身份与组合空间 | Pass |
| 冻结/界面/二硫键/末端约束 | Pass |
| 性质矩阵与幅度档完整性 | Pass |
| 最终CSV/FASTA与manifest | Pass |
| 结构证据完整性 | Pass with AF3-only and no-double-model caveats |
| 预测器对BL21产量的有效性 | Not validated; hypothesis use only |
| 原始远程结果本地归档 | Incomplete：51个manifest绑定中47个本地文件通过哈希，4个远程raw文件未归档 |
| 现有V2报告/PPT/delivery | Expected historical version；不参与V3生成 |
| V3报告DOCX/PDF与终选结果同步 | Pass：manifest已绑定，15＋15身份与AntiFold负向范围检查通过 |
| V3 PPT/压缩交付包 | Out of scope；未制作且不阻断本轮报告finalization |
| V3报告finalization | Pass：report-only release complete |

本审计没有发现需要回滚或重新生成最终30条序列的工程性错误。V3项目报告已完成，后续工作的重点是通过实验验证表达量及结合/功能保持，并保留正向、中性和负向结果的完整溯源。
