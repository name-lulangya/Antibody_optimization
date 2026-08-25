# Codex 项目交接

Last updated: 2026-08-25 20:17:02

Timezone: Asia/Shanghai (UTC+08:00)

> 本文件是当前状态快照，不是追加式历史。现行路线以 Nb252 在 BL21 体系中的表达量优化为唯一设计目标。

## 当前目标与硬约束

- 最终向合作者提供30条Nb252候选，现行最终面板合同为V3：15条单突和15条双突。15条双突只能由最终选定的15条父单突两两组合产生；同位点不同替换互斥，不得组合，也不生成三突或更高阶突变。
- 已释放的`expression_single_mutant_selection_v3_20260825` 30条单突仍是不可改写的V3上游短名单；其后加入用户指定的`T99F`形成31条专家审查池，`v3_parent_single_selection_20260825`冻结15条父单突。102条有效双突已完成共同评分和统一专家审查，15条双突及最终15单突加15双突面板已在`v3_final_15plus15_panel_20260825`释放；该目录manifest是当前候选身份、顺序、理由和gate的唯一权威来源。
- 冻结实验复合物中已复现的 24 个 VHH 界面位点，不再显式优化亲和力，也不在当前路线使用 Rosetta 排序。
- 冻结天然 VHH 邻域中“Nb252亲本残基等于全局/邻域共同优势残基”的高置信保守位点、Cys22/Cys95 和末端 SSGS（reported 125–128）；对高保守但亲本偏离共识的位置只开放共识回变。
- 候选必须保持完整 128-aa 父序列长度、末端 SSGS、两枚原有 Cys，且不得引入新 Cys。
- 现行核心预测工具仅为 NetSolP、NanoMelt 和实验复合物视图 AntiFold。新工具必须先在 47 条可比产量数据上验证，证明有独立且可重复的样本外信息后才能纳入筛选。
- RP3Net 0.0.2 已完成47条正式验证，最终证据等级为`no_supported_use`，不得加入847条候选的生成、筛选或排序。
- PLM_Sol V1.0已完成47条正式验证，最终证据等级为`no_supported_use`，不得加入847条候选的生成、筛选或排序；Nb252单序列smoke分数只证明调用链可运行，不是实测溶解度或产量。

## 最终面板版本边界

- `V3（现行）`：15条单突加15条双突已释放。102条完整128-aa双突全部使用同一规则审查；58条`enhanced`和44条`standard`只代表复核解释深度，不是资格或排序层。终选15条均有两个或三个U/S/Tm中/强改善且无中/强不利，AntiFold仍只保留组成单突负向否决证据。`T99F`组合不设专属触发、配额、奖励或否决，未入选是通用证据与组合约束的结果。
- `V2（历史）`：19条父单突加11条双突，并使用F30/Q1/T27固定配额。其19条父单突、162条历史双突、11条终选双突、代码及结果均保持原样，仅作历史provenance，不参与V3候选生成、评分、排序或交付。
- `V1（历史）`：旧30条单突试选路线，包括旧幅度短名单、受控权衡和稳定词探索例外。其代码和结果保持原样，仅作历史provenance，不参与V3候选生成、评分、排序或交付。
- 版本名指最终面板合同；`expression_single_mutant_selection_v3_20260825`、31条专家审查和`v3_parent_single_selection_20260825`分别是只读上游短名单、父单突专家证据和15条父单突权威来源；最终30条只以`v3_final_15plus15_panel_20260825`为准。

## 权威输入基线

- Nb252 权威设计父本为 collaborator-confirmed 128-aa reported sequence；末端 SSGS 是构建体组成而非 linker。
- 实验结构 VHH 为链 C，NK2R 为链 R；24 个界面残基采用项目已复现的严格 `<4.0 Å` 聚合物重原子距离定义。
- AF3 VHH 是预测结构，只补充结构视图；不得把预测 CDR3 当作实验界面证据。
- 47 条产量记录来自 BL21 表达体系，LTT/WCC 个体近似值可作为数值，LLJ 仍保持分组/删失语义。

## Natural-VHH conservation contract

- 数据源固定为 TNP 论文仓库提交 `a9ba3edc3d967ecf8a2b9b5c2c29bf7495bbc9a0` 的最终 VHH-OAS 描述符表，共 4,059 条非冗余天然 VHH；输入文件 SHA-256 为 `D87A2E66CE0E46D34547D25DF10BF07ABC10B06CCF0D9D0C4A304A36A9D0EBE5`。
- 采用项目已固定的 ANARCII 2.0.8 / IMGT 编号。4,057 条通过 H 链与 framework coverage 审核，2 条因重复编号失败而排除。
- 先按完整 IMGT 域 90% identity 单连接聚类，每簇总权重为 1；再以 framework-only IMGT identity `>=0.80` 且 coverage `>=0.80` 定义 Nb252 邻域。
- 得到 3,784 个去冗余簇；Nb252 邻域含 1,564 条序列、1,532 个有效簇。
- v2将reported 128个位置分为54个`hard_conserved`、1个`conserved_nonconsensus`、33个`cautious`、33个`variable`和7个`insufficient_evidence`。硬保守除要求邻域dominant frequency `>=0.90`、coverage `>=0.80`、有效簇数`>=50`、全局/邻域优势残基一致且全局频率`>=0.80`外，还要求Nb252亲本残基等于该共同优势残基。
- 唯一`conserved_nonconsensus`为reported Q5：邻域/全局优势残基均为V，频率分别为0.9980/0.9976，因此不把Q5亲本状态称为天然硬保守，也不开放任意扫描，只允许`Q5V`共识回变。
- 54个亲本匹配的硬保守位点与界面、Cys22/Cys95、末端SSGS合并后冻结80个reported positions；47个常规可扫描位置产生846条非Cys单突，另加Q5V，共847条。这是待预测的完整约束空间，不是最终30条。
- 权威机器可读合同与结果位于 `docs/result_artifacts/input_baseline/vhh_conservation_consensus_v2_20260819/`；下游必须读取合同，不得从本文复制残基列表。旧`vhh_conservation_20260818`仅保留为历史v1 provenance。
- 下游阶段边界preflight已核对847条候选、48个可变位置、80个硬冻结位置、24个界面位置及CSV/FASTA/亲本一致性；后续三个评分入口复用该通过结果，不再重复同一检查。
- 已完成旧结果的精确回映计划：当前847条中，721条具有可复用NetSolP/NanoMelt及三个AntiFold视图，126条仅缺NetSolP/NanoMelt且AntiFold为AF3-only。缺口严格位于reported 11、14、24、26–29，每个位点18条；Q5V已有三工具旧结果。
- 已生成带IMGT FR/CDR标注的全局天然VHH、Nb252邻域及项目表达序列Logo。项目Logo以47条源序列为审核范围，仅纳入45条编号成功的H链序列；编号失败和非H链各1条保持显式排除，且产量不作为频率权重。

## 当前工具证据

- NetSolP：保留 S（solubility）与 U（usability）原始连续值；在 47 条产量数据中仅显示有限关联，不能单独决定候选。
- NanoMelt：43/47条可评分，其中27条具有LTT/WCC个体数值产量。嵌套逐样本/序列簇留一结果相同：ROC-AUC 0.571、PR-AUC 0.554、MCC 0.408、balanced accuracy 0.701、sensitivity 0.786、specificity 0.615；连续与分类联合证据仍为`no_supported_use`，不得作为BL21产量排序器，只保留为预测稳定性约束。
- AntiFold：47条产量序列中只有Nb252具有匹配实验NK2R复合物，另外46条没有可比实验复合物，而且输出是逐位点结构条件概率而不是统一表达分数；正式分类状态为`not_applicable`，不得报告伪AUC/MCC。V3不给AntiFold改善任何正向选择权重，只在`ΔlogP<=-3`且突变残基位于该位置20种氨基酸最差4名时否决；实验复合物视图优先，实验坐标不可评价位置保留AF3-only来源标签。
- TNP 与 nanoBERT 已完成探索性验证，但不在现行精简筛选工具集中。
- RP3Net：31条数值记录的直接合并Spearman为0.476，但来源内分层Spearman仅0.198且95% bootstrap区间跨0；LLJ有序Kendall为-0.451，与预声明方向相反。嵌套分类ROC-AUC为0.621、PR-AUC为0.620、MCC为0.313，未达到预声明综合门，因此不支持候选使用。
- PLM_Sol：47/47条评分成功。31条数值记录的来源内分层Spearman为0.473，95% bootstrap区间为0.096–0.749，但WCC内部Spearman为-0.096；嵌套分类ROC-AUC为0.638、PR-AUC为0.662、MCC为0.313。PLM_Sol与NetSolP U/S高度重叠，在NetSolP S基础上的序列簇外增量为-0.140，因此gate为`plm_sol_not_supported_for_candidate_use`，不纳入候选使用。固定5 mg结果仍仅供展示。
- 固定5 mg探索图：31条数值记录显示为高产14条、低产17条；RP3Net、NetSolP U和NetSolP S均展示训练折最大MCC阈值及留出指标。该图仅用于直观展示，不作为工具准入、候选筛选或5 mg阈值有效性的正式证据。
- 所有工具验证须同时报告连续关联和预注册方向下的离散分类性能；阈值必须在训练折内选择，外层交叉验证报告 ROC-AUC、PR-AUC、MCC、balanced accuracy、sensitivity、specificity 及阈值稳定性。

## 下一执行路线

1. 三工具全空间补全已经完成，不在输入和工具合同不变时重复运行。重复一致性门共183项、0失败；NetSolP和NanoMelt逐值一致，AntiFold最大绝对差为`5.72e-6`并通过固定`1e-5`容差。
2. 847行完整原始评分矩阵已经释放：721条NetSolP/NanoMelt沿用经重复验证的旧值，126条缺口采用本轮新算值；AntiFold为721条三视图和126条AF3-only证据。实验缺失位置的实验视图保持`not_evaluable`。
3. 1,336条简并稳定词已按固定12符号、大小写敏感、允许重叠的精确连续子串合同加入847行矩阵。Nb252 WT没有稳定词命中；22条单突新增24个命中，825条不变，无减少项。该步骤没有筛选、排名或生成Tier。
4. 47条BL21产量验证不支持“稳定词密度越高、产量越高”：主指标来源分层Spearman为-0.185，95% bootstrap区间为-0.572–0.216，序列簇留一ROC-AUC为0.431、MCC为-0.277，证据等级为`no_supported_use`。稳定词只作为用户指定的可解释软偏好，不能覆盖硬约束、明确性质恶化或AntiFold/NanoMelt的既定用途。
5. 847条单突的四指标全景图已经生成：四张位置×替换残基热图分别展示NetSolP ΔU、NetSolP ΔS、NanoMelt预测ΔTm和AntiFold ΔlogP，并另有ΔU–ΔS及AntiFold–ΔTm两张散点图。AntiFold优先采用721条实验复合物视图结果，126条实验坐标不可评价候选以独立AF3 VHH-only结果补充；来源在逐候选表、图例和散点点形中显式区分。22条稳定词新增候选均用星号标记。
6. V3把NetSolP U、NetSolP S和NanoMelt预测Tm作为三个独立正向指标，仍使用冻结幅度档且不以档内小数排序。常规候选至少有一个指标达到中等/明显改善，三个性质指标不得明显恶化；稳定词通常只作最后平局条件。唯一例外是用户明确指定的`T99F`稳定词假设探索项，其上游不合格事实和弱不利性质证据必须同时保留。
7. AntiFold位置内排名由WT与19种替换的完整20状态分数计算，包括不允许进入候选的新Cys状态。联合否决规则排除151/847条，其中实验复合物来源130条、AF3-only补充来源21条；AntiFold改善不参与分层、排序或平局处理。
8. 性质、AntiFold联合门和硬序列风险后得到61条V3合格单突：A层8条、B层11条、C层42条，无D层。上游30条单突短名单按位置第1轮、第2轮、第3轮依次覆盖，每个位点最多3条；轮内按A>B>C>D、明显改善指标数、改善指标数、较少中等恶化、较少软风险、稳定词平局条件和候选ID排序，不使用连续加权总分。
9. V3上游30条单突短名单由A/B/C层6/7/17条组成，覆盖23个reported位置，实际每个位点最多2条；另有31条合格替补。30条均为唯一128-aa单突、无AntiFold联合否决、无硬序列风险且至少一个U/S/Tm指标达到中等或明显改善。精确审计、合同、CSV/FASTA和图以`expression_single_mutant_selection_v3_20260825`机器制品为唯一权威来源。
10. 31条父单突审查池的结构/VHH专家审查已经完成：27条以实验复合物为主证据，L11Y、L11M、I28Y、F29Q因实验坐标缺失仅作AF3单体低置信度判断；31条均有独立ChimeraX单rotamer视图和逐条结构、溶解度、热稳定性解释。补充的`T99F`判断为`structurally_concerning`，预测溶解度和热稳定性方向均为不利，置信度中等；该结论是审查意见，不是选择或实验结果。
11. 15条父单突已经确定：其中9条至少有一个U/S/Tm强改善档、5条为中等性质证据补充、`T99F`为唯一稳定词探索例外。专家意见只对5条高置信明确物理风险执行硬排；其余中低置信担忧仅作标注或降级。
12. 15条父单突的102条有效双突已完成NetSolP/NanoMelt共同评分和统一专家审查。复核深度按通用三触发并集得到58条`enhanced`和44条`standard`，只控制解释详略；14条含`T99F`组合自然分为2/12，没有专属双突规则。`N76G+F30N`新增reported 30 `NG`脱酰胺基序已用overlay记录，未改写正式源矩阵。
13. 最终15条双突已显式冻结：6条有三个中/强正向指标，9条有两个，均无中/强不利；覆盖13/15个父单突和10/12个reported位置，4条为AF3-only结构距离来源、2条带软序列风险、无局部邻近组合。最终30条CSV/FASTA、102行审计和选择合同均以`v3_final_15plus15_panel_20260825` manifest为准，展示顺序不是效力排名。
14. V1旧30单突试选、V2的19条父单突、162条双突及19+11面板均只保留为历史计算provenance；不修改这些既有结果，也不把其候选、配额、分数或排名用于V3。
15. V3计算链与报告就绪性审计已经完成，机器检查18/18通过，未发现序列错配、冻结/界面违规或最终面板身份错误；结论为`pass_with_material_caveats`。权威审计及证据位于`docs/result_artifacts/weekly_report_result/Nb252_V3_expression_report/`。
16. 下一步在上述目录新建面向导师/合作者的V3项目报告、PPT和交付材料，使其只读取V3终选manifest与最终30条CSV/FASTA并纳入审计披露；现有19＋11的V2材料是预期的只读历史版本，不构成V3审计问题，也不参与V3材料生成。

## 阶段门

- `structure_and_interface_identity=pass`
- `natural_vhh_conservation_contract=pass`
- `expression_single_mutant_constraint_space=pass`：847 条仅表示允许进入预测的单突空间。
- `predictor_continuous_and_classification_validation=pass_with_tool_specific_roles`：NanoMelt正式判为仅稳定性约束；AntiFold因结构覆盖不足正式判为分类不适用、仅实验复合物相容性约束；不得把该门解释为产量预测器已验证。
- `expression_property_completion_plan_v2=pass`：721条精确复用、126条补算和12条抽样复核计划已冻结。
- `expression_property_complete_matrix_v2=pass`：847条NetSolP、NanoMelt和AntiFold证据已按视图可评价范围完成并通过完整性审核；该门不代表候选已筛选。
- `stable_word_feature_evaluation_v1=pass`：847条单突已完成稳定词新增/减少审计，47条产量验证也已完成；`pass`只表示制品完整，不表示稳定词已被验证为产量预测器。
- `expression_single_mutant_four_metric_landscape_v1=pass`：847条、48个位置、721条实验复合物AntiFold值和126条AF3补充值均已准确绘图并保留来源；该门不执行候选选择。
- `expression_single_mutant_trial_selection_v2=pass`：40条幅度短名单、26条严格核心、10条受控权衡、30条试选和7条替补已生成；试选含1条用户指定的`T99F`稳定词探索例外，原始同档小数没有参与排序。
- `expression_single_mutant_selection_v3=pass`：847条完整审计、151条AntiFold联合否决、61条合格池、30条上游单突短名单和31条替补均已释放；该门是现行V3父单突选择的来源，不再表示最终交付面板。
- `final_panel_contract_v3=pass`：现行组成已冻结为15条单突加15条双突，V1/V2明确为只读历史版本。
- `v3_parent_single_expert_review=pass`：31条审查池均已完成来源分层的结构量化、独立单rotamer可视化和VHH专家判断；原始30条上游短名单未改写，补充`T99F`仍保持上游不合格/未选状态，该门不执行父单突选择。
- `v3_parent_single_selection=pass`：15条父单突已冻结，覆盖12个位点；31条均有详细选择或淘汰理由，`T99F`保持上游不合格事实并以下游探索例外入选。未来105个理论无序对中排除3个同位点对，得到102个有效双突。
- `v3_double_mutant_plan=pass`：105个理论无序对已完整审计，3个同位点互斥对被显式排除，102条唯一完整双突及103条共同评分样本已释放；没有在评分前按性质删减候选。
- `v3_double_mutant_property_matrix=pass`：NetSolP与NanoMelt均为103/103覆盖，102条双突性质矩阵和幅度档已释放；该门不执行候选选择。
- `v3_double_expert_review=pass`：102条均有共同审查字段；58/44只表示复核深度，未作为资格、排序或淘汰规则，`T99F`没有双突专属处理。
- `final_panel_v3_double_selection=pass`：15条双突已冻结且通过性质、组成AntiFold、序列风险、结构来源和组合多样性核验。
- `final_panel_v3_release=pass`：最终15条单突加15条双突的CSV/FASTA、102行决策审计、manifest和图已释放；报告/PPT同步仍为`not_performed`。
- `v3_report_audit=pass_with_material_caveats`：18/18机器完整性检查通过；当前面板可进入V3报告书写。现有V2报告/PPT/delivery按预期保留为历史材料；4个远程raw评分文件未本地归档，且T99F、AF3-only、软责任与未做双突侧链建模等边界必须在新报告披露。

## 本轮验证状态

- V3已回读847条完整审计、61条合格池、30条上游单突短名单和31条替补；AntiFold联合否决精确为151条，全部位置排名均基于20种氨基酸状态。上游30条覆盖23个位点，A/B/C为6/7/17条，实际每个位点最多2条，无稳定词单独驱动的入选项。
- V3本地真实数据运行低于1分钟，无需Slurm、checkpoint或resume；专项测试覆盖真实计数、联合否决双条件、20状态完整性、CLI输出和上游字节不变。选择图已人工检查，漏斗、位置覆盖和三独立指标热图均可辨识。
- V3专家审查覆盖31条/23个位点；27条实验复合物主证据、4条AF3-only低置信度证据。结构制品共63张图：31张主候选、7张AF3敏感性、23张WT位点近景和2张来源总览。该上游评价制品仍如实保持`parent_single_selection=not_performed`；下游`v3_parent_single_selection_20260825`已独立完成15条选择，不回写专家审查原件。
- V3父单突选择制品已回读为31条决策、15条父单突和15条FASTA序列；12个位点、3个同位点非法对和102个有效双突计数一致。选择总览图使用冻结幅度档，弱变化统一显示为中性，已人工检查无标题、图例或色标遮挡。
- V3双突结果已回读为102条唯一128-aa双突；两工具均对WT加102条共同表103/103通过，NanoMelt全部评分前126 aa并仅裁末端`GS`。终选制品已回读为102条完整审计、15条双突和30条唯一128-aa最终序列；组成/位置覆盖、AF3-only、软风险及局部邻近上限均通过，责任基序勘误以overlay应用且未改写源矩阵。
- 现有项目报告和PPT描述V2 19+11历史面板，按预期原样保留；新的V3 15+15报告、PPT和交付材料将在V3报告目录中单独生成。
- V3审计独立回接847→30→31→15→102→15→最终30，确认CSV/FASTA、突变重建、硬约束和47个本地manifest哈希绑定均一致；审计报告已标记Complete。当前无须回滚候选，V3报告需继续披露远程raw缺档和科学适用边界。
- RP3Net正式运行覆盖47/47条序列；连续、分类、逐折预测、结果图和gate均已生成并完成schema及计数核验。
- 4,059条输入、4,057条合格序列、1,564条邻域序列、128个位置和847条单突均已回读核对；v2相对v1只新增Q5V，未删除或意外开放其他候选。
- v2下游preflight为pass：847条候选中846条来自常规非Cys扫描、1条为Q5V共识回变，无多突、新Cys、冻结或界面突变。
- v2性质完整矩阵为pass：847条ID和序列均唯一且均为合法128-aa单突；721条复用值与126条新算值来源明确，全部保留末端SSGS，无冻结位点突变或新增Cys。NetSolP、NanoMelt和AF3 AntiFold均覆盖847条；实验单体/复合物AntiFold仅覆盖坐标可评价的721条，其余126条明确为`not_evaluable`。
- 稳定词输入共1,336条且均符合固定12符号合同；847条单突与完整性质矩阵逐ID联接，24条长表变化均覆盖实际突变位点。47条验证保留31条个体数值和16条LLJ有序/删失语义，未训练高容量模型，也未执行候选选择。
- 四指标景观图的数据表和gate已回读，847条候选、48个位置、22条稳定词新增、721条实验复合物AntiFold值与126条AF3补充值计数一致；热图和散点PNG均已人工检查，位置标签、FR/CDR标记、独立色标、星号及AntiFold来源点形均可辨识。
- 试选流程逐条回读847条审计、40条短名单、30条试选和7条替补；30条均为唯一128-aa单突且无明显不利档位和声明的硬序列风险。29条至少有一个中等/明显有利工具家族；唯一例外`T99F`新增一个长度5的稳定词且四指标仅为中性/微弱负向。试选包含24条实验复合物AntiFold和6条AF3补充来源，覆盖13个位置；PNG已人工检查，`T99F`以金色星号和边框单独标识。
- NanoMelt分类制品覆盖27条数值样本并保留4条未评分状态；正式图已人工检查。AntiFold分类不适用合同确认1条匹配结构、46条结构缺失且不生成任何分类指标。
- RP3Net gate固定为`rp3net_not_supported_for_candidate_use`；其模型分数不能解释为mg/L或通用表达阈值。
- PLM_Sol正式运行覆盖47/47条序列；31条数值记录和16条LLJ有序/删失记录语义保持不变，结果图、gate和run summary已回读核对。其固定5 mg展示的序列簇留一ROC-AUC/PR-AUC/MCC为0.761/0.685/0.411，但该展示不参与工具准入。
