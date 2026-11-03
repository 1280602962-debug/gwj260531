# 英译前逐处文字与引用修改清单

基准提交：`d1ca17512eabe5c9d102ebe75c3f9c59668452dc`。五个既有 Markdown 文件有文字变动。以下逐项记录修改操作；同一文件中的来源 M 显示名最终统一为 S1–S5，复合内部字段保持原样。科学数据、图表文件及数值未改。

## submission/Supplementary_Methods_CN.md

| 原文 | 修改后 | 理由 |
| --- | --- | --- |
| 主要读数为 seed 42 的 mode 1，并取原始评分相反数。 | S1 为 seed 42 的 Vina mode 1 评分，并取原始评分相反数。 | 明确 S1 定义 |
| 固定姿态重评分使用 score_only、vina 经验评分与 CNNscore；S2 为 Vina mode 1，S4 为同次保存姿态中的最高 CNNscore。 | 固定姿态重评分使用 score_only 和 vina 经验评分设置；S2 为对 S1 所用 Vina mode 1 姿态计算的 GNINA CNNscore，S4 为同次 Vina 对接全部保存姿态中的最高 GNINA CNNscore。 | 纠正编号与姿态/分数歧义 |
| 独立 GNINA 对接使用 seed 42、exhaustiveness=16、num_modes=9，以 CNNscore 排序，读取第一姿态。 | S5 为 GNINA 独立对接，使用 seed 42、exhaustiveness=16、num_modes=9，以 CNNscore 排序，读取第一姿态的 CNNscore。 | 完整写明 S5 读数 |
| RTMScore 使用 `rtmscore_model1.pth`，固定 Vina mode 1，口袋截取距离为 10 Å，不进行模型集成。 | S3 为 RTMScore，使用 `rtmscore_model1.pth` 对 S1 所用 Vina mode 1 姿态评分，口袋截取距离为 10 Å，不进行模型集成。 | 完整写明 S3 读数 |
| ## 5. 统计、人群与联合排序 | GNINA 原始记录中另保存了 CNNaffinity，但未将其纳入本研究的正式评价。S2、S4 和 S5 的评价读数均为 CNNscore。 /  / ## 5. 统计、人群与联合排序 | 原始记录及正式评分读取函数支持该说明 |
| summary_min 在每次重采样中取较小方向。 | D_vs_A 使用靶点 B 的评分比较 dual 与 A-only，D_vs_B 使用靶点 A 的评分比较 dual 与 B-only。summary_min 的点估计取两个方向 AUROC 的较小值；在每次重采样内部也取两个方向 AUROC 的较小值，再由该重采样分布给出区间。 | 明确方向记号和汇总指标，未改统计 |
| 本次表格整理不拟合模型。 | 本次文字与 SI 修订不拟合模型。配体化学组成文件中的 `max_ecfp4_tanimoto_other_class` 为同一靶点对、同一方向内相对于另一活性类别的最大 ECFP4 Tanimoto 相似性（Figure S6；Supplementary Data）。该字段描述标签层的跨类别相似性，不表示测试折配体相对于训练折的最近邻相似性。 | 区分已有相似性与缺失的跨折分析 |
| S2、S4 和 S5 的评价读数均为 CNNscore。 | S2、S4 和 S5 的评价读数均为 CNNscore。五种评分方案的关系见 Scheme S1，正式编号见 `SCORING_SCHEMES.csv`。 | 评分工作流引用 |

## manuscript/Results_CN.md

| 原文 | 修改后 | 理由 |
| --- | --- | --- |
| 阈值、合格性和 alias 规则见 Scheme 1 与补充方法。 | 阈值、合格性和 alias 规则见补充方法第 1 节。 | 旧研究设计 Scheme 与现行评分 Scheme 不同，改引对应文字 |
| 每臂人数和缺失名单见 Table S2。 | 每臂人数见 Table S4A，缺失名单见 Supplementary Data。 | 现行 S2 为受体表 |
| 完整来源记录见 Tables S19–S21。 | 来源覆盖汇总见 Table S1B，完整来源对照及逐条活性记录见 Supplementary Data。 | 旧来源长表移入数据包 |
| Figure 1a 显示 | Figure 2a,b 显示 | 方向与汇总分面对应现行 Figure 2 |
| F2/F10 的两个方向都低，分别为 | F2/F10 的两个方向 AUROC 点估计均低于 0.5，分别为 | 明确点估计，不暗示两方向区间均排除 0.5 |
| Figure 1b 使用 | Figure 3 使用 | 评分方法配对差现为 Figure 3 |
| 五方法的自身结果见 Table S3 和 Figure S4，配对差值见 Table S4。 | 五方法的自身结果见 Table S5A 和 Figure S3，配对差值见 Table S5B。 | 按投稿目录映射 |
| Figure 2 在同一 dual 集合上 | Figure 4 在同一 dual 集合上 | 负类图改号 |
| 完整人数和区间见 Table S5。 | 完整人数和差值区间见 Table S6A，各条件 AUROC 区间见 Supplementary Data。 | 当前表含差值区间，单条件区间在长表 |
| Figure 3 的三个化学模型 | Figure 5 的三个化学模型 | 化学图改号 |
| Figure 3 没有区间 | Figure 5 没有区间 | 化学图改号 |
| 描述符和异类相似度只描述第 1 层组成，见 Figures S1–S2、Figure S3 和 Tables S7–S9。 | 描述符和异类相似度只描述第 1 层组成，见 Figures S4–S6 及 Supplementary Data 中的 `module_D1_composition.csv`；模型表现与折外预测分别见 Table S7 和 `module_D_chemistry_oof_ligands.csv`。 | 区分模型表、组成图和 OOF 长表 |
| Figure 4 在同一组成员上比较规定口袋和对侧口袋。 | Figure 6 在同一组成员上比较对应口袋和非对应口袋。 | 更新图号及与 Methods 一致的口袋用语 |
| 完整成员见 Table S6。 | 结果汇总见 Table S6B，完整成员见 Supplementary Data。 | 成员名单在长表，不在汇总面板 |
| Figure 5a 只显示 | Figure 7a 只显示 | 种子图改号 |
| 两个方向的完整结果在 Figure S5，汇总在 Tables S10–S11。 | 两个方向的完整结果在 Figure S7，汇总在 Table S8A，逐种子记录见 Supplementary Data。 | 种子图表改号 |
| 见 Figure 5b、Figure S6 和 Table S12。 | 见 Figure 7b 和 Table S8B。 | 旧替代受体补充图已并入正文图与表 |
| Figure 6 的每个面板 | Figure 9 的每个面板 | 联合排序图改号 |
| 方法自身的 M0 名单和四态计数见 Figure S8 与 Table S17。 | 方法自身的 M0 名单和四态计数见 Table S10A,C；方法间配对变化见 Table S10B。 | 当前 Figure S9 是共同名单，不能替代方法自身组成 |
| 完整点估计、区间、人数及名单变化见 Tables S23–S31。 | 人数、Vina 汇总及化学模型见 Table S9A,B,D；完整方法结果、区间、排除成员及名单变化见 Supplementary Data 中的 CAP5 文件。 | CAP5 旧长表已移入数据包 |
| 显示名 M0/M1/M1b/M2/M3 | 显示名 S1/S2/S4/S3/S5 | 与 SCORING_SCHEMES.csv 对齐；复合内部字段不变 |
| neither 111。dual 可以进入两个方向 | neither 111（Figure 1a）。dual 可以进入两个方向 | 为类别构成补图引用 |
| 每臂人数见 Table S4A，缺失名单见 Supplementary Data。 | 每臂人数见 Table S4A，评分可用性见 Table S4B 和 Figure S2，缺失名单见 Supplementary Data。 | 缺失图与表引用 |
| 来源覆盖汇总见 Table S1B， | 来源覆盖汇总见 Figure 1b 和 Table S1B， | 来源图引用 |
| ## R2. 双方向表现与评分工作流差异 /  | ## R2. 双方向表现与评分工作流差异 /  / 共晶配体重对接的 RMSD 结果见 Table S3 和 Figure S1。 /  | 指向已有 RMSD 证据 |
| D_vs_A 的正确口袋是 B，D_vs_B 的正确口袋是 A。差值为 correct−wrong。 | D_vs_A 的对应口袋是 B，D_vs_B 的对应口袋是 A。差值为对应口袋−非对应口袋。 | 统一口袋术语 |
| 16 个臂的错口袋额外减员都是 0。 | 16 个臂的非对应口袋额外减员都是 0。 | 统一口袋术语 |
| 原因是 S1 超时，不是对侧口袋再删人。 | 原因是 S1 超时，不是非对应口袋再删人。 | 统一口袋术语 |
| 各比较的 N、k 和 dual 数见 Table 2； | 各比较的 N、k 和 dual 数见 Table 2，共同名单的类别组成另见 Figure S9； | 共同名单组成图引用 |
| 保留 861 条、排除 6 条。排除发生在 | 保留 861 条、排除 6 条（Figure 8a）。排除发生在 | CAP5 成员图引用 |
| 人数、Vina 汇总及化学模型见 Table S9A,B,D； | 人数、Vina 汇总及化学模型见 Table S9A,B,D，全部评分的方向性比较另见 Figure S8； | CAP5 全评分图引用 |
| 完整方法结果、区间、排除成员及名单变化见 Supplementary Data 中的 CAP5 文件。 | 完整方法结果、区间、排除成员及名单变化见 Supplementary Data 中的 CAP5 文件。PIK3CA/mTOR 的严格活性子集结果见 Figure 8b 和 Table S9C。 | 为已有严格子集表图补引用 |

## manuscript/Supporting_Information_CN.md

| 原文 | 修改后 | 理由 |
| --- | --- | --- |
| 当前来源总表见 Table S20。101 个 Q110 分子的 238 条 ChEMBL 37 原始活性见 Table S21； | 来源覆盖汇总见 Table S1B，当前来源总表见 Supplementary Data 中的 `CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv`。101 个 Q110 分子的 238 条 ChEMBL 37 原始活性见数据文件 `Q110_CHEMBL37_ACTIVITY_ROWS_LONG.csv`； | 旧来源长表已移入数据包 |
| Table S2 给出每臂人数和名单。Table S13 中 expected = valid + missing。 | Table S4A 给出每臂人数，名单见 Supplementary Data。Table S4B 及其类别明细中 expected = valid + missing。 | 现行人群及缺失表位置 |
| 最大异类 ECFP4 Tanimoto 只在同一 pair×arm 内相对另一标签类计算。 | 最大异类 ECFP4 Tanimoto 只在同一 pair×arm 内相对另一标签类计算，见 Figure S6 及 Supplementary Data 中的 `module_D1_composition.csv`。该结果不表示测试折配体相对于训练折的最近邻相似性。 | 避免将现有跨类别分析当成跨折分析 |
| Table S5 和 Figure 2 给出 16 个臂的单靶负类、neither、差值和区间。Table S6 和 Figure 4 给出正确口袋、对侧口袋、差值和额外减员。 | Table S6A 和 Figure 4 给出 16 个臂的单靶负类、neither、差值和差值区间，各条件 AUROC 区间见补充数据。Table S6B 和 Figure 6 给出对应口袋、非对应口袋、差值和额外减员。 | 现行图表位置及区间范围 |
| Table S10 为 40 个 pair×seed 原点，Table S11 为 24 个汇总。Figure S5 同时显示 auc_B、auc_A 和 summary_min。Figure 5a 只显示 summary_min。 | Supplementary Data 中的 `module_E_five_seed.csv` 保留 40 个 pair×seed 原点，Table S8A 为 24 个汇总。Figure S7 同时显示 auc_B、auc_A 和 summary_min。Figure 7a 只显示 summary_min。 | 种子原点、汇总和图分别引用 |
| Table S12 和 Figure S6 覆盖 | Table S8B 和 Figure 7b 覆盖 | 替代受体补充图已并入现行正文图 |
| Table S16 保留全部 9734 行，不在本 PDF 中逐行排印。 | Supplementary Data 中的 `module_I_ligand_ranks.csv` 保留全部 9734 行，不在正文中逐行排印。 | 旧排名长表不再单独编号 |
| Table S17 的四行与 Table S15 的 32 行计数一致。 | Table S10C 的四行与 Table S10B 的 32 行计数一致。 | 四态及方法配对表面板映射 |
| Supplementary_Source_Data.zip 收录当前八对的主结果、人群、分折、来源对照、逐条活性证据，以及 CAP5 的保留名单、统计、折外预测和排名。 | 投稿数据包 `submission/03_Supplementary_Data/Supplementary_Data.zip` 收录当前八对的主结果、人群、分折、来源对照、逐条活性证据，以及 CAP5 的保留名单、统计、折外预测和排名。 | 明确现行投稿数据包 |
| 来源与支持分析的逐表对应见 Source_Register.csv。 | 投稿表格与来源对应见 `submission/TABLE_SOURCES.csv`，数据包内另有 `Data_Index.csv` 和 `Numbered_tables/`。 | 旧来源登记仍为历史编号，改引现行索引 |
| 完整结果见 Table S18。 | 完整结果见 Table S9C。 | 严格子集表改号 |
| 七对靶点的沉积名单、脚本记录及当前身份核对见 Table S22。 | 七对靶点的选样证据汇总见 Table S1A，沉积名单、脚本记录及当前身份核对见 Supplementary Data 中的选样来源文件。 | 选样明细不再独立编号 |
| 见 Tables S23–S24。 | 见 Table S9A，以及 Supplementary Data 中的 `CAP5_CLASS_COUNTS.csv` 和 `CAP5_EXCLUDED_LIGANDS.csv`。 | CAP5 人数汇总和明细分别引用 |
| 完整统计和与主分析的对照见 Tables S25–S31。 | Vina 汇总与化学模型见 Table S9B,D；全部评分的方向性比较见 Figure S8，完整统计和与主分析的对照见 Supplementary Data 中的 CAP5 文件。 | 现行面板不涵盖全部旧长表 |
| ## 5. 负类与口袋的完整结果 | GNINA 原始记录中另保存了 CNNaffinity，但未将其纳入正式评价；正式的固定姿态、保存姿态重排及独立对接方案均使用 CNNscore，见补充方法第 4 节及 `submission/SCORING_SCHEMES.csv`。 /  / ## 5. 负类与口袋的完整结果 | 同步已核实的评分输出范围 |
| 显示名 M0/M1/M1b/M2/M3 | 显示名 S1/S2/S4/S3/S5 | 与 SCORING_SCHEMES.csv 对齐；复合内部字段不变 |
| 新方法使用 pairwise_common_vs_M0，S1 使用 pairwise_common_vs_该方法， | 新方法使用 `pairwise_common_vs_M0`，S1 使用 `pairwise_common_vs_该方法`（该方法的内部标识仍采用来源 M 名称）， | 保留真实源字段，说明内部映射 |
| 完整结果见 Table S9C。 | 完整结果见 Figure 8b 和 Table S9C。 | 严格子集图引用 |
| 包内 SHA256 清单只用于校验文件字节，包内不含 ChEMBL 数据库原库。 | 包内 `Data_Index.csv` 列出文件路径、字节数和数据行数，包内不含 ChEMBL 数据库原库。 | 现行 Data_Index.csv 表头为 path_in_archive,bytes,data_rows；不误称为 SHA256 清单 |

## manuscript/README.md

| 原文 | 修改后 | 理由 |
| --- | --- | --- |
| 本目录中的既有文稿、图版、审阅 PDF 与工作簿保留来源版本。其中旧工作簿使用 S1–S31，图版使用 M0、M1、M1b、M2、M3；这些工作标签尚未在完整文稿和图版中替换。准备最终稿时按投稿目录中的编号对照和评分名称映射统一。 | 本目录的 `Results_CN.md` 和 `Supporting_Information_CN.md` 已于 2026-11-03 按投稿目录更新图表引用和评分显示名。旧工作簿、`Figure_Captions_CN.md`、`Source_Register.csv` 和审阅 PDF 保留来源版本，仍使用旧表号及 M 名称，不能与现行投稿图表混用。当前图注、表注分别使用 `submission/FIGURE_CAPTIONS.md` 和 `submission/TABLE_CAPTIONS.md`。 | 明确活跃文本与历史来源的边界 |
| 本目录的完整 Results 工作稿仍由作者随后修改。 | 本目录 Results 工作稿的引用已对齐；本轮未重写科学结果或替换审阅 PDF。 | 更新已执行状态 |

## submission/README.md

| 原文 | 修改后 | 理由 |
| --- | --- | --- |
| 原 `results/jcim_stage2_figures/` 与 `manuscript/` 图文保留为工作来源，投稿图号以本目录为准。Results 原稿由作者后续修改，本次提供逐节引用对照与核对意见，未将原稿替换为新的完整 Results。 | 原 `results/jcim_stage2_figures/`、旧工作簿、旧图注及审阅 PDF 保留为来源版本，投稿图号以本目录为准。2026-11-03 已将 `manuscript/Results_CN.md` 和 `Supporting_Information_CN.md` 的图表引用及评分显示名对齐本目录，科学结果未作重写。 /  / 英译前文字核查见 [PRETRANSLATION_METHODS_SI_AUDIT_20261103_CN.md](04_Editorial_Checks/PRETRANSLATION_METHODS_SI_AUDIT_20261103_CN.md)。补充方法已修正 S1–S5 并说明 CNNaffinity 未纳入评价。跨折最近邻原始结果和脚本未取得，因此 Table S7 未新增 Panel B。当前材料不含可核查的 Discussion 正文，该部分仍未完成逐句核查。 | 更新执行状态，明确不能声称全部核查通过 |

## 未执行的修改

- Table S7 未新增 Panel B：缺少跨折最近邻逐成员结果、分布表与原脚本。
- 未改 Discussion：当前正文未取得。已读中文稿没有 Discussion 节，不能以大纲或历史文字代替。
- 共晶二维相似性未计算，也未写入新结论；取得的活跃文本中未发现待删除句。共晶重对接 RMSD 保持原样。
- 原始评分、标签、成员、分折、OOF、AUROC、CI、表格单元格及图片全部保持原样。
