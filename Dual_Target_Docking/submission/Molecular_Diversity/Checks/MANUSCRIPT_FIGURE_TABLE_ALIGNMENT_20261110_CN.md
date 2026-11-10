# 正文与补充图表最终排序及待翻译前交叉核查（2026-11-10）


**Table Sx关闭更新（2026-11-10）**：探索性单靶对照已补充为Table S11A/B。下文原审核时记录的该占位缺口已关闭；原始bootstrap输出未恢复，本次以冻结S1评分复现并归档区间。详见`../../../analysis/single_target_controls_20261110/README.md`及`../Tables/TableS11.md`。

## 核查范围和证据边界

基线：cursor/v4-2-final-audited-results-20260930 中 Molecular_Diversity 投稿包；正文：用户提供的《双靶对接第二版中文6.docx》。此次只进行编号、文本来源和现有文件清点，不修改冻结评分、置信区间、原始CSV、图形或表格数值。已有 Release_QA.json 报告1566项数值检查、3764项表格源单元格检查、159个科学输入文件未改动、共晶重对接70条（首位42、保存最佳66）；这些属于仓库既有QA记录，不代表此次重新独立计算或逐页视觉复核。

## 正文图表最终顺序

| 顺序 | 正式图表 | 内容 | 正文对应位置 | 图表来源 |
| --- | --- | --- | --- | --- |
| 1 | Figure 1 | 八对靶点四类人数及活性来源覆盖 | Results 3.1 | Main_Figures/Fig1 |
| 2 | Table 1 | 八对靶点、主要受体PDB及四类人数 | Results 3.1 | Tables/Table1.md |
| 3 | Figure 2 | Vina两方向AUROC及summary_min区间 | Results 3.2 | Main_Figures/Fig2 |
| 4 | Figure 3 | S2−S1、S3−S1、S4−S2、S5−S1配对差值矩阵 | Results 3.2 | Main_Figures/Fig3 |
| 5 | Figure 4 | 单靶负类与neither负类比较 | Results 3.3 | Main_Figures/Fig4 |
| 6 | Figure 5 | 理化、ECFP4、联合Vina及增量 | Results 3.4 | Main_Figures/Fig5 |
| 7 | Figure 6 | 对应/非对应口袋配对AUROC | Results 3.5 | Main_Figures/Fig6 |
| 8 | Figure 7 | 随机种子及替代受体 | Results 3.6 | Main_Figures/Fig7 |
| 9 | Figure 8 | 方向性AUROC变化与前10%单靶配体比例变化 | Results 3.7 | Main_Figures/Fig8 |
| 10 | Table 2 | 32组配对Top-10% dual入选数及召回率 | Results 3.7 | Tables/Table2.md |

正文八图与两表合计10项；Molecular Diversity Research Article当前指南上限为10个Tables plus Figures。此处的编号依据当前投稿目录，不使用旧提交中的九图方案。

## 补充图及示意图的正式顺序

| 正式编号 | 内容 | 推荐首次引用 |
| --- | --- | --- |
| Figure S1 | 14受体×5种子重对接RMSD | Methods 2.3 / Results 3.2 |
| Figure S2 | 各类评分有效数/应有数 | Results 3.1 |
| Figure S3 | 五种评分完整AUROC概览 | Results 3.2 |
| Figure S4 | dual与A-only的八项描述符分布 | Results 3.4 |
| Figure S5 | dual与B-only的八项描述符分布 | Results 3.4 |
| Figure S6 | 正负活性类别之间的最大ECFP4相似性 | Results 3.4 |
| Figure S7 | 五种子的方向性AUROC绝对值 | Results 3.6 |
| Figure S8 | CAP5各评分方案影响 | Results 3.6 |
| Figure S9 | 配对Top-10%三类候选组成 | Results 3.7 |
| Figure S10 | 原Figure 8：CAP5及严格子集敏感性 | Results 3.6 |
| Figure S11 | 原评分方法差值的全部置信区间 | Results 3.2 |
| Figure S12 | 原口袋比较的原始估计及全部区间 | Results 3.5 |
| Scheme S1 | S1–S5评分与姿态关系 | Methods 2.4 |

补充表依次为 Table S1：选样与活性来源；S2：受体与盒子；S3：重对接；S4：人群与可用性；S5：方向性AUROC和方法差值；S6：负类与口袋；S7：化学模型；S8：种子与受体；S9：CAP5及严格子集；S10：Top-k名单。补充表与正文Table 1、Table 2分开编号。

## 与中文6版Word正文的明确不一致

1. Results 3.6原“Figure 8a”“Figure 8b”应改为“Figure S10a”“Figure S10b”，前提是按现有S10图注准确核对面板标签。不能继续称为正文Figure 8。
2. Results 3.7出现的两处“Figure 9”均须改为“Figure 8”；目前投稿目录没有正文Figure 9。
3. Results 3.3首段“Table Sx”没有已发布的正式表格承接。现有Table S6A只包含dual与两类负类的比较，不能作为探索性单靶对照（单靶正类）的直接证据。建议取得原始16方向人数、AUROC、95% CI及脚本后新增Table S11，并同步更新SI索引、图注、表格工作簿；在证据归档前，不能把Table Sx简单替换为Table S6。
4. Results 3.4“训练–测试最近邻相似性”不是Figure S6所展示的“跨活性类别最近邻相似性”。核对补充数据的专门诊断文件位置，不能将二者混为同一指标；Table S7现有内容为化学模型结果。
5. Results 3.7前5%、前20%分析应在Supplementary Data提供对应成员和截取比例文件。现有Table S10主表为Top-10%；不得仅以Top-10%表支撑更改名单长度后的数值。
6. 第6版Methods 2.3生产配体ETKDGv3种子20260727需恢复，以与补充方法及冻结配置一致；2.9的数据和代码位置应恢复为可检索的正式数据可用性说明。
7. Figure 3/6中的框线、符号仅用于报告已有未校正95%配对区间的状态；不能误解为区间宽度或多重比较调整后的统计显著性。图S11/S12保留完整区间。

## 逐项存在性与QA结果

- Main_Figures包含Fig1–Fig8四种导出格式（PDF/EPS/SVG/PNG）。
- Supplementary_Figures包含FigS1–FigS12及SchemeS1四种格式。
- Source_Data包含每幅现有图对应CSV；Tables包含Table1、Table2及TableS1–S10的可阅读表格与Excel工作簿。
- Figure_Captions_CN.md、Figure_Captions_EN.md、Figure_Captions_EN_Submission.md和Figure_Index.csv已存在。投稿与审阅英文图注末尾标点规则不同，应以期刊提交版为准。
- 仓库既有Release_QA.json记录图宽174 mm、至少8 pt字体、600 dpi RGB PNG、PDF/EPS嵌入字体、29页表格审阅版、21页图审阅版。此次核查采用记录与文本文件，未重新打开全部PDF逐页检查；视觉细节仍应在最终Word排版后复测。
- Table 1核对主受体JAK2=8BXH、F10=2JKH，避免沿用历史误置的PDB信息。
- 正式共晶重对接：70次、top-1小于2Å 42次、best-saved小于2Å 66次；EGFR/HER2 S1两个方向AUROC分别0.6274和0.2467；替换负类EGFR通道0.2467→0.7470；化学模型EGFR D_vs_B共同记录59条ECFP4 0.8790、S1 0.2408；这些与已发布Table S3/S5/S6/S7一致。
- 数值不一致时应修改图表/文本而不能改写原始冻结表。未计算或缺来源的数据不得自行补全。

## 参考标准

Molecular Diversity官方：https://link.springer.com/journal/11030/submission-guidelines（Research Article正文图表合计≤10；References≤80；摘要150–250英文单词；关键词4–6个）。
Zhou等（2013）已报道双靶激酶筛选的单靶抑制剂假阳性；Chaput等（2016）讨论评价集与对照组成；Truchon和Bayly（2007）讨论ROC与早期识别。因而新增图表应服务于方向性评价、负类影响及Top-k，而不是重复展示已有AUROC。

## 期刊额外提交文件核查

Molecular Diversity指南明确要求Graphical Abstract图片（至少531×1328像素，高×宽，或同比例更高尺寸）及少于20个英文单词的说明。当前Molecular_Diversity投稿目录的文件清单中未发现明确命名的Graphical Abstract成品，应在正式投稿前另行准备并检查版权/画面科学准确性；本次没有自行生成。期刊还要求图注在稿件文本文件中、Fig.及数字加粗、按先后顺序引用，图注结尾不加标点。第6版尚未整体应用图表最终编号，需同步正文后再做一次完整交叉引用扫描。
