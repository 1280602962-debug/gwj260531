# 投稿表格与补充数据

本目录采用统一投稿编号：正文 Tables 1–2、补充 Tables S1–S10。一个补充表的 A、B、C 面板在同一个工作表中，面板不另计表号。此前工作材料中的 S1–S31 不作为本目录的投稿编号。

整理来源为提交 `1090ed86d8cb553ee365e52bc51c34e48c0e5f6c` 的已验收结果。本次只整理展示、合并表格和复制数据，未重新对接、重评分、拟合模型或计算 bootstrap。主分析仍为 867 条靶点对观察，CAP5 支持分析保留 861 条，PIK3CA/mTOR 采用固定 110 分子名单。

| 分类 | 文件 | 用途 |
| --- | --- | --- |
| 正文表格 | [Main_Tables.xlsx](01_Main_Tables/Main_Tables.xlsx) | Table 1：靶点对及四类人数；Table 2：共同成员上的双靶配体检出 |
| SI 汇总表 | [Supplementary_Tables.xlsx](02_SI_Tables/Supplementary_Tables.xlsx) | Tables S1–S10，每个表号对应一个工作表 |
| 完整补充数据 | [Supplementary_Data.zip](03_Supplementary_Data/Supplementary_Data.zip) | 原始来源 CSV、当前完整评分、活性来源、化学预测、完整排名及 CAP5 输出 |
| 正文图 | [05_Main_Figures](05_Main_Figures/) | Figures 1–9，按 Results 首次出现顺序编号 |
| 补充图 | [06_SI_Figures](06_SI_Figures/) | Figures S1–S9 与 Scheme S1 |
| 图像审阅 | [Figure_Review.pdf](Figure_Review.pdf) | 19 页矢量图册，书签对应图号；完整图注另见 FIGURE_CAPTIONS.md |
| 图片目录 | [FIGURE_INDEX.csv](FIGURE_INDEX.csv) | 图片、引用位置、数据来源和表格的对应 |
| Results 核对 | [RESULTS_AND_FIGURE_AUDIT_CN.md](04_Editorial_Checks/RESULTS_AND_FIGURE_AUDIT_CN.md) | 205 项数值核对及需要修改的具体句子 |
| Results 引用 | [RESULTS_CITATION_MAP_CN.md](RESULTS_CITATION_MAP_CN.md) | 本轮 Results 原引用与新图表编号对照 |
| 表格目录 | [TABLE_INDEX.csv](TABLE_INDEX.csv) | 表号、面板、首次引用位置、工作簿和数据包内 CSV 位置 |
| 编号对照 | [RENUMBERING.csv](RENUMBERING.csv) | 原工作材料 S1–S31 的完整去向 |
| 来源对应 | [TABLE_SOURCES.csv](TABLE_SOURCES.csv) | 每个表及面板对应的已验收源文件 |
| 表题及表注 | [TABLE_CAPTIONS.md](TABLE_CAPTIONS.md) | 英文表题与表注，并附中文标题 |
| 方法引用 | [METHODS_TABLE_CITATIONS_CN.md](METHODS_TABLE_CITATIONS_CN.md) | 按当前 Methods 2.1–2.9 顺序使用的引用编号 |
| 参数说明 | [Supplementary_Methods_CN.md](Supplementary_Methods_CN.md) | 配体准备、评分方案和数据范围说明 |
| 投稿前需处理的事项 | [Outstanding_Evidence_CN.md](04_Editorial_Checks/Outstanding_Evidence_CN.md) | 已知证据边界及全文图表编号衔接事项 |

## 表格安排

| 编号 | 内容 |
| --- | --- |
| Table 1 | 靶点对、正式受体与 dual、A-only、B-only、neither 人数 |
| Table 2 | 每种方法与 Vina 的共同成员、入选人数、dual 检出数和召回率 |
| Table S1 | A，评价集选择条件；B，活性来源覆盖与标签核查 |
| Table S2 | 14 个正式受体及当前 6 个替代受体的链、共晶配体与对接盒 |
| Table S3 | 14 个正式受体在五个种子下的重对接汇总 |
| Table S4 | A，方向性与化学人群；B，双端评分可用性 |
| Table S5 | A，方向性指标及区间；B，共同成员上的评分方法差值 |
| Table S6 | A，负类替换；B，对应与非对应口袋 |
| Table S7 | 描述符、ECFP4、ECFP4＋Vina 的折外表现及评分增量 |
| Table S8 | A，五种子分布；B，六种单端受体替换 |
| Table S9 | A，CAP5 人数；B，Vina 主指标变化；C，PIK3CA/mTOR 严格子集；D，CAP5 化学模型与评分增量 |
| Table S10 | A，方法自身候选组成；B，共同人群上的方向性与入选组成变化；C，四状态计数 |

表格同时提供 Markdown 预览。Markdown 中小数显示到四位，Excel 保留源数值；完整 CSV 保留来源精度。人数和比例的分母见表注。原始长表移至补充数据，不按逐分子文件另设正式 Table S 编号。

## 评分方案

| 投稿名称 | 来源名称 | 定义 |
| --- | --- | --- |
| S1 | M0 | Vina seed 42 mode 1，评分取相反数 |
| S2 | M1 | 相同 Vina mode 1 姿态的 GNINA CNNscore |
| S3 | M2 | 相同 Vina mode 1 姿态的 RTMScore model 1 |
| S4 | M1b | Vina seed 42 全部保存姿态中的最高 CNNscore |
| S5 | M3 | GNINA 独立对接、按 CNNscore 排名第一姿态的 CNNscore |

工作簿按 S1、S2、S3、S4、S5 排列。原始 CSV 保留既有 M 名称，并随包附带名称映射。S4−S2 是姿态重新排序比较；S4−S1 是候选名单分析中的支持性比较，两者不能互换。

## 生成与核对

生成脚本为 `../scripts/rerun_v4_2/build_submission_tables.py` 和 `build_submission_workbooks.mjs`。前者只读取既有统计、整理视图，后者使用 `@oai/artifact-tool` 写入 Excel。Python 的 openpyxl 仅用于只读核对导出单元格。

运行 `python scripts/rerun_v4_2/build_submission_tables.py --project-root .`，工作目录为 `Dual_Target_Docking/`，并提供已安装 artifact-tool 的 Node 环境。`CODEX_PRIMARY_RUNTIME_NODE` 和 `CODEX_PRIMARY_RUNTIME_NODE_MODULES` 可指定该环境。生成程序不运行任何统计或对接模块。

数值、人数、评分映射、工作簿单元格和压缩包来源字节的核对结果见 [Validation_Report.json](04_Editorial_Checks/Validation_Report.json)。已有工作文稿、图版和历史核查文件继续保留其来源身份；投稿表格编号以本目录为准。

## 投稿图片

`build_submission_figures.py` 使用 Matplotlib 从既有结果 CSV 绘图，输出 7 in 宽的 PDF、SVG 和 600 dpi PNG；`07_Figure_Source_Data/` 中提供每图读入的数值。图内评分名称统一为 S1–S5，图注见 [FIGURE_CAPTIONS.md](FIGURE_CAPTIONS.md)，旧图去向见 [FIGURE_RENUMBERING.csv](FIGURE_RENUMBERING.csv)。

原 `results/jcim_stage2_figures/` 与 `manuscript/` 图文保留为工作来源，投稿图号以本目录为准。Results 原稿由作者后续修改，本次提供逐节引用对照与核对意见，未将原稿替换为新的完整 Results。
