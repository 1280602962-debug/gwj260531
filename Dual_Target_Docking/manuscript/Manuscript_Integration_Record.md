# 2026-10-08 文稿整合记录

本文件记录本轮展示整合。它不是新的科学计算，也不评估投稿就绪。

## 基准

- 科学实现基准 `AUDITED_IMPLEMENTATION_BASE_SHA=3f3e0d26052178c38ccb43dec637683a60f97a9d`
- 已验收结果与本轮写作起点 `ACCEPTED_RESULTS_REFERENCE_SHA=WRITING_BASE_SHA=319375f10cb468d3a36e9cd0cf3e420bb5aa9632`
- `SECOND_STAGE_COMPUTE_EXECUTED=YES` 是既有事实
- 本轮 `NEW_SCIENTIFIC_COMPUTE_EXECUTED=NO`
- 本轮 `G_J_K_EXECUTED=NO`
- `SUBMISSION_READINESS=NOT_ASSESSED`

## 展示操作

排序、宽表转列、读取既有四态计数、把 Top-k 单靶占比差乘以 100 后以百分点显示、D1 箱线用 NumPy linear 百分位、以及只改变类别轴位置的确定性抖动。这些操作不产生新的 AUROC、区间、化学预测或名次。

`jcim_stage2_figures.py` 不再调用 `write_concordance`。Figure S8 和 Table S17 读取现有 `results/jcim_stage2_figures/module_I_concordance_counts.csv`。核对脚本从 `module_I_method_compare.csv` 的 32 行独立计数，结果写入本记录和 `manuscript/source_views/display_audit_report.txt`，不覆盖原计数 CSV。独立计数与文件一致：M1 为 2/2/2/2，M1b 为 3/3/1/1 且 supporting 为 YES，M2 为 4/3/1/0，M3 为 3/2/2/1。M2 的 NO_CHANGE 为 0。

## 主张到展示位置

| 主张 | 来源 | 区间 | 成员 | 单位 | 位置 |
|---|---|---|---|---|---|
| 三层人数与身份 | Table S1、Table S2、人群表 | 无 | 808/807/805；938/934/928 | 记录或配体，按表内字段 | 结果 R1；Table 1；Table S2 |
| PRIMARY 方向与配对差 | Table S3、Table S4 | 已存 paired bootstrap | 各方法自身方向任务；差值按 PRIMARY 定义 | AUROC | Figure 1；Tables S3–S4 |
| 负类对照 | Table S5 | 已存 paired bootstrap | 同一 dual，两种负类 | neither−single | Figure 2；Table S5 |
| 化学与对接增量 | Tables S7–S9 | 增量无区间 | 第 3 层 928 | AUROC 点差 | Figure 3；Figures S1–S3 |
| 口袋对照 | Table S6 | 已存 paired bootstrap | 与方向任务同成员 | correct−wrong | Figure 4；Table S6 |
| 五种子 | Tables S10–S11 | IQR 不是置信区间 | 五个已存种子 | AUROC | Figure 5a 只显示 summary_min；两方向完整结果在 Figure S5 |
| 替代受体与未替换端 | Table S12 | 受影响端与 summary_min 使用已存区间 | 八个已存条件 | AUROC 差 | Figure 5b；Figure S6；Table S12 |
| 方向与 Top-k 四态 | Tables S14–S17 | 不新增区间 | 共同比较，不用 method_specific M0 做统一基线 | 方向差为 AUROC；纵轴为百分点 | Figure 6；Figure S8；Table S17 |

## 环境

正式统计环境为 Python 3.11.15、NumPy 2.4.6、RDKit 2025.09.5、scikit-learn 1.9.0。本轮绘图使用同一 Python 3.11.15 与 NumPy 2.4.6，另加 Matplotlib 3.10.9。绘图脚本不导入 RDKit 或 scikit-learn。两份审阅 PDF 由 Matplotlib 嵌入从 `msyh.ttc` 抽出的微软雅黑生成。ReportLab 曾被试用，第二个字体子集在页面上变成缺字，因此没有采用。

## 来源缺口

`Source_Register.csv`、补充信息第 1 节和可采用的 Methods 文字都写明：805 条标签可以回溯到正式表，但没有逐条核对原始实验，也没有同一测定平台或逐条文献端点。`ACTIVITY_LABEL_SEMANTICS_AUDIT.csv` 只覆盖少量特殊行。SHA256 只核对数据包字节。

## 核对与未验证项

核心数字由独立脚本对照正式 CSV 与正文，记录为 `display_audit_report.txt` 的 38/38 PASS。正式图 PDF 页宽为 504.0 pt，即 177.8 mm。最终图和两份审阅 PDF 已实际打开查看：中文可读，图题与对应图留在同一页，Figure 1 横轴为 1–8。

专用几何对齐、PDF 文本和碰撞工具没有运行。该项为 NOT_VERIFIED。替代检查是页宽测量、独立数字核对和人工查看。Figure S4 格内文字为 5.5 pt，Figure S7 格内文字为 4.8 pt，低于本轮 7 pt 工作目标；S7 仍高于 4.5 pt 的最低文字规格。本记录不把该替代检查写成专用碰撞检查已经通过。

本轮没有克隆或完整阅读 scientific-writing 与 nature-figure 的远端技能树。本地表格和 PDF 读取约定只用于生成和检查文件，不作为科学证据。
