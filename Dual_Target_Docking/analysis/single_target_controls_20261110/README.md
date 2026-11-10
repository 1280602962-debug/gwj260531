# Table S11 单靶探索性对照补充归档

## 结果与引用

中文第6稿3.3节的占位引用已对应至 **Table S11**。A部分比较本靶单靶活性配体与neither；B部分比较本靶单靶活性配体与对侧单靶活性配体。每部分包含八对靶点的16个评分通道，共32条AUROC记录。

| 核查项 | 结果 |
| --- | --- |
| Table S11A：AUROC点估计高于0.5 | 16/16 |
| Table S11A：95% CI完全高于0.5 | 4/16 |
| Table S11B：AUROC点估计高于0.5 | 11/16 |
| Table S11B：95% CI完全低于0.5 | 2/16 |
| EGFR评分，EGFR-only为正类、HER2-only为负类 | n=37/32；0.1807 [0.0845, 0.2914] |
| JAK1评分，JAK1-only为正类、TYK2-only为负类 | n=32/32；0.3125 [0.1865, 0.4512] |

两处区间完全低于0.5的说法获得本次复现支持。这里的“两个”指上述两个评分通道，不是同一靶点对的两个dual方向。Table S11中的正类为单靶活性类别，不能与Table S6A中以dual为正类的负类替换混用。Table S11和Figure S11的编号相互独立。

## 来源与计算状态

输入为提交`2017a5d7f9815715f5203a7b4e413c9719e9e032`中的`results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv`。869行中867行符合活性条件，同时包含`main`和`q110_primary`。分析排除dual及activity_eligible=0，按评分靶点检查S1是否有效。AChE/BChE的A-only有效人数为25，B-only为28，neither为15；不另要求配体在未评分靶点也有有效分数。所有类别仍按正式6.0阈值定义。

此前文稿已报告0.1807和0.3125的点估计。已核查当前分支中正式表、Source_Data、results/current及分析脚本；未恢复这项探索性单靶对照的原始bootstrap输出。因此，本次归档状态为**从冻结评分复现并补算区间**，不能标记为找回原始bootstrap文件。既有主要方向性指标、Table S6A及其置信区间保持原值。

使用已取负的`M0_score_A/B`，不再次反转分数。AUROC按Mann-Whitney定义，评分并列贡献0.5。采用按活性类别分层的10,000次有放回重采样；每个靶点对×评分靶点重置NumPy default_rng(271828)。保留输入行序；每次依次抽取正类、对侧单靶类、neither的索引，两个比较共用本次正类抽样。区间取linear 2.5/97.5百分位。区间未作多重比较校正，也不属于同时置信区间。

本次只计算表S11所需的探索性统计，不重新对接、不训练化学模型、不修改冻结评分或正式评价集。输入前后SHA-256一致。`Table_S11_config.json`记录精确环境、协议与来源指纹。

## 文件与复现

- `Table_S11_metrics.csv`：32条全精度点估计、样本数、区间及区间位置。
- `Table_S11_membership.csv`：逐分子纳入、排除原因、实际评分及类内抽样顺序。
- `Table_S11_bootstrap_replicates.npz`：16个通道各10,000次重采样的两个AUROC。键01–08对应固定靶点对顺序，A/B表示评分端；第0列为对侧单靶负类、第1列为neither负类。
- `Table_S11_QA.csv`：32项AUROC与scikit-learn独立实现的核对、正负互换关系与有限重采样数。
- `TableS11.md`、`Table_S11_sections.html`：可编辑出版表格。
- 投稿目录`submission/Molecular_Diversity/Tables/`新增独立`Table_S11.xlsx`和`Table_S11.pdf`；既有`Supplementary_Tables.xlsx`继续保留S1–S10，S11以独立补充工作簿提供。完整PDF审阅稿与HTML已追加S11。

从仓库根目录执行：

```bash
python Dual_Target_Docking/scripts/rerun_v4_2/supplement_single_target_controls.py
python Dual_Target_Docking/scripts/rerun_v4_2/render_table_s11.py
```

统计脚本需要Python、NumPy、scikit-learn；排版脚本另需PyMuPDF、reportlab和Liberation Sans字体。脚本不写入原始评分表。

## 中文稿的精确修改

3.3节探索性单靶对照段落中，第一处“其中4个方向的95%置信区间完全高于0.5”后增加`（Table S11A）`；将“两个方向的置信区间均完全低于0.5（Table Sx）”替换为“上述两个评分通道的95%置信区间均完全低于0.5（Table S11B）”。

已将这两处修改应用于`双靶对接中文6_投稿前修订建议版.docx`同一段落，移除已关闭占位引用的黄色标记，未改写其余段落。
