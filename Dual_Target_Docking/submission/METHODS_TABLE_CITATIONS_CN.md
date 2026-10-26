# Methods 表格引用与首次出现顺序

正文 Table 与补充 Table S 分别按各自首次出现顺序编号。Table 2 是结果表，应在联合排序结果段首次引用，不必为了在 Methods 中出现而提前插入。

| 文章位置 | 建议引用 | 对应内容 |
| --- | --- | --- |
| 2.2 活性来源核查 | Table S1B；Supplementary Data | 来源覆盖与标签一致性；逐分子对应及原始活性另附数据 |
| 2.2 八对靶点及类别人数 | Table 1 | 四类评价成员人数 |
| 2.2 评价集构建与选择条件 | Table S1A；Supplementary Data | 各对实际选择条件及冻结成员，逐分子 mapping 不单独编号 |
| 2.3 受体结构、链和对接盒 | Table S2 | 正式及当前替代结构信息 |
| 2.3 共晶配体重对接 | Table S3 | 五种子汇总；逐种子记录见 Supplementary Data |
| 2.3–2.4 软件、配体准备与评分参数 | Supplementary Methods | 参数用文字说明，评分关系保留工作流图 |
| 2.5 人群与评分缺失 | Table S4 | 有效人数与缺失；名单见 Supplementary Data |
| 2.5 方向性指标及方法差值 | Table S5 | 两个方向、summary_min 及配对方法差值 |
| 2.6.1 负类组成 | Table S6A | 单靶选择性与 neither 负类的配对对照 |
| 2.6.2 化学模型 | Table S7；Supplementary Data | 折外模型表现；描述符、骨架、相似性和预测明细另附 |
| 2.6.3 对应与非对应口袋 | Table S6B | 复用已首次出现的 Table S6，不新增编号 |
| 2.7 随机种子 | Table S8A | 五种子分布 |
| 2.7 替代受体 | Tables S2、S8B | 结构参数与六种替换结果 |
| 2.7 CAP5 | Tables S9A,B；Supplementary Data | 人数、Vina 汇总变化及全部支持分析 |
| 2.7 CAP5 化学模型 | Table S9D | 化学模型重训练后的 ECFP4 表现及加入 Vina 的点估计增量 |
| 2.7 严格活性子集 | Table S9C | PIK3CA/mTOR 的 56 分子严格子集 |
| 2.8 联合排序 | Table S10；Supplementary Data | 候选组成、配对变化及全排名 |
| 联合排序 Results | Table 2 | 共同成员的 dual 检出与召回 |
| 2.9 软件与数据可用性 | Supplementary Methods；Supplementary Data | 软件记录及完整数据索引 |

当前粘贴的 Methods 没有单列严格子集分析。如果保留 Table S9C，可在 2.7 末尾补充：

“此外，在 PIK3CA/mTOR 评价集中按 6.5/5.5 条件形成严格子集，使用已有评分进行支持性评价（Table S9C）。”

评分方案名称按 `SCORING_SCHEMES.csv` 使用。正式差值为 S2−S1、S3−S1、S4−S2、S5−S1；联合排序另有 S4−S1 支持性点估计。方法差值采用各自共同成员，不能直接用 Table S5A 中不同人群的两个 summary_min 相减。

图片引用：2.3 的重对接可引 Figure S1，2.5 的评分可用性可引 Figure S2；评分工作流改用补充 Scheme S1。正文结果 Figure 1 从 Results 3.1 开始，不为流程图额外占用一个正文 Figure 编号。
