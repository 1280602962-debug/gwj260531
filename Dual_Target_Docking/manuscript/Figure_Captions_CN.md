# 图注与表注

编号 1–8 对应 EGFR/HER2、JAK1/JAK2、JAK1/TYK2、PIK3CA/mTOR、AChE/BChE、F2/F10、PPARG/PPARA、PPARA/PPARD。D_vs_A 使用 B 端分数，D_vs_B 使用 A 端分数。

## Scheme 1

研究设计与评分关系。左侧为四类实验标签，primary 阈值为 6.0；activity_eligible=0 不是 neither。右侧为两个方向任务。中部说明 summary_min 是两个方向 AUROC 的较小值。下方区分 M0 固定姿势、M1/M2 对同一姿势重评分、M1b 在已保存姿势上重排，以及 M3 独立对接。该图是定义示意，不是实测活性散点，也不画结合残基。八对 PDB 见 Table 1。三层人数见 Table S2。

## Figure 1

双侧方向判别与配对评分对照。a，M0 在八对上的 D_vs_A、D_vs_B 和 summary_min。点为已存 AUROC，竖线为该行自己的 95% 区间；summary_min 的区间来自 contrast=summary_min，不复制较弱臂的区间。参考线为 0.5。b，M1−M0、M1b−M1、M2−M0 和 M3−M0 的 summary_min 配对差值及各自区间。M1b−M1 是正式重排对照，本图不给 M1b−M0 加区间。区间来自分层或配对 bootstrap，B=10000，seed=271828，linear 2.5/97.5 百分位。样本量见 Table S3 和 Table S4。

## Figure 2

单靶与 neither 负类的方向判别。a，同一 pair×方向上，dual 对单靶负类和 dual 对 neither 的 AUROC，连线只连接同一单元。b，neither−single 的配对差值及 95% 区间。纵轴注明 n_dual/n_single/n_neither。PIK3CA/mTOR 的 neither 每臂为 4。重采样共享 dual 索引，两种负类分别抽取。完整数值见 Table S5。

## Figure 3

固定分折的化学预测与 M0 增量。a，理化描述符、ECFP4、ECFP4+M0 和同成员 chemistry-matched M0 的折外 AUROC。b，ECFP4+M0−ECFP4 的点增量。四个点估计使用第 3 层同一成员，没有误差线、折间标准差或区间。matched M0 不替换 PRIMARY M0。每臂人数见 Table S2 和 Table S7，OOF 记录见 Table S9。

## Figure 4

规定口袋与对侧口袋的评分对照。a，correct 与 wrong 的方向 AUROC，连线只连接同一臂。b，correct−wrong 的配对差值及 95% 区间。D_vs_A 的正确口袋为 B，D_vs_B 的正确口袋为 A。两口袋使用同一单元成员。16 个臂的额外减员均为 0，见 Table S6。该图描述评分差异，不是机制验证。

## Figure 5

五种子与替代受体的核心敏感性。a 只显示 summary_min 的五个已存种子点。竖线为该行已存最小值到最大值，短横线为已存中位数，轴下 IQR 属于 summary_min，不是置信区间。两个方向的完整结果在 Figure S5 和 Tables S10–S11。AChE/BChE 的 seed 17 比其他种子少 1 个 A_only。b，四对八个替代条件的受影响端差值和 summary_min 差值及已有配对区间，顺序为固定靶点对顺序和 A 端然后 B 端。未替换端见 Figure S6 和 Table S12。summary_min 先在各条件内取两方向较小值，再相减。

## Figure 6

同成员的方向增量与名单变化。四个面板分别为 M1、M1b、M2、M3 对 M0 的八个比较。横轴为已存 summary_min 差值。纵轴为 Top-k 单靶占比变化，单位是百分点，由原比例乘 100 得到；四态仍按原比例点值读取。点号与 Table 1 的 1–8 对应。实心或空心符号同时用颜色区分一致改善、一致变差、不一致和无变化。M1b 使用空心符号，表示 supporting：横轴是三组逐类成员一致后的两个 PRIMARY 点估计之和，没有新区间。无变化表示至少一项恰为 0。N、k 和 dual 数见 Table 2。名单来自实验标签构成的挑战集。

## Figure S1

D_vs_A 的八个理化描述符。每个面板是八对的 dual 与 A_only，成员为 Table S8 的第 1 层，不包含 neither。箱体为 linear 的 P25–P75，中线为 P50，须延伸到 1.5×IQR 范围内最远实测值。全部原始点另层绘出，离群点不重复成第二套点。抖动只改变类别轴位置。TPSA 单位为 Å²，MolLogP 为计算描述符。FormalCharge 的 0 保留。

## Figure S2

D_vs_B 的同八个描述符，比较 dual 与 B_only。变量、单位和纵轴范围与 Figure S1 对应。成员仍为第 1 层。

## Figure S3

每臂最大异类 ECFP4 Tanimoto。左为 D_vs_A，右为 D_vs_B。相似度在同一 pair×arm 内相对另一标签类计算，不是对训练集的最近邻，也不是适用域验证。分布来自 Table S8。

## Figure S4

五方法的 D_vs_A、D_vs_B 和 summary_min 点估计。格子为 PRIMARY 的 120 行 AUROC，颜色范围 0–1。该图不把两格相减当作配对差值。区间和方法自身有效人数见 Table S3。

## Figure S5

五种子的 auc_B、auc_A 和 summary_min。符号在三个面板中一致。竖线为已存 min–max，短横线为已存中位数，标注为该面板指标自己的 IQR。每个汇总有五个有限点。IQR 不是置信区间。缺失名单见 Table S10 和 Table S11。

## Figure S6

四对、八个替代受体条件的受影响端、未替换端和 summary_min 差值及已有 95% 区间。未替换端八行均为 0，区间也为 0，这是只替换指定一端的构造核对。共同成员和替换侧见 Table S12。

## Figure S7

四类标签和两个目标端上，各方法的 valid/expected。颜色表示该比值，格内文字为人数。expected、valid、missing 的定义和原因见 Table S13。本图不说明缺失是否随机。

## Figure S8

左图为 M0 在方法自身人群上的 Top-k 三类人数，不是四个共同比较共用的基线。右图为已存四态计数，M1b 标为 supporting。计数与 Table S17 相同，并与 Table S15 的 32 行一致。四态不是独立性检验。

## Table 1

八对靶点及实验标签构成。人数来自 activity_eligible=1 且类别为 dual、A_only、B_only 或 neither 的人群记录。合计 805 条合格分类观察。807 行人群表另有 EH120_059 和 AB_087 两条不合格、无类别记录。原 mapping 因 AB_040 alias 而为 808 行。这些人数早于 M0 有效性筛选。dual 进入两个方向，跨靶点对可以共享实体。938、934 和 928 的关系见 Table S2。

## Table 2

各方法与 M0 的共同成员及 Top-k dual 变化。每一列使用该方法自己的共同名单：新方法行的 population_kind 为 pairwise_common_vs_M0，对应 M0 行的 population_kind 为 pairwise_common_vs_该方法。M1b 列为 supporting。N 和 k 在同一比较的两侧相同。本表不重复 Figure 6 的差值，也不把 dual 数写成天然库命中数。三类组成和全排名见 Tables S14–S16。

## Tables S1–S17

S1 为 mapping 的 808 行，保留 alias、无类别和 eligibility。S2 为 16 臂三层人数及缺失名单。S3、S4 为 PRIMARY 的 120 行和 32 行。S5–S7 为负类、口袋和化学单元。S8 为 938 条第 1 层描述符，Murcko 空字符串表示无环。S9 为 928 条 OOF。S10、S11 为五种子原点和汇总。S12 为八个替代条件及共同 ID。S13 为 320 行缺失计数。S14–S16 为排名描述、32 个比较和 9058 条名次，in_topk 为 YES/NO。S17 为四行既有四态计数。工作表首行保留原字段名，完整精度以数据包中的原 CSV 为准。
