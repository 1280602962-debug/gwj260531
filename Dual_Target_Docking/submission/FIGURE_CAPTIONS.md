## Figure 1

Experimental classes and activity-source coverage. (a) Numbers of dual, A-only, B-only and neither observations in each target pair. (b) Counts of observations with raw records on both targets, complete specified document identifiers, mixed sources or endpoints, and threshold-discordant repeated measurements. Cell color represents the fraction within each pair; text gives counts. Audit flags overlap. Complete identifiers require DOI, PubMed or patent identifiers for all tied activity maxima at both endpoints. The 867 observations represent 845 distinct global ligand entities.

评价集的四类组成及活性来源覆盖。颜色表示各靶点对内比例，格内数字为观察数；来源标记可以重叠。

## Figure 2

Vina directional discrimination and summary_min. Points and horizontal bars show estimates and stored 95% class-stratified percentile bootstrap intervals. (a) Dual versus A-only uses target B; dual versus B-only uses target A. (b) summary_min is the smaller directional AUROC, calculated within each resample. The vertical line marks AUROC 0.5 as a reference. The summary interval is not a simultaneous test that both directions exceed or fall below 0.5.

Vina 两个方向的 AUROC 和 summary_min。点与横线分别为点估计及已有的 95% bootstrap 区间，虚线为 0.5 参考值。

## Figure 3

Paired changes in summary_min between scoring workflows. (a) S2-S1, fixed-pose GNINA versus Vina. (b) S3-S1, fixed-pose RTMScore versus Vina. (c) S4-S2, selection of the highest CNNscore among saved Vina poses versus the original mode-1 pose. (d) S5-S1, independent GNINA docking versus Vina. Every contrast uses its own complete-case paired members. Bars are stored 95% paired bootstrap intervals; no across-pair pooled estimate is calculated.

四类评分工作流的配对差值。每个比较使用自身的共同有效成员，横线为已有的 95% 配对 bootstrap 区间。

## Figure 4

Effect of negative-class composition. (a) AUROC using single-target-selective or neither negatives, with the dual set and score target held fixed within each directional comparison. (b) Neither minus single-target AUROC with stored 95% bootstrap intervals. All scores are Vina S1. Directional counts are listed in Table S6A; neither and single-target negatives are sampled separately while sharing dual resamples.

单靶选择性与 neither 负类的 AUROC 对照及差值，固定每个比较中的 dual 成员与评分靶点。

## Figure 5

Ligand chemistry and the incremental contribution of Vina. (a) Descriptor, ECFP4 and ECFP4-plus-S1 models use existing out-of-fold predictions. Vina S1 is a direct score evaluated on the same members and is not an out-of-fold learned model. (b) AUROC(ECFP4 plus S1) minus AUROC(ECFP4). All 16 directional tasks use matched chemistry members. Points are estimates; no confidence interval for the chemistry-model increment was calculated.

配体化学模型与 Vina 的共同成员比较。三个化学模型使用折外预测，单独 Vina 为直接评分；图中增量为点估计，没有新增置信区间。

## Figure 6

Corresponding and noncorresponding pocket controls. (a) AUROC values on the same members using corresponding and noncorresponding target scores. The corresponding target is B for dual versus A-only and A for dual versus B-only. (b) Corresponding minus noncorresponding AUROC with stored 95% paired bootstrap intervals. All scores are Vina S1. No additional member is lost relative to the corresponding directional S1 population.

对应与非对应口袋的同成员对照及配对差值。

## Figure 7

Sensitivity to search seeds and receptor structures. (a) Five seed-specific Vina summary_min estimates. Short horizontal segments show the stored minimum-to-maximum range; black ticks show medians. These segments are not confidence intervals. Each seed uses its valid members; AChE/BChE seed 17 loses one A-only ligand relative to seed 42. (b) Paired changes in affected-direction AUROC and summary_min for six single-target receptor replacements. Parentheses identify replaced side A or B; the other receptor is held fixed. Bars in (b) are stored 95% paired bootstrap intervals.

随机种子与单端受体替换敏感性。a 中短线表示五种子的范围，非置信区间；b 中横线为受体比较已有的 95% 配对 bootstrap 区间。

## Figure 8

Sensitivity to evaluation-set membership. (a) Vina summary_min for the primary and CAP5 populations; labels give original and retained pair sizes. CAP5 retains 861 of 867 observations. (b) summary_min for the five scoring schemes in the 110-member PIK3CA/mTOR set and its 56-member strict-activity subset. Population-size labels include all four experimental classes; directional sample sizes are reported in Tables S5A and S9C. Bars are stored 95% bootstrap intervals for each population; no paired interval for the between-population change is shown. Full directional CAP5 comparisons, including the PPARA/PPARD S5 dual-versus-A-only interval-boundary change, are shown in Figure S8. Chemistry sensitivity is tabulated in Table S9D.

评价成员组成敏感性。a 为 CAP5 前后的 Vina 指标，b 为 PIK3CA/mTOR 主评价集与严格子集的五种评分。不同人群的原有区间分别展示，未新增人群差值区间。

## Figure 9

Relationship between directional discrimination and candidate-list composition. Each point is one target pair. Both axes show the new scheme minus S1. The horizontal axis is the paired directional summary_min change; the vertical axis is the change in single-target-selective fraction among the top 10% of joint ranks, expressed in percentage points. Ranking and directional changes use their respective complete-case populations, as defined in Methods. Colors and symbols identify the four change categories. Open symbols identify the supporting S4-S1 comparison. No new confidence interval is added. Matched top-list counts are shown in Table 2 and Figure S9.

方向性指标变化与联合名单组成变化的关系。纵轴为百分点，不是相对百分比；每个方法比较分别使用自己的共同成员。空心点表示 S4 对 S1 的支持性比较。

## Figure S1

Cocrystal-ligand redocking for 14 primary receptors and five search seeds. (a) Top-ranked-pose RMSD. (b) Minimum RMSD among all saved poses. RMSDs are symmetry-aware heavy-atom distances computed in the receptor frame. The dashed line marks 2 angstrom. Each point is a recorded seed result; no coordinate or RMSD is recomputed for this figure. Axes differ between panels to preserve visibility of large top-ranked-pose errors.

14 个正式受体的五种子重对接，分别展示首位姿态及保存姿态中的最佳 RMSD。两个面板的横轴范围不同。

## Figure S2

Availability of scores by experimental class, target and scoring scheme. Each cell gives valid/expected counts; its color represents that fraction. Pair numbers are identified in the figure. Missing and invalid scores retain their original status, and are not imputed. The S1-S5 order follows the Methods scoring definitions.

按活性类别、评分靶点与评分方案展示有效评分人数；每格为有效数/应有数。

## Figure S3

Directional AUROC and summary_min for all five scoring schemes. Cell annotations show estimates rounded to two decimals; colors use a common 0-1 scale. Intervals and sample sizes are given in Table S5A. The heat map does not replace paired method comparisons, which use their own shared members.

五种评分方案的两个方向 AUROC 与 summary_min；图中仅为点估计，完整区间见 Table S5A。

## Figure S4

Physicochemical descriptor distributions for dual and A-only ligands. Points are labeled-population observations before score and fold filtering. Boxes show the median and interquartile range, with whiskers to the most extreme values within 1.5 IQR. All observations are retained as points, including those beyond whiskers. Pair numbers are identified in the figure. These distributions are descriptive and do not imply independently validated classification.

dual 与 A-only 配体的八项理化描述符分布。展示标签层全部观察；点为配体，箱线为中位数、四分位区间及 1.5 IQR 范围。

## Figure S5

Physicochemical descriptor distributions for dual and B-only ligands. Points are labeled-population observations before score and fold filtering. Boxes show the median and interquartile range, with whiskers to the most extreme values within 1.5 IQR. All observations are retained as points, including those beyond whiskers. Pair numbers are identified in the figure. These distributions are descriptive and do not imply independently validated classification.

dual 与 B-only 配体的八项理化描述符分布。展示标签层全部观察；点为配体，箱线为中位数、四分位区间及 1.5 IQR 范围。

## Figure S6

Structural similarity between the two experimental classes within each directional task. Each point is the maximum ECFP4 Tanimoto similarity to any ligand in the other class of the same pair and direction. Fingerprints use radius 2, 2048 bits, useChirality=False and useFeatures=False. Pooled boxes summarize the two classes together; class identity is shown by point color and shape.

每个配体与同一靶点对、同一方向另一活性类别配体的最大 ECFP4 Tanimoto 相似性。箱线汇总两类配体，点的颜色与形状区分类别。

## Figure S7

Full seed sensitivity for directional AUROC and summary_min. Each panel shows five recorded seed-specific estimates, stored minimum-to-maximum ranges and medians. Range segments are not confidence intervals. Effective populations follow each seed's valid scores; seed 17 loses one A-only AChE/BChE observation relative to seed 42.

五种子下两个方向 AUROC 和 summary_min 的完整变化。短线为种子间范围，非置信区间。

## Figure S8

Full directional CAP5 sensitivity in the two pairs with excluded members. Each scheme has AUROC_B (dual versus A-only), AUROC_A (dual versus B-only) and summary_min. Primary and CAP5 estimates retain their respective stored 95% intervals; no interval for the difference is added. The asterisk identifies PPARA/PPARD S5 AUROC_B, whose interval changes from [0.4921875, 0.767578125] to [0.5030241935, 0.7772177419]. This boundary change is reported explicitly without a multiplicity-adjusted significance claim. The other six pairs have unchanged members.

CAP5 实际改变成员的两个靶点对中，五种评分的两个方向及 summary_min 的完整比较；标明 PPARA/PPARD S5 的区间边界变化。

## Figure S9

Class composition of matched top-10% joint-ranking lists. Each target pair has an upper hatched S1 bar and a lower solid new-scheme bar, evaluated on that comparison's common two-target complete-case members. Segment labels give counts. Population sizes and k may differ between method comparisons; the same comparison has identical N and k for its two bars. S4-S1 is supporting. No method-specific S1 list is substituted for these paired baselines.

各方法与 Vina 在各自共同成员上的前 10% 候选名单。上方斜线为 S1，下方实心为补充评分；每段数字为类别人数。

## Scheme S1

Scoring workflows S1-S5. S2 and S3 rescore the same Vina mode-1 pose used by S1. S4 chooses the highest CNNscore among all poses saved from that Vina search. S5 uses an independent GNINA search and CNNscore pose ordering. S4-S2 tests reranking, whereas S5-S1 compares complete docking workflows.

S1–S5 的评分和姿态选择关系。S4−S2 为保存姿态重新排序比较，S5−S1 为独立工作流比较。
