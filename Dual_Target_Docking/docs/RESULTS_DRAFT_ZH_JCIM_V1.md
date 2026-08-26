# Results（中文工作稿 · JCIM Articles）

## 3. 结果

### 3.1 实验数据供给限制了严格双靶基准的构建

为确定公开生物活性数据是否能够支持严格的双靶点识别评测，我们首先对 49 对有 ChEMBL 缓存的候选靶标进行供给审计（Figure 2）。一端达到活性阈值、对端明确低活性的配体定义为方向性选择性硬负样本。

在严格标签规则下（dual：两端 pChEMBL ≥ 6.5；选择性类：活性端 ≥ 6.5 且对端 ≤ 5.5），能够同时提供足量 A-only 与 B-only 硬负样本的靶对十分有限。两端严格硬负均不少于 50 的厚面板条件仅有 4 对满足。排除金属依赖 HDAC1/HDAC6 后，PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 构成三个规模相对充足的靶对；EGFR/HER2 仅有 7 个严格 B-only 配体，因此被保留为供给受限案例（Table 1）。BindingDB 与 PubChem 的零对接计数核对支持同一供给稀缺结论（Table S12）。

因此，最终基准的规模主要由公开实验数据中方向性选择性硬负样本的可获得性所约束。严格 6.5/5.5 规则用于量化供给并记录面板构建，而 θ = 6.0 定义全部主 AUROC 的实验状态标签（Methods 2.1）。

### 3.2 基准设定改变了表观双靶判别

在冻结的四对靶标上，采用统一 θ = 6.0 标签规则和口袋匹配方向 AUROC 对 Vina 对接分数进行评价（Figure 1B；Methods 2.4）。分数定义为 \(S=-E_{\mathrm{Vina}}\)（越大越好），dual 为正类。预先指定的最差方向判别摘要为 `summary_min`。算术、几何与调和平均下四对排序不变（Table S26）。AChE/BChE 与 PIK3CA/PIK3CB 上，θ = 6.0 与更严格建造门槛给出相同分类与 AUROC（Table S4）；EGFR/HER2 与 PIK3CA/mTOR 在严格规则下 B_only 过少。整张阈值网格内排序趋势保持一致（Figure S1A）。

EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 的方向性 summary_min 分别为 0.430、0.606、0.500 和 0.692（Table 2；Figure 4A）。EGFR/HER2 的 dual-versus-B-only AUROC 为 0.430，PIK3CA/PIK3CB 为 0.500；PIK3CA/mTOR 两个方向分别为 0.714 和 0.692。相对池化协议，口袋匹配抬高了点估计但未改变排序（Table S6）。

同一套冻结分数再按 Dual versus neither 以及 Dual versus all non-duals 计分（Table 3；Figure 3）。Dual versus neither 是使用实验 inactive 的描述性对照（`vina_mean`）。EGFR/HER2 提供了最清晰的设定对照。Dual versus neither 的 AUROC 为 0.756 [0.562, 0.920]（n_neg = 12），而方向性 summary_min 仍为 0.430 [0.282, 0.578]。Dual versus all non-duals 降至 0.551 [0.443, 0.666]，表明引入选择性配体后判别任务明显变难。在 110 个 EGFR/HER2 配体的混合库中按 `vina_mean` 取 Top-10：1 个 dual、5 个 A-only、4 个 B-only、0 个 neither（EF10 = 0.393；硬负比例 = 0.90）；EF5 也低于随机（Table S25）。AChE/BChE 与 PIK3CA/PIK3CB 的 Dual-versus-neither 增量很小（0.649 与 0.559），区间与方向性臂重叠。PIK3CA/mTOR Dual versus neither 因 neither n = 4 而效能不足；该对 Dual versus all non-duals 为 0.674，接近 summary_min 0.692。

**Table 2.** 冻结 K = 4 评价集上的口袋匹配方向 AUROC（Vina，统一 θ = 6.0），并列出四个预先指定描述符的 `summary_min`。表中类别样本量为 n_scored（dual / A-only / B-only），即对应方向分析中具有有效 Vina 分数的配体数，因而可以低于冻结面板定额；neither 不进入主终点。n_panel 与两端对接覆盖见 Table 1 和 Table S27。最高描述符是最佳单一描述符参考。错口袋与配体效率见 Table S6；描述符双臂见 Table S28。

| 靶对 | n_scored (dual / A-only / B-only) | dual 对 A_only（口袋 B） | dual 对 B_only（口袋 A） | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|------:|---:|------:|-----:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.350, 0.650] | 0.622 | 0.620 | 0.595 | 0.418 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] | 0.463 | 0.448 | 0.310 | 0.260 |

**Table 3.** 同一套 Vina 分数在 Dual-versus-neither 与方向性设定下的 AUROC（统一 θ = 6.0）。Dual-versus-neither 使用实验 inactive（`vina_mean`）；Dual versus all non-duals 把 A-only、B-only 与 neither 都计为负类。方向性 CI 来自 Table 2。负样本集合不同。PIK3CA/mTOR Dual versus neither 效能不足（n_neg = 4）。

| 靶对 | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|------|--------------------------------:|------------------------------:|----------:|----------------------:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/PIK3CB | 0.500 [0.350, 0.650] | 0.559 [0.373, 0.746] | 16 | 0.556 [0.437, 0.672] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |

四个靶对的 `summary_min` 95% bootstrap CI 均包含 0.5；因此在本研究的样本量下，没有一个靶对获得排除随机水平的明确证据。按观察得的类别样本量做可分辨效应模拟表明，当前样本更容易分辨较大的方向性效应，而对中等效应则较弱。当两臂真实 AUROC 均为 0.70 时，`summary_min` CI 排除 0.5 的概率在 EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 上分别为 0.62、0.50、0.56 和 0.22；真实 AUROC 为 0.60 时，相应概率为 0.03–0.07（Table S31；Figure S6）。因此，CI 未能排除 0.5 并不能建立与随机等价。

主面板两端均得分：EGFR/HER2 110/110，AChE/BChE 95/100，PIK3CA/PIK3CB 99/100，PIK3CA/mTOR 48/48（Table S27）。一个 A-only 配体因计算超时而持续无法完成 PIK3CA 对接，因此从需要该分数的分析中剔除（Tables S27、S30）。AUROC 因此以 AutoDock Vina 能够处理的化合物为条件。同一组 Vina 姿态上的替代打分器未改变总体排序（Tables S14–S15；Figure S1B）。独立姿态生成见 Results 3.7。

### 3.3 配体性质与化学型解释了相当一部分表观信号

口袋匹配对接首先与四种预先定义的物化性质进行比较（Figure 4B；Table 2）。相对于每个靶对的最佳单一描述符参考，docking summary_min 的配对差值在 EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 中分别为 −0.052、−0.128、−0.122 和 +0.229；四个 95% 置信区间均包含 0（Table S19；Figure S3C）。即使 PIK3CA/mTOR 的点估计表现出最大的正向差异，现有样本仍不足以将其与配体性质参考明确区分。

AChE/BChE 提供了一个较为直接的混淆案例。dual 配体平均 TPSA 约为 75，而选择性硬负配体约为 51（Figure 4C）；TPSA 单独获得约 0.769 的 AUROC，高于相同比较下的 Vina（约 0.56）。进一步加入重原子数和 TPSA 后，dual-versus-B-only AUROC 从 0.606 增至 0.807，而对接分数的优势比仅约为 1.18（Figure 7C）。该结果表明，该方向上的部分表观对接判别可由配体物化信息解释。

PIK3CA/mTOR 的情况有所不同。加入重原子数和 TPSA 后，AUROC 的变化约为 +0.07 至 +0.11，对接分数的优势比约为 2.19 和 3.08。与描述符的配对差值置信区间仍包含 0。配体效率归一后，仅 PIK3CA/mTOR 仍高于重原子数基线（0.657 对 0.463）。

二维化学结构基线进一步说明了这一问题（Figure 7A）。ECFP4 逻辑回归在 Bemis–Murcko 支架 GroupKFold 下多个方向获得约 0.78–0.91 的 fold AUROC，明显高于部分对应对接对照，例如 EGFR/HER2 dual-versus-B-only 中 ECFP4 AUROC 为 0.89，而对接 AUROC 仅为 0.43。实验标签因此包含可由配体二维结构表征捕获的信息，即使不使用受体信息也可以获得较强判别。PIK3CA/mTOR 上 \(n_{\mathrm{scaffolds}} \approx n\)，该折接近 leave-one-scaffold。同一设定下随机 `StratifiedKFold` 相对支架折的平均差为 +0.026（Table S20；Figure S3D）。

在当前支架分组基准下，加入口袋匹配对接分数后的最大绝对变化为 0.020（未四舍五入值：PIK3CA/mTOR dual versus A-only 的 −0.0198），且若干方向为负（Table S24）。T ≥ 0.3 时，PIK3CA/PIK3CB dual versus A-only 从 0.691 降至 0.503（n_neg = 11），远缘硬负（T < 0.3）为 0.819；T ≥ 0.4/0.5 的格子常为 n_neg ≤ 7，T ≥ 0.7 为空（Table S23）。效价匹配或尺寸匹配子集上，EGFR/HER2 与 PIK3CA/PIK3CB 的 dual 对 B_only 仍偏弱或接近随机（约 0.45–0.52），各臂 n 常低于 15（Table S5；Figure 7D）。

### 3.4 活性聚合、配体面板与受体实现影响评价结果

主标签使用可用 pChEMBL 的最大值。将最大值换成重复测定中位数后，改变了 7/110、1/95、1/99 和 0/48 个配体状态分配，标签一致率分别为 93.6%、98.9%、99.0% 和 100%。API 重拉标签上的靶对 `summary_min` 仅适度变化（0.417→0.424、0.606→0.629、0.500→0.500 和 0.692→0.692；Table S29）。

将 exhaustiveness 从 16 降至 8 后，PIK3CA/mTOR 的 summary_min 从 0.692 降至 0.660（Figure S1D）。在 PM110 面板中（n_scored = 115；dual / A_only / B_only 各 30），Vina summary_min 为 0.648 [0.51, 0.76]，相比 PM48 下降约 0.04，排序趋势保持一致（Figure S1C）。在未使用配体池留出集中（每对 20 / 20 / 20，种子 20260731；EGFR/HER2 不具备同等配额），PIK3CA/mTOR 的 summary_min 为 0.765 [0.603, 0.891]，AChE/BChE 为 0.618 [0.422, 0.759]，PIK3CA/PIK3CB 则下降至 0.425 [0.241, 0.618]（Tables S8、S16）。该留出集共享同一 ChEMBL 抓取批次。

我们进一步测试受体实现：一端受体冻结，只替换另一端（Figure 5；Tables S9、S30）。三个替代晶体均通过共晶重对接 QC（best-of-9 RMSD 分别为 4JPS 0.607 Å、5DXT 0.624 Å、4JSX 0.515 Å）。在 PIK3CA/mTOR 上，当 PIK3CA 4L23 替换为 4JPS 或 5DXT、mTOR 4JT6 保持不变时，PM48 的 summary_min 分别由 0.692 降至 0.486 [0.259, 0.692] 和 0.505 [0.292, 0.696]（Figure 5A）。变化主要发生在依赖替代 PIK3CA 结构的 D/B 方向，D/A 保持 0.714。将 mTOR 4JT6 替换为 4JSX 后 summary_min 为 0.639 [0.418, 0.776]。

同一套 PIK3CA 晶体再用于 PIK3CA/PIK3CB 面板、2WXF 保持冻结时，summary_min **上升**：由 0.500 至 0.691 [0.516, 0.779]（4JPS）和 0.685 [0.506, 0.768]（5DXT）（Figure 5B）。dual versus A-only 保持 0.691；dual versus B-only 由 0.500 升至 0.707 和 0.685。三种 PIK3CA 条件使用相同的 99 个已打分配体（Table S30）。因此，替换同一 PIK3CA 受体在两个靶对上使表观判别向相反方向变化。

### 3.5 错口袋对照揭示面板外尚未解决的失败模式

在主面板中，口袋匹配 summary_min 均高于错口袋对照，四对的 matched-minus-wrong 差值分别为 0.170、0.161、0.151 和 0.090；其中 EGFR/HER2 和 AChE/BChE 的差异置信区间排除 0，PIK3CA/PIK3CB 与 PIK3CA/mTOR 的区间包含 0（Tables S6、S17；Figure 6A；Figure S3A）。错口袋对照的 summary_min 分别为 0.260、0.444、0.349 与 0.602。

未使用配体池留出集中的点估计关系发生反转（Figure 6B）。PIK3CA/mTOR、AChE/BChE 和 PIK3CA/PIK3CB 的错口袋 summary_min 分别为 0.788、0.643 和 0.520，而匹配口袋分别为 0.765、0.618 和 0.425。相应的 matched-minus-wrong 点差值均为负（−0.023 / −0.025 / −0.095），但其 95% 置信区间均包含 0（Table S17；Figure S3B）。效价或尺寸匹配后，错口袋仍 ≥ 匹配口袋（Table S13）。不依赖打分的接触计数在 B 方向获得 0.698–0.714 的 AUROC，但这一粗粒度替代不能解释观察到的 Vina 错口袋判别幅度（Figure 6D；Table S11）。因此，留出集中的错口袋反转是当前基准暴露出的尚未解决的失败模式。

### 3.6 探索性结构背景

结构相似性并不能预测筛选判别：5DXT 与 4L23 的局部口袋 Cα RMSD 仅 0.343 Å，但二者的基准表现仍有明显差异（Table S10）。在 PM48 姿态上，以重原子距离 ≤ 4.5 Å 对 20 个冻结 PIK3CA 口袋残基做几何占有率快照，占有率变化最大的残基为 Met772、Leu807、Gln859、Thr856、Cys838、Glu849、Phe930 和 Asp933（Table S33）。表观判别位移与这些接触模式变化同时出现，只提供受体敏感性的结构假说。占有率变化不能证明某一残基导致了 AUROC 改变，也不能解释 PIK3CA/PIK3CB 上的相反位移。姿态类型学保留在 Supporting Information（Table S7；Note S1）。

### 3.7 独立姿态生成仍保留 EGFR/HER2 的设定差距

Results 3.2 中的 GNINA CNN 与 RTMScore 是对冻结 Vina 姿态的重打分。为检验设定效应是否依赖于该姿态生成引擎，在同一套冻结 EGFR/HER2 与 PIK3CA/mTOR 配体、受体、对接盒和 exhaustiveness 上以 GNINA 1.3.2 对接搜索模式独立生成姿态（Methods 2.3；Table S32）。分析仅使用两端均有分数的配体（EGFR/HER2 在 neither 配体 EH120_109 两端失败后 n = 109；PIK3CA/mTOR 在 A-only 配体 PM48_19 两端失败后 n = 47）。方向性 n_scored 在 EGFR/HER2 上仍为 28 / 38 / 32，在 PIK3CA/mTOR 上变为 18 / 13 / 12。

EGFR/HER2 上 Dual versus neither 仍高（AUROC 0.783 [0.610, 0.922]；n_neg = 11），而方向性 `summary_min` 为 0.220 [0.109, 0.343]。按口袋均分取混合库 Top-10 仍含 1 个 dual 与 9 个实验选择性配体（4 个 A-only、5 个 B-only；EF10 = 0.389）。因此，在该独立姿态生成协议下，设定效应不是 Vina 特有的。PIK3CA/mTOR 的方向性 `summary_min` 符号不变（0.633 对 Vina 0.692）。该比较只检验设定差距是否仍在，不是主张 GNINA 优于或劣于 Vina。
