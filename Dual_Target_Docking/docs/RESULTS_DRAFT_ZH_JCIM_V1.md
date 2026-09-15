## 3. 结果

### 3.1 成对实验数据供给与对接评价集

能够支持四状态对接评价的成对实验标签随样本要求提高而迅速减少。在人源单组分 SINGLE PROTEIN 靶点的无序对中，2,164,618 对至少有 1 个两端均被测定的配体，63,790 对具有不少于 10 个此类配体。若再要求在 \(\theta=6.0\) 下 dual、A-only 和 B-only 均不少于 10 个，仅剩 5,253 对；采用严格 6.5/5.5 活性标准后减至 86 对。Figure 1C 概括这些计数。后续结构和对接位点门槛见 2.3 节。

经过供给、结构和对接相容性检查后，主要评价包括 EGFR/HER2、JAK1/JAK2、JAK1/TYK2、PIK3CA/mTOR、AChE/BChE、F2/F10、PPARG/PPARA 和 PPARA/PPARD（Table 1）。EGFR/HER2 与 PIK3CA/mTOR 从 \(\theta=6.0\) 候选池抽样，其余六对从严格 6.5/5.5 池抽出。四状态分类要求同一配体在两个靶点上均有实验测量，因此成对覆盖共同限制了对接评价集规模。

![Figure 1](../figures/jcim_article/Fig1_four_state_and_supply.png)

**Figure 1.** 四状态双靶评价与数据供给。(A) 按阈值 \(\theta\) 定义的四种实验状态；(B) 两个方向性评价任务：dual vs A-only 使用靶点 B 评分，dual vs B-only 使用靶点 A 评分，\(\mathrm{summary}_{\min}=\min[\mathrm{AUROC}_{D/A}(B),\;\mathrm{AUROC}_{D/B}(A)]\)；(C) 成对实验覆盖在四状态要求下迅速减少。主要评价的八个靶对按 Table 1 的面板构建与结构条件纳入，不是普查计数的直接延续。

### 3.2 不同虚筛对照下的方向性对接表现

固定同一靶点对接评分后，更换实验状态对照所产生的 AUROC 差异因靶对和方向而异。EGFR/HER2 与 JAK1/TYK2 的固定评分差最大。EGFR/HER2 使用 EGFR 口袋评分时，dual–B-only 的 AUROC 为 0.430，将对照换成 neither 后升至 0.808，差值 0.378 [0.205, 0.547]。JAK1/TYK2 使用 JAK1 口袋评分时，对应差值为 0.444 [0.263, 0.620]。其余靶对的差异更小，或其置信区间包含 0（Figure 2A；Table S4）。簇重采样见 3.5 节。

在统一 \(\theta=6.0\) 标签下，两个方向性 AUROC 介于 0.345 至 0.728，描述性较弱臂 \(\mathrm{summary}_{\min}\) 介于 0.345 至 0.692。PPARG/PPARA 的 \(\mathrm{summary}_{\min}\) 区间完全高于 0.5（0.649 [0.504, 0.751]），F2/F10 完全低于 0.5（0.345 [0.211, 0.477]），其余六对穿过 0.5（Figure 2B；Table 2）。这些数值是靶对特异的，不是八对排行。

**Table 2.** 八个主评价靶对上的口袋匹配方向 AUROC（Vina，统一 \(\theta=6.0\)）。表中类别样本量为 n_scored（dual / A-only / B-only）。物化描述符基线见 Table S5。

| 靶对 | n_scored (dual / A-only / B-only) | Dual vs A-only（口袋 B）[95% CI] | Dual vs B-only（口袋 A）[95% CI] | summary_min [95% CI] |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 [0.524, 0.793] | 0.430 [0.282, 0.579] | 0.430 [0.282, 0.578] |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 [0.444, 0.729] | 0.728 [0.595, 0.848] | 0.588 [0.444, 0.725] |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 [0.434, 0.725] | 0.365 [0.231, 0.505] | 0.365 [0.231, 0.503] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 [0.506, 0.899] | 0.692 [0.495, 0.874] | 0.692 [0.470, 0.813] |
| AChE/BChE | 27 / 25 / 28 | 0.650 [0.483, 0.801] | 0.606 [0.442, 0.751] | 0.606 [0.437, 0.730] |
| F2/F10 | 31 / 32 / 32 | 0.413 [0.259, 0.562] | 0.345 [0.214, 0.486] | 0.345 [0.211, 0.477] |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 [0.507, 0.778] | 0.706 [0.569, 0.833] | 0.649 [0.504, 0.751] |
| PPARA/PPARD | 32 / 32 / 32 | 0.646 [0.504, 0.776] | 0.446 [0.296, 0.586] | 0.446 [0.296, 0.584] |

采用双靶虚拟筛选常用的双口袋平均排序时，dual–neither 比较 AUROC 介于 0.514 至 0.770。EGFR/HER2 为 0.756 [0.562, 0.920]，JAK1/TYK2 为 0.770 [0.597, 0.906]，而其方向性 \(\mathrm{summary}_{\min}\) 分别为 0.430 和 0.365（Figure 2C；Table 3）。PIK3CA/mTOR 的 neither 样本为 \(n=4\)。Table 3 同时改变了评分聚合和对照类别，对照类别的影响以 Figure 2A 的固定口袋分析为准。

![Figure 2](../figures/jcim_article/Fig2_negative_class_formulation.png)

**Figure 2.** 对接表现取决于实验状态比较。(A) 固定评分通道后更换对照类别的 \(\Delta\)AUROC（dual–neither 减去 dual–单靶活性），误差棒为配体水平 bootstrap 95% 置信区间；(B) 两个方向性 AUROC；(C) 方向性 \(\mathrm{summary}_{\min}\) 与 dual–neither 比较。菱形标出 PIK3CA/mTOR 的 neither n = 4。面板 C 不是固定评分通道比较。(D) JAK1/TYK2 在双口袋平均 Vina 排序下的 Top-10 类别组成（n = 109）。

**Table 3.** 同一套 Vina 对接分数在方向性与 dual–neither 比较设定下的 AUROC（统一 \(\theta=6.0\)）。dual–neither 比较采用双口袋平均评分 \(S_{\mathrm{mean}}\)。PIK3CA/mTOR 的 neither 样本量（n = 4）较少。

| 靶对 | directional summary_min [95% CI] | Dual vs neither (vina_mean) | n_neither |
|------|--------------------------------:|------------------------------:|----------:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 |
| JAK1/JAK2 | 0.588 [0.444, 0.725] | 0.730 [0.547, 0.875] | 14 |
| JAK1/TYK2 | 0.365 [0.231, 0.503] | 0.770 [0.597, 0.906] | 14 |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 |
| F2/F10 | 0.345 [0.211, 0.477] | 0.519 [0.350, 0.688] | 12 |
| PPARG/PPARA | 0.649 [0.504, 0.751] | 0.685 [0.493, 0.848] | 14 |
| PPARA/PPARD | 0.446 [0.296, 0.584] | 0.565 [0.368, 0.766] | 14 |

该传统排序在 JAK1/TYK2 评价集上按操作点检查。按双口袋平均评分对全部 109 个配体排序时，Top-10 含 1 个 dual、2 个 A-only 和 7 个 B-only，没有 neither 配体（Figure 2D）。较高的 dual–neither AUROC 并未对应较少的高排名单靶活性配体。按 dual 的中位 \(S_{\mathrm{worst}}\) 对 Dual+A-only+B-only（\(n=95\)）作联合过滤、排除 neither 后，保留 16 个 dual，同时保留 12 个 A-only 和 22 个 B-only，dual precision 为 0.32（Figure S1）。

同一套双口袋平均排序和 AND 过滤应用于全部八个主评价集（Table S10）。Top-10 是固定条数，不是固定筛选比例：EGFR/HER2 为 10/110，PIK3CA/mTOR 为 10/48。Top-10 中 dual 计数从 1（EGFR/HER2、JAK1/TYK2）到 6（PIK3CA/mTOR、PPARG/PPARA）。四个靶对的 Top-10 不含 neither，另外四对各占 1 个位置。每个靶对的 Top-10 中都至少保留 1 个单靶活性配体。这些计数是回顾性评价面板上的描述性操作点，不是对接质量的八对排行。

### 3.3 配体化学基线与对接增量判别

若双靶对接所排序的实验类别本身已有化学差异，则不依赖受体的基线也能分开同一套筛选状态。AChE/BChE 仅用 TPSA 时，dual–A-only 和 dual–B-only 的 AUROC 分别为 0.733 和 0.801（Figure 3C；Table S5）。PIK3CA/mTOR 中最佳单一描述符（重原子数）的 \(\mathrm{summary}_{\min}\) 为 0.463（Table S5）。Vina 与最佳单一描述符的差值因靶对而异：八对中六对的 95% 区间包含 0，F2/F10 与 JAK1/TYK2 不包含 0（Table S5）。Vina \(\mathrm{summary}_{\min}\) 区间与最佳描述符点估计见 Figure S2。

在 Bemis–Murcko 骨架分组交叉验证下，仅含配体的 ECFP4 模型已能分开若干方向性比较（Figure 3A）。将对应方向的对接评分加入后，16 个方向的 AUROC 最多变化 0.023，且没有一致的变化方向（Figure 3B；Table S5）。在该回顾性设定下，对接评分未在配体模型之上提供稳定的排序增量。

![Figure 3](../figures/jcim_article/Fig3_ligand_chemistry.png)

**Figure 3.** 配体化学作为竞争解释，不是公平的预测竞赛。(A) 原始分数排序的 Vina AUROC 与骨架分组交叉验证下 ECFP4 折外预测 AUROC；蓝为 Vina，橙为 ECFP4；(B) 在相同划分下将对应方向 Vina 评分加入 ECFP4 后的 AUROC 变化，连线连接同一方向的两个点，不是置信区间；(C) AChE/BChE 的 TPSA 分布。

### 3.4 口袋对应性与对接实现敏感性

若方向性 AUROC 来自基于结构的对接，则在不重新对接的前提下交换两条已计算的口袋评分，较弱臂汇总应变差。主评价集中仅 EGFR/HER2 和 AChE/BChE 的 matched−mismatched \(\mathrm{summary}_{\min}\) 95% 区间排除 0，差值分别为 0.170 [0.060, 0.280] 和 0.161 [0.037, 0.269]。其余六对主集区间包含 0。七个可用内部留出集的差值区间亦均包含 0（−0.079 至 +0.150）（Figure 4A；Table S6）。EGFR/HER2 无未使用池留出集。

![Figure 4](../figures/jcim_article/Fig4_mismatched_pocket.png)

**Figure 4.** 对应口袋与非对应口袋评分对照。(A) 主评价集与留出集的 matched−mismatched \(\Delta\mathrm{summary}_{\min}\)，点与横线为估计值和配体水平 bootstrap 95% 置信区间；(B) 主评价集与留出集的 \(\mathrm{summary}_{\min}\)。\(\dagger\) 表示 EGFR/HER2 无未使用池留出集。底部图例区分主评价集与留出集。

独立 GNINA 1.3.2 姿态生成使用与主分析相同的受体、配体和对接盒，覆盖三对（Figure 5A；Table S7）。EGFR/HER2 的 dual–neither AUROC 为 0.783 [0.610, 0.922]（\(n_{\mathrm{neither}}=11\)），dual–B-only 为 0.220 [0.109, 0.343]。JAK1/TYK2 的 dual–neither 为 0.705 [0.517, 0.876]，方向性 \(\mathrm{summary}_{\min}\) 为 0.317 [0.183, 0.463]。PIK3CA/mTOR 的 \(\mathrm{summary}_{\min}\) 为 0.633，较弱臂 dual–A-only 为 0.633 [0.427, 0.825]，dual–neither 为 0.569 [0.222, 0.889]（\(n=18/13/12/4\)）。这些运行是另一套姿态生成流程，不是对 Vina 口袋对调结果的确认。

在 PIK3CA/mTOR 中，将 PIK3CA 受体由 4L23 替换为 4JPS 后，\(\mathrm{summary}_{\min}\) 从 0.692 [0.470, 0.813] 降至 0.486 [0.259, 0.692]；替换为 5DXT 后为 0.505 [0.292, 0.696]；将 mTOR 4JT6 替换为 4JSX 后为 0.639 [0.418, 0.776]（Figure 5B；Table S7）。

五个固定 Vina 随机种子产生的数值波动相对任务和受体效应更为有限（Figure 5C）。EGFR/HER2 的固定评分任务差在五个种子上均为正。

PPARG/PPARA 是唯一主分析 Vina \(\mathrm{summary}_{\min}\) 区间完全高于 0.5 的靶对（0.649 [0.504, 0.751]）。对同一套 Vina 姿态作 RTMScore 重评分后降至 0.369 [0.233, 0.475]，GNINA CNN 重评分降至 0.500 [0.356, 0.623]（Table S7）。

PIK3CA/mTOR 的面板规模和 exhaustiveness 敏感性结果见 Figure S3。PM48 为主评价集（配额 n = 48，exhaustiveness = 16）；PM110 为同一靶对的更大协议敏感性面板。

14 个主受体在统一化学对应 CalcRMS 表中均至少有一个保存姿态低于 2.0 Å。该结果只是搜索覆盖检查。AChE 4EY7 与 TYK2 3LXP 保存 8 个姿态，其余主受体保存 9 个姿态。EGFR 3POZ 的 top-1 为 9.505 Å，最低保存姿态为 0.760 Å。BChE、mTOR、JAK2、PPARG 和 PPARA 的 top-1 也未回到 2 Å 以内；PPARA 6LXA 的 top-3 最低为 7.848 Å（Figure S4；Table S2）。

![Figure 5](../figures/jcim_article/Fig5_computational_realization.png)

**Figure 5.** 计算实现敏感性。(A) 独立 GNINA 与 Vina：实心为 Vina，空心为 GNINA；蓝为方向性 \(\mathrm{summary}_{\min}\)，橙为 dual–neither；灰线连接同一引擎的两项任务，不是置信区间。(B) PIK3CA/mTOR 受体替换，误差棒为配体水平 bootstrap 95% 置信区间。(C) 八个靶对在五个 Vina 随机种子上的 \(\mathrm{summary}_{\min}\) 范围、中位数与生产种子；该面板不展示设定差距。

### 3.5 方向性评价的标签与样本组成敏感性

从严格 6.5/5.5 候选池构建的评价集，在多个阈值下保留了相同的主要类别组成，因此相应 AUROC 变化较小。EGFR/HER2 和 PIK3CA/mTOR 的类别组成随阈值变化更明显，其方向性估计也发生变化（Figure 6A；Table S3）。在同一套 ChEMBL 37 记录上将最大 pChEMBL 改为中位数后，EGFR/HER2 有 6/110 个类别翻转（\(\mathrm{summary}_{\min}\) 0.430 变为 0.424），AChE/BChE 有 1/95 个翻转（CHEMBL659；0.606 变为 0.629），PPARA/PPARD 有 1/110 个翻转（CHEMBL121；\(\mathrm{summary}_{\min}\) 仍为 0.446）。其余五对保持类别组成和 \(\mathrm{summary}_{\min}\) 点估计（Table S3）。

在排除主评价集分子后，基于剩余候选分子构建的未使用池留出集显示，结果存在一定的样本组成依赖：AChE/BChE、PIK3CA/mTOR 与 JAK1/JAK2 与主评价接近；JAK1/TYK2 有所上升；F2/F10 与 PPARA/PPARD 仍处于较低水平；PPARG/PPARA 则由 0.649 降至 0.535 [0.350, 0.717]（Figure 4B；Table S6）。EGFR/HER2 无同等留出集。固定评分通道差值的簇重采样中，EGFR/HER2 骨架簇区间为 [0.168, 0.562]、文献簇区间为 [0.083, 0.529]，均排除 0；JAK1/TYK2 骨架簇区间为 [0.234, 0.633]，排除 0，文献簇区间为 [−0.034, 0.682]，包含 0（Figure 6B；Table S4）。

### 3.6 独立外部对接集的可用性

对全部八个靶对检索了 BindingDB[16] 与 PubChem。按 \(\theta=6.0\) 独立性过滤剔除共享文献、重复结构以及与开发集 ECFP4 Tanimoto 相似度 \(\geq 0.70\) 的分子后，没有靶对满足独立外部评价准入标准（dual、A-only 和 B-only 各类至少 20 个配体且每类至少 3 个独立来源）。因此未形成外部对接集；该检索不作为外部验证（Figure 6C,D；Table S8）。

![Figure 6](../figures/jcim_article/Fig6_evidence_boundary.png)

**Figure 6.** 证据边界。(A) 活性阈值；(B) 配体、骨架簇和文献簇重采样下，靶点 A 评分的 dual–neither 与 dual–B-only 差值；(C) BindingDB 过滤后各类分子数，颜色饱和于每类 n = 20；(D) BindingDB 过滤后独立来源数，颜色饱和于每类 3 个来源。\(\dagger\) 表示该类 n < 10。
