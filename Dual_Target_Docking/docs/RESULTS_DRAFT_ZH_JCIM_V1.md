# Results（中文工作稿 · JCIM Articles）

## 3. 结果

### 3.1 双端实验数据供给与四状态评价集构建

ChEMBL 中能够支持四状态评价的双端测量数据随样本要求提高而迅速减少。在所有至少存在双端测量的靶对中，2,164,618 对至少有 1 个双端测量配体，63,790 对具有不少于 10 个双端测量配体。若再要求在 \(\theta=6.0\) 下 dual、A-only 和 B-only 均不少于 10 个，仅剩 5,253 对；采用严格 6.5/5.5 双向选择性供给标准后，减至 86 对。这一数据供给变化见 Figure 1C。

上述普查用于描述能够支持双方向评价的数据供给情况。按照 2.3 节的评价面板与结构准入条件，最终主要评价包括 EGFR/HER2、JAK1/JAK2、JAK1/TYK2、PIK3CA/mTOR、AChE/BChE、F2/F10、PPARG/PPARA 和 PPARA/PPARD（Table 1）。各靶对的候选池规则、配额和骨架上限按样本供给规模与结构可行性分别设定（Table 1）。EGFR/HER2 与 PIK3CA/mTOR 的面板候选池采用 \(\theta=6.0\)，其余六对采用严格 6.5/5.5 规则。四状态分类要求同一配体在两个靶点上均具有实验测量，双端实验测量覆盖和双向选择性样本数量共同限制了严格四状态评价集的规模。

![Figure 1](../figures/jcim_article/Fig1_four_state_and_supply.png)

**Figure 1.** 四状态双靶评价与数据供给。(A) 按阈值 \(\theta\) 定义的四种实验状态；(B) 两个方向性评价任务；(C) ChEMBL 数据供给由 2,164,618 对减至 86 对，主要评价保留八个靶对。

### 3.2 实验状态定义与方向性对接评价

固定相同评分通道后，改变实验状态比较所产生的 AUROC 差异在不同靶对和评价方向之间并不一致。其中 EGFR/HER2 和 JAK1/TYK2 出现了较明显的差异。EGFR/HER2 使用靶点 A（EGFR）评分时，dual–B-only 比较的 AUROC 为 0.430，将对照类别替换为双端低活性的 neither 后升至 0.808，差值 0.378 [0.205, 0.547]；另一方向上的差异较小。JAK1/TYK2 使用靶点 A（JAK1）评分时，两项比较的差值为 0.444 [0.263, 0.620]。EGFR/HER2 在骨架簇和文献簇重采样下区间均排除 0；JAK1/TYK2 在骨架簇重采样下保持，而文献簇区间包含 0（Figure 2A；Figure 6B；Tables S4 and S10）。其余靶对和评价方向的差异总体较小或不确定，多数配体水平 95% 置信区间包含 0（Figure 2A；Table S4）。

在统一 \(\theta=6.0\) 标签下，八个靶对方向性判别的较弱方向描述性汇总指标 \(\mathrm{summary}_{\min}\) 介于 0.345 至 0.692。仅 PPARG/PPARA 的配体水平 bootstrap 95% 置信区间完全高于 0.5（0.649 [0.504, 0.751]）；其余靶对的区间包含 0.5，或整体位于 0.5 以下（Figure 2B；Table 2）。

**Table 2.** 八个主评价靶对上的口袋匹配方向 AUROC（Vina，统一 \(\theta=6.0\)）。表中类别样本量为 n_scored（dual / A-only / B-only）。物化描述符基线见 Table S5。

| 靶对 | n_scored (dual / A-only / B-only) | Dual vs A-only（口袋 B） | Dual vs B-only（口袋 A） | summary_min [95% CI] |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 | 0.728 | 0.588 [0.444, 0.725] |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 | 0.365 | 0.365 [0.231, 0.503] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] |
| F2/F10 | 31 / 32 / 32 | 0.413 | 0.345 | 0.345 [0.211, 0.477] |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 | 0.706 | 0.649 [0.504, 0.751] |
| PPARA/PPARD | 32 / 32 / 32 | 0.646 | 0.446 | 0.446 [0.296, 0.584] |

采用双口袋平均评分作描述性比较时，EGFR/HER2 和 JAK1/TYK2 的 dual–neither 比较 AUROC 分别达到 0.756 [0.562, 0.920] 和 0.770 [0.597, 0.906]，而其方向性 \(\mathrm{summary}_{\min}\) 分别为 0.430 和 0.365（Figure 2C；Table 3）。该比较同时改变了评分聚合形式和对照定义，对照组成的单独影响以固定评分通道分析为准。

![Figure 2](../figures/jcim_article/Fig2_negative_class_formulation.png)

**Figure 2.** 对接表现取决于实验状态比较。(A) 固定评分通道后更换对照类别的 \(\Delta\)AUROC；(B) 两个方向性 AUROC；(C) 方向性 \(\mathrm{summary}_{\min}\) 与 dual–neither 比较。菱形标出 PIK3CA/mTOR 的 neither n = 4。

**Table 3.** 同一套 Vina 对接分数在方向性与 dual–neither 比较设定下的 AUROC（统一 \(\theta=6.0\)）。dual–neither 比较采用双口袋平均评分 \(S_{\mathrm{mean}}\)。PIK3CA/mTOR 的 neither 样本量（n = 4）较少。Dual versus all non-duals 见 Table S4。

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

EGFR/HER2 评价集共包含 110 个配体，其中 28 个为 dual。按双口袋平均评分排序时，Top-10 中包含 1 个 dual、5 个 A-only 和 4 个 B-only，没有 neither 配体（Figure S1A；Table S13）。较高的 dual–neither 比较 AUROC 并未对应较少的高排名单靶选择性配体。按 dual 的中位 \(S_{\mathrm{worst}}\) 进行联合过滤时，保留 14 个 dual，同时保留 9 个 A-only 和 24 个 B-only，dual precision 为 0.298（Figure S1B；Table S13）。

### 3.3 配体化学基线与方向性判别

部分靶对中，配体自身的物理化学特征即可对 dual 和单靶选择性配体产生一定区分。AChE/BChE 仅使用 TPSA 时，dual–A-only 比较和 dual–B-only 比较的 AUROC 分别为 0.733 和 0.801（Figure 3C；Table S5）。单一物化性质的判别能力在不同靶对之间差异明显，PIK3CA/mTOR 中最佳单一描述符（重原子数）的 \(\mathrm{summary}_{\min}\) 为 0.463（Table S5）。在八个靶对中，Vina 与最佳单一描述符之间的差异方向和幅度因靶对而异，Vina 并未在各靶对上一致优于单一物化描述符；多数差值的 95% 置信区间包含 0（Figure S2；Table S5）。

在 Bemis–Murcko 骨架分组交叉验证下，不使用受体结构的 ECFP4 分子指纹在多个方向上获得与对接相当甚至更高的 AUROC（Figure 3A）。将对应方向的对接评分加入 ECFP4 逻辑回归模型后，16 个方向的 AUROC 最大绝对变化为 0.023（Figure 3B；Table S5）。在当前评价集和交叉验证设置下，未观察到加入对接评分后跨方向一致的 AUROC 提升。

![Figure 3](../figures/jcim_article/Fig3_ligand_chemistry.png)

**Figure 3.** 配体化学作为竞争解释。(A) Vina 与 ECFP4；(B) 加入对应方向 Vina 评分后的 AUROC 变化；(C) AChE/BChE 的 TPSA 分布。

### 3.4 口袋对应性与计算敏感性

为考察方向性判别是否与实验活性差异所在的靶点口袋相对应，对应口袋与非对应口袋的 \(\mathrm{summary}_{\min}\) 差值见 Figure 4A 与 Tables S6 and S7。主评价集中仅 EGFR/HER2（0.170 [0.060, 0.280]）和 AChE/BChE（0.161 [0.037, 0.269]）两对的差值区间排除 0；其余六对的 95% 置信区间均包含 0。在七个未使用的留出集中，对应口袋与非对应口袋的 \(\mathrm{summary}_{\min}\) 差值区间均包含 0（差值介于 −0.079 至 +0.150；Figure 4A；Table S7）。对应口袋优势未在留出集中稳定重现。

![Figure 4](../figures/jcim_article/Fig4_mismatched_pocket.png)

**Figure 4.** 对应口袋与非对应口袋评分对照。(A) 主评价集与留出集的 matched−mismatched \(\Delta\)；(B) 主评价集与留出集的 \(\mathrm{summary}_{\min}\)。

在 PIK3CA/mTOR 中，将 PIK3CA 受体由 4L23 替换为 4JPS 后，\(\mathrm{summary}_{\min}\) 从 0.692 [0.470, 0.813] 降至 0.486 [0.259, 0.692]；替换为 5DXT 后为 0.505 [0.292, 0.696]；将 mTOR 4JT6 替换为 4JSX 后为 0.639 [0.418, 0.776]（Figure 5B；Table S8）。

采用 GNINA 1.3.2 独立生成姿态和评分后，EGFR/HER2 的 dual–neither 比较 AUROC 为 0.783 [0.610, 0.922]，n_neither = 11；方向性较弱臂 dual–B-only 比较为 0.220 [0.109, 0.343]（n_dual = 28，n_B = 32）。JAK1/TYK2 同样呈现该模式，dual–neither 比较为 0.705，方向性 \(\mathrm{summary}_{\min}\) 为 0.317 [0.183, 0.463]（Figure 5A；Table S9）。在独立姿态生成流程下，类似差异仍可观察到。

PPARG/PPARA 在主要 Vina 评价中是唯一 \(\mathrm{summary}_{\min}\) 置信区间完全高于 0.5 的靶对（0.649 [0.504, 0.751]），但同姿态 RTMScore 重评分降至 0.369 [0.233, 0.475]，GNINA CNN 重评分降至 0.500，未使用池留出集降至 0.535 [0.350, 0.717]（Table S7；Table S9）。五个固定 Vina 随机种子产生的数值波动相对有限，EGFR/HER2 的设定差距在五个 Vina 种子上均为正（Figure 5C；Table S9）。PIK3CA/mTOR 的面板规模和 exhaustiveness 敏感性结果见 Figure S3。共晶重对接中，全部主受体均可产生重原子 RMSD < 2.0 Å 的近天然姿态；其中 EGFR 3POZ 的 top-1 RMSD 为 9.505 Å，而 best-of-9 为 0.760 Å，说明该体系能够生成近天然姿态，但未将其排在首位（Figure S4；Table S2）。

![Figure 5](../figures/jcim_article/Fig5_computational_realization.png)

**Figure 5.** 计算实现敏感性。(A) 独立 GNINA 与 Vina；(B) PIK3CA/mTOR 受体替换；(C) 五个 Vina 随机种子。

### 3.5 方向性评价的标签与样本组成敏感性

改变活性阈值（\(\theta=5.5\)、6.0、6.5 及严格 6.5/5.5）后，部分样本较充足的靶对变化较小，例如 AChE/BChE 稳定在 0.606；当单侧选择性样本量减少时估计波动增大（Figure 6A；Table S3）。将重复活性记录由最大值改为中位数后，主要结果总体稳定；在可进行高置信记录复核的评价集中，更严格筛选也未改变相应方向性结果（Table S3）。

在排除主评价集分子后，基于剩余候选分子构建的未使用池留出集显示，结果存在一定的样本组成依赖：AChE/BChE、PIK3CA/mTOR 与 JAK1/JAK2 与主评价接近；JAK1/TYK2 有所上升；F2/F10 与 PPARA/PPARD 仍处于较低水平；PPARG/PPARA 则由 0.649 降至 0.535 [0.350, 0.717]（Figure 4B；Table S7）。七个可构建留出集的靶对结果见 Figure S5 和 Table S7。对于 EGFR/HER2，配体、骨架簇和文献簇重采样下的固定评分差值区间均排除 0；JAK1/TYK2 的文献簇重采样区间则包含 0（Figure 6B；Table S10）。

### 3.6 外部评价数据的可用性

BindingDB 与 PubChem 按相同四状态规则清点八个靶对的数据供给。随后对 BindingDB 排除共享来源、重复结构和高相似分子，并按外部评价准入条件进行独立来源筛选。经过独立性过滤后，没有靶对同时满足 dual、A-only 和 B-only 每类至少 20 个分子且至少来自 3 个独立来源的准入条件，因此未进行外部对接（Figure 6C,D；Table S11）。按配体发表年份进行时间切分同样受到双向选择性样本供给的限制，未能构建出可同时独立评价两个方向的时间切分测试集（Table S12）。未使用池留出集与年份切分都只是内部敏感性，不作为外部验证。在本研究采用的数据来源、独立性过滤规则和样本准入门槛下，未获得可同时独立评价两个选择性方向的外部数据集。

![Figure 6](../figures/jcim_article/Fig6_evidence_boundary.png)

**Figure 6.** 证据边界。(A) 活性阈值；(B) 配体、骨架簇和文献簇重采样；(C) BindingDB 过滤后各类分子数；(D) BindingDB 过滤后独立来源数。
