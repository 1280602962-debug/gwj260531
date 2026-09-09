# 基于分子对接的双靶识别：多靶对评价设定与判别来源分析

## 摘要

分子对接常用于双靶点虚拟筛选中候选分子的排序与优先选择。然而，在回顾性评价中，对照化合物既可以是在两个靶点上均低于活性阈值的化合物，也可以是仍对其中一个靶点保持活性的单靶选择性配体。本研究在八个人源靶对上建立包含 dual、A-only、B-only 和 neither 的四状态评价体系。在固定同一靶点评分通道的前提下，仅改变对照类别即可改变表观对接判别，而且这种影响在不同靶对和评价方向之间并不一致：EGFR/HER2 中 EGFR 口袋评分对 dual–B-only 比较的 AUROC 为 0.430，而在相同评分通道下与 neither 比较则升至 0.808（差值 0.378 [0.205, 0.547]）；JAK1/TYK2 同样出现类似差值（0.444 [0.263, 0.620]），且独立 GNINA 姿态生成再现了这一任务差异。基于骨架分组的仅配体 ECFP4 基线在多个方向上达到与对接相当的判别力，引入对应口袋对接评分后的最大 AUROC 增量仅为 0.023。在结构归因方面，仅 EGFR/HER2 和 AChE/BChE 主面板的对应口袋较弱方向评分优于非对应口袋，但在未使用的留出集分子中，对应口袋优势未能稳定重现。上述结果表明，在具有双端实验活性的化合物集中，dual 与双端低活性配体之间较好的判别，并不必然对应对单靶选择性配体的良好区分；回顾性双靶评价应分别考察两个方向的选择性对照，并结合配体化学基线与口袋对应性分析以厘清判别来源。

**关键词：** 双靶对接；实验状态；选择性；AutoDock Vina；虚拟筛选评价

## 1. 引言

多靶点药物设计用单一分子同时调控两个或多个疾病相关靶点。与多药联合相比，单一多靶配体有望减少不同药物之间的药代动力学不匹配及制剂复杂性。[1,2] 生成式设计与超大规模化合物库对接已经能够前瞻性得到经实验验证的双靶候选物：POLYGON 对 32 个 MEK1/mTOR 候选物进行了实验测试；[20] 超大规模库对接则从数亿级化学空间中发现具有多靶活性的配体，其中部分预测结合模式获得了实验结构支持。[19] 这些工作展示的是计算方法从候选库中找到双靶活性分子的能力，而不是在已有双端实验标签的化合物中，对接评分如何随实验状态比较而改变含义。

分子对接通过搜索配体在结合位点中的可能构象并对其评分，是虚拟筛选中常用的排序方法。[3,4] 许多基准把已知活性分子与人工诱饵构成二分类任务，例如 DUD 和 DUD-E；LIT-PCBA 则使用高通量实验确定的活性与非活性化合物。[5–7] 对照化合物的来源和定义可以改变虚拟筛选的评价结果，甚至影响不同方法之间的性能排序。[12,13,21] 因此，对照化合物的构成也是虚拟筛选性能评价的一部分。

对同时具有两端实验活性信息的化合物，双靶筛选可以形成 dual、A-only、B-only 和 neither 四种实验状态。dual 与 A-only 或 B-only 仅在其中一个靶点的实验状态上不同，而 dual 与 neither 则在两个靶点上均存在差异。因此，dual-versus-selective 与 dual-versus-neither 得到的对接判别并不具有完全相同的含义。双靶结构筛选通常将候选分子分别对接至两个靶点，并根据两端评分优先选择同时表现有利的分子。Zhou 等人曾在多个激酶靶对上评价基于分子对接的双靶虚拟筛选，并观察到单靶抑制剂可构成预测双靶抑制剂时的重要假阳性来源。[9] 该研究分析了双靶筛选中的不同类型假阳性，但没有将 dual 与 A-only、dual 与 B-only 分别作为两个与对应靶点相关的方向性判别任务来评价。Kinase-Bench 进一步表明，实验选择性配体可用于评价相近靶点之间的结构区分能力，其目标是识别对某一激酶具有选择性的抑制剂；这与在已测两端活性的化合物中区分 dual 与 A-only、B-only 不是同一问题。[22]

现有研究已经表明单靶选择性配体会影响双靶筛选结果，但在已有双端实验活性的数据中，评价任务的定义如何改变对接性能的解释，仍缺乏系统考察。基于此，本研究重点回答三个问题：第一，固定评分通道时，改变所比较配体的实验状态是否改变对双靶对接性能的判断；第二，表观观察到的判别能在多大程度上由不使用受体结构的配体化学特征解释；第三，剩余判别是否与相应靶点的口袋结构信息一致，并具有一定稳定性。本研究据此考察双靶对接性能评价的含义及其判别来源。

## 2. 方法

### 2.1 研究设计

配体按两个靶点上的实验活性划分为 dual、A-only、B-only 和 neither 四类。方向性评价分别比较 dual 与两类单靶选择性配体，并以 dual 为正类。dual 与 A-only 在靶点 A 上均达到活性阈值，实验状态差异位于靶点 B，因此使用靶点 B 的对接评分；dual 与 B-only 的实验状态差异位于靶点 A，因此使用靶点 A 的对接评分。

在保持同一靶点评分通道不变时，分别比较 dual 与单靶选择性配体、dual 与 neither。随后用不使用受体结构的配体化学基线考察该判别能否由配体自身化学特征解释，再用口袋对应、受体替换和独立对接实现检验剩余判别是否与相应靶点结构信息一致。数据处理、样本组成和外部数据可用性相关分析用于评估稳定性和适用范围。

### 2.2 生物活性数据处理与活性状态定义

实验活性数据来自 ChEMBL。IC50、Ki、Kd、EC50 和 Potency 等标准化定量记录均采用 pChEMBL，统一到负对数摩尔尺度。同一配体–靶点组合存在多条有效记录时，主要分析取最大 pChEMBL 作为代表值。盐、混合物等含多个组分的记录先按连通片段拆分，再保留重原子数最多的有机片段。仅将在两个靶点上均有可用实验活性数据的配体纳入四状态分类；缺少任一端实验记录的配体不视为该靶点低活性。

对任一靶对 A/B，以活性阈值 \(\theta\) 定义四种活性类别：dual，\(p_{\mathrm{A}}\geq\theta\) 且 \(p_{\mathrm{B}}\geq\theta\)；A-only，\(p_{\mathrm{A}}\geq\theta\) 且 \(p_{\mathrm{B}}<\theta\)；B-only，\(p_{\mathrm{A}}<\theta\) 且 \(p_{\mathrm{B}}\geq\theta\)；neither，两端均 \(<\theta\)。主要分析统一采用 \(\theta=6.0\)。

候选靶对的双向选择性供给另用严格 6.5/5.5 活性分离规则评估：活性端 pChEMBL \(\geq 6.5\)，低活性端 \(\leq 5.5\)，位于 5.5–6.5 灰区的配体不计入该严格分类。严格规则用于候选靶对供给评估；主要分析统一使用 \(\theta=6.0\)。

为考察活性数据处理方式的影响，分别以中位 pChEMBL 替代最大值，并在具有同日 API 快照的已打分配体上按更严格的人源单蛋白高置信记录重新确定实验状态；两项分析均保持评价集成员和对接评分不变（Table S3）。

### 2.3 靶对筛选与评价集构建

**靶对供给筛选。** 候选靶对来自同一批 ChEMBL 数据，其中靶点限定为人源、单组分 SINGLE PROTEIN、两端均有可用定量活性、以最大 pChEMBL 为代表值，并按 2.2 节划分为四状态。双向选择性供给另用严格 6.5/5.5 规则清点（活性端 \(\geq 6.5\)，低活性端 \(\leq 5.5\)，灰区不入该候选池）。最终是否纳入还取决于蛋白类别、结合位点性质、人源实验结构是否适合统一非共价对接，以及该靶对在上述规则下能否抽出可分析的四状态面板。据此纳入 PIK3CA/mTOR、AChE/BChE、EGFR/HER2、F2/F10、JAK1/TYK2、JAK1/JAK2、PPARG/PPARA 和 PPARA/PPARD。

**评价面板构建。** 因样本供给规模和结构可行性不同，部分靶对采用不同的候选池规则、类别配额或骨架限制；具体规则见 Table 1。EGFR/HER2 与 PIK3CA/mTOR 从 \(\theta=6.0\) 候选池抽样，其余六对从严格 6.5/5.5 候选池按配额抽出。

**统一主要分析。** 所有主要 AUROC 均按 \(\theta=6.0\) 重新确定实验状态，并仅纳入相应方向具有有效对接评分的配体，因此 n_scored 可以低于 n_panel。

**Table 1.** 双靶评价集的组成与主要对接设置。配额为构建目标（dual / A-only / B-only / neither）。n_panel 为面板成员数（含 neither）；n_scored 为两端均有有效 Vina 分数、进入方向性主 AUROC 的 dual / A-only / B-only 计数。PPARG 口袋 9V8H 中保留与受体复合的 PG08-NL 肽。

| 靶对 | 候选池 | 配额 (D / A / B / N) | 骨架上限 | 受体 PDB (A / B) | 分辨率 (Å) | n_panel | n_scored (dual / A-only / B-only) | Vina exhaustiveness |
|------|--------|---------------------:|:--------:|------------------|------------:|-------:|----------------------------------:|--------------------:|
| PIK3CA/mTOR | \(\theta=6.0\) | 18 / 14 / 12 / 4 | \(\leq 2\) | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | 严格 6.5/5.5 | 28 / 28 / 28 / 16 | 无 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| EGFR/HER2 | \(\theta=6.0\) | 28 / 38 / 32 / 12 | \(\leq 5\) | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |
| F2/F10 | 严格 6.5/5.5 | 32 / 32 / 32 / 14 | 无 | 4UDW / 2JKH | 1.16 / 1.25 | 110 | 31 / 32 / 32 | 8 |
| JAK1/TYK2 | 严格 6.5/5.5 | 32 / 32 / 32 / 14 | 无 | 6N7A / 3LXP | 1.33 / 1.65 | 110 | 31 / 32 / 32 | 8 |
| JAK1/JAK2 | 严格 6.5/5.5 | 32 / 32 / 32 / 14 | 无 | 6N7A / 8BXH | 1.33 / 1.30 | 110 | 32 / 32 / 32 | 8 |
| PPARG/PPARA | 严格 6.5/5.5 | 32 / 32 / 32 / 14 | 无 | 9V8H / 6LXA | 1.39 / 1.23 | 110 | 32 / 31 / 32 | 8 |
| PPARA/PPARD | 严格 6.5/5.5 | 32 / 32 / 32 / 14 | 无 | 6LXA / 5U3Q | 1.23 / 1.50 | 110 | 32 / 32 / 32 | 8 |

### 2.4 受体与配体准备及分子对接

#### 2.4.1 受体准备与对接区域定义

各靶点使用的人源实验结构及其共晶配体见 Table 1。相应共晶配体分别用于确定各受体的对接区域。JAK1 6N7A 同时用于 JAK1/TYK2 和 JAK1/JAK2，PPARA 6LXA 同时用于 PPARG/PPARA 和 PPARA/PPARD。TYK2 3LXP 使用 JH1 ATP 结合位点。PPARG 9V8H 中保留与受体复合的 PG08-NL 肽，仅移除共晶小分子和结晶水。

所有受体在使用前均核对蛋白身份、物种、结合位点和共晶配体位置。初始对接盒由共晶配体重原子的笛卡尔坐标范围确定，并沿 \(x\)、\(y\) 和 \(z\) 三个方向分别向两侧扩展 5 Å。若扩展后任一方向的边长仍小于 20 Å，则继续扩展至 20 Å。去除结晶水和用于确定对接区域的共晶配体后，使用 Meeko 准备受体并转换为 PDBQT。存在替代原子位置时，保留无 altLoc 标记或标记为 A 的原子。对接盒坐标见 Table S2。

#### 2.4.2 配体准备

以 2.2 节处理后的 ChEMBL SMILES 作为配体输入。使用 RDKit 添加氢原子，并采用 ETKDGv3 生成一个三维初始构象。随后使用 MMFF 力场进行局部几何优化，最大迭代次数为 200。优化后的配体使用 Meeko 转换为 PDBQT。三维构象生成使用固定随机种子（Table S1）。未进行系统的质子化状态和互变异构体枚举；未指定立体化学按所给 SMILES 由 RDKit 默认处理。

#### 2.4.3 AutoDock Vina 对接与评分

主要对接采用 AutoDock Vina 1.2.7 和默认 Vina scoring function。每个配体–受体组合最多输出 9 个姿态，`energy_range` 设为 3 kcal mol\(^{-1}\)。PIK3CA/mTOR 使用 exhaustiveness = 16，以获得满足重对接近天然构象覆盖要求的搜索设置；其余靶对使用 8。所有主要 Vina 分析均以排名第 1 姿态的 affinity 作为配体–受体评分。

#### 2.4.4 共晶配体重对接

在评价集对接前，对各主受体的共晶配体进行重对接，检查所定义的结合区域和搜索设置能否生成近天然构象。每个共晶配体最多输出 9 个姿态，并计算各保存姿态与实验共晶构象之间的重原子 RMSD。采用全部保存姿态中的最低 RMSD 评估对接搜索对近天然构象的覆盖能力。当该值低于 2.0 Å 时，认为搜索能够产生近天然构象。另报告排名第 1 姿态的 RMSD 和前 3 个姿态中的最低 RMSD，以区分姿态生成能力和评分排序能力。对接盒与共晶重对接结果同列于 Table S2。

#### 2.4.5 替代评分与独立对接

使用 RTMScore 和 GNINA 1.3.2 CNN 对 Vina 生成的全部可用姿态进行补充评分。RTMScore 取全部保存姿态中的最高评分作为配体级结果；GNINA 以 CNN affinity 作为主要 CNN 重评分读数，CNNscore 作为补充分析。上述结果均来自 Vina 已生成姿态的重新评分，不属于独立姿态生成。

此外，在 EGFR/HER2、PIK3CA/mTOR 和 JAK1/TYK2 上使用 GNINA 1.3.2 独立生成姿态。该分析使用与主要分析相同的受体、配体和对接区域，三个靶对的 exhaustiveness 分别为 8、16 和 8，每个配体最多输出 9 个姿态。配体级评分取排名第 1 姿态的 minimizedAffinity，并使用其相反数进行统计分析，使数值越高表示预测结合越有利。独立 GNINA 并未对每个配体都返回双端评分：EGFR/HER2 的 dual–neither 比较使用 n_neither = 11（缺 EH120_109），而 Vina Table 3 为 12；PIK3CA/mTOR 的 A-only 为 n = 13 而非 14；JAK1/TYK2 的 Dual/A-only/B-only/neither 计数为 30/32/29/14，而非 31/32/32/14。该分析用于考察主要观察在另一套姿态生成和评分流程下是否仍然存在，而不用于比较 GNINA 与 Vina 的总体性能（Table S9）。

### 2.5 评价指标与统计分析

#### 2.5.1 主要方向性判别指标

比较 dual 与 A-only 时，两类配体在靶点 A 上均达到实验活性阈值，因此使用靶点 B 的评分：\(\mathrm{AUC}_{D/A}(B)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{A\text{-}only};\;S_{B})\)。比较 dual 与 B-only 时，使用靶点 A 的评分：\(\mathrm{AUC}_{D/B}(A)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{B\text{-}only};\;S_{A})\)。所有方向性分析均以 dual 为正类。

AutoDock Vina 输出的 affinity 越低表示预测结合越有利。为统一 AUROC 的评分方向，将 Vina affinity 转换为 \(S_{\mathrm{Vina}}=-E_{\mathrm{Vina}}\)，因此较高的 \(S_{\mathrm{Vina}}\) 表示更有利的预测结合。两个方向的 AUROC 分别报告，并将两者中的较低值定义为 \(\mathrm{summary}_{\min}\)，仅用于描述较弱方向的总体表现，不作为新的配体评分函数：\(\mathrm{summary}_{\min}=\min(\mathrm{AUC}_{D/A}(B),\;\mathrm{AUC}_{D/B}(A))\)。

#### 2.5.2 不同对照类别与双口袋过滤的描述性比较

另以两个靶点均低活性的 neither 配体作为对照类别，评价对接评分区分 dual 与 neither 的能力。对于两个靶点均获得评分的配体，计算平均评分 \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\)，随后以 dual 为正类、neither 为对照计算 AUROC。将 A-only、B-only 和 neither 合并为对照得到的 dual 与所有非 dual 配体比较 AUROC 仅作为混合候选库的描述性参考，并在 Table 3 中报告。

为单独考察对照类别变化对判别结果的影响，保持所用靶点评分不变，仅改变对照组成。使用靶点 A 的评分时，分别计算 dual–B-only 比较和 dual–neither 比较；使用靶点 B 的评分时，分别计算 dual–A-only 比较和 dual–neither 比较（Table S4）。另取两个靶点评分中的较低值作为双口袋联合过滤评分 \(S_{\mathrm{worst}}=\min(S_{A},S_{B})\)，用于评价配体是否在两个靶点上均获得有利评分。以 dual 配体 \(S_{\mathrm{worst}}\) 分布的中位数作为双口袋过滤阈值，并统计不低于该阈值的 dual 配体数、dual recall、dual precision 以及保留的单靶选择性配体数（Table S13）。由于双口袋平均评分同时改变了对照组成和评分形式，对照类别变化的单独影响以固定同一口袋评分的比较为准。

#### 2.5.3 置信区间与重采样分析

Table 2 中 \(\mathrm{summary}_{\min}\) 的 95% 置信区间基于 2000 次配体水平的非分层百分位 bootstrap 估计。每次从参与该靶对方向性分析的 dual、A-only 和 B-only 联合配体库中有放回抽取与原集合相同数量的配体，重新计算 dual–A-only 比较和 dual–B-only 比较的 AUROC，并取两者较小值作为该次重采样的 \(\mathrm{summary}_{\min}\)。95% 置信区间由所有有效 \(\mathrm{summary}_{\min}\) 的第 2.5 和第 97.5 百分位数确定。表中点估计均由完整分析样本直接计算，而非 bootstrap 均值。dual–neither 比较的 AUROC 置信区间采用按类别重采样的百分位 bootstrap。

比较不同评分方法或对应与非对应口袋评分时，各方案在每次 bootstrap 中使用相同的配体重采样结果以保持配对。为评估化学骨架相关性和文献来源相关性的影响，另以 Bemis–Murcko 骨架簇和文献连通簇为单位进行 cluster bootstrap，算法细节见 Supporting Information（Table S4；Table S10；Figure S6）。Table S10 另报告一项按设计保持类别样本量的独立情景模拟；该模拟不是 Table 2 的配体层区间。

### 2.6 基线、对照与敏感性分析

#### 2.6.1 配体化学对照

为考察在不使用受体结构信息的情况下，配体自身化学特征对实验活性类别的区分能力，计算 ECFP4 指纹（radius = 2，2048 bits）以及分子量、重原子数、cLogP 和 TPSA。对四个单一物化描述符分别计算两个方向的 AUROC 及其 \(\mathrm{summary}_{\min}\)，并以每个靶对 \(\mathrm{summary}_{\min}\) 最高的单一描述符作为物化性质基线。由于最佳单一描述符在同一评价集中确定，其与 Vina 的差值仅作为描述性比较（Table 2；Table S5）。

ECFP4、docking-only 和 ECFP4+docking 模型均采用未做特征缩放的逻辑回归。以 Bemis–Murcko 骨架为分组变量，在相同的 GroupKFold 划分下获得 out-of-fold 预测，并据此计算 AUROC。docking-only 模型仅以相应方向的对接评分作为输入。该 AUROC 由 out-of-fold 预测值计算，而主要分析直接按原始对接评分排序，两者计算方式不同，数值并不完全一致。ECFP4 与 ECFP4+docking 的 AUROC 差值用于描述加入对接评分后的增量判别信息。敏感性分析在相同划分的各训练折上拟合 StandardScaler（Table S5）。

#### 2.6.2 结构归因与评分对应性检验

主要口袋对应性分析比较 matched 与 mismatched 条件下的 \(\mathrm{summary}_{\min}\)，并以配对 bootstrap 估计差值区间；单个方向的 AUROC 差异作为补充分析（Figure 5；Table S6）。差值 \(\Delta=\mathrm{summary}_{\min}^{\mathrm{matched}}-\mathrm{summary}_{\min}^{\mathrm{mismatched}}\) 为正值表示对应口袋的较弱方向更高。

受体结构敏感性主要在 PIK3CA/mTOR 中评价。在保持 mTOR 4JT6 不变时，分别以 4JPS 和 5DXT 替换 PIK3CA 4L23；另以 4JSX 替换 mTOR 4JT6。除被替换的受体外，其余配体集合、实验状态、评分定义和统计方法保持不变（Table S8）。

#### 2.6.3 稳健性与外部数据可用性分析

通过改变活性阈值及重复活性记录的汇总方式重新确定实验状态，并在保持评价集成员和对接评分不变的条件下重复主要分析（Table S3）。对于具有足够剩余候选配体的靶对，在排除主评价集成员后构建未使用池留出集。留出配体使用固定抽样规则，并限制单一 Bemis–Murcko 骨架的过度重复。该分析使用不同配体重新计算主要方向性 AUROC，用于评价结果对评价集成员组成的敏感性。由于这些配体仍来自同一数据来源，该分析属于内部稳健性检验（Table S7）。另通过多个固定随机种子重复主要 Vina 对接，并在 PIK3CA/mTOR 中比较不同 exhaustiveness 设置，以评价搜索随机性和搜索强度对结果的影响（Table S9）。备选种子按完整病例计：AChE/BChE 生产种子 n_complete 为 95，其余四种子为 89–90。

另按配体最早来源文献的发表年份进行时间切分，以考察结果对文献时间分布的敏感性（Table S12）。BindingDB 和 PubChem 按与 ChEMBL 相同的四状态规则清点八个靶对的双向选择性供给（Table S11a）。进入外部对接还要求对同一八个靶对做独立来源剩余清点：去掉与主评价共享的文献来源、结构重复和 ECFP4 Tanimoto \(\geq 0.70\) 的分子，并满足 dual、A-only、B-only 各 \(n\geq 20\) 且每类至少 3 个来源（Table S11b）。只有同时满足上述条件的靶对才进入外部评价并进行对接。相关清点、准入标准和时间界点见 Table S11 与 Table S12。

### 2.7 软件与可重复性

分子结构处理、指纹和描述符计算及统计分析均在 Python 环境中完成，主要使用 RDKit、NumPy、pandas、SciPy 和 scikit-learn。受体和配体 PDBQT 文件使用 Meeko 准备。主要分子对接采用 AutoDock Vina 1.2.7；GNINA 1.3.2 用于补充评分和部分靶对的独立姿态生成；RTMScore 用于 Vina 姿态的另一种补充重评分。

软件版本、随机种子和主要计算参数汇总于 Table S1。各受体的对接区域和共晶配体重对接结果见 Table S2。

## 3. 结果

### 3.1 双端实验数据供给与四状态评价集构建

ChEMBL 中能够支持四状态评价的双端活性数据随样本要求提高而迅速减少。在所有至少存在双端测量的靶对中，2,164,618 对至少有 1 个双端测量配体，63,790 对具有不少于 10 个双端测量配体。进一步要求在 \(\theta=6.0\) 下 dual、A-only 和 B-only 均不少于 10 个后，仅剩 5,253 对；采用严格 6.5/5.5 双向选择性供给标准后，进一步减少至 86 对。

对满足数据供给条件的靶对进一步检查后，部分体系因共同高通量测定背景、靶点类型、结合位点性质或缺少合适的人源实验结构而不适合纳入统一非共价对接评价。最终纳入主要评价的靶对为 PIK3CA/mTOR、AChE/BChE、EGFR/HER2、F2/F10、JAK1/TYK2、JAK1/JAK2、PPARG/PPARA 和 PPARA/PPARD（Table 1）。各靶对的候选池规则、配额和骨架上限按样本供给规模与结构可行性分别设定（Table 1）。四状态分类要求同一配体在两个靶点上均具有实验测量，双端实验测量覆盖和双向选择性样本数量共同限制了严格四状态评价集的规模。

### 3.2 实验状态定义与方向性对接评价

固定相同评分通道后，改变实验状态比较所产生的 AUROC 差异在不同靶对和评价方向之间并不一致。其中 EGFR/HER2 和 JAK1/TYK2 出现了较明显的差异。EGFR/HER2 使用靶点 A（EGFR）评分时，dual–B-only 比较的 AUROC 为 0.430，将对照类别替换为双端低活性的 neither 后升至 0.808，差值 0.378 [0.205, 0.547]；另一方向上的差异较小。EGFR/HER2 的差异在骨架簇 [0.168, 0.562] 和文献簇 [0.083, 0.529] 重采样下仍保持；JAK1/TYK2 使用靶点 A（JAK1）评分时，两项比较的差值为 0.444 [0.263, 0.620]，其在骨架簇重采样 [0.234, 0.633] 下保持，而文献簇重采样下区间扩大并包含 0（[−0.034, 0.682]；Tables S4 and S10）。AChE/BChE、JAK1/JAK2、PPARG/PPARA 和 PPARA/PPARD 在相同评分通道下的差值较小，配体水平 95% 置信区间均包含 0（Figure 2C；Table S4）。

在统一 \(\theta=6.0\) 标签下，八个靶对方向性判别的较弱方向描述性汇总指标 \(\mathrm{summary}_{\min}\) 介于 0.345 至 0.692。仅 PPARG/PPARA 的配体水平 bootstrap 95% 置信区间完全高于 0.5（0.649 [0.504, 0.751]）；其余靶对的区间包含 0.5，或整体位于 0.5 以下（Table 2）。

采用双口袋平均评分作描述性比较时，EGFR/HER2 和 JAK1/TYK2 的 dual–neither 比较 AUROC 分别达到 0.756 [0.562, 0.920] 和 0.770 [0.597, 0.906]，而其方向性 \(\mathrm{summary}_{\min}\) 分别为 0.430 和 0.365（Table 3；Figure 2B）。该比较同时改变了评分聚合形式和对照定义，对照组成的单独影响以固定评分通道分析为准。

评价集的实际排序反映了这一设定的筛选后果。在包含 110 个配体的 EGFR/HER2 评价集按双口袋平均评分排序时，Top-10 中仅有 1 个 dual，其余 9 个均为单靶选择性配体，而没有任何 neither 配体（Table S13）。较高的 dual–neither 比较 AUROC 并未对应较少的高排名单靶选择性配体。若以 dual 配体的中位 \(S_{\mathrm{worst}}\) 进行双口袋联合过滤，在保留 14 个 dual 的同时仍会保留 33 个单靶选择性配体，precision 为 0.298（Table S13）。

**Table 2.** 八个主评价靶对上的口袋匹配方向 AUROC（Vina，统一 \(\theta=6.0\)）及四个物化描述符的 \(\mathrm{summary}_{\min}\) 参考。表中类别样本量为 n_scored（dual / A-only / B-only）。最高描述符代表最佳单一描述符基线。

| 靶对 | n_scored (dual / A-only / B-only) | dual 对 A_only（口袋 B） | dual 对 B_only（口袋 A） | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|------:|---:|------:|-----:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] | 0.369 | 0.416 | 0.482 | 0.428 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] | 0.463 | 0.448 | 0.310 | 0.260 |
| F2/F10 | 31 / 32 / 32 | 0.413 | 0.345 | 0.345 [0.211, 0.477] | 0.432 | 0.482 | 0.515 | 0.345 |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 | 0.365 | 0.365 [0.231, 0.503] | 0.369 | 0.389 | 0.580 | 0.425 |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 | 0.728 | 0.588 [0.444, 0.725] | 0.578 | 0.565 | 0.480 | 0.570 |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 | 0.706 | 0.649 [0.504, 0.751] | 0.507 | 0.478 | 0.485 | 0.627 |
| PPARA/PPARD | 32 / 32 / 32 | 0.646 | 0.446 | 0.446 [0.296, 0.584] | 0.490 | 0.436 | 0.564 | 0.351 |

**Table 3.** 同一套 Vina 对接分数在方向性与 Dual-versus-neither 设定下的 AUROC 比较（统一 \(\theta=6.0\)）。Dual-versus-neither 采用双口袋平均评分 \(S_{\mathrm{mean}}\)。PIK3CA/mTOR 的 neither 样本量（n = 4）较少。

| 靶对 | directional summary_min [95% CI] | Dual vs neither (vina_mean) | n_neither | Dual vs all non-duals |
|------|--------------------------------:|------------------------------:|----------:|----------------------:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |
| F2/F10 | 0.345 [0.211, 0.477] | 0.519 [0.350, 0.688] | 12 | 0.405 [0.273, 0.534] |
| JAK1/TYK2 | 0.365 [0.231, 0.503] | 0.770 [0.597, 0.906] | 14 | 0.527 [0.411, 0.638] |
| JAK1/JAK2 | 0.588 [0.444, 0.725] | 0.730 [0.547, 0.875] | 14 | 0.668 [0.561, 0.770] |
| PPARG/PPARA | 0.649 [0.504, 0.751] | 0.685 [0.493, 0.848] | 14 | 0.675 [0.571, 0.780] |
| PPARA/PPARD | 0.446 [0.296, 0.584] | 0.565 [0.368, 0.766] | 14 | 0.522 [0.406, 0.640] |

### 3.3 配体化学基线与方向性判别

部分靶对中，配体自身的物理化学特征即可对 dual 和单靶选择性配体产生一定区分。AChE/BChE 仅使用 TPSA 时，dual–A-only 比较和 dual–B-only 比较的 AUROC 分别为 0.733 和 0.801（Figure 3C；Table 2）。单一物化性质的判别能力在不同靶对之间差异明显，PIK3CA/mTOR 中最佳单一描述符（重原子数）的 \(\mathrm{summary}_{\min}\) 为 0.463（Table 2）。在八个靶对中，Vina 与最佳单一描述符的 \(\mathrm{summary}_{\min}\) 差异整体较小，且多数差值的 95% 置信区间包含 0（Figure S4；Table S5）。

在 Bemis–Murcko 骨架分组交叉验证下，不使用受体结构的 ECFP4 分子指纹在多个方向上获得与对接相当甚至更高的 AUROC（Figure 3A）。将对应方向的对接评分作为输入特征加入 ECFP4 逻辑回归模型后，全部 16 个方向上 AUROC 的最大绝对变化为 0.023（Table S5）。在当前评价集和交叉验证设置下，未观察到稳定的增量判别。

### 3.4 口袋对应性与计算敏感性

为考察方向性判别是否与实验活性差异所在的靶点口袋相对应，Figure 5 与 Table S6 比较对应口袋与非对应口袋的 \(\mathrm{summary}_{\min}\) 差值。主评价集中仅 EGFR/HER2（0.170 [0.060, 0.280]）和 AChE/BChE（0.161 [0.037, 0.269]）两对的差值区间排除 0；其余六对的 95% 置信区间均包含 0。在七个未使用的留出集中，对应口袋与非对应口袋的 \(\mathrm{summary}_{\min}\) 差值区间均包含 0（差值介于 −0.079 至 +0.150；Figure 5B；Table S7）。对应口袋优势未在留出集中稳定重现。

在 PIK3CA/mTOR 中，将 PIK3CA 受体由 4L23 替换为 4JPS 后，\(\mathrm{summary}_{\min}\) 从 0.692 [0.470, 0.813] 降至 0.486 [0.259, 0.692]；替换为 5DXT 后为 0.505 [0.292, 0.696]；将 mTOR 4JT6 替换为 4JSX 后为 0.639 [0.418, 0.776]（Figure 4B；Table S8）。

采用 GNINA 1.3.2 独立生成姿态和评分后，EGFR/HER2 的 dual–neither 比较 AUROC 为 0.783 [0.610, 0.922]，n_neither = 11；方向性较弱臂 dual–B-only 比较为 0.220 [0.109, 0.343]（n_dual = 28，n_B = 32）。JAK1/TYK2 同样呈现该模式，dual–neither 比较为 0.705，方向性 \(\mathrm{summary}_{\min}\) 为 0.317 [0.183, 0.463]（Figure 4A；Table S9）。在独立姿态生成流程下，类似差异仍可观察到。

PPARG/PPARA 在主要 Vina 评价中是唯一 \(\mathrm{summary}_{\min}\) 置信区间完全高于 0.5 的靶对（0.649 [0.504, 0.751]），但在同姿态 RTMScore 重评分下降至 0.369 [0.233, 0.475]，GNINA CNN 重评分降至 0.500，且在未使用池留出集中降至 0.535 [0.350, 0.717]（Table S7；Table S9）。五个固定 Vina 随机种子及 PIK3CA/mTOR 的 exhaustiveness 设置产生的数值波动相对有限（Figure 4C；Figure 6C；Table S9）。EGFR/HER2 的设定差距在五个 Vina 种子上均为正。共晶配体重对接作为协议质控：所有主受体的保存姿态中均存在重原子 RMSD < 2.0 Å 的近天然构象；EGFR 3POZ 的 top-1 RMSD 为 9.505 Å，全部保存姿态中的最低 RMSD 为 0.760 Å（Table S2）。

### 3.5 方向性评价的标签与样本组成敏感性

改变活性阈值（\(\theta=5.5\)、6.0、6.5 及严格 6.5/5.5）后，样本较充足靶对的 \(\mathrm{summary}_{\min}\) 数值变化较小（AChE/BChE 稳定在 0.606）；当单侧选择性样本量减少时估计波动增大（Figure 6A；Table S3）。将重复活性记录由最大值改为中位数后，主要结果总体稳定；在可进行高置信记录复核的评价集中，进一步筛选也未改变相应方向性结果（Table S3）。

在排除主评价集分子后，基于剩余候选分子构建的未使用池留出集显示，结果存在一定的样本组成依赖：AChE/BChE、PIK3CA/mTOR 与 JAK1/JAK2 与主评价接近；JAK1/TYK2 有所上升；F2/F10 与 PPARA/PPARD 仍处于较低水平；PPARG/PPARA 则由 0.649 降至 0.535 [0.350, 0.717]（Figure 5C；Table S7）。改用文献连通簇重采样表明，文献来源集中度增大了统计估计的不确定性（Table S10；Figure S6）。Table S10 另报告一项保持类别样本量的情景模拟，不用于解释 Table 2 区间。

### 3.6 外部评价数据的可用性

BindingDB 与 PubChem 按同一四状态规则清点八个靶对的双向供给，并在排除共享文献、重复结构和高相似分子后进行独立来源评估。没有任何靶对满足独立外部评价准入条件，因此未进行外部对接（Figure 6D；Table S11）。按配体发表年份进行时间切分同样受到双向选择性样本供给的限制，未能构建出可同时独立评价两个方向的时间切分测试集（Table S12）。未使用池留出集与年份切分都只是内部敏感性，不作为外部验证。当前公开数据尚不足以为这些靶对构建独立的双向四状态外部评价集。

## 4. 讨论

### 4.1 实验状态定义与评价任务

双靶对接评价所回答的问题取决于对照配体及其实验活性状态。dual–neither 比较检验双靶活性配体能否与两个靶点均低活性的化合物区分；dual–A-only 比较和 dual–B-only 比较则分别检验对接评分能否区分 dual 与仍保留一个靶点活性的选择性配体。这两类比较回答的是不同问题，而不能理解为同一评价任务中负样本难度不同。

固定评分通道排除了评分汇总方式改变这个因素，但在从 B-only 换成 neither 的同时更换了具体分子，其骨架、大小、电荷和实验来源也可能发生变化。因此，该分析反映的是评价类别及其对应化合物集合整体变化的影响，而不能将观察到的差异单独归因于实验活性状态。EGFR/HER2 与 JAK1/TYK2 在保持同一靶点评分不变、仅改变对照类别及其对应化合物集合后仍表现出明显差异，说明这一现象不能仅由双口袋评分的聚合方式解释。其他靶对和评价方向的变化幅度并不一致，说明选择性配体不能简单视为一类统一的更难区分的负类，对照类别选择的影响还取决于具体靶对和评价方向。在 Zhou 等工作的基础上，[9] 本研究进一步将两类单靶选择性状态分别构建为对应的方向性判别任务，并在固定评分通道下考察评价设定变化，再结合配体基线和口袋对应性检验判别信息可以被解释到什么程度。

经典虚拟筛选基准采用了不同类型的负类，例如 DUD 和 DUD-E 的人工诱饵，以及 LIT-PCBA 中由高通量实验确定的非活性化合物。[5–7] 近期比较研究进一步表明，改变对照化合物的来源和定义可以改变虚拟筛选性能，甚至影响不同方法之间的相对排序。[12,13,21] 双靶评价还需要考虑两个选择性方向。dual–neither 比较中较高的 AUROC 只能说明对双端低活性配体的区分较好，不能自动作为双靶选择性识别的证据。对于以双靶选择性为目标的回顾性评价，两个方向的选择性对照仍需要分别考察。Kinase-Bench 检验的是相对激酶 decoy 的选择性富集，[22] 与在已测两端活性的化合物中区分 dual 与选择性配体不是同一问题。

### 4.2 配体化学与表观判别来源

回顾性评价中的活性类别本身可能带有化学结构差异。AChE/BChE 中简单物化描述符已能区分部分实验状态，ECFP4 在骨架分组交叉验证下也在多个方向获得较高 AUROC。进一步加入相应方向的对接评分后，16 个方向的 AUROC 最大绝对变化仅为 0.023。这一结果并不意味着对接评分完全不含受体相关信息，但说明较高的回顾性 AUROC 不能仅凭自身数值归因于三维结构互补。判断对接是否增加了判别能力，需要与不使用受体结构的基线放在相同划分下比较。

这种类别间化学差异可能与公开生物活性数据的来源有关。数据库中的化合物往往来自围绕有限先导骨架开展的药物化学研究，这可能使不同实验活性类别同时具有不同的骨架、取代模式和物化性质分布。Caba 等在结构基础虚拟筛选中同时考察了蛋白–配体表征和仅配体 Morgan 指纹：加入 Morgan 指纹改善了多个模型的表现，但没有进一步改善表现最好的 PLEC 模型。[23] 该结果同样提示，在回顾性结构基础虚拟筛选中，仅配体信息可能携带较强的分类信号。简单物化性质匹配可以减少部分低维性质差异，但不能完全控制骨架和药物化学系列之间的高维差异，并且严格匹配会进一步减少可用样本。

### 4.3 结构归因与口袋对应性

已有激酶交叉对接研究表明，受体构象选择可以显著影响姿态恢复和虚拟筛选性能。[14] 本研究中，PIK3CA/mTOR 在替换 PIK3CA 或 mTOR 受体结构后，方向性 AUROC 也发生明显变化。由于这一分析仅系统覆盖一个主评价靶对，现有结果只能说明该体系存在受体结构敏感性，而不能据此推断所有双靶体系都会表现出相同程度的受体依赖。

共晶配体重对接还显示，能够搜索到近天然姿态并不意味着该姿态一定获得最高评分。因此，最低 RMSD 主要反映搜索是否能够产生近天然姿态，并不能评价这些姿态能否被评分函数正确排序，也不能单独说明其虚拟筛选判别能力。更换姿态生成和评分实现后，具体 AUROC 会发生变化，但 EGFR/HER2 和 JAK1/TYK2 中 dual 与单靶选择性配体比较、dual 与 neither 比较之间的明显差异在独立 GNINA 流程下仍可观察到。这说明评价设定带来的差异并非完全依赖于主要 Vina 实现。与此同时，RTMScore、GNINA CNN 重评分和受体替换又显示，具体方向性 AUROC 对评分和结构实现并不稳定。独立 GNINA 检验的是评价设定差异能否在另一套姿态生成流程下再现，而不是对应口袋优势。

评分通道交换提供了更直接的结构对应检验。Figure 5 与 Table S6 报告的是对应与非对应口袋的 \(\mathrm{summary}_{\min}\) 差值。主评价集中，EGFR/HER2 和 AChE/BChE 显示对应口袋评分优势，但其余六个靶对的差值区间均包含 0；在七个可评价的未使用配体留出集中，该差值区间均包含 0。留出集结果未能提供与主面板一致的对应口袋优势；这应理解为结构对应证据不足，而不是证明不存在结构信息。在同源靶点上，非对应口袋也可能产生相关且具有生物学意义的评分。因此，仅凭某一评价集中较高的方向性 AUROC，仍不足以把该判别归因于相应靶点口袋。

### 4.4 双靶对接的实际作用与证据边界

在本文考察的回顾性评价条件下，对接评分应视为候选排序或筛选信号，而不能单独作为真实双靶活性或选择性的证据。两个口袋中均获得有利对接评分，只表示该分子在当前计算条件下对两个靶点均获得了有利预测评分，不能据此确认真实双靶活性或选择性。

已有前瞻性研究表明，多靶对接和生成式设计可以产生经实验验证的多靶候选物。[19,20] 这类研究评价的是从候选库中找到活性分子的能力；本文评价的是在已有双端实验标签的配体中，不同实验状态之间的回顾性判别能力。因此，回顾性 AUROC 既不能替代前瞻实验命中率，也不应被用来否定对接在候选库压缩中的用途。同样，较高的回顾性 AUROC 也不能单独证明方法能稳定找到具有预期活性平衡的双靶分子。

对于回顾性双靶对接评价，至少应分别报告 dual–A-only 比较和 dual–B-only 比较两个方向，并结合仅配体基线、对应/非对应口袋比较以及必要的受体或评分敏感性分析，以判断表观 AUROC 是否具有稳定的结构解释。最终的双靶活性、活性平衡和作用机制仍需要直接实验测定。

### 4.5 适用范围与数据限制

本研究的结论范围主要受到数据层、评价面板构建层以及计算实现层三个维度的限制。

在数据层面上，双端实验活性测量的相对稀缺和异质性构成了根本制约。四状态分类要求同一配体在两个靶点上均有实验测量，而公开数据库中的可用记录在活性指标、实验条件和蛋白构建形式上并不完全一致，且往往高度集中于少数研究较充分的药物化学系列。虽然活性阈值、重复记录汇总方式及高置信筛选未改变主要观察，但仍无法完全消除不同实验来源和化学系列固有的数据异质性。在排除来源重叠、重复结构和高相似分子后，BindingDB 未形成满足本研究准入条件的独立外部评价集，时间切分也未能提供可同时评价两个选择性方向的测试集。

在评价面板构建层面上，由于不同靶对的双向选择性数据供给差异明显，评价面板在候选池和抽样约束上存在实际差异。虽然主要统计分析统一采用 \(\theta=6.0\) 标签并进行了留出集和敏感性分析，但构建方式差异对跨靶对可比性的潜在影响仍不可完全排除。此外，纳入的八个靶对用于考察不同体系中的一致性，不应视为相互独立的统计重复：JAK1 与 PPARA 各被使用两次，部分靶对属于相近激酶或核受体家族。除 PPARG/PPARA 外，主要 Vina \(\mathrm{summary}_{\min}\) 的 bootstrap 95% 置信区间均未完全高于 0.5，反映出较宽的统计不确定性。

在计算实现层面上，本研究主要依赖于每个靶点单一的代表性晶体结构和 AutoDock Vina 默认参数。尽管共晶重对接显示能够覆盖近天然构象，但共晶重对接主要反映近天然姿态的搜索覆盖，不能替代对评分排序能力和筛选性能的评价。独立 GNINA 姿态生成、替代受体构象和其他评分函数的分析仅覆盖了部分体系。同时，配体准备未系统枚举质子化状态、互变异构体和构象系综，这些近似对预测精度的影响仍待进一步系统评估。因此，不宜将本文观察到的具体数值外推为双靶对接在其他蛋白家族和更广阔化学空间中的普适表现。

提高此类评价的证据强度，需要获得同一批化合物在两个靶点上的成对活性数据，并尽量统一实验条件。这样可以减少不同文献和测定条件带来的异质性；如果再按化学空间取样，还能降低结果对少数药物化学系列的依赖。这类数据将更适合用于检验四状态评价能否推广至更多蛋白家族，并支持时间切分、跨数据库和前瞻性实验验证。

## 5. 结论

本研究在八个具有双端实验活性的人源靶对上系统考察了双靶分子对接评价设定及其判别来源，主要得出三项结论：

第一，固定同一靶点评分通道时，更换对照类别及其对应的化合物集合可以改变对双靶对接性能的判断，而且这种影响具有靶对和方向依赖性。能够区分双靶活性与双端低活性配体，并不保证能够区分双靶活性与单靶选择性配体，二者反映了不同的判别问题。

第二，部分回顾性判别可以由不使用受体结构的配体化学信息获得，加入对接评分后的判别增量较为有限；同时，对应口袋评分优势在未使用的留出集中未能稳定重现，表明回顾性判别不能简单直接归因于口袋特异性三维结构互补。

第三，建议回顾性双靶对接评价分别明确报告两个方向的选择性判别任务，并结合配体化学基线与口袋对应性等对照分析以辅助解释。此类回顾性评价指标反映的是在已知实验标签集合中的区分能力，不能直接等同于前瞻性虚拟筛选中的真实命中能力。

## 数据与软件可用性

评价面板成员、实验状态标签、受体与对接盒定义、逐配体对接分数、分析表，以及重建本文统计与图件所需的全部脚本，均可在公开仓库 https://github.com/1280602962-debug/gwj260531 的 `Dual_Target_Docking` 目录中获取。

排入 Supporting Information 的核实表为 **Tables S1–S13**（软件与种子、对接盒子与共晶 RMSD、阈值与 pChEMBL 聚合、固定评分负类、ECFP4 增量、对应/非对应口袋、未使用池留出集、PIK3CA/mTOR 晶体替换、独立 GNINA 与五种子 Vina、文献簇与样本量情景、BindingDB 外部准入、2018 年切分计数、EGFR/HER2 操作点）。历史工作表 S1–S54（含原 Table S54 多种子表）的合并对照见仓库内 `Dual_Target_Docking/data/manuscript_lock/SI_TABLE_MERGE_MAP_v1.csv`。

更长的逐配体分数、holdout 成员、多种子长表、物化 caliper、骨架相近的选择性对照、聚合均值、完整病例与 assay-context 底表、J0 候选对普查、BindingDB REST 历史计数、leave-cognate-out、PIK3CA 占有率快照、接触计数、全链序列一致性，以及 MCL1/Bcl-xL 适用性压力测试（不进入 Table 2），见 **Note S14**，随代码与 SHA-256 清单归档至 GitHub Release；Zenodo DOI 将从打标签快照签发，而不是从当前仍可能变化的分支签发。

`data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` 索引主要数值结果及其来源 CSV。面向稿件的表 SHA-256 校验和见 `REVISION_CHECKSUM_MANIFEST_v1.csv`。原生切片合约为 `protocol/external_slice_contract.yaml`；评价合约为 `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`。ChEMBL 供给审计冻结于 2026-07-23；高置信 activity 视图抓取于 2026-08-26；BindingDB 原生归档锁定为 202609。分析环境与零新对接的复现命令见仓库 README。BindingDB TSV 归档本身不随仓库分发；CI 只核对已提交的 CSV。

## 参考文献

(1) Anighoro, A.; Bajorath, J.; Rastelli, G. Polypharmacology: Challenges and Opportunities in Drug Discovery. *J. Med. Chem.* **2014**, *57*, 7874–7887. DOI: 10.1021/jm5006463.

(2) Proschak, E.; Stark, H.; Merk, D. Polypharmacology by Design: A Medicinal Chemist’s Perspective on Multitargeting Compounds. *J. Med. Chem.* **2019**, *62*, 420–444. DOI: 10.1021/acs.jmedchem.8b00760.

(3) Kitchen, D. B.; Decornez, H.; Furr, J. R.; Bajorath, J. Docking and Scoring in Virtual Screening for Drug Discovery: Methods and Applications. *Nat. Rev. Drug Discov.* **2004**, *3*, 935–949. DOI: 10.1038/nrd1549.

(4) Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. Model.* **2021**, *61*, 3891–3898. DOI: 10.1021/acs.jcim.1c00203.

(5) Huang, N.; Shoichet, B. K.; Irwin, J. J. Benchmarking Sets for Molecular Docking. *J. Med. Chem.* **2006**, *49*, 6789–6801. DOI: 10.1021/jm0608356.

(6) Mysinger, M. M.; Carchia, M.; Irwin, J. J.; Shoichet, B. K. Directory of Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking. *J. Med. Chem.* **2012**, *55*, 6582–6594. DOI: 10.1021/jm300687e.

(7) Tran-Nguyen, V.-K.; Jacquemard, C.; Rognan, D. LIT-PCBA: An Unbiased Data Set for Machine Learning and Virtual Screening. *J. Chem. Inf. Model.* **2020**, *60*, 4263–4273. DOI: 10.1021/acs.jcim.0c00155.

(8) Su, M.; Yang, Q.; Du, Y.; Feng, G.; Liu, Z.; Li, Y.; Wang, R. Comparative Assessment of Scoring Functions: The CASF-2016 Update. *J. Chem. Inf. Model.* **2019**, *59*, 895–913. DOI: 10.1021/acs.jcim.8b00545.

(9) Zhou, S.; Li, Y.; Hou, T. Feasibility of Using Molecular Docking-Based Virtual Screening for Searching Dual Target Kinase Inhibitors. *J. Chem. Inf. Model.* **2013**, *53*, 982–996. DOI: 10.1021/ci400065e.

(10) Zhou, X.; Guan, J.; Zhang, Y.; Peng, X.; Wang, L.; Ma, J. Reprogramming Pretrained Target-Specific Diffusion Models for Dual-Target Drug Design. In *The Thirty-eighth Annual Conference on Neural Information Processing Systems (NeurIPS 2024)*; 2024. arXiv:2410.20688.

(11) Wu, J.; Qiao, A.; Wang, Z.; Wei, Z.; Chen, S. FuseDiff: Symmetry-Preserving Joint Diffusion for Dual-Target Structure-Based Drug Design. In *Proceedings of the 32nd ACM SIGKDD Conference on Knowledge Discovery and Data Mining, Vol. 2*; ACM: New York, 2026; pp 12432–12443. DOI: 10.1145/3770855.3819050.

(12) Tran-Nguyen, V.-K.; Ballester, P. J. Beware of Simple Methods for Structure-Based Virtual Screening: The Critical Importance of Broader Comparisons. *J. Chem. Inf. Model.* **2023**, *63*, 1401–1405. DOI: 10.1021/acs.jcim.3c00218.

(13) Ahmed, F.; Soellner, M. B.; Brooks, C. L., III. Real-World Assessment of Machine-Learned Docking Using Bioassay-Derived Benchmarks. *J. Chem. Inf. Model.* **2026**, *66*, 8752–8759. DOI: 10.1021/acs.jcim.5c03020.

(14) Schaller, D. A.; Christ, C. D.; Chodera, J. D.; Volkamer, A. Benchmarking Cross-Docking Strategies in Kinase Drug Discovery. *J. Chem. Inf. Model.* **2024**, *64*, 8848–8858. DOI: 10.1021/acs.jcim.4c00905.

(15) Sindt, F.; Bret, G.; Rognan, D. On the Difficulty to Rescore Hits from Ultralarge Docking Screens. *J. Chem. Inf. Model.* **2025**, *65*, 5553–5566. DOI: 10.1021/acs.jcim.5c00730.

(16) Liu, T.; Hwang, L.; Burley, S. K.; Nitsche, C. I.; Southan, C.; Walters, W. P.; Gilson, M. K. BindingDB in 2024: a FAIR Knowledgebase of Protein-Small Molecule Binding Data. *Nucleic Acids Res.* **2025**, *53*, D1633–D1644. DOI: 10.1093/nar/gkae1075.

(17) Tanaka, Y.; Aikawa, K.; Nishida, G.; Homma, M.; Sogabe, S.; Igaki, S.; Hayano, Y.; Sameshima, T.; Miyahisa, I.; Kawamoto, T.; Tawada, M.; Imai, Y.; Inazuka, M.; Cho, N.; Imaeda, Y.; Ishikawa, T. Discovery of Potent Mcl-1/Bcl-xL Dual Inhibitors by Using a Hybridization Strategy Based on Structural Analysis of Target Proteins. *J. Med. Chem.* **2013**, *56*, 9635–9645. DOI: 10.1021/jm401170c.

(18) García-Ortegón, M.; Simm, G. N. C.; Tripp, A. J.; Hernández-Lobato, J. M.; Bender, A.; Bacallado, S. DOCKSTRING: Easy Molecular Docking Yields Better Benchmarks for Ligand Design. *J. Chem. Inf. Model.* **2022**, *62*, 3486–3502. DOI: 10.1021/acs.jcim.1c01334.

(19) Wu, Y.; Vigneron, S.; Braz, J.; Srinivasan, K.; Fink, E. A.; Huang, X.-P.; Xu, X.; Huebner, H.; Kim, J. Y.; Wang, J.; Pfeiffer, T.; Sakamoto, K.; Moroz, Y. S.; Radchenko, D. S.; Rodriguiz, R. M.; Irwin, J. J.; Gmeiner, P.; Billesboelle, C.; Roth, B. L.; Basbaum, A. I.; Manglik, A.; Wetsel, W. C.; Shoichet, B. K. Large Library Docking for Polypharmacology. *J. Med. Chem.* **2026**, *69*, 6210–6229. DOI: 10.1021/acs.jmedchem.5c03810.

(20) Munson, B. P.; Chen, M.; Bogosian, A.; Kreisberg, J. F.; Licon, K.; Kuenzi, B. M.; Ideker, T. De novo generation of multi-target compounds using deep generative chemistry. *Nat. Commun.* **2024**, *15*, 3636. DOI: 10.1038/s41467-024-47120-y.

(21) Gu, S.; Shen, C.; Zhang, X.; Sun, H.; Cai, H.; Luo, H.; Zhao, H.; Liu, B.; Du, H.; Zhao, Y.; Fu, C.; Zhai, S.; Deng, Y.; Liu, H.; Hou, T.; Kang, Y. Benchmarking AI-powered docking methods from the perspective of virtual screening. *Nat. Mach. Intell.* **2025**, *7*, 509–520. DOI: 10.1038/s42256-025-00993-0.

(22) Wei, T.-H.; Zhou, S.-S.; Jing, X.-L.; Liu, J.-C.; Sun, M.; Zhao, Z.-H.; Li, Q.-Q.; Wang, Z.-X.; Yang, J.; Zhou, Y.; Wang, X.; Ling, C.-X.; Ding, N.; Xue, X.; Yu, Y.-C.; Wang, X.-L.; Yin, X.-Y.; Sun, S.-L.; Cao, P.; Li, N.-G.; Shi, Z.-H. Kinase-Bench: Comprehensive Benchmarking Tools and Guidance for Achieving Selectivity in Kinase Drug Discovery. *J. Chem. Inf. Model.* **2024**, *64*, 9528–9550. DOI: 10.1021/acs.jcim.4c01830.

(23) Caba, K.; Tran-Nguyen, V.-K.; Rahman, T.; Ballester, P. J. Comprehensive machine learning boosts structure-based virtual screening for PARP1 inhibitors. *J. Cheminform.* **2024**, *16*, 40. DOI: 10.1186/s13321-024-00832-1.
