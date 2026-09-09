# Methods（中文工作稿 · JCIM Articles）

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

此外，在 EGFR/HER2、PIK3CA/mTOR 和 JAK1/TYK2 上使用 GNINA 1.3.2 独立生成姿态。该分析使用与主要分析相同的受体、配体和对接区域，三个靶对的 exhaustiveness 分别为 8、16 和 8，每个配体最多输出 9 个姿态。配体级评分取排名第 1 姿态的 minimizedAffinity，并使用其相反数进行统计分析，使数值越高表示预测结合越有利。该分析用于考察主要观察在另一套姿态生成和评分流程下是否仍然存在，而不用于比较 GNINA 与 Vina 的总体性能（Table S9）。

### 2.5 评价指标与统计分析

#### 2.5.1 主要方向性判别指标

比较 dual 与 A-only 时，两类配体在靶点 A 上均达到实验活性阈值，因此使用靶点 B 的评分：\(\mathrm{AUC}_{D/A}(B)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{A\text{-}only};\;S_{B})\)。比较 dual 与 B-only 时，使用靶点 A 的评分：\(\mathrm{AUC}_{D/B}(A)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{B\text{-}only};\;S_{A})\)。所有方向性分析均以 dual 为正类。

AutoDock Vina 输出的 affinity 越低表示预测结合越有利。为统一 AUROC 的评分方向，将 Vina affinity 转换为 \(S_{\mathrm{Vina}}=-E_{\mathrm{Vina}}\)，因此较高的 \(S_{\mathrm{Vina}}\) 表示更有利的预测结合。两个方向的 AUROC 分别报告，并将两者中的较低值定义为 \(\mathrm{summary}_{\min}\)，仅用于描述较弱方向的总体表现，不作为新的配体评分函数：\(\mathrm{summary}_{\min}=\min(\mathrm{AUC}_{D/A}(B),\;\mathrm{AUC}_{D/B}(A))\)。

#### 2.5.2 不同对照类别与双口袋过滤的描述性比较

另以两个靶点均低活性的 neither 配体作为对照类别，评价对接评分区分 dual 与 neither 的能力。对于两个靶点均获得评分的配体，计算平均评分 \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\)，随后以 dual 为正类、neither 为对照计算 AUROC。将 A-only、B-only 和 neither 合并为对照得到的 dual 与所有非 dual 配体比较 AUROC 仅作为混合候选库的描述性参考，并在 Table 3 中报告。

为单独考察对照类别变化对判别结果的影响，保持所用靶点评分不变，仅改变对照组成。使用靶点 A 的评分时，分别计算 dual–B-only 比较和 dual–neither 比较；使用靶点 B 的评分时，分别计算 dual–A-only 比较和 dual–neither 比较（Table S4）。另取两个靶点评分中的较低值作为双口袋联合过滤评分 \(S_{\mathrm{worst}}=\min(S_{A},S_{B})\)，用于评价配体是否在两个靶点上均获得有利评分。以 dual 配体 \(S_{\mathrm{worst}}\) 分布的中位数作为双口袋过滤阈值，并统计不低于该阈值的 dual 配体数、dual recall、dual precision 以及保留的单靶选择性配体数（Table S13）。由于双口袋平均评分同时改变了对照组成和评分形式，对照类别变化的单独影响以固定同一口袋评分的比较为准。

#### 2.5.3 置信区间与重采样分析

Table 2 中 \(\mathrm{summary}_{\min}\) 的 95% 置信区间基于 2000 次配体水平的非分层百分位 bootstrap 估计。每次从参与该靶对方向性分析的 dual、A-only 和 B-only 联合配体库中有放回抽取与原集合相同数量的配体，重新计算 dual–A-only 比较和 dual–B-only 比较的 AUROC，并取两者较小值作为该次重采样的 \(\mathrm{summary}_{\min}\)。95% 置信区间由所有有效 \(\mathrm{summary}_{\min}\) 的第 2.5 和第 97.5 百分位数确定。表中点估计均由完整分析样本直接计算，而非 bootstrap 均值。dual–neither 比较的 AUROC 置信区间采用按类别重采样的百分位 bootstrap。

比较不同评分方法或对应与非对应口袋评分时，各方案在每次 bootstrap 中使用相同的配体重采样结果以保持配对。为评估化学骨架相关性和文献来源相关性的影响，另以 Bemis–Murcko 骨架簇和文献连通簇为单位进行 cluster bootstrap，具体算法细节及样本量可检测效应模拟见 Supporting Information（Table S4；Table S10；Figure S6）。

### 2.6 基线、对照与敏感性分析

#### 2.6.1 配体化学对照

为考察在不使用受体结构信息的情况下，配体自身化学特征对实验活性类别的区分能力，计算 ECFP4 指纹（radius = 2，2048 bits）以及分子量、重原子数、cLogP 和 TPSA。对四个单一物化描述符分别计算两个方向的 AUROC 及其 \(\mathrm{summary}_{\min}\)，并以每个靶对 \(\mathrm{summary}_{\min}\) 最高的单一描述符作为物化性质基线。由于最佳单一描述符在同一评价集中确定，其与 Vina 的差值仅作为描述性比较（Table 2；Table S5）。

ECFP4、docking-only 和 ECFP4+docking 模型均采用逻辑回归。以 Bemis–Murcko 骨架为分组变量，在相同的 GroupKFold 划分下获得 out-of-fold 预测，并据此计算 AUROC。docking-only 模型仅以相应方向的对接评分作为输入。该 AUROC 由 out-of-fold 预测值计算，而主要分析直接按原始对接评分排序，两者计算方式不同，数值并不完全一致。ECFP4 与 ECFP4+docking 的 AUROC 差值用于描述加入对接评分后的增量判别信息（Table S5）。

#### 2.6.2 结构归因与评分对应性检验

主要口袋对应性分析比较 matched 与 mismatched 条件下的 \(\mathrm{summary}_{\min}\)，并以配对 bootstrap 估计差值区间；单个方向的 AUROC 差异作为补充分析（Figure 5；Table S6）。差值 \(\Delta=\mathrm{summary}_{\min}^{\mathrm{matched}}-\mathrm{summary}_{\min}^{\mathrm{mismatched}}\) 为正值表示对应口袋的较弱方向更高。

受体结构敏感性主要在 PIK3CA/mTOR 中评价。在保持 mTOR 4JT6 不变时，分别以 4JPS 和 5DXT 替换 PIK3CA 4L23；另以 4JSX 替换 mTOR 4JT6。除被替换的受体外，其余配体集合、实验状态、评分定义和统计方法保持不变（Table S8）。

#### 2.6.3 稳健性与外部数据可用性分析

通过改变活性阈值及重复活性记录的汇总方式重新确定实验状态，并在保持评价集成员和对接评分不变的条件下重复主要分析（Table S3）。对于具有足够剩余候选配体的靶对，在排除主评价集成员后构建未使用池留出集。留出配体使用固定抽样规则，并限制单一 Bemis–Murcko 骨架的过度重复。该分析使用不同配体重新计算主要方向性 AUROC，用于评价结果对评价集成员组成的敏感性。由于这些配体仍来自同一数据来源，该分析属于内部稳健性检验（Table S7）。另通过多个固定随机种子重复主要 Vina 对接，并在 PIK3CA/mTOR 中比较不同 exhaustiveness 设置，以评价搜索随机性和搜索强度对结果的影响（Table S9）。

另按配体最早来源文献的发表年份进行时间切分，以考察结果对文献时间分布的敏感性（Table S12）。BindingDB 和 PubChem 按与 ChEMBL 相同的四状态规则清点八个靶对的双向选择性供给（Table S11a）。进入外部对接还要求对同一八个靶对做独立来源剩余清点：去掉与主评价共享的文献来源、结构重复和 ECFP4 Tanimoto \(\geq 0.70\) 的分子，并满足 dual、A-only、B-only 各 \(n\geq 20\) 且每类至少 3 个来源（Table S11b）。只有同时满足上述条件的靶对才进入外部评价并进行对接。相关清点、准入标准和时间界点见 Table S11 与 Table S12。

### 2.7 软件与可重复性

分子结构处理、指纹和描述符计算及统计分析均在 Python 环境中完成，主要使用 RDKit、NumPy、pandas、SciPy 和 scikit-learn。受体和配体 PDBQT 文件使用 Meeko 准备。主要分子对接采用 AutoDock Vina 1.2.7；GNINA 1.3.2 用于补充评分和部分靶对的独立姿态生成；RTMScore 用于 Vina 姿态的另一种补充重评分。

软件版本、随机种子和主要计算参数汇总于 Table S1。各受体的对接区域和共晶配体重对接结果见 Table S2。
