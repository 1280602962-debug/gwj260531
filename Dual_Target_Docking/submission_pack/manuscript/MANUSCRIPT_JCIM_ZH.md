# 双靶分子对接评价：多靶对实验状态比较与判别来源分析

## 摘要

分子对接常用于双靶点虚拟筛选中的候选排序。回顾性评价中的对照既可以是两端均低于活性阈值的化合物，也可以是仍对其中一个靶点保持活性的单靶选择性配体。我们在八个人源靶对上，按 dual、A-only、B-only 和 neither 四种实验状态建立评价集，并在固定同一靶点评分通道后比较不同对照。更换对照类别及其对应化合物集合改变了表观对接判别，且该变化因靶对和方向而异：EGFR/HER2 中 EGFR 口袋评分对 dual–B-only 的 AUROC 为 0.334，对 neither 升至 0.780（差值 0.446 [0.259, 0.632]）；JAK1/TYK2 出现类似差值（0.444 [0.261, 0.630]）。其余靶对的效应更小或更不确定。不依赖受体的 ECFP4 在多个方向上携带了明显的类别信息。在相同的骨架分组逻辑回归设定下，加入对应口袋评分后，AUROC 最多变化 0.023。八个主评价集中 AChE/BChE 与 EGFR/HER2 的 matched−mismatched 95% 区间排除 0，七个可用留出集均未排除 0；这些区间未进行多重比较校正。这些结果支持更窄的解释：在两端均有实验测量的化合物中，dual 与双端低活性配体之间较好的判别，并不等于对单靶选择性配体的良好区分。回顾性双靶评价应分别报告两个方向的选择性对照，并结合配体化学与口袋对应性判断表观判别的来源。

**关键词：** 双靶对接；实验状态；选择性；AutoDock Vina；虚拟筛选评价

## 1. 引言

多靶点药物设计旨在利用单一分子同时调控两个或多个疾病相关靶点。与多药联合相比，单一多靶配体有望减少不同药物之间的药代动力学不匹配及制剂复杂性。[1,2] 生成式设计与超大规模化合物库对接已经能够前瞻性得到经实验验证的双靶候选物：POLYGON 对 32 个 MEK1/mTOR 候选物进行了实验测试；[19] 超大规模库对接则从数亿级化学空间中发现具有多靶活性的配体，其中部分预测结合模式获得了实验结构支持。[18] 这些工作展示的是计算方法从候选库中找到双靶活性分子的能力。结构导向的双靶生成模型亦可同时针对两个口袋进行设计。[10,11] 当基准配体在两个靶点上均有实验测量时，对接评分的解释如何改变，仍需单独评价。

分子对接通过搜索配体在结合位点中的可能构象并对其评分，是虚拟筛选中常用的排序方法。[3,4] 许多基准把已知活性分子与人工诱饵构成二分类任务，例如 DUD 和 DUD-E；LIT-PCBA 则使用高通量实验确定的活性与非活性化合物；DOCKSTRING 等基准则提供可重复的对接流程以便比较配体设计方法。[5–7,17] 同一计算流程在不同基准构成下可能得到明显不同的性能评价。[12,13,20] 超大规模对接命中的再评分也并不总是能够稳定保留原先排序。[15] 由于不同基准往往同时改变靶点、配体组成和实验背景，这类比较并不能单独确定造成性能差异的具体因素。在固定靶点和评分通道的条件下比较不同实验状态，可以减少靶点和评分方法差异对结果解释的影响。

当同一配体在两个靶点上都有可用实验测定时，可以按活性阈值把它分成四类：两端均达到阈值的双活性配体，仅在靶点 A 或仅在靶点 B 达到阈值的单靶选择性配体，以及两端均低于阈值的双端低活性配体。下文分别记为 dual、A-only、B-only 和 neither。dual 与 A-only 或 B-only 只在其中一个靶点上不同，与 neither 则在两个靶点上均不同，因此 dual–A-only、dual–B-only 与 dual–neither 回答的不是同一个问题。双靶结构筛选通常将候选分子分别对接至两个口袋，并优先选择两端评分均有利的分子。Zhou 等人在多个激酶靶对上评价了对接双靶虚拟筛选，并观察到单靶抑制剂是预测双靶抑制剂时的重要假阳性来源；[9] 该工作分析了不同类型假阳性，但没有把两个选择性状态明确表述为分别与对应靶点对齐的方向性任务。Kinase-Bench 用实验选择性配体检验相近激酶之间的结构区分，目标是识别对某一激酶具有选择性的抑制剂，这与在两端均有实验测量的化合物中区分 dual 与 A-only、B-only 不是同一问题。[21]

已有工作表明单靶选择性配体会影响双靶筛选。当两端实验测量已经可得时，评价任务本身如何改变对接性能的解释，仍需要专门考察。本文据此提出三个问题。第一，固定评分通道后，更换所比较配体的实验状态是否改变对双靶对接性能的判断。第二，不使用受体结构的配体化学特征能够获得怎样的判别。第三，观察到的对接判别是否与对应口袋一致，并在已有对照中是否稳定。

## 2. 方法

### 2.1 研究设计

双靶虚拟筛选通常将每个候选分子分别对接至两个结合位点，并优先排序两端评分均有利的分子。对该排序作回顾性评价，需要两端都有实验活性，因为被对接的分子可能是双靶活性、仅对其中一个靶点有活性，或两端均低于阈值。配体按两个靶点上的实验活性划分为 dual、A-only、B-only 和 neither 四类。

方向性对接评价分别比较 dual 与两类单靶选择性配体，并以 dual 为正类。这两项比较对应双靶筛选中的难负样本问题：单靶活性配体仍可能在一个口袋获得有利评分，并被排成双靶候选。dual 与 A-only 在靶点 A 上均达到活性阈值，实验状态差异位于靶点 B，因此使用口袋 B 的对接评分；dual 与 B-only 的实验状态差异位于靶点 A，因此使用口袋 A 的对接评分。四种实验状态及两个方向性对接任务见 Figure 1A,B。

固定同一靶点评分后的两个方向性 AUROC 是主要对接终点。两者较低值 \(\mathrm{summary}_{\min}=\min\{\mathrm{AUROC}_{D/A}(B),\mathrm{AUROC}_{D/B}(A)\}\) 只是较弱方向的描述性汇总，不是单独的总体性能终点。固定口袋的 dual–neither 比较考察对照中不再含单靶活性配体时 AUROC 如何变化，相当于传统虚筛中的双端低活性对照（Table S4）。双口袋平均分 \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\) 仅作为传统双靶排序读数，用于 dual–neither 比较（Table 3）。配体化学基线检验不依赖受体的性质能否分开同一套筛选类别；对应与非对应口袋评分检验判别是否跟随所指定的对接位点。受体替换、独立姿态生成、标签敏感性和未使用池留出用于考察对接结果的敏感性。检索 BindingDB 与 PubChem，以判断能否形成独立的外部对接集。

### 2.2 生物活性数据处理与活性状态定义

实验活性数据来自 ChEMBL 37。Table 2 主标签取各配体–靶点组合的最大 pChEMBL。具有 IC50、Ki、Kd、EC50、Potency、IC50app 或 Ki app 且有 pChEMBL 的定量记录予以保留，并统一到负对数摩尔尺度。同一配体–靶点组合存在多条有效记录时，主要分析取最大 pChEMBL 作为代表值。该选择偏向最强记录，可受单次高值、测定条件或终点类型影响。盐、混合物等含多个组分的记录先按连通片段拆分，再保留重原子数最多的有机片段。仅将在两个靶点上均有可用实验测量的配体纳入四状态分类；缺少任一端实验记录的配体不视为该靶点低活性。靶对宇宙普查使用同一套 ChEMBL 37 记录。

对任一靶对 A/B，以活性阈值 \(\theta\) 定义四种活性类别：dual，\(p_{\mathrm{A}}\geq\theta\) 且 \(p_{\mathrm{B}}\geq\theta\)；A-only，\(p_{\mathrm{A}}\geq\theta\) 且 \(p_{\mathrm{B}}<\theta\)；B-only，\(p_{\mathrm{A}}<\theta\) 且 \(p_{\mathrm{B}}\geq\theta\)；neither，两端均 \(<\theta\)。主要分析统一采用 \(\theta=6.0\)。A-only 和 B-only 仅表示相对于当前靶对及给定活性阈值的实验状态，不表示该配体对其他蛋白具有广泛选择性；dual 仅表示两个研究靶点上的实验活性均达到阈值，不涉及两端活性平衡、细胞效应或治疗价值。后文所称单靶选择性配体均采用这一相对定义。

候选靶对的双向选择性供给另用严格 6.5/5.5 活性分离规则评估：活性端 pChEMBL \(\geq 6.5\)，低活性端 \(\leq 5.5\)，位于 5.5–6.5 灰区的配体不计入该严格分类。严格规则用于候选靶对供给评估和六对的面板构建，不是 Table 2 的第二套主定义。主分析标签是 \(\theta=6.0\) 四状态。活性来源不统一：EGFR/HER2、PIK3CA/mTOR 和多数 AChE/BChE 为 `chembl_assay_adjudicated`；Track-B 五对为 `chembl37_dump_panel`；五个 AChE 配体为 `panel_pchembl_no_audit_rows`（仅 AB_056 进入主集）。不得声称八对均为 assay-level adjudicated。

为考察活性数据处理方式的影响，在保持评价集成员和对接评分不变的条件下，对全部八对用同一套 ChEMBL 37 记录比较最大值与中位数聚合（Table S3）。该检查不替代 Table 2。

### 2.3 靶对筛选与评价集构建

**靶对供给筛选。** 候选靶对须同时支持四状态实验评价和非共价小分子对接。ChEMBL 37 中共有人源、单组分 SINGLE PROTEIN 靶点 5,869 个，其中 4,672 个至少有一个合格分子。候选范围为这 4,672 个靶点构成的 10,911,456 个无序对，不含自身对。至少有一个配体在两端均有定量活性的靶对进入普查（2,164,618 对）。后续门槛使用同一套 ChEMBL 37 记录。两端均至少 10 个配体时剩 63,790 对；\(\theta=6.0\) 下 dual、A-only、B-only 各不少于 10 时剩 5,253 对；严格 6.5/5.5 双向选择性规则要求 A-only 和 B-only 各不少于 50 个配体时剩 86 对。剔除 qHTS 枢纽蛋白、CYP ADME 面板和锌依赖金属酶后剩 26 对。两端各至少 5 个人源全配体结构（\(\leq 3.5\) Å，至少一个非聚合物配体）后剩 19 对。该 3.5 Å 门槛是靶对供给筛选。实际对接的 mTOR 结构 4JT6 分辨率为 3.60 Å，不改写为 3.5 Å 结构。分子身份和大小筛选要求具有有效分子结构、游离形式分子量为 150–750 Da、重原子数为 10–60，且不含规则中列出的金属元素。筛选后两个严格选择性类别仍各不少于 50 个配体的靶对共 17 对，覆盖 12 个靶点体系。最后要求具有与对接流程相容的非共价小分子位点。膜蛋白 GPCR、转运体、所选可逆共价复合物及结构域归属不明确的体系未纳入；F2/PRSS1 则因 PRSS1 在本研究的药理学设定中作为非期望作用靶点而排除。

经上述研究范围和对接流程检查后剩 7 对：JAK1/JAK2、JAK1/TYK2、PIK3CA/mTOR、AChE/BChE、F2/F10、PPARG/PPARA 和 PPARA/PPARD。CTSK/CTSS 通过了更早的供给与结构门槛，但因两端共晶均为可逆共价、需不同对接处理而被排除。EGFR/HER2 作为供给受限的对接案例保留。它未满足严格 6.5/5.5 选择性供给标准（最小选择性类别计数为 7）。它具有合适的人源全配体结构、可由共晶定义的对接位点，并在 \(\theta=6.0\) 下有足够 dual、A-only 和 B-only 配体做方向性评价。这八对即主对接评价集（Table 1）。可执行筛选规则、人工判断及各靶对的纳入和排除理由见 Table S9。

**评价面板构建。** 因样本供给规模和结构可行性不同，部分靶对采用不同的候选池规则、类别配额或骨架限制；具体规则见 Table 1。EGFR/HER2 与 PIK3CA/mTOR 从 \(\theta=6.0\) 候选池抽样，其余六对从严格 6.5/5.5 候选池按配额抽出。EGFR/HER2 和 PIK3CA/mTOR 分别限制同一类别中单一 Bemis–Murcko 骨架的最大重复数为 5 和 2，其余靶对未另设骨架上限。

**统一主要分析。** 所有主要 AUROC 均按 \(\theta=6.0\) 重新确定实验状态，并仅纳入两端均有有效 Vina 分数的配体进入主方向性表，因此 n_scored 可以低于 n_panel。各方向 AUROC 使用对应口袋评分。

**Table 1.** 双靶评价集的组成与主要对接设置。配额为构建目标（dual / A-only / B-only / neither）。n_scored 为两端均有有效 Vina 分数、进入方向性主 AUROC 的 dual / A-only / B-only 计数。受体分辨率见 Table S2。PPARG 口袋 9V8H 中保留与受体复合的 PG08-NL 肽。

| 靶对 | 候选池 | 配额 (D / A / B / N) | 受体 PDB (A / B) | n_scored (dual / A-only / B-only) | Vina exhaustiveness |
|------|--------|---------------------:|------------------|----------------------------------:|--------------------:|
| EGFR/HER2 | \(\theta=6.0\) | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 37 / 31 | 8 |
| JAK1/JAK2 | 严格 6.5/5.5 | 32 / 32 / 32 / 14 | 6N7A / 8BXH | 32 / 32 / 32 | 8 |
| JAK1/TYK2 | 严格 6.5/5.5 | 32 / 32 / 32 / 14 | 6N7A / 3LXP | 31 / 32 / 32 | 8 |
| PIK3CA/mTOR | \(\theta=6.0\) | 18 / 14 / 12 / 4 | 4L23 / 4JT6 | 18 / 14 / 12 | 16 |
| AChE/BChE | 严格 6.5/5.5 | 28 / 28 / 28 / 16 | 4EY7 / 4BDS | 27 / 26 / 28 | 8 |
| F2/F10 | 严格 6.5/5.5 | 32 / 32 / 32 / 14 | 4UDW / 2JKH | 31 / 32 / 32 | 8 |
| PPARG/PPARA | 严格 6.5/5.5 | 32 / 32 / 32 / 14 | 9V8H / 6LXA | 32 / 31 / 32 | 8 |
| PPARA/PPARD | 严格 6.5/5.5 | 32 / 32 / 32 / 14 | 6LXA / 5U3Q | 32 / 32 / 32 | 8 |

### 2.4 受体与配体准备及分子对接

#### 2.4.1 受体准备与对接区域定义

准备脚本记录了各受体的链选择、替代原子位置处理和残基模板设置，详见 Table S2。归档的准备记录未规定统一的 pH 依赖质子化或缺失环区重建流程。

各靶点使用的人源实验结构 PDB 见 Table 1，共晶配体和对接盒见 Table S2。相应共晶配体分别用于确定各受体的对接区域。JAK1 6N7A 同时用于 JAK1/TYK2 和 JAK1/JAK2，PPARA 6LXA 同时用于 PPARG/PPARA 和 PPARA/PPARD。TYK2 3LXP 使用 JH1 ATP 结合位点。PPARG 9V8H 中保留与受体复合的 PG08-NL 肽，仅移除共晶小分子和结晶水。

纳入受体前均核对蛋白身份、物种、结合位点和共晶配体位置。初始对接盒由共晶配体重原子的笛卡尔坐标范围确定，并沿 \(x\)、\(y\) 和 \(z\) 三个方向分别向两侧扩展 5 Å。若扩展后任一方向的边长仍小于 20 Å，则继续扩展至 20 Å。去除用于确定对接区域的结晶水和共晶配体后，生产分析使用已沉积的受体 PDBQT（`data/provenance/receptor_input_registry.csv`）。14 个槽位并非都有可恢复的 Meeko 或 Protein Preparation Wizard 制备记录，本文不声称存在统一受体准备协议。存在替代原子位置时，沉积文件保留无 altLoc 标记或标记为 A 的原子。对接盒坐标见 Table S2。

#### 2.4.2 配体准备

以 2.2 节处理后的 ChEMBL SMILES 作为配体输入。Track-B、PIK3CA/mTOR 与 AChE/BChE 的已恢复脚本使用 RDKit 加氢、ETKDGv3（seed 20260727）、MMFF ≤200 和 Meeko PDBQT。这些包的精确生产配体 PDBQT 未入库；零对接统计使用已沉积分数表。EGFR/HER2 原始配体 PDBQT 不可恢复；uniform RDKit/Meeko 重建在提升后作为当前 EGFR 配体准备。未进行系统的质子化状态和互变异构体枚举。本仓库不能从 SMILES 完整重跑全部原始对接。

#### 2.4.3 AutoDock Vina 对接与评分

主要对接采用 AutoDock Vina 1.2.7 和默认 Vina 评分函数，作为双靶排序评价的主引擎。每个配体–受体组合最多输出 9 个姿态，`energy_range` 设为 3 kcal mol\(^{-1}\)。AChE 4EY7 与 TYK2 3LXP 保存 8 个姿态。PIK3CA/mTOR 使用 exhaustiveness = 16，以获得满足重对接近天然构象覆盖要求的搜索设置；其余靶对使用 8。所有主要 Vina 分析均以第一保存模型中 `REMARK VINA RESULT` 报告的排名第 1 姿态 affinity 作为配体–受体评分。该字段即 Table 2 所用生产 `*_affinity` 列。主要 Vina 设置见 Table S1。

#### 2.4.4 共晶配体重对接

共晶配体重对接用于评价实验结合构象的恢复情况。分别报告排名第 1 姿态的重原子 RMSD、前 3 个排名姿态中的最低 RMSD（top-3）以及全部保存姿态中的最低 RMSD。最后一项以 2.0 Å 为界评价搜索覆盖，并不表示排名第 1 姿态即为近天然构象；top-1 用于评价姿态排序。14 个主受体的 RMSD 由保存姿态通过化学对应的 RDKit CalcRMS 重新计算，未重新对接（Table S2）；AChE 4EY7 与 TYK2 3LXP 保存 8 个姿态，其余受体保存 9 个；PIK3CA/mTOR 使用图自同构映射，EGFR 3POZ 使用重建 QC。[8]

#### 2.4.5 替代评分与独立对接

独立 GNINA 分析是针对 EGFR/HER2、PIK3CA/mTOR 和 JAK1/TYK2 的敏感性检查，用于对照类别差异明显的案例及 PIK3CA/mTOR 流程敏感性案例，并非预设的八个靶对引擎比较。

使用 RTMScore 和 GNINA 1.3.2 CNN 对 Vina 生成的全部可用姿态进行补充评分。RTMScore 取全部保存姿态中的最高评分作为配体级结果；GNINA 以 CNN affinity 作为主要 CNN 重评分读数，CNNscore 作为补充分析。上述结果均来自 Vina 已生成姿态的重新评分，不属于独立姿态生成。

在 EGFR/HER2、PIK3CA/mTOR 和 JAK1/TYK2 上另用 GNINA 1.3.2 独立生成姿态。该分析使用与主要分析相同的受体、配体和对接区域，三个靶对的 exhaustiveness 分别为 8、16 和 8，每个配体最多输出 9 个姿态。配体级评分取排名第 1 姿态的 minimizedAffinity，并使用其相反数进行统计分析，使数值越高表示预测结合越有利。独立 GNINA 并未对每个配体都返回双端评分：EGFR/HER2 的 dual–neither 比较使用 n_neither = 11（缺 EH120_109），而主要 Vina dual–neither 比较为 12；PIK3CA/mTOR 的 A-only 为 n = 13 而非 14；JAK1/TYK2 的 Dual/A-only/B-only/neither 计数为 30/32/29/14，而非 31/32/32/14。该分析用于观察主要结果在另一套姿态生成和评分流程下是否仍然存在，而不用于比较 GNINA 与 Vina 的总体性能（Table S7）。

### 2.5 评价指标与统计分析

#### 2.5.1 主要方向性判别指标

比较 dual 与 A-only 时，两类配体在靶点 A 上均达到实验活性阈值，因此使用靶点 B 的评分：\(\mathrm{AUROC}_{D/A}(B)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{A\text{-}only};\;S_{B})\)。比较 dual 与 B-only 时，使用靶点 A 的评分：\(\mathrm{AUROC}_{D/B}(A)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{B\text{-}only};\;S_{A})\)。所有方向性分析均以 dual 为正类。

AutoDock Vina 输出的 affinity 越低表示预测结合越有利。为统一 AUROC 的评分方向，将 Vina affinity 转换为 \(S_{\mathrm{Vina}}=-E_{\mathrm{Vina}}\)，因此较高的 \(S_{\mathrm{Vina}}\) 表示更有利的预测结合。两个方向的 AUROC 是一级终点。两者中的较低值定义为 \(\mathrm{summary}_{\min}\)，仅用于描述较弱方向，不是总体性能指标，也不作为新的配体评分函数：\(\mathrm{summary}_{\min}=\min[\mathrm{AUROC}_{D/A}(B),\;\mathrm{AUROC}_{D/B}(A)]\)。两个 AUROC 取最小具有选择性下偏，因此不把它当作普通独立终点。

#### 2.5.2 传统双口袋排序与对照类别比较

双靶虚拟筛选通常合并两个口袋的评分来排序候选分子。对两端均有评分的配体，传统读数为平均分 \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\)，以 dual 为正类、neither 为对照（Table 3）。另将 A-only、B-only 和 neither 合并为非 dual 类的 AUROC 随分数表归档。

为把虚筛对照类别的影响与评分聚合分开，保持同一口袋评分不变，比较不同实验状态对照下的 AUROC。靶点 A 评分用于 dual–B-only 及对应的 dual–neither 比较；靶点 B 评分用于 dual–A-only 及对应的 dual–neither 比较（Table S4）。更换对照类别也同时更换了具体化合物。另取两端评分中的较低值 \(S_{\mathrm{worst}}=\min(S_{A},S_{B})\) 作为 AND 型双口袋过滤。以 dual 的 \(S_{\mathrm{worst}}\) 中位数为阈值，统计 dual 保留数、dual recall、dual precision 以及保留的单靶活性配体数（Table S10）。八对双口袋平均排序下的 Top 10% 类别组成见 Figure 2D。截断为完整四状态面板上的 \(k=\lceil 0.10\,n\rceil\)，以便不同规模的靶对按同一筛选比例比较。Top 10% dual 比例为 \(\mathrm{dual}/k\)。相对全面板 dual 基线，\(\mathrm{EF}_{\mathrm{dual},10\%}=(\mathrm{dual}/k)/(n_{\mathrm{dual}}/n)\)。八个主评价靶对的同一套排序与过滤规则见 Table S10。

#### 2.5.3 置信区间与重采样分析

Table 2 同时报告两条方向性 AUROC 及其描述性最小值的逐项 95% 置信区间。主分析采用类别分层非参数百分位 bootstrap，B = 2000，种子 20260729。方向性 AUROC 在对应的两个实验状态类别内分别有放回抽样，并保持原类别样本量。对于 \(\mathrm{summary}_{\min}\)，每次重采样只抽取一次 dual 配体，并将同一次 dual 抽取同时用于两个口袋评分；A-only 与 B-only 在各自类别内独立重采样；在该次重采样内重新计算两个方向 AUROC 后再取较小值。固定评分 \(\Delta\)AUROC 中，dual–selective 与 dual–neither 共用同一次 dual 抽取，两类阴性对照各自独立重采样。matched 与 mismatched 口袋比较采用配体水平配对重采样，使同一配体的两条口袋评分同时被抽中或同时不被抽中。点估计由完整分析样本直接计算。配体水平非分层 bootstrap 仅作为敏感性分析，不替换类别分层主区间。按 Bemis–Murcko 骨架簇和文献连通簇进行的簇 bootstrap 用于两个最大固定评分差的来源依赖性敏感性（Table S4），不替换配体水平主区间。二项正态 detectable-effect 仿真复用同一类别分层、共享 dual 的 bootstrap、当前八对类别样本量和 B = 2000、种子 20260729；它估计在指定真 AUROC 下 CI 排除 0.5 的概率，不是观测功效。区间未做多重比较校正。

比较不同评分方法或对应与非对应口袋评分时，各方案在每次 bootstrap 中使用相同的配体重采样结果以保持配对。为评估化学骨架相关性和文献来源相关性的影响，另以 Bemis–Murcko 骨架簇和文献连通簇为单位进行簇水平 bootstrap。两个最大固定评分差值的骨架簇与文献簇区间见 Table S4。

### 2.6 基线、对照与敏感性分析

#### 2.6.1 配体化学对照

每个方向性模型的 GroupKFold 折数取 5、Bemis–Murcko 骨架数及两类样本量中的最小值。逻辑回归采用 C = 1.0、最多 4000 次迭代，GroupKFold 使用默认的不打乱划分；完整设置见 Table S1。

为检验不经对接的配体排序能否分开同一套实验类别，计算 ECFP4 指纹（radius = 2，2048 bits）以及分子量、重原子数、cLogP 和 TPSA。对四个单一物化描述符分别计算两个方向的 AUROC 及其 \(\mathrm{summary}_{\min}\)。全面板上 \(\mathrm{summary}_{\min}\) 最高的单一描述符只作为描述性单变量筛选，不是经过选择校正的预测估计。化学对照的预测读数来自嵌套骨架 GroupKFold、仅在训练折上选择描述符后的折外 AUROC（Table S5）。ECFP4 骨架 GroupKFold 流程未改。

ECFP4、仅对接评分（docking-only）和 ECFP4+对接评分模型均采用未做特征缩放的逻辑回归。以 Bemis–Murcko 骨架为分组变量，在相同的 GroupKFold 划分下获得折外（out-of-fold）预测，并据此计算 AUROC。仅对接评分模型仅以相应方向的对接评分作为输入。该 AUROC 由折外预测值计算，而主要分析直接按原始对接评分排序，两者计算方式不同，数值并不完全一致。ECFP4 与 ECFP4+对接评分模型的 AUROC 差值，用于衡量对接评分对配体模型的增量。Figure 3A 比较结构对接排序与配体排序，不是公平的预测竞赛。敏感性分析在相同划分的各训练折上拟合 StandardScaler（Table S5）。

#### 2.6.2 结构归因与评分对应性检验

由于 \(\Delta\) 比较的是两条方向性结果中的较低值，正值不表示两个方向都具有优势；更换评分通道后，较弱方向的身份可能发生变化。

双靶对接把每个评分指派给一个口袋。对应评分使用口袋 B 作 dual–A-only 比较、口袋 A 作 dual–B-only 比较。非对应对照在不重新对接的前提下交换这两条已计算的评分通道。分析比较两种指派下的 \(\mathrm{summary}_{\min}\)，并以配对 bootstrap 估计 \(\Delta=\mathrm{summary}_{\min}^{\mathrm{matched}}-\mathrm{summary}_{\min}^{\mathrm{mismatched}}\)。正值表示对应口袋的较弱方向更高。较弱方向是否切换见 Table S6。

受体结构敏感性主要在 PIK3CA/mTOR 中评价。在保持 mTOR 4JT6 不变时，分别以 4JPS 和 5DXT 替换 PIK3CA 4L23；另以 4JSX 替换 mTOR 4JT6。除被替换的受体外，其余配体集合、实验状态、评分定义和统计方法保持不变（Table S7）。

#### 2.6.3 稳健性与外部数据可用性分析

通过改变活性阈值及重复活性记录的汇总方式重新确定实验状态，并在保持评价集成员和对接评分不变的条件下重复主要分析（Table S3）。对于具有足够剩余候选配体的靶对，在排除主评价集成员后构建未使用池留出集。留出配体使用固定抽样规则，并限制单一 Bemis–Murcko 骨架的过度重复。该分析使用不同配体重新计算主要方向性 AUROC，用于评价对接结果对评价集成员组成的敏感性。由于这些配体仍来自同一数据来源，该分析属于内部稳健性检验（Table S6）。另通过五个固定随机种子重复主要 Vina 对接，所用对接盒与主分析相同（EGFR/HER2 使用校正后的 cognate 重原子盒子），并在 PIK3CA/mTOR 中比较不同 exhaustiveness 设置。完整病例计数和固定成员交集留在仓库。二项正态 detectable-effect 仿真使用当前八对类别样本量和与 Table 2 相同的类别分层、共享 dual bootstrap。

BindingDB[16] 和 PubChem 对全部八个靶对检索独立外部对接集的候选数据。外部对接采用 \(\theta=6.0\) 独立性过滤（Table S8）。开发分子集合包括主评价集、扩大的 PIK3CA/mTOR PM110 面板和内部留出集。剔除共享文献来源、结构重复和 ECFP4 Tanimoto \(\geq 0.70\) 的分子后，要求 dual、A-only、B-only 各 \(n\geq 20\) 且每类至少 3 个来源。只有同时满足上述条件的靶对才进入外部对接评价。

### 2.7 软件与可重复性

分子结构处理、指纹和描述符计算及统计分析均在 Python 环境中完成，主要使用 RDKit、NumPy、pandas、SciPy 和 scikit-learn（分析冻结：Python 3.12.3；Table S1）。Track-B、PIK3CA/mTOR 与 AChE/BChE 配体 PDBQT 遵循已恢复的 RDKit/Meeko 脚本，精确生产文件未入库。受体使用已沉积 PDBQT；14 个槽位没有统一的 Meeko 或 Protein Preparation Wizard 记录（`docs/RECEPTOR_PREP_PROVENANCE.md`）。主要分子对接采用 AutoDock Vina 1.2.7；GNINA 1.3.2 用于补充评分和部分靶对的独立姿态生成；RTMScore 用于 Vina 姿态的另一种补充重评分。

软件版本、随机种子和主要计算参数汇总于 Table S1。各受体的对接区域和共晶配体重对接结果见 Table S2。

## 3. 结果

### 3.1 成对实验数据供给与对接评价集

能够支持四状态对接评价的成对实验标签随样本要求提高而迅速减少。在人源单组分 SINGLE PROTEIN 靶点的无序对中，2,164,618 对至少有 1 个两端均被测定的配体，63,790 对具有不少于 10 个此类配体。若再要求在 \(\theta=6.0\) 下 dual、A-only 和 B-only 均不少于 10 个，仅剩 5,253 对；采用严格 6.5/5.5 活性标准后减至 86 对。Figure 1C 概括这些计数。后续结构和对接位点门槛见 2.3 节。

经过供给、结构和对接相容性检查后，主要评价包括 EGFR/HER2、JAK1/JAK2、JAK1/TYK2、PIK3CA/mTOR、AChE/BChE、F2/F10、PPARG/PPARA 和 PPARA/PPARD（Table 1）。EGFR/HER2 与 PIK3CA/mTOR 从 \(\theta=6.0\) 候选池抽样，其余六对从严格 6.5/5.5 池抽出。四状态分类要求同一配体在两个靶点上均有实验测量，因此成对覆盖共同限制了对接评价集规模。

![Figure 1](../figures/jcim_article/Fig1_four_state_and_supply.png)

**Figure 1.** 四状态双靶评价与数据供给。(A) 按阈值 \(\theta\) 定义的四种实验状态；(B) 两个方向性评价任务：dual vs A-only 使用靶点 B 评分，dual vs B-only 使用靶点 A 评分，\(\mathrm{summary}_{\min}=\min[\mathrm{AUROC}_{D/A}(B),\;\mathrm{AUROC}_{D/B}(A)]\)；(C) 成对实验覆盖在四状态要求下迅速减少。主要评价的八个靶对按 Table 1 的面板构建与结构条件纳入，不是普查计数的直接延续。

### 3.2 方向性对接表现与候选排序后果

固定同一靶点对接评分后，更换实验状态对照所产生的 AUROC 差异因靶对和方向而异。EGFR/HER2 与 JAK1/TYK2 的固定评分差最大。EGFR/HER2 使用 EGFR 口袋评分时，dual–B-only 的 AUROC 为 0.334，将对照换成 neither 后升至 0.780，差值 0.446 [0.259, 0.632]。JAK1/TYK2 使用 JAK1 口袋评分时，对应差值为 0.444 [0.261, 0.630]。其余靶对的差异更小，或其置信区间包含 0（Figure 2A；Table S4）。簇重采样见 3.6 节。

在统一 \(\theta=6.0\) 标签下，两个方向性 AUROC 介于 0.334 至 0.728，描述性较弱臂 \(\mathrm{summary}_{\min}\) 介于 0.334 至 0.692。PPARG/PPARA 的 \(\mathrm{summary}_{\min}\) 区间完全高于 0.5（0.649 [0.511, 0.746]），EGFR/HER2 与 F2/F10 完全低于 0.5（0.334 [0.197, 0.471] 和 0.345 [0.216, 0.482]），其余 5 对穿过 0.5（Figure 2B；Table 2）。这些数值是靶对特异的，不是八对排行。在复用当前类别样本量和 Table 2 同类 bootstrap 的二项正态仿真中，当真较弱臂 AUROC 为 0.60 时，`summary_min` CI 排除 0.5 的概率最高为 0.071；为 0.75 时，除 PIK3CA/mTOR（0.466）外均不低于 0.813。该仿真不是观测功效（Table S4）。

**Table 2.** 八个主评价靶对上的口袋匹配方向 AUROC（Vina，统一 \(\theta=6.0\)）。表中类别样本量为 n_scored（dual / A-only / B-only）。物化描述符基线见 Table S5。

| 靶对 | n_scored (dual / A-only / B-only) | Dual vs A-only（口袋 B）[95% CI] | Dual vs B-only（口袋 A）[95% CI] | summary_min [95% CI] |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|
| EGFR/HER2 | 28 / 37 / 31 | 0.663 [0.521, 0.792] | 0.334 [0.197, 0.471] | 0.334 [0.197, 0.471] |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 [0.448, 0.723] | 0.728 [0.602, 0.842] | 0.588 [0.448, 0.716] |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 [0.434, 0.717] | 0.365 [0.233, 0.508] | 0.365 [0.233, 0.505] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 [0.528, 0.889] | 0.692 [0.491, 0.868] | 0.692 [0.480, 0.802] |
| AChE/BChE | 27 / 26 / 28 | 0.652 [0.492, 0.806] | 0.606 [0.443, 0.749] | 0.606 [0.439, 0.735] |
| F2/F10 | 31 / 32 / 32 | 0.413 [0.271, 0.571] | 0.345 [0.218, 0.491] | 0.345 [0.216, 0.482] |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 [0.512, 0.781] | 0.706 [0.573, 0.830] | 0.649 [0.511, 0.746] |
| PPARA/PPARD | 32 / 32 / 32 | 0.647 [0.509, 0.776] | 0.446 [0.301, 0.592] | 0.446 [0.301, 0.590] |

采用双靶虚拟筛选常用的双口袋平均排序时，dual–neither 比较 AUROC 介于 0.514 至 0.770。EGFR/HER2 为 0.759 [0.551, 0.926]，JAK1/TYK2 为 0.770 [0.613, 0.906]，而其方向性 \(\mathrm{summary}_{\min}\) 分别为 0.334 和 0.365（Figure 2C；Table 3）。PIK3CA/mTOR 的 neither 样本为 \(n=4\)。Table 3 同时改变了评分聚合和对照类别，对照类别的影响以 Figure 2A 的固定口袋分析为准。

![Figure 2](../figures/jcim_article/Fig2_negative_class_formulation.png)

**Figure 2.** 核心双靶候选排序结果。(A) 固定评分 \(\Delta\)AUROC（dual–neither 减去 dual–匹配单靶活性），八对两个评分方向，误差棒为配体水平 bootstrap 95% 置信区间；(B) 两个方向性 AUROC：dual vs A-only 使用口袋 B，dual vs B-only 使用口袋 A；(C) 方向性 \(\mathrm{summary}_{\min}\) 与双口袋平均 dual–neither AUROC。面板 C 不是固定评分比较。n = 4 标注 PIK3CA/mTOR 的 neither 样本。(D) 八对在双口袋平均 Vina 排序下的 Top 10% 类别组成，\(k=\lceil 0.10\,n\rceil\)。条形为 dual / A-only / B-only / neither 的 100% 堆叠；右侧标注 \(k\)。\(\mathrm{EF}_{\mathrm{dual},10\%}\) 见 Table 3。

**Table 3.** 双口袋平均 Vina 排序的候选排序后果（统一 \(\theta=6.0\)）。dual–neither 使用 \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\)。Top 10% 组成在完整四状态面板上取 \(k=\lceil 0.10\,n\rceil\)。\(\mathrm{EF}_{\mathrm{dual},10\%}=(\mathrm{dual}/k)/(n_{\mathrm{dual}}/n)\)。PIK3CA/mTOR 的 neither 样本较小（n = 4）。

| 靶对 | 双口袋平均 D-vs-neither AUROC [95% CI] | n_neither | Top 10% D/A/B/N | \(\mathrm{EF}_{\mathrm{dual},10\%}\) |
|------|--------------------------------------------:|----------:|----------------:|-------------------------------------:|
| EGFR/HER2 | 0.759 [0.551, 0.926] | 12 | 1 / 4 / 6 / 0 | 0.351 |
| JAK1/JAK2 | 0.730 [0.558, 0.882] | 14 | 6 / 5 / 0 / 0 | 1.875 |
| JAK1/TYK2 | 0.770 [0.613, 0.906] | 14 | 1 / 3 / 7 / 0 | 0.320 |
| PIK3CA/mTOR | 0.514 [0.222, 0.806] | 4 | 4 / 1 / 0 / 0 | 2.133 |
| AChE/BChE | 0.659 [0.484, 0.820] | 14 | 5 / 3 / 1 / 1 | 1.759 |
| F2/F10 | 0.519 [0.341, 0.691] | 12 | 4 / 1 / 6 / 0 | 1.255 |
| PPARG/PPARA | 0.685 [0.498, 0.859] | 14 | 7 / 3 / 0 / 1 | 2.168 |
| PPARA/PPARD | 0.565 [0.373, 0.748] | 14 | 5 / 2 / 3 / 1 | 1.563 |

同一套双口袋平均排序按各面板 Top 10%、\(k=\lceil 0.10\,n\rceil\) 应用于全部八个主评价集，以便不同规模的靶对按同一筛选比例比较（Figure 2D；Table 3；Table S10）。Top 10% 中 dual 计数从 1/11（EGFR/HER2、JAK1/TYK2；dual 比例 0.091）到 4/5（PIK3CA/mTOR；0.800）和 7/11（PPARG/PPARA；0.636）。相对各面板 dual 基线，\(\mathrm{EF}_{\mathrm{dual},10\%}\) 从 0.320（JAK1/TYK2）和 0.351（EGFR/HER2）到 2.133（PIK3CA/mTOR）和 2.168（PPARG/PPARA）。EGFR/HER2 为 1 / 4 / 6 / 0（\(\mathrm{EF}_{\mathrm{dual},10\%}=0.351\)）；AChE/BChE 为 5 / 3 / 1 / 1（\(\mathrm{EF}_{\mathrm{dual},10\%}=1.759\)）。小于 1 表示 Top 10% 中 dual 比例低于全面板基线。五个靶对的 Top 10% 不含 neither；AChE/BChE、PPARG/PPARA 和 PPARA/PPARD 各占 1 个位置。每个靶对的 Top 10% 中都至少保留 1 个单靶活性配体。这些计数是回顾性评价面板上的描述性操作点，不是对接质量的八对排行，也不是前瞻命中率。Dual+A-only+B-only 上的 AND 过滤见 Table S10。

### 3.3 配体化学基线与对接增量判别

若双靶对接所排序的实验类别本身已有化学差异，则不依赖受体的基线也能分开同一套筛选状态。AChE/BChE 仅用 TPSA 时，dual–A-only 和 dual–B-only 的 AUROC 分别为 0.742 和 0.801（Figure S1；Table S5）。PIK3CA/mTOR 中最佳单一描述符（重原子数）的 \(\mathrm{summary}_{\min}\) 为 0.463（Table S5）。Vina 与最佳单一描述符的差值因靶对而异：八对中 5 对的 95% 区间包含 0，JAK1/TYK2、PIK3CA/mTOR、F2/F10 不包含 0（Figure S1；Table S5）。

在 Bemis–Murcko 骨架分组交叉验证下，仅含配体的 ECFP4 模型已能分开若干方向性比较（Figure 3A）。将对应方向的对接评分加入后，16 个方向的 AUROC 最多变化 0.023，且没有一致的变化方向（Figure 3B；Table S5）。这些面板比较的是同一套实验类别上的仅配体排序与结构基础排序，不是正式算法排行。在该回顾性设定下，对接评分未在配体模型之上提供稳定的排序增量。

![Figure 3](../figures/jcim_article/Fig3_ligand_chemistry.png)

**Figure 3.** 配体化学与对接增量信息，不是公平的预测竞赛。(A) 原始分数排序的 Vina AUROC 与骨架分组交叉验证下 ECFP4 折外预测 AUROC，八对 16 个方向；(B) 在相同划分和相同配体上将对应方向 Vina 评分加入 ECFP4 后的 \(\Delta\)AUROC，连线连接同一方向的两个点，不是置信区间。AChE/BChE 的 TPSA 分布与单一描述符森林见 Figure S1。

### 3.4 口袋对应性

若方向性 AUROC 来自基于结构的对接，则在不重新对接的前提下交换两条已计算的口袋评分，较弱臂汇总应变差。主评价集中 AChE/BChE（0.177 [0.053, 0.291]）与 EGFR/HER2（0.107 [0.006, 0.220]）的 matched−mismatched \(\mathrm{summary}_{\min}\) 95% 区间排除 0。其余六对主集区间包含 0。七个可用未使用池内部留出集的差值区间亦均包含 0（−0.079 至 +0.150）（Figure 4；Table S6）。EGFR/HER2 无未使用池留出集。区间包含 0 并不证明受体没有信息，只说明该对照不支持将该对解释为 matched-pocket 优势。

![Figure 4](../figures/jcim_article/Fig4_mismatched_pocket.png)

**Figure 4.** 口袋对应性。(A) 八对主评价集 matched−mismatched \(\Delta\mathrm{summary}_{\min}\) 及配体水平 bootstrap 95% 置信区间；(B) 七对未使用池内部留出集的 matched−mismatched \(\Delta\mathrm{summary}_{\min}\)。颜色区分主集与留出集，不表示区间是否排除 0。AChE/BChE 与 EGFR/HER2 主集区间排除 0；七个留出集区间均包含 0。\(\dagger\) 表示 EGFR/HER2 无未使用池留出集。主集与留出集的 \(\mathrm{summary}_{\min}\) 见 Table S6。

### 3.5 计算实现稳健性

独立 GNINA 1.3.2 姿态生成使用与主分析相同的受体、配体和对接盒，覆盖三对（Figure 5A；Table S7）。EGFR/HER2 的 dual–neither AUROC 为 0.705 [0.465, 0.910]（\(n_{\mathrm{neither}}=10\)），dual–B-only 为 0.227 [0.104, 0.373]。JAK1/TYK2 的 dual–neither 为 0.705 [0.524, 0.872]，方向性 \(\mathrm{summary}_{\min}\) 为 0.317 [0.187, 0.455]。PIK3CA/mTOR 的 \(\mathrm{summary}_{\min}\) 为 0.633，较弱臂 dual–A-only 为 0.633 [0.410, 0.769]，dual–neither 为 0.569 [0.236, 0.889]（\(n=18/13/12/4\)）。这些运行用于观察主要观察能否在另一套姿态生成流程下保留，不是 GNINA 与 Vina 的总体性能比较。

在 PIK3CA/mTOR 中，将 PIK3CA 受体由 4L23 替换为 4JPS 后，\(\mathrm{summary}_{\min}\) 从 0.692 [0.480, 0.802] 降至 0.486 [0.264, 0.694]；替换为 5DXT 后为 0.505 [0.296, 0.713]；将 mTOR 4JT6 替换为 4JSX 后为 0.639 [0.435, 0.783]（Figure 5B；Table S7）。

五个固定 Vina 随机种子产生的数值波动相对任务和受体效应更为有限（Figure 5C；Table S7）。EGFR/HER2 使用与 Table 2 相同的校正 cognate 重原子盒子；生产种子与主分析 \(\mathrm{summary}_{\min}\) 0.334 一致，五种子范围为 0.317–0.349。EGFR/HER2 的固定评分任务差在五个种子上均为正。PPARG/PPARA 的同姿态 RTMScore 与 GNINA CNN 重评分见 Table S7。PPARG/PPARA 是唯一主分析 Vina \(\mathrm{summary}_{\min}\) 区间完全高于 0.5 的靶对（0.649 [0.511, 0.746]）；RTMScore 重评分后降至 0.369 [0.237, 0.478]，GNINA CNN 重评分降至 0.500 [0.347, 0.632]。

![Figure 5](../figures/jcim_article/Fig5_computational_realization.png)

**Figure 5.** 单一对接实现的计算稳健性。(A) EGFR/HER2、PIK3CA/mTOR 与 JAK1/TYK2 上独立 GNINA 与 Vina：实心为 Vina，空心为 GNINA；蓝为方向性 \(\mathrm{summary}_{\min}\)，橙为 dual–neither。灰线连接同一引擎的两项任务，不是置信区间。(B) PIK3CA/mTOR 受体替换：主结构 4L23 / 4JT6 以及 4JPS、5DXT、4JSX，误差棒为配体水平 bootstrap 95% 置信区间。(C) 八个靶对在五个 Vina 随机种子上的 \(\mathrm{summary}_{\min}\) 范围、中位数与生产种子（`same_protocol_as_primary=1`）。AChE/BChE 的种子实现不使用主分析成员集（`same_membership_as_primary=0`；主分析 n=95，available-case n=94，交集 n=88）；固定成员敏感性不改变定性结论。该图只回答结果是否依赖单一对接实现。

### 3.6 敏感性与外部数据范围

从严格 6.5/5.5 候选池构建的评价集，在多个阈值下保留了相同的主要类别组成，因此相应 AUROC 变化较小。EGFR/HER2 和 PIK3CA/mTOR 的类别组成随阈值变化更明显，其方向性估计也发生变化（Figure S4；Table S3）。最大与中位数聚合使用同一套合格原始记录、同一分子交集和当前对接分数，而不是把全部已评分配体与 dump 缺失记录混为一谈。EGFR/HER2 有 5/109 个类别翻转（主分析 \(\mathrm{summary}_{\min}\) 0.334）；AChE/BChE 有 1/94 个翻转（CHEMBL659；max 0.606 变为 0.629）；PPARA/PPARD 有 1/110 个翻转（CHEMBL121；\(\mathrm{summary}_{\min}\) 仍为 0.446）。其余靶对保持类别组成和 \(\mathrm{summary}_{\min}\) 点估计（Table S3）。

在排除主评价集分子后，基于剩余候选分子构建的未使用池留出集显示，结果存在一定的样本组成依赖：AChE/BChE、PIK3CA/mTOR 与 JAK1/JAK2 与主评价接近；JAK1/TYK2 有所上升；F2/F10 与 PPARA/PPARD 仍处于较低水平；PPARG/PPARA 则由 0.649 降至 0.535 [0.360, 0.705]（Table S6）。EGFR/HER2 无同等留出集。EGFR/HER2 口袋 A 的配体层固定评分差值为 0.446 [0.259, 0.632]。将校正盒评分代入原冻结的骨架簇与文献簇后，EGFR/HER2 簇区间为 [0.232, 0.637]（骨架）和 [0.117, 0.636]（文献），均排除 0。JAK1/TYK2 骨架簇区间为 [0.220, 0.631]，排除 0；文献簇无法在本冻结中重算（配体–文献分组映射与 ChEMBL 37 sqlite 均不可用）；此前存档区间包含 0，不作为当前计算结果（Figure S4；Table S4）。

PIK3CA/mTOR 的面板规模和 exhaustiveness 敏感性结果见 Figure S2。PM48 为主评价集（配额 n = 48，exhaustiveness = 16）；PM110 为同一靶对的更大协议敏感性面板。14 个主受体在统一化学对应 CalcRMS 表中均至少有一个保存姿态低于 2.0 Å（Figure S3；Table S2）。该结果只是搜索覆盖检查。AChE 4EY7 与 TYK2 3LXP 保存 8 个姿态，其余主受体保存 9 个姿态。EGFR 3POZ 的 top-1 为 1.019 Å，HER2 3RCD 的 top-1 为 1.947 Å，均低于 2 Å。BChE、mTOR、JAK2、PPARG 和 PPARA 的 top-1 未回到 2 Å 以内；PPARA 6LXA 的 top-3 最低为 7.848 Å。

对全部八个靶对检索了 BindingDB[16] 与 PubChem。按 \(\theta=6.0\) 独立性过滤剔除共享文献、重复结构以及与开发集 ECFP4 Tanimoto 相似度 \(\geq 0.70\) 的分子后，没有靶对满足独立外部评价准入标准（dual、A-only 和 B-only 各类至少 20 个配体且每类至少 3 个独立来源）。因此未形成外部对接集；该检索是准入筛选，不是外部验证（Figure S5；Table S8）。0/8 对通过完整外部对接门槛。


## 4. 讨论

### 4.1 候选排序问题

双靶对接评价所回答的问题取决于对照配体及其实验活性状态。dual–neither 比较检验双靶活性配体能否与两个靶点均低活性的化合物区分；dual–A-only 比较和 dual–B-only 比较则分别检验对接评分能否区分 dual 与仍保留一个靶点活性的选择性配体。这两类比较回答的是不同问题，而不是同一评价任务中负样本难度不同。

固定评分通道后，从 B-only 换成 neither 仍会更换具体分子，其骨架、大小、电荷和实验来源也可能不同。EGFR/HER2 与 JAK1/TYK2 在同一靶点评分下更换对照后仍出现明显差异，不能仅由双口袋评分的聚合方式解释。

其他靶对和评价方向的变化幅度并不一致。较大的固定评分差值既出现在供给受限的 EGFR/HER2 面板中，也出现在通过标准供给筛选的 JAK1/TYK2 中。F2/F10 两个方向均低于随机，进一步说明对接信号的方向也可以因靶对而异。对照类别选择的影响还取决于具体靶对和评价方向。

在 Zhou 等工作的基础上，[9] 本研究将两类单靶选择性状态分别构建为对应的方向性判别任务，并在固定评分通道下考察评价设定变化。

跨 DUD-E、LIT-PCBA 及相关基准设定的比较已经表明，对照组成可以改变虚拟筛选评价；[5–7,12,13,20] Kinase-Bench 则面向激酶选择性识别。[21] 本文回答的是另一问题：在两端均有实验测量的化合物中，更换实验状态比较如何改变对接判别的解释。主要推断来自靶对内固定评分比较，不是八对排行。

### 4.2 配体化学与表观判别来源

回顾性评价中的活性类别本身可能带有化学结构差异。AChE/BChE 中简单物化描述符已能区分部分实验状态。不依赖受体的 ECFP4 在骨架分组交叉验证下于多个方向携带了明显的类别信息。在 ECFP4 模型中加入对应口袋对接评分后，16 个方向的 AUROC 最大绝对变化为 0.023。较高的回顾性 AUROC 因此不能仅凭自身数值归因于三维结构互补；判断对接是否增加了判别能力，需要与不使用受体结构的基线放在相同划分下比较。

公开数据库中的化合物往往来自围绕有限先导骨架开展的药物化学研究，不同实验活性类别可能同时具有不同的骨架、取代模式和物化性质分布。部分 dual 配体也可能来自连接、融合或其他组合药效团设计，[1,2] 从而相对单靶系列系统改变分子大小、极性、柔性和局部结构基序。本文未对配体架构作显式标注，因此仅配体信号不能归到某一种设计机制。Figure 3 因此支持存在不依赖受体的类别信息，并不证明 ECFP4 在形式上优于对接。Bemis–Murcko 分组降低了相同骨架在训练和测试间直接重复的影响，但不能完全消除同一药物化学系列或文献来源带来的结构相关性。Caba 等也曾观察到，仅配体 Morgan 指纹可在回顾性结构基础虚拟筛选中携带较强分类信号。[22] 简单物化性质匹配可以减少部分低维性质差异，但不能完全控制骨架和药物化学系列之间的高维差异。

### 4.3 口袋对应性

不重新对接的评分通道交换提供了更直接的结构对应检验。主评价集中 AChE/BChE 与 EGFR/HER2 的对应口袋正差值区间排除 0。七个内部留出集未一致再现对应口袋优势。在同源靶点上，非对应口袋仍可能给出相关评分，因此该对照是不完美的特异性对照。仅凭某一评价集中较高的方向性 AUROC，仍不足以把该判别归因于相应靶点口袋。区间包含 0 并不表示受体没有信息。

### 4.4 计算实现稳健性

已有激酶交叉对接研究表明，受体构象选择可以显著影响姿态恢复和虚拟筛选性能。[14] 本研究中，PIK3CA/mTOR 在替换 PIK3CA 或 mTOR 受体结构后，方向性 AUROC 也发生明显变化。这一分析仅覆盖一个主评价靶对，只能说明该体系存在受体结构敏感性。

共晶重对接能够搜索到近天然姿态，并不表示该姿态一定获得最高评分。最低 RMSD 主要反映搜索覆盖，不能评价这些姿态能否被正确排序。更换姿态生成和评分方法后，具体 AUROC 会发生变化。对 EGFR/HER2 和 JAK1/TYK2，独立 GNINA 姿态生成下仍观察到类似的任务差异，表明这两项观察并不特异于主要 Vina 姿态生成流程。RTMScore、GNINA CNN 重评分和受体替换则显示，具体方向性 AUROC 对所用评分方法和结构实现较为敏感。独立 GNINA 用于观察评价设定差异能否在另一套姿态生成流程下出现，而不是检验对应口袋优势。

### 4.5 适用范围、证据边界与前瞻验证

在本文考察的回顾性评价条件下，对接评分应视为候选排序或筛选信号，而不是真实双靶活性或选择性的直接证据。两个口袋中均获得有利对接评分，只表示该分子在当前计算条件下对两个靶点均获得了有利预测评分。

已有前瞻性研究表明，多靶对接和生成式设计可以产生经实验验证的多靶候选物。[18,19] 这类研究评价的是从候选库中找到活性分子的能力；本文评价的是在两端均有实验测量的配体中，不同实验状态之间的回顾性判别能力。回顾性 AUROC 不能替代前瞻实验命中率；本文结果也不评价对接在候选库压缩中的总体价值。较高的回顾性 AUROC 也不能单独证明方法能稳定找到具有预期活性平衡的双靶分子。JAK1/TYK2 Top 10% 的类别组成来自固定构建的回顾性评价面板，其中 109 个配体中包含 31 个 dual，不能直接解释为商业库或前瞻筛选中的命中率。八对 Top 10% 组成见 Figure 2D 与 Table 3。

对于此类回顾性双靶基准，我们建议分别报告 dual–A-only 比较和 dual–B-only 比较两个方向，并结合仅配体基线、对应/非对应口袋比较以及必要的受体或评分敏感性分析，以判断表观 AUROC 是否具有稳定的结构解释。最终的双靶活性、活性平衡和作用机制仍需要直接实验测定。

主要限制是异质的活性数据、不完全统一的面板构建，以及对替代对接设置的部分覆盖。两端实验测量仍然少且来源异质。四状态分类要求同一配体在两个靶点上均有测定。尽管统一到 pChEMBL，IC50、Ki、Kd、EC50 和 Potency 在生物学上并不等价。取最大 pChEMBL 进一步偏向最强记录。在同一套 ChEMBL 37 记录上改用中位数聚合，产生了 Table S3 中的类别翻转，但未改变主要方向性判断。该检查只部分缓解取最强记录的偏向。所得类别是数据库来源的回顾性状态，不是同一实验体系下的生物学真值。外部评价数据的可用性见 Results 3.6。

评价面板构建在各靶对间并不完全统一。EGFR/HER2 与 PIK3CA/mTOR 从 \(\theta=6.0\) 池抽样，其余六对从严格 6.5/5.5 池抽样；骨架上限也因面板而异。因此绝对 AUROC 不用于靶对间定量排序。主要推断来自固定评分通道的靶对内比较。纳入的八个靶对用于考察不同体系中的一致性，不应视为相互独立的统计重复：JAK1 与 PPARA 各被使用两次，部分靶对属于相近激酶或核受体家族。

计算上，本研究主要依赖于每个靶点单一的代表性晶体结构，并以 AutoDock Vina 作为主要对接和评分方法。共晶重对接显示能够覆盖近天然构象，但这主要反映搜索覆盖，不能替代对评分排序能力和筛选性能的评价。独立 GNINA 姿态生成、替代受体构象和其他评分函数的分析仅覆盖了部分体系。配体准备未系统枚举质子化状态、互变异构体和构象系综。因此，不宜将本文观察到的具体数值外推为双靶对接在其他蛋白家族和更广阔化学空间中的普适表现。

提高此类评价的证据强度，需要获得同一批化合物在两个靶点上、尽量统一实验条件下的成对活性测量，并开展独立或前瞻性验证。

## 5. 结论

本研究考察实验状态定义如何影响回顾性双靶对接评价，对象为八个人源靶对，且配体在两个靶点上均有实验测量。固定同一靶点评分通道后，更换对照类别及其对应化合物集合可以改变性能判断，且该影响因靶对和方向而异。能够区分 dual 与双端低活性配体，并不保证能够区分 dual 与单靶选择性配体。

部分回顾性判别可以由不使用受体结构的配体化学获得；在当前 ECFP4 模型和骨架分组划分下，加入对接评分未表现出跨方向一致的 AUROC 提升。七个内部留出集中，对应口袋优势未一致再现。因此，仅凭较高的回顾性 AUROC，并不能证明所观察到的判别来自相应口袋特有的三维信息。对于此类回顾性双靶基准，应分别报告两个方向的选择性任务，并结合配体化学与口袋对应性解释结果。这些指标描述的是已知实验状态集合中的区分能力，不能直接等同于前瞻筛选命中率。最终的双靶活性和选择性仍需实验测定。

## 数据与软件可用性

评价面板定义、实验状态标签、受体与对接盒规格、逐配体分数表、分析脚本以及图件生成代码，均可在公开仓库 https://github.com/1280602962-debug/gwj260531 的 `Dual_Target_Docking` 目录中获取。所报告的统计分析可由已存放的分数与元数据表复现，无需重跑全部对接。排版补充材料为 Tables S1–S10。其余质控表保留在公开仓库。

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

(17) García-Ortegón, M.; Simm, G. N. C.; Tripp, A. J.; Hernández-Lobato, J. M.; Bender, A.; Bacallado, S. DOCKSTRING: Easy Molecular Docking Yields Better Benchmarks for Ligand Design. *J. Chem. Inf. Model.* **2022**, *62*, 3486–3502. DOI: 10.1021/acs.jcim.1c01334.

(18) Wu, Y.; Vigneron, S.; Braz, J.; Srinivasan, K.; Fink, E. A.; Huang, X.-P.; Xu, X.; Huebner, H.; Kim, J. Y.; Wang, J.; Pfeiffer, T.; Sakamoto, K.; Moroz, Y. S.; Radchenko, D. S.; Rodriguiz, R. M.; Irwin, J. J.; Gmeiner, P.; Billesboelle, C.; Roth, B. L.; Basbaum, A. I.; Manglik, A.; Wetsel, W. C.; Shoichet, B. K. Large Library Docking for Polypharmacology. *J. Med. Chem.* **2026**, *69*, 6210–6229. DOI: 10.1021/acs.jmedchem.5c03810.

(19) Munson, B. P.; Chen, M.; Bogosian, A.; Kreisberg, J. F.; Licon, K.; Kuenzi, B. M.; Ideker, T. De novo generation of multi-target compounds using deep generative chemistry. *Nat. Commun.* **2024**, *15*, 3636. DOI: 10.1038/s41467-024-47120-y.

(20) Gu, S.; Shen, C.; Zhang, X.; Sun, H.; Cai, H.; Luo, H.; Zhao, H.; Liu, B.; Du, H.; Zhao, Y.; Fu, C.; Zhai, S.; Deng, Y.; Liu, H.; Hou, T.; Kang, Y. Benchmarking AI-powered docking methods from the perspective of virtual screening. *Nat. Mach. Intell.* **2025**, *7*, 509–520. DOI: 10.1038/s42256-025-00993-0.

(21) Wei, T.-H.; Zhou, S.-S.; Jing, X.-L.; Liu, J.-C.; Sun, M.; Zhao, Z.-H.; Li, Q.-Q.; Wang, Z.-X.; Yang, J.; Zhou, Y.; Wang, X.; Ling, C.-X.; Ding, N.; Xue, X.; Yu, Y.-C.; Wang, X.-L.; Yin, X.-Y.; Sun, S.-L.; Cao, P.; Li, N.-G.; Shi, Z.-H. Kinase-Bench: Comprehensive Benchmarking Tools and Guidance for Achieving Selectivity in Kinase Drug Discovery. *J. Chem. Inf. Model.* **2024**, *64*, 9528–9550. DOI: 10.1021/acs.jcim.4c01830.

(22) Caba, K.; Tran-Nguyen, V.-K.; Rahman, T.; Ballester, P. J. Comprehensive machine learning boosts structure-based virtual screening for PARP1 inhibitors. *J. Cheminform.* **2024**, *16*, 40. DOI: 10.1186/s13321-024-00832-1.
