# 对接双靶识别中的基准设定与化学混淆

## 摘要

两端有利的对接分数是否构成双靶识别证据，尚未在实验定义的单靶选择性配体上得到充分检验。我们构建 DualFourClass-Bench，这是一套经策展的四对、四状态面板，含两条方向性主任务以及保守的最差方向判别摘要（`summary_min`）。在同一套冻结 AutoDock Vina 分数上，EGFR/HER2 显示出明显的 formulation 对照：Dual versus neither 的 AUROC 为 0.756，而方向性 `summary_min` 为 0.430。其他靶对没有显示同样的差距，PIK3CA/mTOR 的 Dual-versus-neither 对照则效能不足。在支架分组模型中，把 docking 加到 ECFP4 后 AUROC 的最大绝对变化为 0.020；将最大 pChEMBL 换成重复测定中位数只产生很小的靶对层变化。替代受体使 PIK3CA/mTOR 从 0.692 变为 0.486/0.505，并使 PIK3CA/PIK3CB 向相反方向变化。这些结果支持把实验定义的选择性硬负样本与混淆感知对照作为双靶对接评价的互补要求，并表明表观判别仍具有明显的靶对与受体实现依赖性。

**关键词：** 双靶对接；基准设定；选择性硬负样本；化学混淆；受体实现；虚拟筛选

## 1. 引言

### 1.1 双靶设计及传统对接基准的不足

多靶点药物设计（multitarget drug design）旨在通过单一小分子同时调控两个或多个生物学靶点，以应对复杂疾病中的通路冗余、代偿性信号以及药物耐药等问题。与传统单靶点药物相比，合理设计的多靶点配体有望通过协同调节相互关联的生物学过程获得更充分的药理效应，因此已成为多靶点药物发现和多药理学（polypharmacology）研究的重要方向。[1] 近年来，多靶点小分子的理性设计逐渐由经验性的多药理筛选，转向结合结构生物学、计算化学与生成式模型的结构导向设计。[2]

在这一过程中，分子对接（molecular docking）仍是结构基础虚拟筛选（structure-based virtual screening, SBVS）中最常用的计算工具之一：先预测配体在蛋白结合口袋中的结合构象，再用打分函数对配体–受体相互作用排序，从而在大规模化合物库中给出结构互补性的近似评价。[3][4] 因此，在双靶点药物发现中，一个自然的计算策略是分别将候选分子对接至两个靶点，并据此判断其是否具有潜在的双靶结合能力。

既有虚拟筛选研究已经表明，对接结果的解释高度依赖数据集构建、负样本定义、化学偏倚和评价指标。DUD 与 DUD-E 使用物化性质匹配 decoy，以避免表观 enrichment 退化为粗粒度配体性质分离。[5][6] LIT-PCBA 进一步采用实验 assay 标签，并系统控制已知 decoy/chemical biases，以提高虚拟筛选评价的现实性。[7] CASF-2016 则评价复合物上的 scoring、ranking、docking 与 screening power，仍然属于单复合物问题。[8] 这些资源都没有在实验标注的四状态配体空间中定义双靶方向判别。

然而，**将单靶点 docking 的评价逻辑直接扩展到双靶点任务并不充分。**

对于双靶点配体而言，任务结构发生了改变。一个严格的双靶点 benchmark 至少需要区分四种具有不同生物学含义的配体状态（四状态数据集，而不是四分类器）：同时作用于两个靶点的 **dual-active** 配体、仅作用于靶点 A 的 **A-selective** 配体、仅作用于靶点 B 的 **B-selective** 配体，以及两个靶点均缺乏足够活性的 **neither** 配体（Figure 1A）：

|  | *B*<sup>+</sup> | *B*<sup>−</sup> |
|--|:--:|:--:|
| *A*<sup>+</sup> | Dual | A-only |
| *A*<sup>−</sup> | B-only | Neither |

其中，A-only 和 B-only 不是普通负样本，而是**选择性硬负样本（selectivity hard negatives）**。它们在一个靶点上已有较强活性，却在另一靶点缺乏相应活性。计算终点因而检验 docking 能否在两个方向上将 dual-active 与对应单靶选择性配体区分开；这一判别本身不证明独立的 pocket-specific recognition 或生物学双靶活性。

Zhou、Li 与 Hou 曾在四对激酶上评价相对 noninhibitor 的 dual-target docking，并报告结构依赖性和预测 dual 中的 false positives。[9] 本文在该评价设定基础上进一步引入实验定义的 A-only/B-only 方向硬负，并直接比较不同 formulation 下的表观判别。Dual-versus-neither 是 **nonselectivity-controlled comparator**，不是 “the conventional dual-target benchmark”。

池化两个口袋分数可能掩盖较弱方向，而两端都优于参考配体只定义计算成功，并不等价于实验双靶活性。因此 benchmark 与读出必须匹配四状态生物学空间。

从公开数据构建此类面板要求同一化合物在两个靶点上均有可比较测量，并要求两个方向都有足量选择性硬负。由于 assay 类型、条件和覆盖不同，能够支持平衡四状态评价的靶对数量本身就是数据供给问题，而不是预先存在的资源。

双靶面板还继承化学混淆：若 dual-active 与选择性配体的分子性质或支架分布不同，AUROC 可以反映 label-associated ligand distributions，而不只是口袋互补性。因此需要显式的物化与 ligand-only controls。[7]

### 1.2 Docking-based 双靶设计使严格评价成为实际问题

DualDiff 与 FuseDiff 说明了这一差异的实际意义：两者均使用较差口袋分数和“生成分子在两个靶点均优于参考配体”的比例评价双靶设计。[10][11] 这些指标衡量相对参考配体的计算成功，而不是相对实验选择性配体的判别。因此生成式 docking metrics 与本文 hard-negative endpoint 是互补而非竞争基准；本文不重对接其生成分子，也不把其 reported metrics 重新解释为本文 primary endpoint 的竞争者。

### 1.3 研究目的与贡献

**本文要问的是：benchmark formulation 是否改变双靶识别的表观证据。** 我们构建实验定义的四状态面板，以针对 A-selective 与 B-selective 硬负的口袋匹配方向判别为主任务，并与 nonselectivity-controlled Dual-versus-neither comparator 比较；随后用 ligand-only、物化性质、wrong-pocket、配体池、活性聚合与受体实现对观察到的判别进行压力测试。

DualFourClass-Bench 是具有两条方向主任务的 curated benchmark：dual 对 A-only 在口袋 B 打分，dual 对 B-only 在口袋 A 打分（Figure 1B）。neither 描述完整实验空间，但不进入 primary AUROC。`summary_min` 是保守的最差方向判别摘要，不是新 docking score。公开数据审计决定多少候选靶对能够支持该构建；评价集规模是审计结果，而非预设目标。

本文贡献是评价协议与资源，而不是新 docking algorithm。它把任务设定、混淆感知评价和评价条件敏感性连接起来，检验表观判别究竟是 docking 的固定属性，还是 benchmark 条件下的条件性结果。该协议用于双靶虚拟筛选与生成设计的下游校验，不是与现有生成模型竞赛，也不是 comprehensive dual-target suite。

## 2. 方法

### 2.1 数据来源与活性数据整理

双靶评价所需的配体活性作为 **experimentally derived activity labels**，通过 ChEMBL Web API 的公开 activity 端点获取。靶对供给审计于 2026-07-23 冻结。pChEMBL 将若干经标准化的定量效力或亲和力测量（如 IC50、EC50、Ki、Kd 和 Potency）转换为近似 −log10 活性尺度，便于大规模公开数据整合。不同 assay 类型、实验条件与测定体系并不等价；本文将 pChEMBL 作为策展中的统一近似，而不解释为同一条件下可直接比较的绝对结合亲和力。

同一配体–靶标若有多条可用 pChEMBL 记录，主策展采用**最大 pChEMBL** 作为一对一代表值。由于 assay 类型、条件与实验体系并不等价，另将活性聚合敏感性作为预先指定的分析：从 ChEMBL activity 端点重拉 assay 级记录，并分别采用该批记录的最大值和中位数进行重新聚合（Table S29）。该替代聚合覆盖全部冻结评价面板，但不改变面板成员、对接参数或 Vina 分数；类别比较仍使用 θ = 6.0。由于一次缓存与 API 数据不一致，EGFR/HER2 的 API-refetched maximum aggregation（0.417）与冻结主面板（0.430）存在轻微差异，因此 A4 数值不直接替代 Table 2 的冻结主结果。任一端缺少有效 pChEMBL 的配体不进入需要双端标签的分析。

ChEMBL 结构常含盐、溶剂化物或多组分形式。对接前按连通片段拆分，并保留重原子数最多的有机片段作为计算母体。

为检验 ChEMBL 供给门槛是否仅为单一库的覆盖假象，对进入冻结评价集的靶对另做 **BindingDB / PubChem 计数核对**（零对接、不重建面板）。BindingDB 使用 REST `getLigandsByUniprots`（cutoff = 1 mM，以免截掉弱端测定）；PubChem 使用 PUG REST `protein/accession/…/concise`。类型限于 IC50/Ki/Kd/EC50；代表值取最大转换 p 活性；分类规则与 2.2 的严格供给门槛相同。配体身份分别用 BindingDB monomerid 与 PubChem CID，**不做**跨库结构合并。主比较采用去掉 `>`/`<` 截尾的**等式测定**；将不等式当作点估计只作敏感性。该核对只报告计数（Supporting Information Table S12）。

### 2.2 靶对供给审计与实验配体状态定义

为判断公开数据能否支持严格的双靶 benchmark，先对候选靶对做数据供给审计。对每一对靶标 A/B，按两端实验活性将配体定义为四种**实验状态**：

- **dual**：两端均有较强活性；
- **A-only**：仅 A 端较强，B 端较弱；
- **B-only**：仅 B 端较强，A 端较弱；
- **neither**：两端均缺乏足够活性。

A-only 与 B-only 是选择性硬负样本，不是 DUD/DUD-E 式假定 decoy。

**严格供给审计规则（construction gate，不是全部最终比较的唯一标签）。** Dual：两端 pChEMBL ≥ 6.5。A-only：A ≥ 6.5 且 B ≤ 5.5。B-only 对称。Neither：两端 ≤ 5.5。介于 5.5 与 6.5 的灰区不进入严格审计。该规则用于判断某一靶对在两个方向上是否具有足够的选择性硬负，以支持规模较均衡的面板。审计通过门槛与最终进入评价集的靶对名单见 Results 3.1；金属依赖体系（如 HDAC）按预先声明排除，不作为常规非共价对接主对象。

**正文主比较采用预先统一的 θ = 6.0 标签。** Dual：两端 ≥ θ；A-only：A ≥ θ 且 B < θ；B-only 对称；neither：两端 < θ。严格 6.5/5.5 规则用于 target-pair supply audit 与面板建造；primary analysis 的实验状态标签在 panel construction 前预先统一定义为 θ = 6.0。建造规则在抽样前按供给审计冻结，并写入 Table 1。阈值选择服务于可分析面板，不根据 docking 结果调整。作为支持性敏感性分析，在 θ ∈ {5.5, 6.5} 与严格 6.5/5.5 规则下重标四种状态并重算口袋匹配 summary_min（Table S4）。该网格不是与 Table 2 竞争的第二套主标准。样本量过小的格子在 Results 中标记 underpowered，Methods 不预判其数值。

### 2.3 DualFourClass-Bench 面板构建

**该资源是保留四种实验状态、但以两条方向判别为主任务的策展基准。** Dual、A-only、B-only 与 neither 均保留以描述实验空间；预先指定的主终点是 dual versus A-only 与 dual versus B-only。neither 不进入 primary directional AUROC。这不是四分类器 benchmark。

候选靶对按 2.2 的严格供给审计筛选。最终冻结评价集包含 PIK3CA/mTOR、AChE/BChE、PIK3CA/PIK3CB 与 EGFR/HER2。EGFR/HER2 保留为**供给受限案例**，其组成不与其余靶对按同一厚面板供给条件等价。

每个靶对从符合相应实验标签的候选池中按预先冻结的类别配额抽样，随机种子为 20260729。在抽样时结构可用的面板上施加 Bemis–Murcko 支架封顶，以降低同一化学系列过度代表：PIK3CA/mTOR（PM48）同一类别内同一支架最多 2 个分子；EGFR/HER2 最多 5 个。AChE/BChE 与 PIK3CA/PIK3CB 抽样时结构不可用，因此只采用类别配额和确定性随机顺序，不施加额外多样性约束。事后计算的 Murcko 支架随数据表报告。各面板的最终成员、状态标签、ChEMBL identifier、SMILES 与抽样脚本随冻结数据包提供；观察对接分数后不再重抽面板。

四对的建造规则并不相同。AChE/BChE 与 PIK3CA/PIK3CB 在严格 6.5/5.5 规则下抽样；EGFR/HER2 与 PIK3CA/mTOR 因严格规则下 B_only 过少而改用 θ = 6.0。因此跨对 AUROC 同时混合靶对生物学与面板构建差异（样本量、阈值、化学系列、受体），不能读成纯粹的 intrinsic docking performance。

配额与建造标签如下。AChE/BChE 与 PIK3CA/PIK3CB：严格 6.5/5.5，目标 dual / A_only / B_only / neither = 28 / 28 / 28 / 16（面板 n = 100）。EGFR/HER2：θ = 6.0 建造规则（n = 110）。PIK3CA/mTOR：θ = 6.0，主比较面板 PM48（n = 48；建造 dual / A_only / B_only / neither = 18 / 14 / 12 / 4），并在其上冻结受体与对接协议。

对接失败的配体–受体组合从该受体分数中剔除；任一端缺少可用分数的配体不进入需要两端分数的口袋匹配 AUROC，故分析用计数可低于建造定额（Table 1）。AUROC 因此是**以对接引擎能够处理的化合物为条件**的。尝试 / 成功 / 失败计数（含 AutoDock 原子类型 `B` 等化学覆盖失败）见 Table S27。

PIK3CA/mTOR 另构建扩面面板 PM110：保留 PM48 全部 48 个配体，并按严格规则追加分子，目标配额 dual / A_only / B_only / neither = 30 / 30 / 30 / 25。PM110 是 PM48 的超集，用于评价面板规模增加后点估计是否同向，不是与其他靶对独立等价的 primary benchmark，也不是独立重复实验。主文跨对比较以 PM48 为准。

配体侧外推采用一次 unused-pool holdout（2.11），因为剩余硬负供给不足以支持互不重叠平衡面板的分布。配体层有放回 bootstrap（2.8）描述固定面板内的不确定度，不是供给池重抽。

**Table 1.** DualFourClass-Bench 评价集组成与对接设置。construction labels 记录各靶对的供给/建面规则；Tables 2–3 的全部 primary AUROC 均使用统一 θ = 6.0 实验状态标签。严格 6.5/5.5 是供给与建面门槛，重标敏感性见 Table S4。

| 靶对 | 建造标签规则 | 受体 PDB (A / B) | 分辨率 (Å) | 面板 n | 分析用 n (dual / A_only / B_only) | Vina exhaustiveness |
|------|--------------|------------------|------------:|-------:|----------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | 严格 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | 严格 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

### 2.4 蛋白结构与结合位点定义

受体取自 Protein Data Bank 中含实验确定结构与小分子共晶配体的条目。冻结主分析使用：PIK3CA/mTOR，4L23 / 4JT6（共晶配体 X6K / PI-103）；AChE/BChE，4EY7 / 4BDS（E20 / THA）；PIK3CA/PIK3CB，4L23 / 2WXF（X6K / 039）；EGFR/HER2，3POZ / 3RCD（03P / TAK-285）。分辨率见表 1。

结合位点由各结构的共晶配体定义。以共晶配体重原子坐标计算轴对齐包围盒（AABB），三方向各外扩 5 Å；任一边若小于 20 Å，则将该边设为至少 20 Å。盒子中心与边长冻结于 JSON，并汇总于 Supporting Information Table S2。

受体准备时去除水分子与共晶配体，再用 Meeko 生成 PDBQT。PIK3CA、mTOR、EGFR 与 HER2：使用冻结目录中已含氢的蛋白坐标，经 `mk_prepare_receptor.py --read_pdb` 转换。AChE、BChE 与 PIK3CB：从沉积 PDB 的 ATOM/TER 记录提取蛋白（去除水与异源原子），以 `mk_prepare_receptor`（默认 alternate location A）转换。未额外运行 PDBFixer 补全缺失原子，也未用 Reduce 做独立的 pH 依赖质子化或组氨酸互变异构枚举；质子化属于冻结准备协议的一部分。主分析均为非共价小分子对接，不在盒子中把金属离子或其他辅因子当作额外可对接组分。完整命令与输入文件随公开数据包提供。

### 2.5 共晶配体重对接质量控制

正式对接前，对每个冻结受体做共晶配体重对接，以检验对接盒子、受体准备与搜索参数能否在**保留的姿态集合**中生成近似共晶构象。

每个共晶配体生成 9 个 docking poses，计算其与实验共晶构象的重原子 RMSD（对接坐标系，不做蛋白叠合）。PIK3CA/mTOR 与 EGFR/HER2 使用 meeko `REMARK SMILES IDX` 映射后，在图自同构上取最小 CalcRMS；AChE/BChE 与 PIK3CB 使用重原子匈牙利匹配（`linear_sum_assignment`）。定义

\[
\mathrm{RMSD}_{\mathrm{best9}} = \min_{i=1,\ldots,9} \mathrm{RMSD}_i.
\]

预先通过标准为 \(\mathrm{RMSD}_{\mathrm{best9}} < 2.0\) Å。

该 QC 检验的是 **pose-generation capability**：协议能否在保留的 pose ensemble 中产生近晶构象。它**不等于**要求 Vina 排名第一的 pose（mode 1）必须为近晶构象。best-of-9 QC 与 mode-1 scoring 是不同层面的评价。

若默认 exhaustiveness 下未满足预设 QC，则在不改变盒子、受体与随机种子的条件下，将搜索强度提高至预先规定的备用水平并重新 QC。主分析因此采用受体特异的冻结 exhaustiveness：PIK3CA/mTOR 为 16，其余主面板为 8。各受体的 QC 数值、mode-1 与 best-of-9 对照见 Supporting Information Table S3。

### 2.6 配体准备与分子对接

配体从冻结 ChEMBL SMILES 统一准备：去盐并保留最大有机片段，RDKit 加显式氢，ETKDGv3 生成三维构象（随机种子 20260727），MMFF 局部优化最多 200 步，再经 Meeko 默认参数转为 PDBQT。不进行系统性的质子化状态、互变异构体或构象枚举；各靶对使用同一 ligand-preparation protocol。不使用 Schrödinger LigPrep。

分子对接采用 AutoDock Vina 1.2.7，默认 `vina` 打分函数。每个配体–受体组合生成 9 个 poses，`energy_range = 3` kcal mol\(^{-1}\)，随机种子 20260727（与 ETKDG 相同）。exhaustiveness 按 Table 1 的受体特异冻结值。配体准备、盒子生成规则与打分函数在各主面板上相同；仅受体坐标、盒子数值与预先定义的 exhaustiveness 按靶标变化。完整参数见 Supporting Information Table S1。

### 2.7 替代打分通道

为检验主观察是否依赖单一打分函数，在**同一组 Vina-generated poses** 上另用 RTMScore 与 GNINA CNN 重打分。

RTMScore 使用公开权重 `rtmscore_model1`，对每个配体–受体组合的 9 个 Vina poses 分别打分，取该口袋最高 RTMScore。

GNINA 1.3.2 在 CPU 模式下做 CNN rescoring（`--cnn_scoring rescore --minimize`）。最终协议对全部 9 个 Vina poses 分别转 SDF（Open Babel）并重打分，取每端最高 CNNscore，与 RTM 的姿态覆盖对齐。仅使用 Vina mode 1 的 GNINA 结果保留为历史 sensitivity control，不是最终通道读出。

Vina 主读出是 mode 1 能量；RTM 与 GNINA 是 best-of-9 重打分。三者对 9 个姿态的聚合并不相同，因此 **不作为 head-to-head docking-engine competition**，而作为 scoring-channel sensitivity analysis。Primary endpoint 始终由 Vina 定义。

### 2.8 主终点与统计分析

#### 2.8.1 口袋匹配方向 AUROC

全文中的“dual-target recognition”指这一计算判别任务，不表示已经验证的生物学识别，也不等同于独立的口袋特异结合。

对每个靶对 A/B 计算两条二分类 AUROC。dual 对 A-only 使用**口袋 B** 的分数：

\[
\mathrm{AUC}_{D/A} = \mathrm{AUROC}(\text{dual},\;\text{A-only};\;S_B),
\]

以检验非选择性口袋 B 的分数能否把 dual-active 与已在 A 端强效的 A-only 分开。dual 对 B-only 使用口袋 A 的分数：

\[
\mathrm{AUC}_{D/B} = \mathrm{AUROC}(\text{dual},\;\text{B-only};\;S_A).
\]

dual 始终为正类。neither 不进入上述对比。

Vina 输出结合能 \(E_{\mathrm{Vina}}\)（kcal mol\(^{-1}\)，通常越负表示预测结合越强）。定义

\[
S_{\mathrm{Vina}} = -E_{\mathrm{Vina}},
\]

使所有 primary scores 遵循“越大表示预测结合越强”。RTMScore 与 GNINA CNN 分数本身已是越高越好，不再取负。

#### 2.8.2 最差方向判别摘要（`summary_min`）

最差方向判别摘要定义为：

\[
\mathrm{summary}_{\min} = \min(\mathrm{AUC}_{D/A},\;\mathrm{AUC}_{D/B}).
\]

`summary_min` 是将两条方向 AUROC 压缩为单值的保守**最差方向判别摘要**。它避免较强方向掩盖另一方向的失败，但不是新的 docking score、独立统计检验，也不表示 calibration、sensitivity、specificity 或生物学亲和力。两条方向 AUROC 的算术平均、几何平均与调和平均作为聚合敏感性报告（Table S26）；四种摘要下四对排序不变，EGFR formulation 对照的方向也不变。全文唯一主终点是统一 θ = 6.0 下的口袋匹配 Vina `summary_min`（Table 2；PIK3CA/mTOR 主面板为 PM48）。次级、敏感性、证伪与探索性终点的完整层级见 Table S16。

#### 2.8.3 物化描述符对照

用 RDKit 计算预先指定的描述符面板：重原子数（GetNumHeavyAtoms）、分子量（MolWt）、cLogP（MolLogP）与 TPSA。每个描述符按与对接分数相同的方向 AUROC 流程评价，**正文与 SI 报告全部四个**（Table 2；Table S28）。其中 AUROC 最高者记为 **best single-descriptor reference**，只是该面板上的事后最大值，不是 confirmatory competitor。为避免先选最优描述符再做正式比较的选择偏倚，docking 与该参考的配对 Δ 不以“击败 best descriptor”作为 confirmatory test（Table S19）。

#### 2.8.4 分数聚合对照

作为辅助分析，同时计算两端分数的 pooled mean、wrong-pocket assignment（定义见 2.9.1）以及 worst-pocket aggregation。它们不是 primary endpoint，只用于判断不同聚合是否改变双靶识别结论（Table S6）。

#### 2.8.5 Bootstrap 不确定度

AUROC 与 summary_min 的不确定度用配体层 bootstrap：在保持类别标签结构的条件下对配体有放回重采样，每次重算两条方向 AUROC 与 summary_min。\(B = 2000\)，随机种子 20260729，百分位数 95% CI 为 \([P_{2.5}, P_{97.5}]\)。错口袋与描述符等配对比较在**同一次**重采样上计算 \(\Delta = \mathrm{Metric}_1 - \mathrm{Metric}_2\)，得到 paired bootstrap 区间（Table S17、S19）。另报 Murcko 支架重采样区间作为对照；正文以配体层为准。置信区间作描述性不确定度；除预先定义的主终点外，不对多靶对、多对照做多重比较意义上的 confirmatory testing，也不把“CI 是否跨越 0.5”单独等同于正式显著性。

#### 2.8.6 Benchmark-formulation comparison

在同一套冻结 Vina 分数上，将 **Dual-versus-neither comparator**（实验 inactive；`vina_mean` 与 `vina_worst`）以及 Dual versus all non-duals 作为辅助对照，与方向性主终点并列。Dual-versus-neither 是本面板上的 **nonselectivity-controlled comparator**，不是声称既有双靶基准都以 Dual versus neither 为官方任务。neither 用于该对照，仍不进入 Table 2。PIK3CA/mTOR 的 neither n = 4 标记 underpowered。该比较只问：省略选择性硬负是否会改变对双靶识别的表观证据；不是第二套主终点，也不是配对显著性检验（负样本集合不同；Table 3；Table S22）。单靶式类比——口袋 A 上 (dual + A-only) 对 (B-only + neither)，以及对称的 B 对照——见 Table S22。

### 2.9 混淆、证伪与化学对照

#### 2.9.1 Wrong-pocket falsification control

将靶点 A 与 B 的分数对调，配体、受体与其余分析设置不变，重算方向 AUROC 与 summary_min。该分析是 **falsification control**，不是用来证明口袋特异的阳性对照。固定面板上 matched > wrong **不**作为 pocket-specific signal 的证据。错口袋接近或高于匹配口袋，则视为对 pocket-specific interpretation 的反证。holdout 的 point-estimate reversal 进一步说明：wrong-pocket **不是在面板迁移下可靠的通用负对照**。

#### 2.9.2 配体层混淆对照

在配体效率归一（\(S_{\mathrm{dock}}/N_{\mathrm{heavy}}\)）、效价约束（\(|\Delta\mathrm{pChEMBL}| \leq 0.5\)）和尺寸约束（\(|\Delta N_{\mathrm{heavy}}| \leq 2\)）后重算方向 AUROC。逻辑回归比较 docking alone 与 docking + heavy-atom count + TPSA（\(C=1.0\)，`max_iter=2000`）。这些分析诊断结果对尺寸、效价和极性的敏感性；缩小后的匹配子集不作为独立 confirmatory evidence（Tables S5、S19）。

#### 2.9.3 Ligand-only 与探索性结构对照

Morgan/ECFP4（半径 2，2048 bit）加逻辑回归用于估计当前面板内的 label-associated discrimination，而非建立可迁移的活性预测器。评价采用最多五折的 Bemis–Murcko scaffold `GroupKFold`；随机 `StratifiedKFold` 只作泄漏核对，logistic docking AUROC 与 Table 2 rank AUROC 分开解释（Tables S20、S24）。最近邻 Tanimoto 子集只作诊断：T ≥ 0.7 为空，现有更低阈值子集样本量有限（Table S23）。另以 4.0 Å 内配体重原子数作为粗粒度 scoring-independent contact count，诊断 wrong-pocket 中的 size/burial contribution（Table S11）；全链序列一致性仅作为探索性背景描述符（Table S7）。两者都不是残基级 PLIF 或口袋相似度指标。

### 2.10 支持性单靶富集参照

支持性的 PIK3CA/mTOR 单靶 active-versus-weak 分析使用 pChEMBL ≥ 6.5 对 ≤ 5.5，并按 MW、cLogP 与 TPSA 匹配；配体准备与 docking 同主面板。AUROC、EF1% 和 EF5% 只在 Supporting Information 中作为背景报告，不替代方向性双靶终点。

### 2.11 未使用配体池 holdout

为检验结论是否依赖于冻结面板的具体成员，从严格标签池中排除已用于主面板与 PM110 的 ChEMBL 条目，在剩余 unused pool 中构建 **unused-pool, panel-external holdout**。它不是跨数据库或跨实验体系的 external validation：配体仍来自同一 ChEMBL 抓取批次、同一靶对与同一标签规则。

Holdout 只在 unused-pool 配额足以按 dual / A-only / B-only 各抽 20 个配体的靶对上构建。预先冻结为 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB；EGFR/HER2 不具备同等未用池抽样条件，记为 not eligible，不补做不等价抽样。PIK3CA/mTOR 排除的是 PM110 超集，从而覆盖 PM48。抽样种子 `HOLDOUT_SEED = 20260731`（区别于建面种子），同一 Murcko 支架在每一状态类中最多 3 个成员。清单先冻结、后对接。

Holdout 不参与主面板构建、对接协议调整或 primary endpoint 选择。受体、盒子、配体准备、exhaustiveness、打分与统计与主 benchmark 相同，并使用同一 `summary_min` 与配体层 bootstrap。未能产生 Vina 分数的配体–受体组合按 2.3 从需要该分数的分析中剔除。同一 holdout 配体上并列计算描述符对照；错口袋、效价/尺寸匹配与 contact count 按 2.9 在 holdout 上重算（Table S8、S13）。效价/尺寸匹配诊断不改写 Table S8 的主 holdout 数字。

### 2.12 受体结构敏感性分析

为评价 benchmark 结论对受体结构选择的敏感性，另选满足以下**预先声明**条件的替代晶体：（i）polymer entity 与目标蛋白真实对应，排除嵌合体或非目标同源骨架；（ii）含 ATP 位点或目标结合位点的小分子共晶；（iii）分辨率可接受；（iv）通过与 2.5 相同的 cognate redocking QC。实际进入对接的替代结构为 PIK3CA 4JPS、5DXT 与 mTOR 4JSX。该分析是 **receptor-structure sensitivity analysis**（receptor-realization effect），不是稳健性检验，也不是用来证明某一晶体“更正确”，更不是把 PIK3CA/mTOR 预设为结构不变的 positive case。目的是量化双靶判别终点对受体实现对的敏感性，而不是挑选更优受体结构。

替换采用**单口袋**设计。在 PIK3CA/mTOR（PM48）上，4JPS/5DXT 替换口袋 A，口袋 B 仍用冻结 4JT6 分数；4JSX 替换口袋 B，口袋 A 仍用冻结 4L23 分数；exhaustiveness = 16，与 PM48 主面板一致。在 PIK3CA/PIK3CB 上，同一套已准备的 4JPS/5DXT 替换口袋 A，口袋 B 仍用冻结 2WXF 分数；exhaustiveness = 8，与该主面板一致。新盒子按该替代晶体自身共晶配体、以 2.4 的同一 AABB 规则生成。配体准备、随机种子（20260727）、打分函数与 primary endpoint 与相应主分析相同。未能产生 Vina 分数的作业按 2.3 剔除；attempted / successful / failed 计数与换晶表一并报告。

作为探索性、零新对接的几何对照，在已冻结晶体坐标上做刚体叠合：Biopython `PDBParser` 提取最长链 Cα，按残基编号与残基名精确匹配，`Superimposer` 一次 Kabsch 拟合得全域 RMSD；口袋残基由参考结构共晶配体重原子 ≤5 Å 界定，在**同一变换**下计算口袋局域 RMSD，不做二次局部拟合。再将替代结构共晶配体按同一变换投影，计算与参考共晶配体质心的距离。不同结构匹配的 Cα 数目可以不同，全域 RMSD 因此不是等覆盖比较。本对照仅含有限数目的替代晶体，不预设 Cα RMSD 能够定量解释 AUROC 变化（Table S10）。

### 2.13 软件与数据可用性

计算在 Python 3 环境下完成。主要软件：RDKit 2026.3.1、meeko 0.7.1、AutoDock Vina 1.2.7、GNINA 1.3.2、RTMScore（`rtmscore_model1`）；Vina 姿态转 SDF 使用 Open Babel。刚体叠合与全链序列比对使用 Biopython（`PDBParser`、`Superimposer`、`PairwiseAligner`）。AUROC、逻辑回归与交叉验证使用 NumPy、SciPy、scikit-learn 与 pandas（版本见公开复现环境）。评价面板、对接分数、分析脚本与完整参数表已在公开仓库提供，见 Data and Software Availability。

## 3. 结果

### 3.1 实验数据供给限制了严格双靶基准的构建

为确定公开生物活性数据是否能够支持严格的双靶点识别评测，我们首先对 49 对有 ChEMBL 缓存的候选靶标进行供给审计（Figure 2）。双靶对接评测需要四种实验标签状态：dual、A-only、B-only 与 neither（Figure 1A）。基准是四状态数据集；预先指定的主终点是两条方向 pairwise 判别（dual 对 A-only、dual 对 B-only）。我们将一端达到活性阈值、对端明确低活性的配体定义为方向性选择性硬负样本，用于检验对接分数能否同时压住两条单靶臂。

在严格标签规则下（dual：两端 pChEMBL ≥ 6.5；选择性类：活性端 ≥ 6.5 且对端 ≤ 5.5），尽管候选靶对数量较多，能够同时提供足量 A-only 与 B-only 硬负样本的靶对十分有限。两端严格硬负均不少于 50 的厚面板条件仅有 4 对满足。排除不适合作为常规小分子对接评测对象的金属依赖 HDAC1/HDAC6 后，PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 构成三个规模相对充足的靶对；EGFR/HER2 仅有 7 个严格 B-only 配体，因此被保留为供给受限案例，而不是与前三对等价的厚面板（Table 1）。

这一供给限制并非 ChEMBL 单一数据库特有的计数现象。对最终四对靶标进行 BindingDB 与 PubChem 的零对接计数核对后（Supporting Information Table S12），与 pChEMBL 更接近的 `equal_only` 规则下，前三对的 min HN 分别为 BindingDB 76 / 92 / 58、PubChem 86 / 97 / 61（ChEMBL 缓存为 80 / 78 / 56），仍全部 ≥ 50。EGFR/HER2 虽可在其他数据库中达到约 30 个 B 端硬负样本，仍不足以满足 ≥ 50 的厚面板标准。将不等式活性记录作为点估计（`as_is`）会显著增加 EGFR/HER2 的表观供给（BindingDB min HN 升至 85），但 92 个 as-is B_only 中有 49 个在 EGFR 端只有 `>` 记录；这一处理改变了“两端具有等式定量测定”的标签定义，因此未用于冻结 benchmark。PubChem 与 BindingDB 计数接近，符合沉积重叠，不作两次独立普查。

因此，最终 benchmark 的规模并非根据对接表现事后筛选，而主要由公开实验数据中方向性选择性硬负样本的可获得性所约束。DualFourClass-Bench 是一套受供给约束、但由实验标签锚定的评价集（a constrained but experimentally grounded benchmark），不是覆盖全部双靶任务的完整抽样，也不是 comprehensive benchmark suite。严格 6.5/5.5 规则用于量化供给并记录面板构建，而 θ = 6.0 定义全部 primary AUROC 的实验状态标签（Methods 2.2）。后续结果依次考察 benchmark formulation（3.2）、配体层化学 baseline（3.3）、评价条件敏感性（3.4）与证伪对照（3.5）。

### 3.2 基准 formulation 改变了表观双靶判别

在冻结的四对靶标上，采用统一 θ = 6.0 标签规则和口袋匹配方向 AUROC 对 Vina docking scores 进行评价（Figure 1B；Methods 2.8）。分数定义为 \(S=-E_{\mathrm{Vina}}\)（越大越好），dual 为正类。预先指定的最差方向判别摘要为 `summary_min`（两条方向 AUROC 的较小值），使较强方向不能掩盖较弱方向。算术平均、几何平均与调和平均只作聚合敏感性；四种摘要下四对排序不变（Table S26）。AChE/BChE 与 PIK3CA/PIK3CB 在建造时使用更严格的 6.5/5.5 规则，但在本数据上 θ = 6.0 给出完全相同的配体分类与 AUROC（Table S4）；EGFR/HER2 与 PIK3CA/mTOR 对阈值更敏感，严格规则下 B_only 过少并标记 underpowered，故严格规则只作支持性敏感性分析，不作第二套主标准。整张阈值网格内排序趋势保持一致（Figure S1A）。

这四个 `summary_min` **不是**可互换的 intrinsic docking performance。AChE/BChE 与 PIK3CA/PIK3CB 在严格供给规则下建面；EGFR/HER2 与 PIK3CA/mTOR 使用 θ = 6.0；面板还在 n、化学系列与受体上不同。跨对差异同时混合这些构建因素与靶对生物学。

EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 的方向性 summary_min 分别为 0.430、0.606、0.500 和 0.692（Table 2；Figure 4A；Figure S4）。不同靶对的主要限制来自不同的弱臂：EGFR/HER2 的 dual-versus-B-only AUROC 为 0.430，PIK3CA/PIK3CB 为 0.500；PIK3CA/mTOR 两个方向分别为 0.714 和 0.692（Figure 4A）。相对池化协议，口袋匹配抬高了点估计但未改变排序（Table S6）。

同一套冻结分数再按 Dual-versus-neither comparator 以及 Dual versus all non-duals 计分（Table 3；Figure 3）。Dual-versus-neither 是本面板上的 **nonselectivity-controlled comparator**（实验 inactive；`vina_mean`），不是声称既有双靶基准都以 Dual versus neither 为官方任务。两套 AUROC 使用不同负样本，是 **descriptive formulation contrast**，不是配对显著性检验。

EGFR/HER2 提供了最清晰的 formulation 例子。Dual versus neither 的 AUROC 为 0.756 [0.562, 0.920]（n_neg = 12），而方向性 summary_min 仍为 0.430 [0.284, 0.576]。Dual versus all non-duals 降至 0.551 [0.443, 0.666]，说明额外难度来自选择性配体。在 110 个 EGFR/HER2 配体的混合库中按 `vina_mean` 取 Top-10：1 个 dual、5 个 A-only、4 个 B-only、0 个 neither（EF10 = 0.393；hard-negative fraction = 0.90）；EF5 也低于随机（Table S25）。因此 Dual-versus-neither 读出在该对上会给出有利的双靶判别印象，而方向性任务与筛选向 Top-10 都优先富集选择性配体。

该 formulation gap **依赖靶对**，不是四对 overestimation 定律。AChE/BChE 与 PIK3CA/PIK3CB 的 Dual-versus-neither 增量很小（0.649 与 0.559），区间与方向性臂重叠。PIK3CA/mTOR Dual versus neither 因 neither n = 4 而 underpowered，不作反向效应解释；该对 Dual versus all non-duals 为 0.674，接近 summary_min 0.692。因此 PIK3CA/mTOR 作为 **conditional directional signal**，而不是全文中心成功案例（Results 3.4）。

**Table 2.** 冻结 K = 4 评价集上的口袋匹配方向 AUROC（Vina，统一 θ = 6.0），并列出四个预先指定描述符的 `summary_min`。`n` 是要求两端均有分数后，实际进入 primary AUROC 的 dual / A_only / B_only 类别样本量；neither 不进入这些 AUROC。建造面板规模与两端对接覆盖分别见 Table 1 和 Table S27。最高描述符是 best single-descriptor reference，不是 confirmatory competitor。错口袋与配体效率见 Table S6；描述符双臂见 Table S28。

| 靶对 | n (dual / A_only / B_only) | dual 对 A_only（口袋 B） | dual 对 B_only（口袋 A） | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|------:|---:|------:|-----:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.284, 0.576] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.440, 0.740] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.347, 0.648] | 0.622 | 0.620 | 0.595 | 0.418 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.464, 0.802] | 0.463 | 0.448 | 0.310 | 0.260 |

**Table 3.** 同一套 Vina 分数在 Dual-versus-neither 与方向性 formulation 下的 AUROC（统一 θ = 6.0）。Dual-versus-neither 使用实验 inactive（`vina_mean`）；Dual versus all non-duals 把 A-only、B-only 与 neither 都计为负类。方向性 CI 来自 Table 2。该对比是描述性的；负样本集合不同。PIK3CA/mTOR Dual versus neither underpowered（n_neg = 4）。

| 靶对 | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|------|--------------------------------:|------------------------------:|----------:|----------------------:|
| EGFR/HER2 | 0.430 [0.284, 0.576] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.440, 0.740] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/PIK3CB | 0.500 [0.347, 0.648] | 0.559 [0.373, 0.746] | 16 | 0.556 [0.437, 0.672] |
| PIK3CA/mTOR | 0.692 [0.464, 0.802] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |

四个靶对的 `summary_min` 95% bootstrap CI 均包含 0.5；因此在本研究的样本量和不确定度下，没有一个靶对获得排除随机水平的明确证据。在主分析冻结受体协议下，PIK3CA/mTOR 的点估计最高（0.692；95% CI 0.464–0.802），但相对 best single-descriptor reference 的配对差值 CI 也包含 0（Table S19），且这一点估计优势对受体实现并非不变（Results 3.4）。AChE/BChE（0.606）低于 TPSA（0.733）；EGFR/HER2（0.430）与 PIK3CA/PIK3CB（0.500）也未显示超过相应描述符参考的明确优势。

对接覆盖并不完整。主面板两端均得分：EGFR/HER2 110/110，AChE/BChE 95/100，PIK3CA/PIK3CB 99/100，PIK3CA/mTOR 48/48（Table S27）。一个 A-only 配体因计算超时而持续无法完成 PIK3CA 对接，因此从需要该分数的分析中剔除；详细标识符与超时设置见 Tables S27 和 S30。AUROC 因此以 AutoDock Vina 能够处理的化合物为条件。采用 RTMScore 或 GNINA 作为替代 scoring channel 未改变总体排序。GNINA 在统一 best-of-9 pose coverage 后，EGFR/HER2、AChE/BChE 与 PIK3CA/mTOR 的口袋匹配 summary_min 仍不超过同面板 Vina；PIK3CA/PIK3CB 上 GNINA best-of-9 为 0.533、Vina 为 0.500，二者均近随机且区间重叠（Table S14–S15；Figure S1B）。GNINA 仍只是单一 CNN 通道对照。协议通过了 cognate pose-generation QC；该 QC 不是 screening-performance validation。

### 3.3 配体性质与化学型解释了相当一部分表观信号

为判断 docking discrimination 是否超越简单的 ligand-level signal，我们首先将口袋匹配对接与四种预先定义的物化性质进行比较（Figure 4B；Table 2）。相对于每个靶对的 **best single-descriptor reference**，docking summary_min 的 paired difference 在 EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 中分别为 −0.052、−0.128、−0.122 和 +0.229；四个 95% confidence intervals 均包含 0（Table S19；Figure S3C）。由此可见，即使 PIK3CA/mTOR 的点估计表现出最大的正向差异，现有样本仍不足以将其与 ligand-property reference 明确区分。该对照使用口袋匹配 summary_min，不是 pooled `vina_mean` 门控（EGFR/HER2 的 `vina_mean` 为 0.2824，≠ Table 2 的 0.4297）。

AChE/BChE 提供了一个较为直接的混淆案例。dual 配体平均 TPSA 约为 75，而选择性硬负配体约为 51（Figure 4C）；TPSA 单独获得约 0.769 的 AUROC，高于相同比较下的 Vina（约 0.56）。进一步加入 heavy-atom count 和 TPSA 后，dual-versus-B-only AUROC 从 0.606 增至 0.807，而 docking score 的 OR 仅约为 1.18（Figure 7C）。该结果表明，该方向上的部分表观 docking discrimination 可由配体物化信息解释，因此不能将其直接归因于独立的 pocket-specific recognition。

PIK3CA/mTOR 的情况有所不同。加入 heavy-atom count 和 TPSA 后，AUROC 的变化约为 +0.07 至 +0.11，docking score 的 OR 约为 2.19 和 3.08，提示该靶对可能存在一定 residual pocket-related signal；然而，与 descriptor 的 paired difference 置信区间仍包含 0，因此这一残余信号不能被视为已确证的独立优势。配体效率归一后，仅 PIK3CA/mTOR 仍高于重原子数基线（0.657 对 0.463）。

二维化学结构 baseline 进一步说明了这一问题（Figure 7A）。该模型用于估计当前面板中仅依赖配体结构即可获得的 label-associated discrimination，而不是建立可迁移的双靶活性预测器。ECFP4 + logistic regression 在 Bemis–Murcko scaffold GroupKFold 下多个方向获得约 0.78–0.91 的 fold AUROC，明显高于部分对应 docking contrasts，例如 EGFR/HER2 dual-versus-B-only 中 ECFP4 AUROC 为 0.85，而 docking AUROC 仅为 0.43。该结果只说明同一 Murcko 支架不跨训练/测试折时判别仍可保持，**不是** target-external generalization。PIK3CA/mTOR 上 \(n_{\mathrm{scaffolds}} \approx n\)，该折接近 leave-one-scaffold。同一设定下随机 `StratifiedKFold` 相对支架折的平均差为 +0.011（八个方向对比；Table S20；Figure S3D）。dual/selective 标签与 chemotype 存在系统性关联，因此单独观察 docking score 的 AUROC 并不足以证明其识别来源于 pocket-specific physical interactions。

在当前支架分组基准下，加入口袋匹配对接分数后的最大绝对变化为 0.020（未四舍五入值：PIK3CA/mTOR dual versus A-only 的 −0.0198），且若干方向为负（Table S24）。因此 docking 在该分析中只产生很小的 CV AUROC 增量改善；这不等于 docking 不提供任何额外信息。logistic 模型较简单，K = 4，没有 nested model comparison，且 logistic docking AUROC 不是 Table 2 的 rank AUROC。由于足够相似的硬负样本稀缺，chemical-similarity 分析是诊断而非正式 chemically matched control。T ≥ 0.3 时，PIK3CA/PIK3CB dual versus A-only 从 0.691 降至 0.503（n_neg = 11），远缘硬负（T < 0.3）为 0.819；T ≥ 0.4/0.5 的格子常为 n_neg ≤ 7，T ≥ 0.7 为空（Table S23）。

效价匹配或尺寸匹配子集上，EGFR/HER2 与 PIK3CA/PIK3CB 的 dual 对 B_only 仍偏弱或接近随机（约 0.45–0.52）；PIK3CA/mTOR 的排序趋势保持一致，但各臂 n 常低于 15、区间较宽（Table S5；Figure 7D）。全部四个描述符见图 7B，均不作 confirmatory competitor。

### 3.4 评价条件敏感性：活性聚合、配体面板与受体实现

主标签使用可用 pChEMBL 的最大值。为检验活性聚合规则对结果的敏感性，我们重新从 ChEMBL activity 端点获取冻结面板配体的 assay-level pChEMBL，并分别采用该批记录的最大值和中位数进行重新聚合。该分析用于比较两种聚合规则，而不改变冻结面板成员、对接参数或主分析 Vina 分数。由于一次缓存与 API 数据不一致（`EH120_060` / CHEMBL24828），EGFR/HER2 的 API-refetched maximum aggregation（0.417）与冻结主面板（0.430）存在轻微差异，因此 A4 数值不直接替代 Table 2 的冻结主结果。中位数聚合改变了 7/110、1/95、1/99 和 0/48 个配体状态分配，标签一致率分别为 93.6%、98.9%、99.0% 和 100%（Table S29）。尽管标签变化有限，靶对 `summary_min` 估计值仅适度变化（0.417→0.424、0.606→0.629、0.500→0.500 和 0.692→0.692），表明主要基准解释对该聚合选择大体不敏感。Assay 间异质性仍然存在，因为 pChEMBL 值并非 assay-equivalent。PIK3CA/mTOR 在 API-max 和 API-median 聚合下均保持 0.692 不变。

为判断 PIK3CA/mTOR 的较高 summary_min 是否仅由特定 panel 构成或 docking 搜索参数造成，我们进行了 ligand-panel 和 protocol-level sensitivity analyses（Figure S5）。将 exhaustiveness 从 16 降至 8 后，summary_min 从 0.692 降至 0.660，变化约 0.03，明显小于不同 target pairs 之间的性能差异（Figure S1D）。

在包含 PM48 全部配体并扩展至实际 n = 115 的 PM110 面板中（分析用 dual / A_only / B_only 各 30），Vina summary_min 为 0.648 [0.51, 0.76]，相比 PM48 的 0.692 下降约 0.04，但排序趋势保持一致（Figure S1C）。这一同向结果提示点估计并非完全由 PM48 的具体成员构成所决定，但 PM110 是嵌套 stability check，不是独立验证集。同面板 RTMScore 为 0.576；GNINA best-of-9 为 0.613 [0.46, 0.74]，PM48 同口径为 0.655 [0.43, 0.81]，仍不高于同面板 Vina。

在未参与主面板构建和协议调优的 unused-pool holdout 中（每对 20 / 20 / 20，种子 20260731；EGFR/HER2 不具备同等配额，记为 not eligible），PIK3CA/mTOR 的 summary_min 进一步达到 0.765 [0.603, 0.891]，高于主面板的 0.692；AChE/BChE 为 0.618 [0.422, 0.759]，与主面板接近但 confidence interval 跨越 0.5；PIK3CA/PIK3CB 则下降至 0.425 [0.241, 0.618]（Table S8 / Table S16）。PIK3CA/PIK3CB holdout 尝试 60 个配体，59 个两端得分；HOAP_028 因 AutoDock 原子类型 `B` 不支持（含硼）而两端失败（Table S27）。AChE 与 PIK3CA/mTOR holdout 为 60/60 成功。硼失败是引擎化学覆盖限制，不是 silent missingness；AUROC 以可处理化合物为条件。该 holdout 共享同一 ChEMBL 抓取批次，不能读成跨数据库独立验证；其作用是支持所观察信号在未参与建面配体池中的持续性。

因而，PIK3CA/mTOR 的方向性趋势在同一 ChEMBL 体系的未参与建面配体中仍可观察，而 PIK3CA/PIK3CB 的趋势未能保持。在当前评价中，这些点估计仍依赖 target-pair context，不能视为可跨靶对迁移的属性。

尽管 PIK3CA/mTOR 在 ligand-panel sensitivity analysis 中保持了方向性信号，我们进一步测试方向性判别是否依赖于特定 receptor realization：一端受体冻结，只替换另一端（Figure 5；Table S9；Table S30）。三个替代晶体结构均通过 cognate redocking QC，best-of-9 RMSD 分别为 0.607 Å（4JPS）、0.624 Å（5DXT）和 0.515 Å（4JSX）；嵌合体 3T8M 已排除。

在 PIK3CA/mTOR 上，当 PIK3CA 4L23 替换为 4JPS 或 5DXT、mTOR 4JT6 保持不变时，PM48 的 summary_min 分别由 0.692 降至 0.486 [0.259, 0.692] 和 0.505 [0.292, 0.696]（Figure 5A）。变化主要发生在依赖替代 PIK3CA 结构的 D/B direction，而依赖原始 mTOR 结构的 D/A direction 保持 0.714。将 mTOR 4JT6 替换为 4JSX 后 summary_min 为 0.639 [0.418, 0.776]。mTOR 端换晶后点估计仍高于 0.5，但 95% CI 包含 0.5。

同一套 PIK3CA 晶体再用于 PIK3CA/PIK3CB 面板，2WXF 分数保持冻结（exhaustiveness = 8，与主面板一致；Figure 5B）。替换后 summary_min **上升**：由 0.500 至 0.691 [0.516, 0.779]（4JPS）和 0.685 [0.506, 0.768]（5DXT）。仍使用冻结 2WXF 的 dual versus A-only 保持 0.691；使用替代 PIK3CA 分数的 dual versus B-only 由 0.500 升至 0.707 和 0.685。弱臂由原始 PIK3CA 结构上的 D/B 切换为 4JPS 后的冻结 PIK3CB 臂；5DXT 后两臂接近平衡。三种 PIK3CA 条件使用相同的 99 个已打分配体，因为一个 A-only 配体持续超时（Table S30）。

因此，同一 PIK3CA 扰动在两个靶对上方向相反：受体实现按靶对提高或降低表观判别点估计。这是评价条件效应，不是稳健性证明或单向 collapse。单口袋设计避免同时替换两个结构。两例均共享 PIK3CA，因此该格局不是 K = 4 上的普遍定律。

Cα structural comparison 进一步显示，5DXT 与 4L23 的口袋局域 Cα RMSD 仅为 0.343 Å，但 PIK3CA/mTOR 的 summary_min 仍降至 0.505，说明简单的 backbone similarity 并不足以保持判别（Table S10）。这批 PIK3CA 沉积结构彼此的整链 Cα RMSD（1.44–1.49 Å）大于这批 mTOR 沉积结构彼此的差异（0.45 Å），与 PIK3CA/mTOR 上 PIK3CA 端变动更大的方向一致，但不能定量解释 PIK3CA/PIK3CB 上的相反位移：5DXT 仅匹配 862 个 Cα，少于 4JPS 的 982 个；替代结构各仅 1–2 个。共晶配体质心距离 2.1–2.6 Å 只说明对接的仍是同一大类 ATP 竞争位点。姿态生成 QC 通过，并不等于 screening discrimination 可迁移。协议通过了 cognate pose-generation QC；它不是 virtual-screening validation。

### 3.5 错口袋证伪暴露出面板迁移下的未解决失败模式

在主面板中，pocket-matched summary_min 均高于 wrong-pocket control，四对的 matched-minus-wrong differences 分别为 0.170、0.161、0.151 和 0.090；其中 EGFR/HER2 和 AChE/BChE 的差异置信区间排除 0，PIK3CA/PIK3CB 与 PIK3CA/mTOR 的区间包含 0（Table S6；Table S17；Figure 6A；Figure S3A）。错口袋对照的 summary_min 分别为 0.260、0.444、0.349 与 0.602。主面板上 matched > wrong **不**作为 pocket-specific signal 的证据。

然而，unused-pool holdout 中的点估计关系发生反转（Figure 6B）。PIK3CA/mTOR、AChE/BChE 和 PIK3CA/PIK3CB 的 wrong-pocket summary_min 分别为 0.788、0.643 和 0.520，而 matched-pocket 分别为 0.765、0.618 和 0.425。相应的 matched-minus-wrong point differences 均为负（−0.023 / −0.025 / −0.095），但其 95% confidence intervals 均包含 0（Table S17；Figure S3B）。这一 point-estimate reversal 并不是统计上已经解析的确定性反转。因此 wrong-pocket **不是在面板迁移下可靠的通用负对照**。

为检验这一反转是否可以由 holdout 中的配体效价或分子大小差异解释，进一步进行 potency- and size-matched comparisons（Figure 6C；Table S13）。wrong-pocket ≥ matched-pocket 的关系仍未翻转（效价匹配后：AChE/BChE 0.642 对 0.593，n_min = 18；PIK3CA/PIK3CB 0.562 对 0.363，n_min = 11；PIK3CA/mTOR 0.734 对 0.715，n_min = 12）。holdout 相对主面板确有抽样偏移——最明显的是 PIK3CA/mTOR：holdout dual / A_only 的 pA 均值比主面板低约 1.1–1.3，B_only 的 pB 低约 1.8——但匹配后悖论仍在。

scoring-independent contact_count 在 B direction 获得 0.698–0.714 的 AUROC，提示该方向存在 ligand size/burial contribution，但这一粗粒度 surrogate 不能解释观察到的 Vina wrong-pocket discrimination 的幅度（Figure 6D；Table S11）。例如 PIK3CA/mTOR 中 Vina wrong-pocket summary_min 为 0.788，而 contact_count 的较弱一臂为 0.552。A 臂上 dual 与 A_only 尺寸差较小，contact_count AUROC 更接近随机（0.552–0.622）。

因此，holdout 中 wrong-pocket reversal 应被视为当前 benchmark 暴露出的 unresolved failure mode，而不是可以由单一尺寸或效价因素解释的现象。

### 3.6 结构背景只提供探索性线索

全链序列一致性只提供探索性结构背景（Table S7）。四对靶标的排序不符合简单的“靶点越相似越难区分”规则，但 n = 4，且全链一致性不是口袋相似度指标。PIK3CA/mTOR 的代表性姿态诊断还显示，选择性硬负可在两个口袋形成看似合理的 ATP-like pose，而 alternative rescoring 可能偏好远离共晶位的姿态。详细案例保留在 Supporting Information；它们用于说明 failure mode，而不是残基级机制，pose-generation QC 也不是 screening validation。

## 4. 讨论

### 4.1 基准 formulation 改变了双靶对接的证据标准

本研究的首要发现不是某一种 docking scoring function 在双靶任务上取得了最高性能，而是 **benchmark formulation 可以改变“双靶对接成功”的表观含义**。传统 docking benchmark 通常将活性配体与 decoy 进行区分，而本文要求模型同时区分 dual-active 配体与两个方向上的 single-target selective hard negatives。后者在一个靶点上具有较强实验活性，因此不能被简单视为普通 decoy。同一套分数上的 Dual-versus-neither 是 **nonselectivity-controlled comparator**，不是“the conventional dual-target benchmark”。供给审计显示，公开生物活性数据很难同时提供两个方向上足够数量的这类配体；49 个候选靶对中只有少数能够满足严格厚面板要求（Results 3.1；Figure 2）。

相对于 Zhou 等采用的 dual-target evaluation setting，[9] 本文进一步引入实验定义的 A-only/B-only 方向硬负，并直接比较不同 formulation 下的表观判别。EGFR/HER2 上，Dual versus neither 为 0.756，方向性 `summary_min` 为 0.430，混合库 Top-10 中 9/10 为选择性配体。**这是依赖靶对的清晰例子，不是四对定律或配对显著性检验：** AChE/BChE 与 PIK3CA/PIK3CB 只有小且区间重叠的增量，PIK3CA/mTOR comparator 则效能不足。

这一数据限制本身具有方法学意义。DUD、DUD-E 与 LIT-PCBA 已经表明，decoy construction、chemical bias 和真实 assay 标签会显著改变虚拟筛选的性能判断。[5–7] 简单方法或不恰当的 unbiasing 也可以通过学习配体分布而高估 structure-based virtual screening。[12] 近期基于 bioassay-derived data 的评价则进一步强调，与人工构造的 ligand/decoy 集相比，真实 assay-derived benchmarks 可以揭示模型在更接近实际筛选环境中的局限。[13] DualFourClass-Bench 并未使用这些单靶数据集，也不评价 DiffDock-Pocket；它把同一关切延伸到双靶任务：评价结论取决于硬负样本如何被实验定义，而不是取决于候选靶对清单有多长。

因此，DualFourClass-Bench 的主要价值并不在于提供一个规模很大的数据集，而在于将双靶识别问题转化为一个实验标签驱动的 hard-negative discrimination task，并显式要求两个方向同时成立。该资源是 curated four-pair panel + evaluation protocol，不是 comprehensive dual-target suite。

错口袋对照属于同一证据标准下尚未解决的面板外失败模式。主面板上 pocket-matched 高于 wrong-pocket；unused-pool holdout 上点估计反转，但配对区间仍包含零（Results 3.5；Figure 6）。因此 wrong-pocket 不是在面板迁移下可靠的通用负对照；主面板 matched > wrong 也不作为口袋特异证明。

### 4.2 配体化学基线揭示了表观对接判别中的混淆来源

双靶点 docking 的困难并不能简单理解为两个单靶 docking 任务的性能相加。四个 `summary_min` bootstrap CI 均包含 0.5，因此在当前样本量下没有一个靶对获得排除随机水平的明确证据。在主分析冻结受体协议下，PIK3CA/mTOR 的点估计最高（0.692），但相对 best single-descriptor reference 的配对差值 CI 也包含 0（Results 3.2；Table S19）。因此该对是 **conditional directional signal**，不是可推广的成功案例。

AChE/BChE 上 TPSA 单独即可获得高于 docking 的 discrimination，而 ECFP4 scaffold-grouped baseline 在多个方向上进一步超过 docking（Results 3.3）。这一强 ligand-only baseline 表明，实验标签与无需受体信息即可利用的化学空间差异相关；因此 docking performance 必须相对于 ligand-only baselines 解释，而不能孤立阅读。在当前支架分组基准下，把 docking 加到 ECFP4 后只产生很小的 CV AUROC 增量改善（最大绝对变化 0.020），且若干方向为负。这一面板内结果不等于 docking 不提供任何额外信息或口袋特异信息。

化学型约束硬负给出同一方向的证据。T ≥ 0.7 的匹配子集为空。T ≥ 0.3 时，未匹配时最强的一臂（PIK3CA/PIK3CB dual versus A-only，0.691）降至 0.503（n_neg = 11），而远缘硬负（T < 0.3）升至 0.819。T ≥ 0.3 是 similarity-constrained subset，不是 chemically matched analogue set。

这一结果与近年来对虚拟筛选 benchmark 中化学偏倚的关注一致：简单模型或不恰当构造的 decoys 可以通过学习 ligand distribution 而获得看似优异的 performance。[7,12] 如果不设置 A-only/B-only hard negatives 以及 ligand-property/chemical baselines，一个看似优秀的 dual-target docking result 可能实际上只是识别了与 dual label 相关的分子属性。

### 4.3 受体实现是评价条件的另一维度

PIK3CA/mTOR 提供了本研究中最值得进一步研究、但也最需要谨慎解释的案例。其主面板 summary_min 为 0.692，PM110 扩展面板为 0.648，unused-pool holdout 为 0.765。三者同向提示结果并非完全由 PM48 的具体成员构成所决定，但这些检查仍处于同一 ChEMBL 数据体系，只是稳定性分析而非独立外部验证。与此同时，该信号并非 receptor-invariant：替换 PIK3CA 后 summary_min 降至 0.486 和 0.505，而替换 mTOR 后为 0.639（Results 3.4；Figure 5A）。更准确的结论是在特定评价条件下观察到有限的计算判别信号，而不是可靠的生物学双靶识别。

同一套 PIK3CA 替换在 PIK3CA/PIK3CB 上方向相反：summary_min 由 0.500 升至 0.691 和 0.685，而 B 端受体保持冻结（Figure 5B）。因此受体选择不仅可以改变表观判别的幅度，也可以改变其方向。这一对比不支持把受体替换简单解释为 docking accuracy 的损失：受体实现是评价条件的另一独立维度，并可增强或削弱表观判别。两对、单口袋设计比单一 collapse 轶事更有解释力，但还不是普遍定律（K = 4；两对共享 PIK3CA）。不宣称具体分子机制。

耦合任务读法值得讨论，但仍是假说。`summary_min` 跟踪弱臂，因此替换 \(S_A\) 可以改变哪一臂成为瓶颈。PIK3CA/PIK3CB 原来的弱臂是 4L23 上的 D/B（0.500）；4JPS 后该臂升至 0.707，冻结 2WXF 臂（0.691）成为限制。其他非互斥假说包括局部侧链/口袋几何重排 dual 与选择性分数，以及两套面板面对同一 PIK3CA 晶体时不同的配体化学分布。由于没有面板范围的残基级 PLIF 分析，这些解释仍未解决。

5DXT 与 4L23 的局部口袋 Cα RMSD 仅为 0.343 Å，但 PIK3CA/mTOR 的 summary_min 仍降至 0.505。这说明“结构相似”与“screening discrimination 可迁移”并不是同一个问题。pose-generation QC 通过，也不等于 screening performance 不变。该结果与近期 cross-docking benchmark 中 receptor representation 被视为独立性能变量的观点相一致；那些工作使用的对接引擎与本文不同，不能当作同一协议的外推。[14]

### 4.4 对双靶虚拟筛选与生成式设计的含义

这些结果对双靶点虚拟筛选与生成式设计具有直接意义。同时在两个口袋获得 favorable docking scores 并不能自动等同于 experimentally plausible dual-active ligands。若双靶生成模型把 docking 当作下游过滤器，就不应只报告 two-pocket score，而应在 selective hard-negative 和 ligand-only controls 下评价。即便在单靶超大规模对接之后，hit 的后处理与再打分也已被证明难以稳健地区分已知结合分子与无活性分子；[15] 双靶场景额外要求同时压住两条实验硬负臂，因此更不能把两端有利分数或其简单平均当作充分证据。

本研究并不证明现有双靶生成模型无效，也没有直接评测 DualDiff、FuseDiff 或其他生成模型。[10,11] DualDiff 的 Dual High Affinity 是生成分子在两个靶上均优于各自参考配体，不是均值池化；FuseDiff 的独立测试集为 DualDiff benchmark（DDF）。这些工作回答的是结构生成能否获得有利对接分数。DualFourClass-Bench 可以作为这类方法的 downstream evaluation layer，用于检验生成分子是否真正超过实验定义的 single-target hard negatives，而不是仅仅优化 docking score。

### 4.5 局限性

以下五条限制界定本文结论的解释范围。

第一，评价集仅含四对靶标，因为实验定义的双靶硬负样本稀缺。K = 4 是受数据供给约束的案例面板，不是 comprehensive suite。四个 `summary_min` 还混合了面板构建差异（严格 6.5/5.5 对 θ = 6.0；不等 n）与靶对生物学，不能读成 intrinsic docking performance 的总体排序。

第二，ground truth 来自 ChEMBL。unused-pool holdout 仍属同一抓取批次，不是跨数据库独立验证。BindingDB/PubChem 核对仅为计数。

第三，活性聚合对照之后，assay heterogeneity 仍然存在。主策展使用最大 pChEMBL；换成重复测定中位数后 pair-level 估计变化很小（Results 3.4；Table S29），因此 max 聚合是 controlled limitation，不是未关闭的 fatal ground-truth 风险。pChEMBL 仍非 assay-equivalent。confidence≥8 与 Homo sapiens 过滤未重建。

第四，受体实现可以提高或降低成对判别，但实验并未给出分子起源。两个受体敏感性例子均共享 PIK3CA，因此不能视为两个互不相关替换靶标上的独立证据。口袋局域 Cα RMSD 不能单独解释性能变化，残基级 PLIF/侧链分析未系统进行。PIK3CA/PIK3CB 的一次对接超时在原始与替代 PIK3CA 晶体上均一致报告，排除原因不是其标签（Table S30）。

第五，本研究评价的是计算判别，而不是对新预测双靶化合物做前瞻实验。benchmark 回答的是对接排序的可靠性，不是入选分子的前瞻生物学效力。本文不旨在证明 docking 对双靶发现普遍有效或无效；它问的是评价 formulation 本身是否改变双靶识别的表观证据。

## 5. 结论

本研究建立 DualFourClass-Bench，作为针对 A-selective、B-selective 硬负样本的实验锚定计算判别环境。在四对冻结靶标上，最差方向判别摘要高度依赖靶对（`summary_min` 0.430–0.692），且把 docking 加到 ECFP4 后只产生很小的支架分组 CV AUROC 增量改善。在这四对面板中，主终点对最大 pChEMBL 与重复测定中位数聚合的选择总体不敏感。PIK3CA/mTOR 的点估计最高，并在未使用配体池稳定性检查中保持同向趋势；但其不确定性与受体敏感性排除了可推广的决策规则。

更广泛的分析表明，表观双靶对接判别取决于任务设定、配体化学组成和受体实现。EGFR/HER2 上，Dual-versus-neither comparator 为 0.756，而方向性最差一臂为 0.430；这是依赖靶对的描述性对照，不是配对显著性检验或四对定律。若干靶对上，ligand-only reference 达到或超过 docking；unused-pool holdout 还暴露出配对区间包含零的未解决 wrong-pocket point-estimate reversal。受体实现改变了表观判别的幅度和方向。这些结果支持把实验定义的选择性硬负样本、配体层混淆对照、面板外稳定性检查与受体敏感性作为双靶对接评价的互补要求。DualFourClass-Bench 的主要贡献是界定证据与可靠性边界的系统协议，而不是普适 docking winner。

## 数据与软件可用性

评价面板成员、实验状态标签、受体与对接盒定义、逐配体对接分数、分析表，以及重建本文统计与图件所需的全部脚本，均可在公开仓库 https://github.com/1280602962-debug/gwj260531 的 `Dual_Target_Docking` 目录中获取。`data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` 索引主要数值结果及其来源表。分析环境与零新对接的复现命令见仓库 README。

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
