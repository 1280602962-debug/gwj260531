# Introduction（中文工作稿 · JCIM Articles）

> 按五段连续科学论证重构，不再在原稿上逐句修补。  
> 投稿以英文为准：[`INTRODUCTION_SECTION_JCIM_EN_V1.md`](INTRODUCTION_SECTION_JCIM_EN_V1.md)。  
> 引用编号与核验边界：[`INTRODUCTION_REFS_JCIM_V1.md`](INTRODUCTION_REFS_JCIM_V1.md)。  
> 定位：[`POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`](POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md)。  
> **Introduction 不写死 K = 4**；冻结评价集规模放 Methods / Results。小节编号便于与 Methods、Figure 1 对齐，定稿时可改为连续段落。

---

## 1. 双靶点药物设计的兴起及结构基础虚拟筛选的作用

多靶点药物设计（multitarget drug design）旨在通过单一小分子同时调控两个或多个生物学靶点，以应对复杂疾病中的通路冗余、代偿性信号以及药物耐药等问题。与传统单靶点药物相比，合理设计的多靶点配体有望通过协同调节相互关联的生物学过程获得更充分的药理效应，因此已成为多靶点药物发现和多药理学（polypharmacology）研究的重要方向。[^1] 近年来，多靶点小分子的理性设计逐渐由经验性的多药理筛选，转向结合结构生物学、计算化学与生成式模型的结构导向设计。[^2]

在这一过程中，分子对接（molecular docking）仍是结构基础虚拟筛选（structure-based virtual screening, SBVS）中最常用的计算工具之一：先预测配体在蛋白结合口袋中的结合构象，再用打分函数对配体–受体相互作用排序，从而在大规模化合物库中给出结构互补性的近似评价。[^3][^4] 因此，在双靶点药物发现中，一个自然的计算策略是分别将候选分子对接至两个靶点，并据此判断其是否具有潜在的双靶结合能力。

与此同时，既有虚拟筛选研究已经表明，对接结果的解释高度依赖于数据集构建、负样本定义、化学偏倚以及评价指标。经典的 DUD 与 DUD-E 基准把评价建成“活性配体相对物化匹配 decoy 的富集”，并强调 decoy 若在粗物化性质上与活性物相差过大，表观 enrichment 可能只是在分离分子量或极性等配体层特征。[^5][^6] LIT-PCBA 进一步指出，DUD、DUD-E、MUV 一类人工构造的 active/decoy 数据集可能包含显著乃至隐蔽的化学偏倚，从而高估虚拟筛选方法的真确度；该基准改用实验剂量–响应标签，并对物化性质范围加以控制。[^7] 另一方面，CASF 一类结构打分基准把评分与对接搜索解耦，用 scoring / ranking / docking / screening power 评价单复合物上的打分函数，并不回答“两端实验活性状态空间”中的双靶识别问题。[^8]

然而，**将单靶点 docking 的评价逻辑直接扩展到双靶点任务并不充分。**

## 2. 双靶点识别与传统单靶点虚拟筛选任务存在本质差异

传统结构基础虚拟筛选通常可以表述为区分目标蛋白的活性配体与非活性化合物，即 Active versus Decoy。在这一任务中，负样本主要用于构建与目标活性配体相对的判别边界。

对于双靶点配体而言，任务结构发生了改变。一个严格的双靶点 benchmark 至少需要区分四种具有不同生物学含义的配体状态（四状态数据集，而不是四分类器）：同时作用于两个靶点的 **dual-active** 配体、仅作用于靶点 A 的 **A-selective** 配体、仅作用于靶点 B 的 **B-selective** 配体，以及两个靶点均缺乏足够活性的 **neither** 配体（Figure 1A）：

|  | *B*<sup>+</sup> | *B*<sup>−</sup> |
|--|:--:|:--:|
| *A*<sup>+</sup> | Dual | A-only |
| *A*<sup>−</sup> | B-only | Neither |

其中，A-only 和 B-only 并不是传统意义上的普通负样本，而是针对双靶识别任务最关键的**选择性硬负样本（selectivity hard negatives）**。它们在一个靶点上已经具有较强活性，因此可以产生看似合理的 docking score，但在另一个靶点上缺乏相应活性。对于一个真正具有双靶识别能力的评价方法而言，关键问题因而并不是候选分子是否能够在两个口袋中分别获得较好的 docking score，而是其是否能够在**两个方向上同时区别真正的 dual-active 配体与对应的单靶选择性配体**。

这一问题比已有的双靶对接评价更严格。Zhou、Li 与 Hou 曾在四对激酶上评估对接虚拟筛选：先做单靶 inhibitor 对 noninhibitor，再做 dual-target identification，并报告结构依赖性以及预测 dual 列表中较高的 false-positive rate。[^9] 该工作已经说明双靶对接可以被基准化，并且相对 inactive 的对接并不能给出干净的 dual hit list。它没有把实验标注的 A-only / B-only 当作方向性硬负，也没有问：同一套对接分数上，Dual-versus-neither（或 Dual-versus-inactive）读出是否会改变对 docking 能力的解释。

这一差异也意味着，简单地将两个靶点的 docking score 进行平均、求和或其他池化处理，并不能充分描述双靶识别能力。例如，一个配体可能在靶点 A 上获得非常有利的评分，而在靶点 B 上表现较差；其平均分仍可能较高，但这一结果并不能支持其具有双靶活性。类似地，与单一参考配体进行相对 docking score 比较，可以用于定义某种计算意义上的“双靶成功”，但并不能直接回答一个更严格的实验验证问题：**计算评分能否将真正的双靶活性配体与具有单靶强活性的选择性配体区分开来？**

因此，双靶点 docking 的关键评测问题不是简单寻找一个更优的 docking score，而是首先建立与其生物学状态空间相匹配的 benchmark 和评价任务。

## 3. 双靶点 benchmark 的构建受到实验数据供给和化学混淆的双重限制

建立上述四状态数据集在实际数据中并不容易。传统单靶点虚拟筛选 benchmark 可以依赖相对成熟的 active/decoy 构建策略，[^5][^6] 而严格的双靶点 benchmark 要求同一化合物在两个靶点上均具有可比较的实验活性信息，并进一步能够识别具有明确选择性差异的 A-only 和 B-only 化合物。换言之，一个可用于双靶点评测的实验数据集不仅需要足够数量的 dual-active 配体，还需要两个方向上数量和活性范围均相对充分的选择性硬负样本。DualFourClass-Bench 保留四种实验状态，但预先指定的主评价是两条方向 pairwise 判别（dual 对 A-only、dual 对 B-only），而不是四分类器。

这一要求显著提高了公开数据集构建的门槛。不同数据库之间的 assay 类型、活性指标、实验条件和记录覆盖范围并不完全一致，而双靶点任务还受到两个靶点同时具有可用实验数据这一额外条件的限制。因此，**可用于严格四状态双靶点评测的 target pair 数量本身就是一个需要量化的数据供给问题，而不能简单假设所有具有药理关联的靶点组合都能够形成平衡 benchmark。** 公开活性库中“两端都被定量测到、并拉开选择性间隔”的化合物并不自动充足；能够通过该完整性门槛的靶对有多少，是 benchmark construction 的科学问题，而不是事后用“数据不够”来解释评价集规模。

此外，双靶点 benchmark 还面临与传统虚拟筛选类似、但更容易被忽略的化学混淆问题。分子量、极性、氢键特征、脂溶性以及化学骨架等配体层面的特征可能同时影响实验活性分布和 docking score。如果 dual-active、A-only 和 B-only 化合物在这些性质上存在系统差异，那么一个看似具有较高 AUROC 的 docking 模型可能实际上是在利用配体层面的统计差异，而不是正确捕捉两个蛋白结合口袋中的结构互补性。LIT-PCBA 已经表明，人工构造的 active/decoy 数据中的化学偏倚能够显著抬高方法的表观性能，因此双靶点 docking benchmark 同样需要显式设置化学和物化性质对照。[^7]

## 4. 现有双靶点生成方法进一步提出了对严格 docking 评价的需求

近年来，双靶点结构基础生成方法进一步把“如何判断一个分子算双靶”推到了实际评价流程里。Zhou、Guan 等将双靶药物设计建成生成任务，并提出 DualDiff：在三维空间中对齐两个口袋，通过共享配体节点的 SE(3) 等变消息组合，把在单靶复合物上预训练的扩散模型迁移到双靶场景。[^10] 其对接评价使用 AutoDock Vina 重对接，报告两端的 Vina Dock、两端中较弱一端的 Max Vina Dock，以及 Dual High Affinity——即生成分子在**两个靶上的亲和力均优于各自参考配体**的比例。[^10] 这是相对参考配体的计算双成功，不是对两端分数做均值池化；Max Vina Dock 已经关注较弱一臂。

Wu 等提出的 FuseDiff 则将共享配体分子图与两个靶点特异的结合构象进行联合扩散建模，并以 DualDiff 基准（DDF）作为独立测试集评价生成结果，同样报告 Vina Dock、Max Vina Dock 与 Dual High Affinity。[^11] 这类工作说明，**如何评价一个分子是否真正满足“双靶点”要求，已经成为双靶点计算药物设计中的实际问题。**

但是，生成式双靶点研究中常用的评价方式与实验活性驱动的硬负判别任务并不完全相同。以两个口袋中的 docking score 是否同时超过某一参考配体为标准，可以衡量计算意义上的双靶 docking success，但它并不能直接检验该分子是否能够区别于仅对其中一个靶点具有较强活性的选择性配体。因而，这类生成方法的 docking-based evaluation 与本文关注的**实验标签驱动的 dual-versus-selective discrimination**属于互补的评价问题，而不是相互替代的 benchmark。本文没有对 DualDiff 或 FuseDiff 的生成分子重新对接；DualFourClass-Bench 的预期用途是下游校验，而不是对这些生成方法做实证打分赛。

## 5. 本研究的目的与贡献

基于上述问题，本研究建立了一个面向双靶点分子识别的系统性 docking benchmarking protocol，并构建公开的 **DualFourClass-Bench** 资源，用于评估 docking score 在严格双靶点识别任务中的可靠性及其适用边界。贡献应理解为系统评测协议与基准资源，而不是新的对接算法或打分函数。

首先，我们对公开活性数据中的候选 target pairs 进行系统审计，根据两个靶点上的实验活性以及选择性间隔，构建由 **dual-active、A-only、B-only 和 neither** 四种实验状态组成的评价面板。该数据审计同时用于量化严格四状态 benchmark 的实际数据供给，而不是预先假设所有候选靶对均能够满足相同的数据要求。neither 保留在面板中以描述完整状态空间；预先指定的主终点只使用 dual 对 A-only 与 dual 对 B-only。评价集是在经过严格数据供给审计后形成的多靶点面板上报告的（建造规则、受体与对接参数见 Methods）。

其次，我们采用与双靶点任务相匹配的 **pocket-matched directional AUROC** 作为主要评价指标（Figure 1B）。对于 dual 与 A-only 的比较，使用靶点 B 的 docking score 评价其对非选择性靶点的额外识别能力；对于 dual 与 B-only 的比较，则相应使用靶点 A 的 docking score。进一步以两个方向中较弱的一臂作为 summary measure（summary_min），从而避免一个靶点上的高评分掩盖另一个靶点上的识别失败。与此同时，将 pooled docking score、wrong-pocket control 以及二维化学和物化性质 baseline 作为辅助对照，以区分真正的 pocket-specific signal 与 ligand-level confounding。

最后，我们在经过供给审计后保留的多个 target pairs 上比较 docking-based discrimination 的一致性，并进一步通过 wrong-pocket、化学性质、scaffold-aware chemical baseline、同一数据批次中的 unused-pool holdout 以及 receptor-structure robustness 等分析考察其可靠边界。研究重点并非提出新的 docking scoring function，而是回答一个更基础的问题：

> **现有 docking scores 在多大程度上能够将实验定义的双靶活性配体与单靶选择性硬负配体区分开来，以及这种区分能力在多大程度上依赖于特定靶点、受体结构或配体化学性质？**

同一设计还引出一个嵌套问题：在同一套分数上，传统 Dual-versus-inactive（或 Dual-versus-neither）读出是否会比方向性硬负任务给出更乐观的解释。

用一句话界定问题：**The central question is therefore an experimentally defined dual-target recognition task: not whether a ligand can obtain favorable docking scores at two targets, but whether docking can distinguish dual-active ligands from target-selective hard negatives across both target directions.**

通过这一设计，本研究旨在为双靶点虚拟筛选和生成式双靶点药物设计提供一个更严格的下游评价基准，并明确 docking-based dual-target recognition 可以被可靠解释的范围及其潜在混淆来源。

---

[^1]: Anighoro, Bajorath, Rastelli, *J. Med. Chem.* **2014**, *57*, 7874–7887. DOI: 10.1021/jm5006463.
[^2]: Proschak, Stark, Merk, *J. Med. Chem.* **2019**, *62*, 420–444. DOI: 10.1021/acs.jmedchem.8b00760.
[^3]: Kitchen et al., *Nat. Rev. Drug Discov.* **2004**, *3*, 935–949. DOI: 10.1038/nrd1549.
[^4]: Eberhardt et al., *J. Chem. Inf. Model.* **2021**, *61*, 3891–3898. DOI: 10.1021/acs.jcim.1c00203.
[^5]: Huang, Shoichet, Irwin, *J. Med. Chem.* **2006**, *49*, 6789–6801. DOI: 10.1021/jm0608356. (DUD)
[^6]: Mysinger et al., *J. Med. Chem.* **2012**, *55*, 6582–6594. DOI: 10.1021/jm300687e. (DUD-E)
[^7]: Tran-Nguyen, Jacquemard, Rognan, *J. Chem. Inf. Model.* **2020**, *60*, 4263–4273. DOI: 10.1021/acs.jcim.0c00155. (LIT-PCBA)
[^8]: Su et al., *J. Chem. Inf. Model.* **2019**, *59*, 895–913. DOI: 10.1021/acs.jcim.8b00545. (CASF-2016)
[^9]: Zhou, Li, Hou, *J. Chem. Inf. Model.* **2013**, *53*, 982–996. DOI: 10.1021/ci400065e. 四对激酶 dual-target docking evaluation；inhibitor vs noninhibitor；结构依赖；高 false-positive。
[^10]: Zhou, Guan et al., *The Thirty-eighth Annual Conference on Neural Information Processing Systems (NeurIPS 2024)*; arXiv:2410.20688. DualDiff；Dual High Affinity 定义为两端均优于参考配体。与 [^9] 不是同一篇 Zhou。
[^11]: Wu et al., arXiv:2603.05567, 2026 (preprint). FuseDiff；独立测试集为 DualDiff benchmark (DDF)。
