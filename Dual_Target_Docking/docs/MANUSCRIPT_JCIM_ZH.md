# DualFourClass-Bench — 组装中文工作稿

> 供阅读与内部核对的中文主稿。**投稿以英文 [`MANUSCRIPT_JCIM_EN.md`](MANUSCRIPT_JCIM_EN.md) 为准。** 规范来源是下列中文章节文件；修改那些文件后运行 `python3 docs/assemble_manuscript_zh.py`。不要另开 `_V2` / `_FINAL` 分文件。
>
> 主张边界：[`../data/jcim_bench_v0/CLAIM_CEILING.md`](../data/jcim_bench_v0/CLAIM_CEILING.md)。主图：[`../figures/jcim_article/`](../figures/jcim_article/)。Methods 文末“写法说明”已在组装时去掉。

---

<!-- BEGIN TITLE_AND_ABSTRACT_JCIM_ZH_V1.md -->
# 题名与摘要（JCIM Articles 中文工作稿）

> 投稿以英文为准：[`TITLE_AND_ABSTRACT_JCIM_EN_V1.md`](TITLE_AND_ABSTRACT_JCIM_EN_V1.md)。  
> 主张边界：[`CLAIM_CEILING.md`](../data/jcim_bench_v0/CLAIM_CEILING.md)。  
> Dual-versus-neither 是 **不控制选择性的对照（nonselectivity-controlled comparator）**，不是“传统双靶基准”。  
> Dual-versus-neither 与方向性 AUROC 是 **描述性 formulation 对照**，不是配对显著性检验。  
> 数字与英文摘要一致，取自冻结 CSV / Table 2–3 / Table S29–S30。

---

## 题名（首选）

**用方向性选择性硬负样本评测对接双靶识别**

英文对应：*Benchmarking Docking-Based Dual-Target Recognition with Directional Selectivity Hard Negatives*

备选（更尖锐，略偏 Perspective，不作为首次投 JCIM 的首选）：

*基准 formulation 会改变表观证据：用方向性选择性硬负评测对接双靶识别*

不要使用：“对接能够/不能识别双靶配体”；把 DualFourClass 写成 comprehensive suite；把 Dual-versus-neither 写成“传统基准”。

---

## 摘要

两端有利的对接分数是否构成双靶识别的证据，尚未在方向性单靶硬负任务上得到充分检验。我们构建 DualFourClass-Bench，这是一套经策展的四对、四状态评价面板，含两条方向性主任务：dual-active 对 A-only 选择性配体在口袋 B 上打分，以及 dual-active 对 B-only 选择性配体在口袋 A 上打分。靶对汇总取较弱一臂（`summary_min`）。我们要问的是：这一基准的 formulation 本身是否会改变双靶识别的表观证据。在同一套冻结 AutoDock Vina 分数上，省略选择性配体的 Dual-versus-neither 对照可在部分靶对情境下给出过于有利的双靶识别印象。EGFR/HER2 是原理性案例：Dual versus neither 的 AUROC 为 0.756，而方向性 `summary_min` 为 0.430；混合库排序的 Top-10 中有 9 个是实验选择性配体。AChE/BChE 与 PIK3CA/PIK3CB 仅显示小且区间重叠的增量；PIK3CA/mTOR 的 Dual versus neither 因 neither n = 4 而效能不足，不解释为反向效应。在支架分组交叉验证下，把对接分数加到 ECFP4 后增量 AUROC 很小。将最大 pChEMBL 换成重复测定的中位数后，主终点的靶对估计大体不变；而受体实现可以改变表观判别的幅度甚至方向。双靶判别高度依赖靶对，并且在支架感知评价下，对接相对配体层化学基线只提供有限增量信息。这些结果并不构成四对高估定律，也不证明对接不含口袋特异信息。DualFourClass-Bench 是受数据供给约束的评价协议，不是全面的双靶套件。
<!-- END TITLE_AND_ABSTRACT_JCIM_ZH_V1.md -->

---

<!-- BEGIN INTRODUCTION_DRAFT_ZH_JCIM_V1.md -->
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

这一问题比已有的双靶对接评价更严格。Zhou、Li 与 Hou 曾在四对激酶上评估对接虚拟筛选：先做单靶 inhibitor 对 noninhibitor，再做 dual-target identification，并报告结构依赖性以及预测 dual 列表中较高的 false-positive rate。[^9] 该工作已经说明双靶对接可以被基准化，并且相对 inactive 的对接并不能给出干净的 dual hit list。它没有把实验标注的 A-only / B-only 当作方向性硬负，也没有问：同一套对接分数上，Dual-versus-neither（inactive）comparator 是否会改变对方向性 Dual-versus-selective 的解释。本文中 Dual-versus-neither 是 **nonselectivity-controlled comparator**，不是 “the conventional dual-target benchmark”。

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

现有双靶分子设计评价通常检验一个分子能否在两个靶点上同时获得有利分数，但这一标准并不直接检验其能否区别于实验上仅对其中一个靶点有活性的配体。这一区分很重要：一端的有利分数可以与对端识别失败并存，而选择性配体也可能在两个对接口袋中都看起来合理。

**因此，本文要问的是：benchmark 的 formulation 本身是否会改变双靶识别的表观证据。** 我们构建由实验定义的四状态配体面板，以针对 A-selective 与 B-selective 硬负的方向性口袋匹配判别作为主任务，并与不控制选择性的 Dual-versus-neither comparator 对照；随后检验该信号在化学、物化性质、配体池、活性聚合与受体结构对照下是否仍然成立。

贡献是评价协议与 curated benchmark 资源，而不是新的对接算法或打分函数。DualFourClass-Bench 是 **four-state curated benchmark with two directional primary tasks**：dual 对 A-only 在口袋 B 打分，dual 对 B-only 在口袋 A 打分（Figure 1B）。neither 保留以描述完整实验状态空间，不进入 primary AUROC。靶对汇总为较弱一臂（`summary_min`），使一端高分不能掩盖另一端失败。同一套分数上的 Dual-versus-neither 是 comparator，不是 “the conventional dual-target benchmark”。

公开数据供给审计首先回答有多少候选靶对能够支持这一四状态构建（Methods 2.1–2.3）。评价集规模是该审计的结果，而不是 Introduction 预先冻结的设计目标。pooled docking score、wrong-pocket control 以及二维化学和物化 baseline 作为辅助对照，用以区分 pocket-specific signal 与 ligand-level confounding。

嵌套的实验问题仍然是：现有 docking scores 在多大程度上能够将实验定义的双靶活性配体与单靶选择性硬负配体区分开来，以及这种区分能力在多大程度上依赖于特定靶点、受体结构或配体化学性质。该协议旨在为双靶虚拟筛选和生成式双靶设计提供更严格的下游校验——不是对这些生成方法做实证打分赛，也不是 comprehensive dual-target suite。

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
<!-- END INTRODUCTION_DRAFT_ZH_JCIM_V1.md -->

---

<!-- BEGIN METHODS_DRAFT_ZH_JCIM_V1.md -->
# Methods（中文工作稿 · JCIM Articles）

> 结构：预定义评价协议 + robustness / sensitivity / falsification，而不是实验记录压缩版。  
> 语气对照：Vu et al., *J. Chem. Inf. Model.* **2025**, 65, 4833–4843（写清做什么、参数与软件；数字进 Results / SI）。  
> 配套：[`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md)、[`INTRODUCTION_DRAFT_ZH_JCIM_V1.md`](INTRODUCTION_DRAFT_ZH_JCIM_V1.md)、英文稿 [`METHODS_SECTION_JCIM_EN_V1.md`](METHODS_SECTION_JCIM_EN_V1.md)。  
> **Methods 只写协议**；供给计数、cognate RMSD、AUROC、holdout 点估计等一律在 Results / SI。  
> 不编造 1000 次互不重叠 panel 重抽，也不重建已冻结面板。全面板 median 与第二靶对受体替换已完成，见 Table S29–S30。  
> DualFourClass-Bench 是 **four-state curated benchmark**；primary endpoint 是两条方向 pairwise AUROC，不是四分类器。

---

## 2. 方法

### 2.1 数据来源与活性数据整理

双靶评价所需的配体活性作为 **experimentally derived activity labels**，通过 ChEMBL Web API 的公开 activity 端点获取。靶对供给审计于 2026-07-23 冻结。pChEMBL 将若干摩尔浓度–响应型测定（如 IC50、Ki、Kd、EC50）转换为近似 −log10 活性尺度，便于大规模公开数据整合。不同 assay 类型、实验条件与测定体系并不等价；本文将 pChEMBL 作为策展中的统一近似，而不解释为同一条件下可直接比较的绝对结合亲和力。

同一配体–靶标若有多条可用 pChEMBL 记录，冻结数据包采用**最大 pChEMBL** 作为一对一代表值，用于主策展。assay 类型、条件与实验体系并不等价；取最大可能抬高单次测定读数。因此将活性聚合敏感性作为**预先指定的敏感性分析**：在从 ChEMBL activity 端点重拉 assay 级记录后，用重复测定的**中位数**替换最大值（Table S29）。该分析覆盖全部冻结评测面板，**不用于重定义面板成员或对接参数**。冻结 Vina 分数不重算。类别比较使用同一 θ = 6.0 规则。冻结文件（`mols_*.json`）仍只保存该代表浮点数；中位数标签存在 A4 表中，不作为重建的主面板。任一端缺少有效 pChEMBL 的配体不进入需要双端标签的分析。

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

**正文主比较采用预先统一的 θ = 6.0 标签。** Dual：两端 ≥ θ；A-only：A ≥ θ 且 B < θ；B-only 对称；neither：两端 < θ。建造阶段允许在严格规则下单端选择性过少时改用该单阈值规则凑齐配额；建造规则在抽样前按供给审计冻结，并写入 Table 1。阈值选择服务于可分析配额，不是在观察对接分数后回改标签。作为支持性敏感性分析，在 θ ∈ {5.5, 6.5} 与严格 6.5/5.5 规则下重标四种状态并重算口袋匹配 summary_min（Table S4）。该网格不是与 Table 2 竞争的第二套主标准。样本量过小的格子在 Results 中标记 underpowered，Methods 不预判其数值。

### 2.3 DualFourClass-Bench 面板构建

**The resource is a four-state curated benchmark with two directional primary tasks.** Dual, A-only, B-only 与 neither 均保留以描述实验空间；预先指定的主终点是 dual versus A-only 与 dual versus B-only。neither 不进入 primary directional AUROC。这不是四分类器 benchmark。

候选靶对按 2.2 的严格供给审计筛选。最终冻结评价集包含 PIK3CA/mTOR、AChE/BChE、PIK3CA/PIK3CB 与 EGFR/HER2。EGFR/HER2 按预先批准的角色保留为**供给受限案例**（`PAIR_ROLES_APPROVED_JCIM.yaml`），其组成不与其余靶对按同一厚面板供给条件等价。

每个靶对从符合相应实验标签的候选池中按预先冻结的类别配额抽样。面板抽样使用固定随机种子 20260729。在能够计算 Bemis–Murcko 支架的面板上施加支架封顶，以降低同一化学系列过度代表：PIK3CA/mTOR（PM48）同一类别内同一支架最多 2 个分子；EGFR/HER2 最多 5 个。AChE/BChE 与 PIK3CA/PIK3CB 在建面时 SMILES 尚未并入抽样表，无法施加 Murcko 封顶；实际抽样仅按类别配额与确定性随机顺序进行，**不再施加额外化学多样性约束**。事后可算的 Murcko 支架随冻结表报告。各面板的最终成员、状态标签、ChEMBL identifier、SMILES 与抽样脚本随冻结数据包提供；本文不在观察对接分数后重抽面板。

四对的建造规则并不相同。AChE/BChE 与 PIK3CA/PIK3CB 在严格 6.5/5.5 规则下抽样；EGFR/HER2 与 PIK3CA/mTOR 因严格规则下 B_only 过少而改用 θ = 6.0。因此跨对 AUROC 同时混合靶对生物学与面板构建差异（样本量、阈值、化学系列、受体），不能读成纯粹的 intrinsic docking performance。

配额与建造标签如下。AChE/BChE 与 PIK3CA/PIK3CB：严格 6.5/5.5，目标 dual / A_only / B_only / neither = 28 / 28 / 28 / 16（面板 n = 100）。EGFR/HER2：沿用既有 θ = 6.0 面板（n = 110）。PIK3CA/mTOR：θ = 6.0，主比较面板 PM48（n = 48；建造 dual / A_only / B_only / neither = 18 / 14 / 12 / 4），并在其上冻结受体与对接协议。

对接失败的配体–受体组合从该受体分数中剔除；任一端缺少可用分数的配体不进入需要两端分数的口袋匹配 AUROC，故分析用计数可低于建造定额（Table 1）。AUROC 因此是**以对接引擎能够处理的化合物为条件**的。尝试 / 成功 / 失败计数（含 AutoDock 原子类型 `B` 等化学覆盖失败）见 Table S27。

PIK3CA/mTOR 另构建扩面面板（历史名 PM110）：保留 PM48 全部 48 个配体，并按严格规则追加分子，目标配额 dual / A_only / B_only / neither = 30 / 30 / 30 / 25。PM110 是 PM48 的超集，用于评价面板规模增加后点估计是否同向，不是与其他靶对独立等价的 primary benchmark，也不是独立重复实验。主文跨对比较以 PM48 为准。

本文不以从供给池重复抽取互不重叠平衡面板的分布作为稳健性读出。该路径受硬负供给限制（定量见 Results）；正式的配体侧外推是一次 unused-pool holdout（2.11）。配体层有放回 bootstrap（2.8）描述固定面板内的不确定度，不称作供给池重抽。

**Table 1.** DualFourClass-Bench 评价集组成与对接设置（建造规则）

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

对每个靶对 A/B 计算两条二分类 AUROC。dual 对 A-only 使用**口袋 B** 的分数：

\[
\mathrm{AUC}_{D/A} = \mathrm{AUROC}(\text{dual},\;\text{A-only};\;S_B),
\]

以检验对接能否利用非选择性靶点 B 的结构信息，把 dual-active 与已在 A 端强效的 A-only 分开。dual 对 B-only 使用口袋 A 的分数：

\[
\mathrm{AUC}_{D/B} = \mathrm{AUROC}(\text{dual},\;\text{B-only};\;S_A).
\]

dual 始终为正类。neither 不进入上述对比。

Vina 输出结合能 \(E_{\mathrm{Vina}}\)（kcal mol\(^{-1}\)，通常越负表示预测结合越强）。定义

\[
S_{\mathrm{Vina}} = -E_{\mathrm{Vina}},
\]

使所有 primary scores 遵循“越大表示预测结合越强”。RTMScore 与 GNINA CNN 分数本身已是越高越好，不再取负。

#### 2.8.2 summary_min

靶对汇总为较弱一臂：

\[
\mathrm{summary}_{\min} = \min(\mathrm{AUC}_{D/A},\;\mathrm{AUC}_{D/B}).
\]

该规则是与双靶任务同构的 **worst-arm aggregation**，不是新的 scoring function。选择最小值是为了避免一端较强的 discrimination 掩盖另一端失败；它不是唯一自然的数学聚合。算术平均与调和平均作为敏感性报告（Table S26）。四对排序以及 EGFR Dual-versus-neither 对照在三种聚合下方向不变。全文只有一个主终点：统一 θ = 6.0 下的口袋匹配 Vina `summary_min`（Table 2；PIK3CA/mTOR 主面板为 PM48）。预指定次级终点为两条方向臂、RTMScore 口袋匹配、GNINA CNN best-of-9 口袋匹配，以及 2.8.3 的描述符面板。稳健性 / 证伪终点为 θ 网格、PM110、E = 8、unused-pool holdout、受体替换、错口袋对照（含配对 Δ）。探索性终点为 ECFP4、contact_count（非 PLIF）以及 pooled `vina_mean` 的 Top-10 硬负计数。完整层级见 Supporting Information Table S16。`vina_mean` 池化方向 AUROC **不是** Table 2。

#### 2.8.3 物化描述符对照

用 RDKit 计算预先指定的描述符面板：重原子数（GetNumHeavyAtoms）、分子量（MolWt）、cLogP（MolLogP）与 TPSA。每个描述符按与对接分数相同的方向 AUROC 流程评价，**正文与 SI 报告全部四个**（Table 2；Table S28）。其中 AUROC 最高者记为 **best single-descriptor reference**，只是该面板上的事后最大值，**不是** confirmatory competitor，也不是“trivial baseline”假设检验。为避免先选最优描述符再做正式比较的选择偏倚，docking 与该参考的配对 Δ 不以“击败 best descriptor”作为 confirmatory test（Table S19）。

#### 2.8.4 分数聚合对照

作为辅助分析，同时计算两端分数的 pooled mean、wrong-pocket assignment（定义见 2.9.1）以及 worst-pocket aggregation。它们不是 primary endpoint，只用于判断不同聚合是否改变双靶识别结论（Table S6）。

#### 2.8.5 Bootstrap 不确定度

AUROC 与 summary_min 的不确定度用配体层 bootstrap：在保持类别标签结构的条件下对配体有放回重采样，每次重算两条方向 AUROC 与 summary_min。\(B = 2000\)，随机种子 20260729，百分位数 95% CI 为 \([P_{2.5}, P_{97.5}]\)。错口袋与描述符等配对比较在**同一次**重采样上计算 \(\Delta = \mathrm{Metric}_1 - \mathrm{Metric}_2\)，得到 paired bootstrap 区间（Table S17、S19）。另报 Murcko 支架重采样区间作为对照；正文以配体层为准。置信区间作描述性不确定度；除预先定义的主终点外，不对多靶对、多对照做多重比较意义上的 confirmatory testing，也不把“CI 是否跨越 0.5”单独等同于正式显著性。

#### 2.8.6 Benchmark-formulation comparison

在同一套冻结 Vina 分数上，将 **Dual-versus-neither comparator**（实验 inactive；`vina_mean` 与 `vina_worst`）以及 Dual versus all non-duals 作为辅助对照，与方向性主终点并列。Dual-versus-neither 是本面板上的 **nonselectivity-controlled comparator**，不是声称既有双靶基准都以 Dual versus neither 为官方任务。neither 用于该对照，仍不进入 Table 2。PIK3CA/mTOR 的 neither n = 4 标记 underpowered。该比较只问：省略选择性硬负是否会改变对双靶识别的表观证据；不是第二套主终点，也不是配对显著性检验（负样本集合不同；Table 3；Table S22）。单靶类比——口袋 A 上 (dual + A-only) 对 (B-only + neither)，以及对称的 B 对照——仅作 Zhou 式背景。

### 2.9 混淆、证伪与化学对照

#### 2.9.1 Wrong-pocket falsification control

将靶点 A 与 B 的分数对调，配体、受体与其余分析设置不变，重算方向 AUROC 与 summary_min。该分析是 **falsification control**，不是用来证明口袋特异的阳性对照。固定面板上 matched > wrong **不**作为 pocket-specific signal 的证据。错口袋接近或高于匹配口袋，则视为对 pocket-specific interpretation 的反证。holdout 反转进一步说明：wrong-pocket **不是在面板迁移下可靠的通用负对照**。

#### 2.9.2 配体效率归一

各口袋分数除以重原子数，\(S_{\mathrm{LE}} = S_{\mathrm{dock}} / N_{\mathrm{heavy}}\)，再按 primary 流程计算方向 AUROC 与 summary_min，以检验对接分是否主要反映分子大小。

#### 2.9.3 效价与尺寸匹配子集

分别构建 \(|\Delta\mathrm{pChEMBL}| \leq 0.5\) 的 potency-matched 子集与 \(|\Delta N_{\mathrm{heavy}}| \leq 2\) 的 size-matched 子集，在子集上重算方向 AUROC。匹配会减小样本量；该分析只判断方向是否明显改变，不把低样本量子集的点估计当作独立强证据（Table S5）。

#### 2.9.4 Covariate-adjusted analysis

逻辑回归比较

\[
\mathrm{Model}_1:\ Y \sim S_{\mathrm{dock}}, \qquad
\mathrm{Model}_2:\ Y \sim S_{\mathrm{dock}} + N_{\mathrm{heavy}} + \mathrm{TPSA},
\]

其中 \(Y\) 为 dual 对相应选择性硬负的二分类标签。使用 scikit-learn `LogisticRegression`（\(C = 1.0\)，`max_iter = 2000`）。报告模型 AUROC、对接分数回归系数及其优势比（OR）。该分析询问对接分在控制分子大小与极性后是否仍有 residual discrimination，协变量模型不是 primary predictor。

#### 2.9.5 二维化学基线

Morgan/ECFP4（半径 2，2048 bit）加与 2.9.4 相同的逻辑回归，建立仅依赖二维结构的基线。评价采用 Bemis–Murcko scaffold `GroupKFold`，折数 \(\min(5, N_{+}, N_{-}, N_{\mathrm{scaffold}})\) 且至少两折，使同一骨架不跨训练/测试折。高 CV AUROC 只说明同一 Murcko 支架的分子不在训练/测试折共享时判别仍可保持，**不是** target-external generalization。PIK3CA/mTOR 上 \(n_{\mathrm{scaffolds}} \approx n\)，该折接近 leave-one-scaffold。随机 `StratifiedKFold` 仅作泄漏核对（Table S20）。增量模型（physchem、ECFP4、docking 及其组合）使用同一折；logistic docking AUROC 不是 Table 2 的 rank AUROC（Table S24）。A-only/B-only 相对 dual 的最近邻 ECFP4 Tanimoto 匹配在 T ≥ 0.3 / 0.4 / 0.5 报告，因为这些面板上 T ≥ 0.7 匹配为空（Table S23）。T ≥ 0.3 是 **similarity-constrained subset**，不是 chemically matched analogue set。

#### 2.9.6 Scoring-independent contact count

在已冻结的 Vina **mode-1** 姿态上计算不依赖打分函数的几何量：配体重原子中与受体重原子距离 ≤ 4.0 Å 的原子数

\[
N_{\mathrm{contact}} = \#\{i:\ \min_j d_{ij} \le 4.0\,\text{Å}\}.
\]

该描述符不使用对接能量函数。用 \(N_{\mathrm{contact}}\) 在口袋 A 上比较 dual 对 A-only、在口袋 B 上比较 dual 对 B-only，与错口袋对照的同口袋比较同构，作为 scoring-independent geometric confounder control，检验错口袋判别是否可能只反映更大分子产生更多埋藏接触。4.0 Å 为粗粒度接触阈值，**不是** PLIF。不预设其幅度与 Vina 错口袋一致（Table S11）。

#### 2.9.7 跨对序列一致性（探索性）

从各冻结受体 `*_protein.pdb` 用 Biopython `PDBParser` 提取最长蛋白链一级序列（仅标准氨基酸 ATOM），以 `PairwiseAligner`（BLOSUM62，全局比对，gap open = −11、extend = −1）计算靶对内全链序列一致性（分别以比对长度与较短链归一；Table S7）。该指标是整体相似度的粗粒度代理，不涉及口袋残基对应或结构叠合，不用于口袋 RMSD 或 PLIF 主张。

### 2.10 单靶富集参照

在 PIK3CA 4L23 与 mTOR 4JT6 上分别构建单靶 active–weak-active 集合。活性分子：pChEMBL ≥ 6.5。弱效分子：同靶已测定且 pChEMBL ≤ 5.5，并按分子量 ±50 Da、cLogP ±1.5、TPSA ±25 Å² 与活性分子做性质匹配。分子量与 logP 窗口沿用 property-matched decoy 的常见设定（Mysinger et al., *J. Med. Chem.* **2012**, *55*, 6582–6594）；TPSA 窗口为同一思想下增加的极性匹配。目标规模约 50 个活性分子与 150 个弱效分子。配体准备、受体、盒子与 Vina 协议与 PIK3CA/mTOR 主面板相同（exhaustiveness = 16）。报告 AUROC、EF1% 与 EF5%。该实验只提供单靶 docking enrichment 的背景参照，不替代 dual-target 的 summary_min。

### 2.11 未使用配体池 holdout

为检验结论是否依赖于冻结面板的具体成员，从严格标签池中排除已用于主面板与 PM110 的 ChEMBL 条目，在剩余 unused pool 中构建 **unused-pool, panel-external holdout**。它不是跨数据库或跨实验体系的 external validation：配体仍来自同一 ChEMBL 抓取批次、同一靶对与同一标签规则。

Holdout 只在 unused-pool 配额足以按 dual / A-only / B-only 各抽 20 个配体的靶对上构建。预先冻结为 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB；EGFR/HER2 不具备同等未用池抽样条件，记为 not eligible，不补做不等价抽样。PIK3CA/mTOR 排除的是 PM110 超集，从而覆盖 PM48。抽样种子 `HOLDOUT_SEED = 20260731`（区别于建面种子），同一 Murcko 支架在每一状态类中最多 3 个成员。清单先冻结、后对接。

Holdout 不参与主面板构建、对接协议调整或 primary endpoint 选择。受体、盒子、配体准备、exhaustiveness、打分与统计与主 benchmark 相同，并使用同一 `summary_min` 与配体层 bootstrap。未能产生 Vina 分数的配体–受体组合按 2.3 从需要该分数的分析中剔除。同一 holdout 配体上并列计算描述符对照；错口袋、效价/尺寸匹配与 contact count 按 2.9 在 holdout 上重算（Table S8、S13）。效价/尺寸匹配诊断不改写 Table S8 的主 holdout 数字。

### 2.12 受体结构敏感性分析

为评价 benchmark 结论对受体结构选择的敏感性，另选满足以下**预先声明**条件的替代晶体：（i）polymer entity 与目标蛋白真实对应，排除嵌合体或非目标同源骨架；（ii）含 ATP 位点或目标结合位点的小分子共晶；（iii）分辨率可接受；（iv）通过与 2.5 相同的 cognate redocking QC。实际进入对接的替代结构为 PIK3CA 4JPS、5DXT 与 mTOR 4JSX。该分析是 **receptor-structure sensitivity analysis**（receptor-realization effect），不是稳健性检验，也不是用来证明某一晶体“更正确”，更不是把 PIK3CA/mTOR 预设为结构不变的 positive case。目的是量化双靶判别终点对受体实现对的敏感性，而不是挑选更优受体结构。

替换采用**单口袋**设计。在 PIK3CA/mTOR（PM48）上，4JPS/5DXT 替换口袋 A，口袋 B 仍用冻结 4JT6 分数；4JSX 替换口袋 B，口袋 A 仍用冻结 4L23 分数；exhaustiveness = 16，与 PM48 主面板一致。在 PIK3CA/PIK3CB 上，同一套已准备的 4JPS/5DXT 替换口袋 A，口袋 B 仍用冻结 2WXF 分数；exhaustiveness = 8，与该主面板一致。新盒子按该替代晶体自身共晶配体、以 2.4 的同一 AABB 规则生成。配体准备、随机种子（20260727）、打分函数与 primary endpoint 与相应主分析相同。未能产生 Vina 分数的作业按 2.3 剔除；attempted / successful / failed 计数与换晶表一并报告。

作为探索性、零新对接的几何对照，在已冻结晶体坐标上做刚体叠合：Biopython `PDBParser` 提取最长链 Cα，按残基编号与残基名精确匹配，`Superimposer` 一次 Kabsch 拟合得全域 RMSD；口袋残基由参考结构共晶配体重原子 ≤5 Å 界定，在**同一变换**下计算口袋局域 RMSD，不做二次局部拟合。再将替代结构共晶配体按同一变换投影，计算与参考共晶配体质心的距离。不同结构匹配的 Cα 数目可以不同，全域 RMSD 因此不是等覆盖比较。本对照仅含有限数目的替代晶体，不预设 Cα RMSD 能够定量解释 AUROC 变化（Table S10）。

### 2.13 软件与数据可用性

计算在 Python 3 环境下完成。主要软件：RDKit 2026.3.1、meeko 0.7.1、AutoDock Vina 1.2.7、GNINA 1.3.2、RTMScore（`rtmscore_model1`）；Vina 姿态转 SDF 使用 Open Babel。刚体叠合与全链序列比对使用 Biopython（`PDBParser`、`Superimposer`、`PairwiseAligner`）。AUROC、逻辑回归与交叉验证使用 NumPy、SciPy、scikit-learn 与 pandas（版本见公开复现环境）。评价面板、对接分数、分析脚本与完整参数表随公开数据包提供，见 Data and Software Availability。

---
<!-- END METHODS_DRAFT_ZH_JCIM_V1.md -->

---

<!-- BEGIN RESULTS_DRAFT_ZH_JCIM_V1.md -->
# Results（中文工作稿 · JCIM 式证据链）

> 结构按 *供给 → formulation → 化学 → 聚合/受体敏感性 → 错口袋失败 → 探索性结构* 重排为 6 节。  
> Results 只报告数字、预定义比较与直接可见的模式；解释、机制假说与使用边界放 Discussion。  
> 全部数字可追溯至 `data/jcim_bench_v0/`、`data/jcim_strengthen_t0t1_v0/`、`data/jcim_holdout_v0/`、`data/jcim_structure_robust_v0/` 与 `data/jcim_supply_crossdb_v0/`。未做的全面板残基级 PLIF 不写入。  
> **主图：** Fig 1 任务定义；Fig 2 供给；Fig 3 formulation 对比（Dual-versus-neither vs directional summary_min）；Fig 4 弱臂/描述符；Fig 5 受体实现敏感性（含 4JSX）；Fig 6 错口袋反转；Fig 7 化学型/协变量检验。森林图与 unused-pool holdout 现为 Figure S4 / S5。详见 [`FIGURE_AND_TOC_PLAN_JCIM_V1.md`](FIGURE_AND_TOC_PLAN_JCIM_V1.md)。  
> **定位：** DualFourClass-Bench 是 *constrained but experimentally grounded* 的四状态评测资源，不是“完整双靶对接排行榜”。禁止 D-DRAF；禁止 “docking can/cannot identify dual-target ligands”。见 [`POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`](POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md)。

---

## 3. 结果

### 3.1 实验数据供给限制了严格双靶基准的构建

为确定公开生物活性数据是否能够支持严格的双靶点识别评测，我们首先对 49 对有 ChEMBL 缓存的候选靶标进行供给审计（Figure 2）。双靶对接评测需要四种实验标签状态：dual、A-only、B-only 与 neither（Figure 1A）。基准是四状态数据集；预先指定的主终点是两条方向 pairwise 判别（dual 对 A-only、dual 对 B-only）。我们将一端达到活性阈值、对端明确低活性的配体定义为方向性选择性硬负样本，用于检验对接分数能否同时压住两条单靶臂。

在严格标签规则下（dual：两端 pChEMBL ≥ 6.5；选择性类：活性端 ≥ 6.5 且对端 ≤ 5.5），尽管候选靶对数量较多，能够同时提供足量 A-only 与 B-only 硬负样本的靶对十分有限。两端严格硬负均不少于 50 的厚面板条件仅有 4 对满足。排除不适合作为常规小分子对接评测对象的金属依赖 HDAC1/HDAC6 后，PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 构成三个规模相对充足的靶对；EGFR/HER2 仅有 7 个严格 B-only 配体，因此被保留为供给受限案例，而不是与前三对等价的厚面板（Table 1）。

这一供给限制并非 ChEMBL 单一数据库特有的计数现象。对最终四对靶标进行 BindingDB 与 PubChem 的零对接计数核对后（Supporting Information Table S12），与 pChEMBL 更接近的 `equal_only` 规则下，前三对的 min HN 分别为 BindingDB 76 / 92 / 58、PubChem 86 / 97 / 61（ChEMBL 缓存为 80 / 78 / 56），仍全部 ≥ 50。EGFR/HER2 虽可在其他数据库中达到约 30 个 B 端硬负样本，仍不足以满足 ≥ 50 的厚面板标准。将不等式活性记录作为点估计（`as_is`）会显著增加 EGFR/HER2 的表观供给（BindingDB min HN 升至 85），但 92 个 as-is B_only 中有 49 个在 EGFR 端只有 `>` 记录；这一处理改变了“两端具有等式定量测定”的标签定义，因此未用于冻结 benchmark。PubChem 与 BindingDB 计数接近，符合沉积重叠，不作两次独立普查。

因此，最终 benchmark 的规模并非根据对接表现事后筛选，而主要由公开实验数据中方向性选择性硬负样本的可获得性所约束。DualFourClass-Bench 是一套受供给约束、但由实验标签锚定的评价集（a constrained but experimentally grounded benchmark），不是覆盖全部双靶任务的完整抽样。构建细节见 Methods 2.1–2.3。

### 3.2 基准 formulation 改变了 EGFR/HER2 上的表观双靶识别

在冻结的四对靶标上，采用统一 θ = 6.0 标签规则和口袋匹配方向 AUROC 对 Vina docking scores 进行评价（Figure 1B；Methods 2.8）。分数定义为 \(S=-E_{\mathrm{Vina}}\)（越大越好），dual 为正类。预先指定的靶对汇总为 `summary_min`（两臂较小值），使较强一臂不能掩盖较弱一臂。算术平均与调和平均只作敏感性聚合；三种聚合下四对排序不变（Table S26）。AChE/BChE 与 PIK3CA/PIK3CB 在建造时使用更严格的 6.5/5.5 规则，但在本数据上 θ = 6.0 给出完全相同的配体分类与 AUROC（Table S4）；EGFR/HER2 与 PIK3CA/mTOR 对阈值更敏感，严格规则下 B_only 过少并标记 underpowered，故严格规则只作支持性敏感性分析，不作第二套主标准。整张阈值网格内排序趋势保持一致（Figure S1A）。

这四个 `summary_min` **不是**可互换的 intrinsic docking performance。AChE/BChE 与 PIK3CA/PIK3CB 在严格供给规则下建面；EGFR/HER2 与 PIK3CA/mTOR 使用 θ = 6.0；面板还在 n、化学系列与受体上不同。跨对差异同时混合这些构建因素与靶对生物学。

EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 的方向性 summary_min 分别为 0.430、0.606、0.500 和 0.692（Table 2；Figure 4A；Figure S4）。不同靶对的主要限制来自不同的弱臂：EGFR/HER2 的 dual-versus-B-only AUROC 为 0.430，PIK3CA/PIK3CB 为 0.500；PIK3CA/mTOR 两个方向分别为 0.714 和 0.692（Figure 4A）。相对池化协议，口袋匹配抬高了点估计但未改变排序（Table S6）。

同一套冻结分数再按 Dual-versus-neither comparator 以及 Dual versus all non-duals 计分（Table 3；Figure 3）。Dual-versus-neither 是本面板上的 **nonselectivity-controlled comparator**（实验 inactive；`vina_mean`），不是声称既有双靶基准都以 Dual versus neither 为官方任务。两套 AUROC 使用不同负样本，是 **descriptive formulation contrast**，不是配对显著性检验。

EGFR/HER2 是 proof-of-principle。Dual versus neither 的 AUROC 为 0.756 [0.562, 0.920]（n_neg = 12），而方向性 summary_min 仍为 0.430 [0.284, 0.576]。Dual versus all non-duals 降至 0.551 [0.443, 0.666]，说明额外难度来自选择性配体。在 110 个 EGFR/HER2 配体的混合库中按 `vina_mean` 取 Top-10：1 个 dual、5 个 A-only、4 个 B-only、0 个 neither（EF10 = 0.393；hard-negative fraction = 0.90）；EF5 也低于随机（Table S25）。因此 Dual-versus-neither 读出在该对上会支持对接双靶识别，而方向性任务与筛选向 Top-10 都优先富集选择性配体。

该 formulation gap **依赖靶对**，不是四对 overestimation 定律。AChE/BChE 与 PIK3CA/PIK3CB 的 Dual-versus-neither 增量很小（0.649 与 0.559），区间与方向性臂重叠。PIK3CA/mTOR Dual versus neither 因 neither n = 4 而 underpowered，不作反向效应解释；该对 Dual versus all non-duals 为 0.674，接近 summary_min 0.692。因此 PIK3CA/mTOR 作为 **conditional directional signal**，而不是全文中心成功案例（Results 3.4）。

**Table 2.** 冻结 K = 4 评价集上的口袋匹配方向 AUROC（Vina，统一 θ = 6.0），并列出四个预先指定描述符的 `summary_min`。最高描述符是 best single-descriptor reference，不是 confirmatory competitor。错口袋与配体效率见 Table S6；描述符双臂见 Table S28。

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

因此 docking discrimination 并未表现出一致的跨靶对能力。PIK3CA/mTOR 是唯一一个 summary_min 点估计同时高于 0.5 和其 best single-descriptor reference（重原子数 0.463）的靶对，但其 95% bootstrap CI 仍为 0.464–0.802，未能排除随机。AChE/BChE（0.606）低于 TPSA（0.733）；EGFR/HER2（0.430）与 PIK3CA/PIK3CB（0.500）也未显示超过相应描述符参考的明确优势。

对接覆盖并不完整。主面板两端均得分：EGFR/HER2 110/110，AChE/BChE 95/100，PIK3CA/PIK3CB 99/100，PIK3CA/mTOR 48/48（Table S27）。PIK3CA/PIK3CB 唯一失败是 `PAB_034`（A-only；CHEMBL5089694），4L23 上对接超时（`timeout_900s`，23 个可旋转键），不是标签过滤；PIK3CB 2WXF 成功。AUROC 因此以 AutoDock Vina 能够处理的化合物为条件。采用 RTMScore 或 GNINA 作为替代 scoring channel 未改变总体排序。GNINA 在统一 best-of-9 pose coverage 后，EGFR/HER2、AChE/BChE 与 PIK3CA/mTOR 的口袋匹配 summary_min 仍不超过同面板 Vina；PIK3CA/PIK3CB 上 GNINA best-of-9 为 0.533、Vina 为 0.500，二者均近随机且区间重叠（Table S14–S15；Figure S1B）。GNINA 仍只是单一 CNN 通道对照。协议通过了 cognate pose-generation QC；该 QC 不是 screening-performance validation。

### 3.3 配体性质与化学型解释了相当一部分表观信号

为判断 docking discrimination 是否超越简单的 ligand-level signal，我们首先将口袋匹配对接与四种预先定义的物化性质进行比较（Figure 4B；Table 2）。相对于每个靶对的 **best single-descriptor reference**，docking summary_min 的 paired difference 在 EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 中分别为 −0.052、−0.128、−0.122 和 +0.229；四个 95% confidence intervals 均包含 0（Table S19；Figure S3C）。由此可见，即使 PIK3CA/mTOR 的点估计表现出最大的正向差异，现有样本仍不足以将其与 ligand-property reference 明确区分。该对照使用口袋匹配 summary_min，不是 pooled `vina_mean` 门控（EGFR/HER2 的 `vina_mean` 为 0.2824，≠ Table 2 的 0.4297）。

AChE/BChE 提供了一个较为直接的混淆案例。dual 配体平均 TPSA 约为 75，而选择性硬负配体约为 51（Figure 4C）；TPSA 单独获得约 0.769 的 AUROC，高于相同比较下的 Vina（约 0.56）。进一步加入 heavy-atom count 和 TPSA 后，dual-versus-B-only AUROC 从 0.606 增至 0.807，而 docking score 的 OR 仅约为 1.18（Figure 7C）。该结果表明，该方向上的 docking discrimination 很大程度上依赖与配体物化性质相关的信号，而不能直接解释为独立的 pocket-specific information。

PIK3CA/mTOR 的情况有所不同。加入 heavy-atom count 和 TPSA 后，AUROC 的变化约为 +0.07 至 +0.11，docking score 的 OR 约为 2.19 和 3.08，提示该靶对可能存在一定 residual pocket-related signal；然而，与 descriptor 的 paired difference 置信区间仍包含 0，因此这一残余信号不能被视为已确证的独立优势。配体效率归一后，仅 PIK3CA/mTOR 仍高于重原子数基线（0.657 对 0.463）。

二维化学结构 baseline 进一步说明了这一问题（Figure 7A）。ECFP4 + logistic regression 在 Bemis–Murcko scaffold GroupKFold 下多个方向获得约 0.78–0.91 的 fold AUROC，明显高于部分对应 docking contrasts，例如 EGFR/HER2 dual-versus-B-only 中 ECFP4 AUROC 为 0.85，而 docking AUROC 仅为 0.43。该结果只说明同一 Murcko 支架不跨训练/测试折时判别仍可保持，**不是** target-external generalization。PIK3CA/mTOR 上 \(n_{\mathrm{scaffolds}} \approx n\)，该折接近 leave-one-scaffold。同一设定下随机 `StratifiedKFold` 相对支架折的平均差为 +0.011（八个方向对比；Table S20；Figure S3D），泄漏很小。dual/selective 标签与 chemotype 存在系统性关联，因此单独观察 docking score 的 AUROC 并不足以证明其识别来源于 pocket-specific physical interactions。

在当前支架分组基准下，把口袋匹配对接分数加到 ECFP4 后，CV AUROC 的变化至多约 0.01，若干方向为负（Table S24）。这不是“docking 一般没有结构信息”：logistic 结构简单、K = 4，也没有 nested model comparison。logistic docking AUROC 不是 Table 2 的 rank AUROC，且常常更低。ECFP4 Tanimoto ≥ 0.7 的 chemotype-constrained A-only/B-only 子集为空。T ≥ 0.3 时，未匹配时最强的一臂（PIK3CA/PIK3CB dual versus A-only，0.691）降至 0.503（n_neg = 11），而远缘硬负（T < 0.3）升至 0.819（Table S23）。T ≥ 0.3 是 similarity-constrained subset，不是 chemically matched analogue set。T ≥ 0.4/0.5 的格子常为 n_neg ≤ 7，不作为第二套主结果解释。足够接近的化学匹配受数据集供给限制，这是四状态标签之外的第二层数据瓶颈。

效价匹配或尺寸匹配子集上，EGFR/HER2 与 PIK3CA/PIK3CB 的 dual 对 B_only 仍偏弱或接近随机（约 0.45–0.52）；PIK3CA/mTOR 的排序趋势保持一致，但各臂 n 常低于 15、区间较宽（Table S5；Figure 7D）。全部四个描述符见图 7B，均不作 confirmatory competitor。

### 3.4 表观性能对 pChEMBL 聚合不敏感，但对受体实现敏感

主标签使用可用 pChEMBL 的最大值。在对每个已打分配体重新拉取 assay 级记录后，把该聚合换成重复测定的中位数，θ = 6.0 下四状态类别翻转：EGFR/HER2 7/110（标签一致率 103/110 = 93.6%），AChE/BChE 1/95（94/95 = 98.9%），PIK3CA/PIK3CB 1/99（98/99 = 99.0%），PIK3CA/mTOR 0/48（48/48 = 100%）（Table S29）。数值上 max ≠ median 比类别翻转更常见（40/110、13/95、25/99、27/48）。API 重拉标签上，`summary_min` 由 0.417→0.424（EGFR/HER2）、0.606→0.629（AChE/BChE）、0.500→0.500（PIK3CA/PIK3CB）、0.692→0.692（PIK3CA/mTOR）。冻结 Table 2 的 EGFR/HER2 是 0.430 而非 0.417，因为一处缓存/API 不一致（`EH120_060` / CHEMBL24828）在 API max 下把该配体标成 dual；相对冻结表，中位数聚合仍使 EGFR/HER2 为 0.424。靶对排序与方向性主结论因此对这一聚合选择不敏感。assay 间异质性仍然存在，因为 pChEMBL 并非 assay-equivalent。

为判断 PIK3CA/mTOR 的较高 summary_min 是否仅由特定 panel 构成或 docking 搜索参数造成，我们进行了 ligand-panel 和 protocol-level sensitivity analyses（Figure S5）。将 exhaustiveness 从 16 降至 8 后，summary_min 从 0.692 降至 0.660，变化约 0.03，明显小于不同 target pairs 之间的性能差异（Figure S1D）。

在包含 PM48 全部配体并扩展至实际 n = 115 的 PM110 面板中（分析用 dual / A_only / B_only 各 30），Vina summary_min 为 0.648 [0.51, 0.76]，相比 PM48 的 0.692 下降约 0.04，但排序趋势保持一致（Figure S1C）。该结果支持 PIK3CA/mTOR 的方向性信号并非完全由 PM48 的特定成员驱动，但 PM110 与 PM48 并非独立验证集，因此该结果应解释为 stability check。同面板 RTMScore 为 0.576；GNINA best-of-9 为 0.613 [0.46, 0.74]，PM48 同口径为 0.655 [0.43, 0.81]，仍不高于同面板 Vina。

更重要的是，在未参与主面板构建和协议调优的 unused-pool holdout 中（每对 20 / 20 / 20，种子 20260731；EGFR/HER2 不具备同等配额，记为 not eligible），PIK3CA/mTOR 的 summary_min 进一步达到 0.765 [0.603, 0.891]，高于主面板的 0.692；AChE/BChE 为 0.618 [0.422, 0.759]，与主面板接近但 confidence interval 跨越 0.5；PIK3CA/PIK3CB 则下降至 0.425 [0.241, 0.618]（Table S8 / Table S16）。PIK3CA/PIK3CB holdout 尝试 60 个配体，59 个两端得分；HOAP_028 因 AutoDock 原子类型 `B` 不支持（含硼）而两端失败（Table S27）。AChE 与 PIK3CA/mTOR holdout 为 60/60 成功。硼失败是引擎化学覆盖限制，不是 silent missingness；AUROC 以可处理化合物为条件。该 holdout 共享同一 ChEMBL 抓取批次，不能读成跨机构独立验证；其作用是支持所观察信号在未参与建面配体池中的持续性。

因而，PIK3CA/mTOR 的方向性 signal 在同一 ChEMBL 体系的未参与建面配体中仍然可观察到，而 PIK3CA/PIK3CB 的 signal 则未能保持。这进一步说明 docking performance 主要由 target-pair context 决定，而不是一个可在不同靶对之间稳定迁移的属性。

尽管 PIK3CA/mTOR 在 ligand-panel sensitivity analysis 中保持了方向性信号，我们进一步测试方向性判别是否依赖于特定 receptor realization：一端受体冻结，只替换另一端（Figure 5；Table S9；Table S30）。三个替代晶体结构均通过 cognate redocking QC，best-of-9 RMSD 分别为 0.607 Å（4JPS）、0.624 Å（5DXT）和 0.515 Å（4JSX）；嵌合体 3T8M 已排除。

在 PIK3CA/mTOR 上，当 PIK3CA 4L23 替换为 4JPS 或 5DXT、mTOR 4JT6 保持不变时，PM48 的 summary_min 分别由 0.692 降至 0.486 [0.259, 0.692] 和 0.505 [0.292, 0.696]（Figure 5A）。变化主要发生在依赖替代 PIK3CA 结构的 D/B direction，而依赖原始 mTOR 结构的 D/A direction 保持 0.714。将 mTOR 4JT6 替换为 4JSX 后 summary_min 为 0.639 [0.418, 0.776]。mTOR 端换晶后点估计仍高于 0.5，但 95% CI 包含 0.5。

同一套 PIK3CA 晶体再用于 PIK3CA/PIK3CB 面板，2WXF 分数保持冻结（exhaustiveness = 8，与主面板一致；Figure 5B）。替换后 summary_min **上升**：由 0.500 至 0.691 [0.516, 0.779]（4JPS）和 0.685 [0.506, 0.768]（5DXT）。仍使用冻结 2WXF 的 dual versus A-only 保持 0.691。使用替代 PIK3CA 分数的 dual versus B-only 由 0.500 升至 0.707（4JPS）和 0.685（5DXT）。弱臂因此切换：原来是 4L23 上的 D/B（0.500）；4JPS 后瓶颈变成冻结 2WXF 臂（0.691）；5DXT 后两臂接近平衡。两套替代作业均尝试 100 个配体、成功 99 个；缺失配体仍是 `PAB_034`（4JPS 与 5DXT 上均 600 s 超时）。该配体在原始 4L23 上已经超时，因此 99 配体集合与 Table 2 相同。本协议下任何 PIK3CA 晶体都没有 100 配体方向 AUROC；失败是对接超时，不是实验标签过滤。

因此，同一 PIK3CA 扰动在两个靶对上方向相反。受体选择既可以增强、也可以削弱表观双靶 discrimination。这是 receptor-realization effect，不是稳健性证明，也不是单向 collapse。设计只扰动一端，Δ 可归因于被替换口袋，而不是同时换掉两个结构。两对共享 PIK3CA；该格局不是 K = 4 上的普遍定律。

Cα structural comparison 进一步显示，5DXT 与 4L23 的口袋局域 Cα RMSD 仅为 0.343 Å，但 PIK3CA/mTOR 的 summary_min 仍降至 0.505，说明简单的 backbone similarity 并不足以保持判别（Table S10）。这批 PIK3CA 沉积结构彼此的整链 Cα RMSD（1.44–1.49 Å）大于这批 mTOR 沉积结构彼此的差异（0.45 Å），与 PIK3CA/mTOR 上 PIK3CA 端变动更大的方向一致，但不能定量解释 PIK3CA/PIK3CB 上的相反位移：5DXT 仅匹配 862 个 Cα，少于 4JPS 的 982 个；替代结构各仅 1–2 个。共晶配体质心距离 2.1–2.6 Å 只说明对接的仍是同一大类 ATP 竞争位点。姿态生成 QC 通过，并不等于 screening discrimination 可迁移。协议通过了 cognate pose-generation QC；它不是 virtual-screening validation。

### 3.5 错口袋对照揭示仍未解决的失败模式

在主面板中，pocket-matched summary_min 均高于 wrong-pocket control，四对的 matched-minus-wrong differences 分别为 0.170、0.161、0.151 和 0.090；其中 EGFR/HER2 和 AChE/BChE 的差异置信区间排除 0，PIK3CA/PIK3CB 与 PIK3CA/mTOR 的区间包含 0（Table S6；Table S17；Figure 6A；Figure S3A）。错口袋对照的 summary_min 分别为 0.260、0.444、0.349 与 0.602。主面板上 matched > wrong **不**作为 pocket-specific signal 的证据。

然而，这一关系在 unused-pool holdout 中发生反转（Figure 6B）。PIK3CA/mTOR、AChE/BChE 和 PIK3CA/PIK3CB 的 wrong-pocket summary_min 分别为 0.788、0.643 和 0.520，而 matched-pocket 分别为 0.765、0.618 和 0.425。相应的 matched-minus-wrong point differences 均为负（−0.023 / −0.025 / −0.095），但其 95% confidence intervals 均包含 0（Table S17；Figure S3B）。holdout 反转是点估计格局，不是区间已排除零的统计结论。因此 wrong-pocket **不是在面板迁移下可靠的通用负对照**。

为检验这一反转是否可以由 holdout 中的配体效价或分子大小差异解释，进一步进行 potency- and size-matched comparisons（Figure 6C；Table S13）。wrong-pocket ≥ matched-pocket 的关系仍未翻转（效价匹配后：AChE/BChE 0.642 对 0.593，n_min = 18；PIK3CA/PIK3CB 0.562 对 0.363，n_min = 11；PIK3CA/mTOR 0.734 对 0.715，n_min = 12）。holdout 相对主面板确有抽样偏移——最明显的是 PIK3CA/mTOR：holdout dual / A_only 的 pA 均值比主面板低约 1.1–1.3，B_only 的 pB 低约 1.8——但匹配后悖论仍在。

scoring-independent contact_count 在 B direction 获得 0.698–0.714 的 AUROC，提示 ligand size/burial 对该方向确实存在影响，但其幅度不足以解释全部 Vina wrong-pocket signal（Figure 6D；Table S11）。例如 PIK3CA/mTOR 中 Vina wrong-pocket summary_min 为 0.788，而 contact_count 的较弱一臂仅为 0.552。A 臂上 dual 与 A_only 尺寸差很小，contact_count AUROC 更接近随机（0.552–0.622）。

因此，holdout 中 wrong-pocket reversal 应被视为当前 benchmark 暴露出的 unresolved failure mode，而不是可以由单一尺寸或效价因素解释的现象。

### 3.6 结构背景只提供探索性线索

作为探索性分析，我们进一步比较了靶对内全链序列一致性与 summary_min（Table S7）。四对靶标中，PIK3CA/mTOR 具有最低的全链序列一致性（对齐长度分母 18.1%）但最高的 summary_min，而 EGFR/HER2 恰好相反（71.4%；ErbB 家族激酶域高度同源）。该现象与简单的“靶点越相似越难区分”假设并不一致。然而，该分析仅包含四个靶对，且全链序列一致性并不是 binding-pocket similarity 的直接度量，因此该观察仅作为结构背景线索，而不作为相关性证据。PIK3CA 与 mTOR 同属 PIKK 相关家族，ATP 竞争位点存在已知的局部结构同源性；表中的低全链一致性不应解读为“两口袋不相似”。

在已导出的姿态级诊断中，PIK3CA/mTOR 上可观察到两类代表性 failure typology（非全面板 PLIF）。T2：选择性硬负在两个口袋都形成几何干净、hinge 阳性的 ATP-like pose（例如 amino-triazine / morpholine–ATP 化学型在弱端 mTOR 上仍占用高、hinge 阳性），使两端分数同时偏高。T5：部分经典 dual（如 Torin1、Omipalisib）在 Vina 中较强，但 alternative rescoring 的优选姿态可偏离 PIK3CA hinge/共晶位。共晶配体 PI-103 / X6K 在协议检查中可回收近晶姿态（Table S3）。这些是观察到的姿态模式，不是残基级机制；pose-generation QC 也不是 screening validation。
<!-- END RESULTS_DRAFT_ZH_JCIM_V1.md -->

---

<!-- BEGIN DISCUSSION_DRAFT_ZH_JCIM_V1.md -->
# Discussion（中文工作稿 · JCIM Articles）

> 与 [`DISCUSSION_LIMITATIONS_DRAFT_ZH_JCIM_V1.md`](DISCUSSION_LIMITATIONS_DRAFT_ZH_JCIM_V1.md) 配套：本文件写解释、替代解释与使用边界；完整局限库存仍在 Limitations 稿，正文 4.5 只保留五条。收束段见 [`CONCLUSIONS_DRAFT_ZH_JCIM_V1.md`](CONCLUSIONS_DRAFT_ZH_JCIM_V1.md)。  
> 口径对照 [`DISCUSSION_RESULT_MAP_V1.md`](DISCUSSION_RESULT_MAP_V1.md)；引用核验 [`DISCUSSION_REFS_JCIM_V1.md`](DISCUSSION_REFS_JCIM_V1.md)。  
> 五节：formulation → chemistry → receptor realization → implications → limitations。不把开放问题写成已解决。

---

## 4. 讨论

### 4.1 基准 formulation 改变了双靶对接的证据标准

本研究的首要发现不是某一种 docking scoring function 在双靶任务上取得了最高性能，而是 **benchmark formulation 可以改变“双靶对接成功”的表观含义**。传统 docking benchmark 通常将活性配体与 decoy 进行区分，而本文要求模型同时区分 dual-active 配体与两个方向上的 single-target selective hard negatives。后者在一个靶点上具有较强实验活性，因此不能被简单视为普通 decoy。同一套分数上的 Dual-versus-neither 是 **nonselectivity-controlled comparator**，不是“the conventional dual-target benchmark”。供给审计显示，公开生物活性数据很难同时提供两个方向上足够数量的这类配体；49 个候选靶对中只有少数能够满足严格厚面板要求（Results 3.1；Figure 2）。

Zhou、Li 与 Hou 已经表明，对接用于双激酶筛选时，相对 noninhibitor 可以看起来有用，表现依赖结构，并且预测 dual 列表仍有较高 false-positive rate。^(9) DualFourClass-Bench 在同一套分数上追问更窄的问题：Dual-versus-neither（inactive）读出与方向性 Dual-versus-selective 读出是否一致。它们在 EGFR/HER2 上并不一致（Results 3.2，Table 3，Figure 3）：Dual versus neither 为 0.756，方向性 summary_min 为 0.430，混合库 Top-10 富集选择性配体（9/10）。**EGFR/HER2 是清晰的例子，不是四对定律。** AChE/BChE 与 PIK3CA/PIK3CB 的增量很小且区间重叠；PIK3CA/mTOR Dual versus neither 因 neither n = 4 而 underpowered。该对照是描述性的，不是配对显著性检验。相对 2013 年工作的增量是这一 formulation gap，而不是又一次四对对接普查。

这一数据限制本身具有方法学意义。DUD、DUD-E 与 LIT-PCBA 已经表明，decoy construction、chemical bias 和真实 assay 标签会显著改变虚拟筛选的性能判断。^(5–7) 简单方法或不恰当的 unbiasing 也可以通过学习配体分布而高估 structure-based virtual screening。^(12) 近期基于 bioassay-derived data 的评价则进一步强调，与人工构造的 ligand/decoy 集相比，真实 assay-derived benchmarks 可以揭示模型在更接近实际筛选环境中的局限。^(13) DualFourClass-Bench 并未使用这些单靶数据集，也不评价 DiffDock-Pocket；它把同一关切延伸到双靶任务：评价结论取决于硬负样本如何被实验定义，而不是取决于候选靶对清单有多长。

因此，DualFourClass-Bench 的主要价值并不在于提供一个规模很大的数据集，而在于将双靶识别问题转化为一个实验标签驱动的 hard-negative discrimination task，并显式要求两个方向同时成立。该资源是 curated four-pair panel + evaluation protocol，不是 comprehensive dual-target suite。

错口袋对照属于同一证据标准下尚未解决的面板外失败模式。主面板上 pocket-matched 高于 wrong-pocket；unused-pool holdout 上点估计反转，但配对区间仍包含零（Results 3.5；Figure 6）。因此 wrong-pocket 不是在面板迁移下可靠的通用负对照；主面板 matched > wrong 也不作为口袋特异证明。

### 4.2 化学信息可以替代表观对接信号

双靶点 docking 的困难并不能简单理解为两个单靶 docking 任务的性能相加。即便采用任务对齐的方向性指标，四对靶标仍表现出明显异质性。三个靶对的 summary_min 位于随机水平附近或低于 best single-descriptor reference，只有 PIK3CA/mTOR 表现出较高的点估计，且其置信区间仍与随机相容（Results 3.2）。该对是 **conditional directional signal**，不是可推广的成功案例。

AChE/BChE 上 TPSA 单独即可获得高于 docking 的 discrimination，而 ECFP4 scaffold-grouped baseline 在多个方向上进一步超过 docking（Results 3.3）。dual/selective 标签本身携带的 ligand-level information 可以产生强烈的 apparent signal。在当前支架分组基准下，把对接分数加到 ECFP4 后 CV AUROC 变化至多约 0.01，若干方向为负。**dual-target discrimination was strongly target-pair dependent, and docking provided limited incremental information beyond ligand-level chemical baselines under scaffold-aware evaluation.** 这是对本面板的陈述，不是 docking 没有口袋特异信息的证明。

化学型约束硬负给出同一方向的证据。T ≥ 0.7 的匹配子集为空。T ≥ 0.3 时，未匹配时最强的一臂（PIK3CA/PIK3CB dual versus A-only，0.691）降至 0.503（n_neg = 11），而远缘硬负（T < 0.3）升至 0.819。T ≥ 0.3 是 similarity-constrained subset，不是 chemically matched analogue set。

这一结果与近年来对虚拟筛选 benchmark 中化学偏倚的关注一致：简单模型或不恰当构造的 decoys 可以通过学习 ligand distribution 而获得看似优异的 performance。^(7,12) 如果不设置 A-only/B-only hard negatives 以及 ligand-property/chemical baselines，一个看似优秀的 dual-target docking result 可能实际上只是识别了与 dual label 相关的分子属性。

### 4.3 受体实现是另一独立的不确定来源

PIK3CA/mTOR 提供了本研究中最值得进一步研究、但也最需要谨慎解释的案例。其主面板 summary_min 为 0.692，PM110 扩展面板为 0.648，unused-pool holdout 为 0.765，说明该方向性信号并非完全由 PM48 中少数配体驱动。与此同时，该信号并未表现出 receptor invariance：替换 PIK3CA 受体结构后，summary_min 降至 0.486 和 0.505，而替换 mTOR 结构后仍为 0.639（Results 3.4；Figure 5A）。更准确的说法不是“PIK3CA/mTOR docking 可以可靠识别双靶配体”，而是该靶对在特定 receptor realization 下存在有限的 directional signal：配体面板替换后可以持续，受体实现替换后不能假定不变。

同一套 PIK3CA 替换在 PIK3CA/PIK3CB 上方向相反：summary_min 由 0.500 升至 0.691 和 0.685，而 B 端受体保持冻结（Figure 5B）。因此受体选择不仅可以改变表观判别的幅度，也可以改变其方向。受体实现是独立的方差来源：既可以增强、也可以削弱表观判别，而不仅仅决定某一结果是否“稳健”。两对、单口袋设计比单一 collapse 轶事更有解释力，但还不是普遍定律（K = 4；两对共享 PIK3CA）。不宣称具体分子机制。

耦合任务读法值得讨论，但仍是假说。`summary_min` 跟踪弱臂，因此替换 \(S_A\) 可以改变哪一臂成为瓶颈。PIK3CA/PIK3CB 原来的弱臂是 4L23 上的 D/B（0.500）；4JPS 后该臂升至 0.707，冻结 2WXF 臂（0.691）成为限制。其他非互斥假说包括局部侧链/口袋几何重排 dual 与选择性分数，以及两套面板面对同一 PIK3CA 晶体时不同的配体化学分布。残基级 PLIF 分析未完成。

5DXT 与 4L23 的局部口袋 Cα RMSD 仅为 0.343 Å，但 PIK3CA/mTOR 的 summary_min 仍降至 0.505。这说明“结构相似”与“screening discrimination 可迁移”并不是同一个问题。pose-generation QC 通过，也不等于 screening performance 不变。该结果与近期 cross-docking benchmark 中 receptor representation 被视为独立性能变量的观点相一致；那些工作使用的对接引擎与本文不同，不能当作同一协议的外推。^(14)

### 4.4 对双靶虚拟筛选与生成式设计的含义

这些结果对双靶点虚拟筛选与生成式设计具有直接意义。同时在两个口袋获得 favorable docking scores 并不能自动等同于 experimentally plausible dual-active ligands。若双靶生成模型把 docking 当作下游过滤器，就不应只报告 two-pocket score，而应在 selective hard-negative 和 ligand-only controls 下评价。即便在单靶超大规模对接之后，hit 的后处理与再打分也已被证明难以稳健地区分已知结合分子与无活性分子；^(15) 双靶场景额外要求同时压住两条实验硬负臂，因此更不能把两端有利分数或其简单平均当作充分证据。

本研究并不证明现有双靶生成模型无效，也没有直接评测 DualDiff、FuseDiff 或其他生成模型。^(10,11) DualDiff 的 Dual High Affinity 是生成分子在两个靶上均优于各自参考配体，不是均值池化；FuseDiff 的独立测试集为 DualDiff benchmark（DDF）。这些工作回答的是结构生成能否获得有利对接分数。DualFourClass-Bench 可以作为这类方法的 downstream evaluation layer，用于检验生成分子是否真正超过实验定义的 single-target hard negatives，而不是仅仅优化 docking score。

### 4.5 局限性

正文只强调五条最高优先级限制；完整清单见 Limitations 稿。收束主张见 Conclusions，不在此重复。

第一，评价集仅含四对靶标，因为实验定义的双靶硬负样本稀缺。K = 4 是受数据供给约束的案例面板，不是 comprehensive suite。四个 `summary_min` 还混合了面板构建差异（严格 6.5/5.5 对 θ = 6.0；不等 n）与靶对生物学，不能读成 intrinsic docking performance 的总体排序。

第二，ground truth 来自 ChEMBL。unused-pool holdout 仍属同一抓取批次，不是跨数据库独立验证。BindingDB/PubChem 核对仅为计数。

第三，活性聚合对照之后，assay heterogeneity 仍然存在。主策展使用最大 pChEMBL；换成重复测定中位数后 pair-level 估计变化很小（Results 3.4；Table S29），因此 max 聚合是 controlled limitation，不是未关闭的 fatal ground-truth 风险。pChEMBL 仍非 assay-equivalent。confidence≥8 与 Homo sapiens 过滤未重建。

第四，受体实现可以提高或降低成对判别，但实验并未给出分子起源。口袋局域 Cα RMSD 不能单独解释性能变化，残基级 PLIF/侧链分析未系统进行。PIK3CA/PIK3CB 的一次对接超时（`PAB_034`）在原始 4L23 与两套替代 PIK3CA 晶体上均为 100 attempted / 99 successful / 1 failed；排除原因不是其标签。

第五，本研究评价的是计算判别，而不是对新预测双靶化合物做前瞻实验。benchmark 回答的是对接排序的可靠性，不是入选分子的前瞻生物学效力。本文不旨在证明 docking 对双靶发现普遍有效或无效；它问的是评价 formulation 本身是否改变双靶识别的表观证据。
<!-- END DISCUSSION_DRAFT_ZH_JCIM_V1.md -->

---

<!-- BEGIN CONCLUSIONS_DRAFT_ZH_JCIM_V1.md -->
# Conclusions（中文工作稿 · JCIM Articles）

> 投稿以英文为准：[`CONCLUSIONS_SECTION_JCIM_EN_V1.md`](CONCLUSIONS_SECTION_JCIM_EN_V1.md)。  
> 两段：做了什么 + 得到了什么；意味着什么 + 未来评价标准。不复制 Results，不写 PDB/引擎数字，不写 validated / robust / 通用决策规则。  
> 主张边界：[`CLAIM_CEILING.md`](../data/jcim_bench_v0/CLAIM_CEILING.md)。术语分层见文末。

---

## 5. 结论

本研究建立 DualFourClass-Bench，作为对接双靶识别的实验锚定评价环境，显式检验对接能否将 dual-active 配体与 A-selective、B-selective 硬负样本区分开来。在四对冻结靶标上，双靶判别高度依赖靶对（`summary_min` AUROC 介于 0.430 与 0.692 之间），并且在支架感知评价下，对接相对配体层化学基线只提供有限增量信息。主终点估计对最大 pChEMBL 与中位数聚合大体不敏感。PIK3CA/mTOR 给出最高点估计，并在未参与建面的配体池中保持正向方向信号；但主面板估计的不确定度及其对受体结构的敏感性，排除将其解释为可迁移的双靶决策规则。

更广泛的分析表明，表观双靶对接性能由任务定义、配体化学组成和受体实现方式共同决定。EGFR/HER2 上，Dual-versus-neither comparator（AUROC 0.756）会支持对接双靶识别，而方向性弱臂仍为 0.430；该 formulation gap 依赖靶对，不是四对统一反转，也是描述性对照而非配对显著性检验。若干靶对上，物化或化学型参考达到或超过对接判别；未使用配体池还暴露出错口袋对照不低于口袋匹配的未解决反转，效价或尺寸匹配未能消除，尽管相应配对置信区间仍包含零。受体实现对可以改变表观判别的幅度甚至方向；这是 realization effect，不是稳健性证书。这些发现反对把两个口袋上的有利分数当作双靶活性的充分证据。双靶虚拟筛选应纳入实验定义的单靶硬负样本、配体层混淆对照、面板外配体评价以及受体敏感性分析。因此，DualFourClass-Bench 的主要贡献不是产生一个普适的 docking winner，而是提供系统协议，用以界定 docking-based dual-target recognition 的证据与可靠性边界。

---

## 术语分层（全文只在 Conclusion 用 grounded）

| 章节 | 推荐用语 |
|------|----------|
| Introduction | experimentally defined dual-target recognition task |
| Methods | experimentally derived activity labels |
| Results | experimentally defined hard negatives |
| Conclusion | experimentally grounded evaluation setting（仅一次） |

禁止：validated；robust performance；generalizable dual-target docking strategy；docking is ineffective；docking can identify dual-target ligands。
<!-- END CONCLUSIONS_DRAFT_ZH_JCIM_V1.md -->

---

