# Introduction（中文工作稿 · JCIM Articles）

## 1. 引言

### 1.1 双靶设计及传统对接基准的不足

多靶点药物设计（multitarget drug design）旨在通过单一小分子同时调控两个或多个生物学靶点，以应对复杂疾病中的通路冗余、代偿性信号以及药物耐药等问题。与传统单靶点药物相比，合理设计的多靶点配体有望通过协同调节相互关联的生物学过程获得更充分的药理效应，因此已成为多靶点药物发现和多药理学（polypharmacology）研究的重要方向。[^1] 近年来，多靶点小分子的理性设计逐渐由经验性的多药理筛选，转向结合结构生物学、计算化学与生成式模型的结构导向设计。[^2]

在这一过程中，分子对接（molecular docking）仍是结构基础虚拟筛选（structure-based virtual screening, SBVS）中最常用的计算工具之一：先预测配体在蛋白结合口袋中的结合构象，再用打分函数对配体–受体相互作用排序，从而在大规模化合物库中给出结构互补性的近似评价。[^3][^4] 因此，在双靶点药物发现中，一个自然的计算策略是分别将候选分子对接至两个靶点，并据此判断其是否具有潜在的双靶结合能力。

既有虚拟筛选研究已经表明，对接结果的解释高度依赖数据集构建、负样本定义、化学偏倚和评价指标。DUD 与 DUD-E 使用物化性质匹配 decoy，以避免表观 enrichment 退化为粗粒度配体性质分离。[^5][^6] LIT-PCBA 进一步采用实验 assay 标签，并系统控制已知 decoy/chemical biases，以提高虚拟筛选评价的现实性。[^7] CASF-2016 则评价复合物上的 scoring、ranking、docking 与 screening power，仍然属于单复合物问题。[^8] 这些资源都没有在实验标注的四状态配体空间中定义双靶方向判别。

然而，**将单靶点 docking 的评价逻辑直接扩展到双靶点任务并不充分。**

对于双靶点配体而言，任务结构发生了改变。一个严格的双靶点 benchmark 至少需要区分四种具有不同生物学含义的配体状态（四状态数据集，而不是四分类器）：同时作用于两个靶点的 **dual-active** 配体、仅作用于靶点 A 的 **A-selective** 配体、仅作用于靶点 B 的 **B-selective** 配体，以及两个靶点均缺乏足够活性的 **neither** 配体（Figure 1A）：

|  | *B*<sup>+</sup> | *B*<sup>−</sup> |
|--|:--:|:--:|
| *A*<sup>+</sup> | Dual | A-only |
| *A*<sup>−</sup> | B-only | Neither |

其中，A-only 和 B-only 不是普通负样本，而是**选择性硬负样本（selectivity hard negatives）**。它们在一个靶点上已有较强活性，却在另一靶点缺乏相应活性。计算终点因而检验 docking 能否在两个方向上将 dual-active 与对应单靶选择性配体区分开；这一判别本身不证明独立的 pocket-specific recognition 或生物学双靶活性。

Zhou、Li 与 Hou 曾在四对激酶上评价相对 noninhibitor 的 dual-target docking，并报告结构依赖性和预测 dual 中的 false positives。[^9] 本文在该评价设定基础上进一步引入实验定义的 A-only/B-only 方向硬负，并直接比较不同 formulation 下的表观判别。Dual-versus-neither 是 **nonselectivity-controlled comparator**，不是 “the conventional dual-target benchmark”。

池化两个口袋分数可能掩盖较弱方向，而两端都优于参考配体只定义计算成功，并不等价于实验双靶活性。因此 benchmark 与读出必须匹配四状态生物学空间。

从公开数据构建此类面板要求同一化合物在两个靶点上均有可比较测量，并要求两个方向都有足量选择性硬负。由于 assay 类型、条件和覆盖不同，能够支持平衡四状态评价的靶对数量本身就是数据供给问题，而不是预先存在的资源。

双靶面板还继承化学混淆：若 dual-active 与选择性配体的分子性质或支架分布不同，AUROC 可以反映 label-associated ligand distributions，而不只是口袋互补性。因此需要显式的物化与 ligand-only controls。[^7]

### 1.2 Docking-based 双靶设计使严格评价成为实际问题

DualDiff 与 FuseDiff 说明了这一差异的实际意义：两者均使用较差口袋分数和“生成分子在两个靶点均优于参考配体”的比例评价双靶设计。[^10][^11] 这些指标衡量相对参考配体的计算成功，而不是相对实验选择性配体的判别。因此生成式 docking metrics 与本文 hard-negative endpoint 是互补而非竞争基准；本文不重对接其生成分子，也不把其 reported metrics 重新解释为本文 primary endpoint 的竞争者。

### 1.3 研究目的与贡献

**本文要问的是：benchmark formulation 是否改变双靶识别的表观证据。** 我们构建实验定义的四状态面板，以针对 A-selective 与 B-selective 硬负的口袋匹配方向判别为主任务，并与 nonselectivity-controlled Dual-versus-neither comparator 比较；随后用 ligand-only、物化性质、wrong-pocket、配体池、活性聚合与受体实现对观察到的判别进行压力测试。

DualFourClass-Bench 是具有两条方向主任务的 curated benchmark：dual 对 A-only 在口袋 B 打分，dual 对 B-only 在口袋 A 打分（Figure 1B）。neither 描述完整实验空间，但不进入 primary AUROC。`summary_min` 是保守的最差方向判别摘要，不是新 docking score。公开数据审计决定多少候选靶对能够支持该构建；评价集规模是审计结果，而非预设目标。

本文贡献是评价协议与资源，而不是新 docking algorithm。它把任务设定、混淆感知评价和评价条件敏感性连接起来，检验表观判别究竟是 docking 的固定属性，还是 benchmark 条件下的条件性结果。该协议用于双靶虚拟筛选与生成设计的下游校验，不是与现有生成模型竞赛，也不是 comprehensive dual-target suite。

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
