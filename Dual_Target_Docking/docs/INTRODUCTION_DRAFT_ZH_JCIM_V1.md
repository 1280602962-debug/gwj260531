# Introduction（中文工作稿 · JCIM Articles）

## 1. 引言

多靶点药物设计（multitarget drug design）旨在通过单一小分子同时调控两个或多个生物学靶点，以应对复杂疾病中的通路冗余、代偿性信号以及药物耐药。与传统单靶点药物相比，合理设计的多靶点配体有望通过协同调节相互关联的生物学过程获得更充分的药理效应，因此已成为多药理学（polypharmacology）研究的重要方向。[^1] 近年来，多靶点小分子的理性设计逐渐由经验性筛选转向结合结构生物学、计算化学与生成式模型的结构导向设计。[^2] 分子对接（molecular docking）仍是结构基础虚拟筛选（structure-based virtual screening, SBVS）中最常用的计算工具之一：先预测配体在蛋白结合口袋中的构象，再用打分函数对配体–受体互补性排序。[^3][^4] 因此，在双靶点药物发现中，一个自然策略是分别将候选分子对接至两个靶点，并据此判断其是否具有潜在双靶结合能力。对接结果的解释高度依赖数据集构建。DUD 与 DUD-E 使用物化性质匹配 decoy，以避免表观富集退化为粗粒度配体性质分离。[^5][^6] LIT-PCBA 采用实验 assay 标签并控制已知 decoy 与化学偏倚。[^7] CASF-2016 评价复合物上的 scoring、ranking、docking 与 screening power，仍然属于单复合物问题。[^8] 这些资源都没有在实验标注的四状态配体空间中定义双靶方向判别。

一个严格的双靶评价需要区分四种具有不同生物学含义的配体状态：同时作用于两个靶点的 **dual-active** 配体、仅作用于靶点 A 的 **A-selective** 配体、仅作用于靶点 B 的 **B-selective** 配体，以及两端均缺乏足够活性的 **neither** 配体（Figure 1A）。A-only 与 B-only 是该任务的**选择性硬负样本（selectivity hard negatives）**：它们在一个靶点上已有较强活性，却在另一靶点缺乏相应活性。计算终点因而检验对接能否在两个方向上将 dual-active 与对应单靶选择性配体区分开。Zhou、Li 与 Hou 曾在四对激酶上评价相对非抑制剂的 dual-target docking，并报告结构依赖性和预测 dual 中的假阳性。[^9] 本文在该设定上引入实验定义的方向性硬负，并在同一套分数上比较不同基准设定。Dual versus neither 用作非选择性对照（nonselectivity-controlled comparator）。双靶面板还可继承化学混淆：若 dual-active 与选择性配体的分子性质或支架分布不同，AUROC 可以反映与标签相关的配体分布，而不只是口袋互补性。[^7] 从公开数据构建此类面板还要求同一化合物在两端均有可比较测量，并要求两个方向都有足量选择性硬负，因此平衡的四状态评价本身受数据供给约束。

近期双靶生成方法仍使用相对参考配体的对接成功指标。[^10][^11] 当对接被用作下游过滤器时，本文的评价问题便具有实际意义。

本文要问的是：基准设定是否改变双靶识别的表观证据。我们构建 DualFourClass-Bench，作为具有两条方向主任务的四状态面板：dual 对 A-only 在口袋 B 打分，dual 对 B-only 在口袋 A 打分（Figure 1B），并以保守的最差方向判别摘要（`summary_min`）汇总。随后从任务设定、配体层混淆以及配体、活性标签与受体扰动三个层面评价观察到的判别。本研究评价的是对接排序的可靠性，而不是开发新的对接算法。

---

[^1]: Anighoro, Bajorath, Rastelli, *J. Med. Chem.* **2014**, *57*, 7874–7887. DOI: 10.1021/jm5006463.
[^2]: Proschak, Stark, Merk, *J. Med. Chem.* **2019**, *62*, 420–444. DOI: 10.1021/acs.jmedchem.8b00760.
[^3]: Kitchen et al., *Nat. Rev. Drug Discov.* **2004**, *3*, 935–949. DOI: 10.1038/nrd1549.
[^4]: Eberhardt et al., *J. Chem. Inf. Model.* **2021**, *61*, 3891–3898. DOI: 10.1021/acs.jcim.1c00203.
[^5]: Huang, Shoichet, Irwin, *J. Med. Chem.* **2006**, *49*, 6789–6801. DOI: 10.1021/jm0608356. (DUD)
[^6]: Mysinger et al., *J. Med. Chem.* **2012**, *55*, 6582–6594. DOI: 10.1021/jm300687e. (DUD-E)
[^7]: Tran-Nguyen, Jacquemard, Rognan, *J. Chem. Inf. Model.* **2020**, *60*, 4263–4273. DOI: 10.1021/acs.jcim.0c00155. (LIT-PCBA)
[^8]: Su et al., *J. Chem. Inf. Model.* **2019**, *59*, 895–913. DOI: 10.1021/acs.jcim.8b00545. (CASF-2016)
[^9]: Zhou, Li, Hou, *J. Chem. Inf. Model.* **2013**, *53*, 982–996. DOI: 10.1021/ci400065e.
[^10]: Zhou, Guan et al., *The Thirty-eighth Annual Conference on Neural Information Processing Systems (NeurIPS 2024)*; arXiv:2410.20688. DualDiff.
[^11]: Wu et al., *Proceedings of the 32nd ACM SIGKDD Conference on Knowledge Discovery and Data Mining, Vol. 2*; ACM: New York, 2026; pp 12432–12443. DOI: 10.1145/3770855.3819050. FuseDiff.
