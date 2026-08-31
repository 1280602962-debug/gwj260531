# Introduction（中文工作稿 · JCIM Articles）

## 1. 引言

多靶点药物设计（multitarget drug design）旨在通过单一小分子同时调控两个或多个生物学靶点，以应对复杂疾病中的通路冗余、代偿性信号以及药物耐药。与传统单靶点药物相比，合理设计的多靶点配体有望通过协同调节相互关联的生物学过程获得更充分的药理效应，因此已成为多药理学（polypharmacology）研究的重要方向。[^1] 近年来，多靶点小分子的理性设计逐渐由经验性筛选转向结合结构生物学、计算化学与生成式模型的结构导向设计。[^2] 分子对接（molecular docking）仍是结构基础虚拟筛选（structure-based virtual screening, SBVS）中最常用的计算工具之一：先预测配体在蛋白结合口袋中的构象，再用打分函数对配体–受体互补性排序。[^3][^4] 因此，在双靶点药物发现中，一个自然策略是分别将候选分子对接至两个靶点，并据此判断其是否具有潜在双靶结合能力。对接结果的解释高度依赖数据集构建。DUD 与 DUD-E 使用物化性质匹配 decoy，以避免表观富集退化为粗粒度配体性质分离。[^5][^6] LIT-PCBA 采用实验 assay 标签并控制已知 decoy 与化学偏倚。[^7] CASF-2016 评价复合物上的 scoring、ranking、docking 与 screening power，仍然属于单复合物问题。[^8] 这些资源都没有在实验标注的四状态配体空间中定义双靶方向判别。

一个严格的双靶评价需要区分 **dual-active**、**A-selective**、**B-selective** 与 **neither** 四种实验状态（Figure 1A）。A-only 与 B-only 是该任务的**选择性硬负样本（selectivity hard negatives）**：它们在一个靶点上已有较强活性。计算终点因而检验对接能否在两个方向上将 dual-active 与对应单靶选择性配体区分开。Zhou、Li 与 Hou 曾在四对激酶上评价相对非抑制剂的 dual-target docking。[^9] 本文在该设定上引入实验定义的方向性硬负，并在同一套分数上比较不同基准设定。Dual versus neither 按实验 inactive 计分，作为基准设定对照。平衡的四状态面板还受两端可比较测量与双向选择性硬负供给的约束。

近期工作已经把双靶问题做实，但回答的并不是同一评价问题。Wu 等表明，大规模对接可以前瞻性地得到选定靶对的联合结合配体，同时也报告后续优化仍然困难。[^19] 该研究问的是对接能否找到双靶活性分子；本文问的是：当负类由实验选择性配体而不是非结合配体定义时，回顾性双靶识别证据是否改变。POLYGON 在十万级结合数据上生成双靶化学空间，并合成 32 个 MEK1/mTOR 化合物做实验验证，[^20] 因此不宜把近期生成式多药理学一概写成对接成功指标。生成式双靶方法也使用相对参考配体的对接度量。[^10][^11] Kinase-Bench 汇集 6875 个选择性配体、75 个激酶与 422,799 个 decoy，检验相对激酶特异性 decoy 的选择性富集；[^22] DualFourClass 则在每一对上直接构造实验测定的 dual、A-only、B-only 与 neither 四状态。一项 147 靶的 AI 对接基准进一步表明，方法排序取决于负类是实验低活性 TrueDecoy 还是商业库随机 decoy，[^21] 说明配方本身就是科学主张的一部分。

本文要问的是：基准设定是否改变双靶识别的表观证据。我们构建 DualFourClass-Bench，作为具有两条方向主任务的四状态面板：dual 对 A-only 在口袋 B 打分，dual 对 B-only 在口袋 A 打分（Figure 1B），并以二者较弱一臂汇总为最弱臂 AUROC（`summary_min`）。我们进一步考察该判别是否能够在不同配体、活性聚合方式及受体结构条件下保持。

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
[^19]: Wu, Vigneron, Braz et al., *J. Med. Chem.* **2026**, *69*, 6210–6229. DOI: 10.1021/acs.jmedchem.5c03810.
[^20]: Munson, Chen, Bogosian et al., *Nat. Commun.* **2024**, *15*, 3636. DOI: 10.1038/s41467-024-47120-y.
[^21]: Gu, Shen, Zhang et al., *Nat. Mach. Intell.* **2025**, *7*, 509–520. DOI: 10.1038/s42256-025-00993-0.
[^22]: Wei, Zhou, Jing et al., *J. Chem. Inf. Model.* **2024**, *64*, 9528–9550. DOI: 10.1021/acs.jcim.4c01830.
