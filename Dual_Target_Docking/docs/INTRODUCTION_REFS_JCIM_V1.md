# Introduction 引用核验（JCIM 工作稿共用）

> 中文稿 [`INTRODUCTION_DRAFT_ZH_JCIM_V1.md`](INTRODUCTION_DRAFT_ZH_JCIM_V1.md) 与英文稿 [`INTRODUCTION_SECTION_JCIM_EN_V1.md`](INTRODUCTION_SECTION_JCIM_EN_V1.md) 使用同一套编号。  
> 本文件只收录**已核对题名、作者、年份、卷页与 DOI/arXiv**的条目，并写明允许引用的句子边界。未核验或与本文任务不符的文献不进入 Introduction。

编号按正文首次出现顺序。ACS 期刊用 *J. Med. Chem.* / *J. Chem. Inf. Model.* 的 Cite-this 年卷页。

---

## 著录（投稿用）

1. Anighoro, A.; Bajorath, J.; Rastelli, G. Polypharmacology: Challenges and Opportunities in Drug Discovery. *J. Med. Chem.* **2014**, *57*, 7874–7887. DOI: [10.1021/jm5006463](https://doi.org/10.1021/jm5006463).

2. Proschak, E.; Stark, H.; Merk, D. Polypharmacology by Design: A Medicinal Chemist’s Perspective on Multitargeting Compounds. *J. Med. Chem.* **2019**, *62*, 420–444. DOI: [10.1021/acs.jmedchem.8b00760](https://doi.org/10.1021/acs.jmedchem.8b00760). (ACS Cite this: **2019**, *62*, 420–444; published online 2018-07-23.)

3. Kitchen, D. B.; Decornez, H.; Furr, J. R.; Bajorath, J. Docking and Scoring in Virtual Screening for Drug Discovery: Methods and Applications. *Nat. Rev. Drug Discov.* **2004**, *3*, 935–949. DOI: [10.1038/nrd1549](https://doi.org/10.1038/nrd1549).

4. Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. Model.* **2021**, *61*, 3891–3898. DOI: [10.1021/acs.jcim.1c00203](https://doi.org/10.1021/acs.jcim.1c00203). (Europe PMC: *61*(8), 3891–3898.)

5. Huang, N.; Shoichet, B. K.; Irwin, J. J. Benchmarking Sets for Molecular Docking. *J. Med. Chem.* **2006**, *49*, 6789–6801. DOI: [10.1021/jm0608356](https://doi.org/10.1021/jm0608356). (Directory of Useful Decoys, DUD.)

6. Mysinger, M. M.; Carchia, M.; Irwin, J. J.; Shoichet, B. K. Directory of Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking. *J. Med. Chem.* **2012**, *55*, 6582–6594. DOI: [10.1021/jm300687e](https://doi.org/10.1021/jm300687e).

7. Tran-Nguyen, V.-K.; Jacquemard, C.; Rognan, D. LIT-PCBA: An Unbiased Data Set for Machine Learning and Virtual Screening. *J. Chem. Inf. Model.* **2020**, *60*, 4263–4273. DOI: [10.1021/acs.jcim.0c00155](https://doi.org/10.1021/acs.jcim.0c00155). (ACS title uses American spelling *Data Set*.)

8. Su, M.; Yang, Q.; Du, Y.; Feng, G.; Liu, Z.; Li, Y.; Wang, R. Comparative Assessment of Scoring Functions: The CASF-2016 Update. *J. Chem. Inf. Model.* **2019**, *59*, 895–913. DOI: [10.1021/acs.jcim.8b00545](https://doi.org/10.1021/acs.jcim.8b00545). (Test set: 285 protein–ligand complexes with crystal structures and binding constants, compiled as a scoring benchmark from PDBbind-quality complexes.)

9. Zhou, X.; Guan, J.; Zhang, Y.; Peng, X.; Wang, L.; Ma, J. Reprogramming Pretrained Target-Specific Diffusion Models for Dual-Target Drug Design. In *The Thirty-eighth Annual Conference on Neural Information Processing Systems (NeurIPS 2024)*; 2024. arXiv: [2410.20688](https://arxiv.org/abs/2410.20688). (Author-provided booktitle from the official DualDiff repository; the method name DualDiff appears **inside** this paper; do not invent a separate paper titled “DualDiff”.)

10. Wu, J.; Qiao, A.; Wang, Z.; Wei, Z.; Chen, S. FuseDiff: Symmetry-Preserving Joint Diffusion for Dual-Target Structure-Based Drug Design. arXiv: [2603.05567](https://arxiv.org/abs/2603.05567), **2026**. (Preprint; not a journal article at the time of this draft.)

---

## 每条允许写什么 / 禁止写什么

| # | 正文位置 | 该文实际主张（核验） | Introduction 允许的用法 | 禁止的用法 |
|---|----------|----------------------|-------------------------|------------|
| 1 | §1 多靶/双靶治疗动机 | 复杂疾病可能需要同时调控网络中多个节点；多靶药相对单靶或联用有潜在优势；并讨论计算设计多靶配体。 | 通路冗余/复杂疾病下的多靶药理动机。 | 当作双靶 docking benchmark；堆具体适应证病例。 |
| 2 | §1 理性多靶设计 | 多靶配体识别长期依赖偶然；现可借助计算、晶体学、片段设计等工具进行理性设计。 | 由经验筛选走向结构/计算导向设计。 | 声称该文评价了对接双靶识别。 |
| 3 | §1 docking / SBVS | 对接把小分子放入大分子结合位点，打分函数评价互补性，广泛用于 hit identification 与 lead optimization。 | docking 仍是 SBVS 常用工具。 | 暗示 Kitchen 2004 已覆盖四类双靶任务。 |
| 4 | §1 当代对接引擎 | 作者称 AutoDock Vina 为使用最广泛的开源对接程序之一；1.2.0 扩展力场与 Python 接口。本文主实验用 Vina 1.2.x（Methods）。 | 对接在虚拟筛选中的实际地位；后文实验引擎的文献锚点。 | 把 Vina 1.2.0 论文写成双靶评测。 |
| 5 | §1–§3 单靶 VS 基准 | DUD：decoy 应在物化性质上接近活性配体、拓扑上可区分，否则 enrichment 可能只是在分大小/极性等粗特征；相对未校正库，DUD 上的 enrichment 更低，说明偏倚会抬高表观性能。 | 虚拟筛选评价依赖负样本构建。 | 把 DUD decoy 说成本文的 A_only/B_only。 |
| 6 | §1–§3 物化匹配 decoy | DUD-E：102 个靶、ChEMBL 配体、每个配体 50 个 ZINC 物化匹配 decoy。 | 经典单靶 docking benchmark；负样本定义影响评价。 | 声称 DualFourClass-Bench 使用 DUD-E 分子。 |
| 7 | §1、§3 化学偏倚 | 明确指出 DUD / DUD-E / MUV 等人工 active/decoy 集存在明显与隐蔽化学偏倚，可能**高估**虚拟筛选真确度；LIT-PCBA 改用 PubChem 剂量–响应实验标签，并控制物化性质范围。 | benchmark 构建本身改变评价结论；双靶面板也需要化学/物化对照。 | 声称 LIT-PCBA 是四类双靶基准，或本文使用了 LIT-PCBA 数据。 |
| 8 | §1 结构打分基准 | CASF 把打分与对接过程解耦；四项指标为 scoring / ranking / docking / screening power；CASF-2016 测试集为 285 个高质量复合物。 | 已有结构评价体系针对单复合物打分、姿态与单靶筛选，不是 dual-versus-selective。 | 把 CASF screening power 等同于本文的分臂硬负 AUROC。 |
| 9 | §4 双靶生成背景 | DualDiff（及 CompDiff）为零样本组合单靶扩散模型的双靶生成方法。评测用 AutoDock Vina 重对接：P-1/P-2 Vina Dock、**Max Vina Dock**（两端 Vina Dock 的最大值，即较弱一端）、**Dual High Affinity**＝生成分子在**两个靶上的结合亲和力均优于各自参考配体**的比例。摘要提到双靶策略与克服肿瘤耐药的潜在价值。 | 生成式双靶设计把 docking-based 双成功当作实用评价；与实验硬负判别互补。 | 把 Dual High Affinity 写成均值池化；声称本文重跑了 DualDiff 生成分子；把论文题名写成 DualDiff。 |
| 10 | §4 最新双靶结构生成 | 端到端扩散：联合生成共享配体分子图与两个口袋特异结合构象。训练数据自建；**独立测试集为 DualDiff benchmark（DDF）**。报告 Vina Score/Min/Dock、Max Vina Dock 与 Dual High Affinity。摘要：polypharmacology、疗效与降低耐药。 | 说明“如何评价双靶结构质量”已成为实际问题；FuseDiff 仍用对接口径。 | 当作已发表期刊论文；声称 FuseDiff 做了实验 dual vs A_only/B_only；声称本文重对接了 FuseDiff 分子。 |

---

## 有意不引入 Introduction 的文献

| 候选 | 原因 |
|------|------|
| Hall et al., *J. Chem. Inf. Model.* **2025**（超大规模对接结果库） | 不是双靶识别或四类硬负基准；旧稿用来暗示“虚拟筛选规模”过宽。 |
| Vu et al., *J. Chem. Inf. Model.* **2025**, *65*, 4833–4843 | 单靶 VS（DUDE-Z）上的姿态采样 + 打分评价；Methods 语气对照可用，Intro 引用易被读成“同类双靶评测”。 |
| Trott & Olson, *J. Comput. Chem.* **2010**（原始 Vina） | 本文引擎文献锚点用 Eberhardt 1.2.0 已足够；Methods 可另引。 |
| 未核验的疾病个案双靶综述 | Intro 不堆适应证；Anighoro / Proschak 已覆盖药理动机。 |

---

## 与本文数据/主张的边界（写 Intro 时自检）

- **K = 4** 与靶对名单属于 Methods / Results，不在 Introduction 写死。
- 不把 ChEMBL 供给审计写成“数据不够所以只做了四对”；应写成严格四类 ground-truth 的 completeness 本身是 methodological bottleneck。
- DualDiff / FuseDiff 是应用场景与互补评测口径，不是竞争对象；本文**没有**对其生成分子重新对接。
- 贡献表述：systematic benchmarking / evaluation protocol + 资源名 **DualFourClass-Bench**；禁止 D-DRAF 与 “novel framework named …”。
