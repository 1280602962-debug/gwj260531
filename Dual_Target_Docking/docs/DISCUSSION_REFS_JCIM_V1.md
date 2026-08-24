# Discussion 引用核验（JCIM 工作稿）

> Introduction 编号 1–11 见 [`INTRODUCTION_REFS_JCIM_V1.md`](INTRODUCTION_REFS_JCIM_V1.md)。本文件从 12 起收录 Discussion 新增条目。  
> 只写已核验的题名、作者、年卷页/DOI，并写明允许引用的句子边界。未核到 ACS Cite-this 卷页的条目用 DOI + 发表日期，不编造页码。

---

## 著录

12. Tran-Nguyen, V.-K.; Ballester, P. J. Beware of Simple Methods for Structure-Based Virtual Screening: The Critical Importance of Broader Comparisons. *J. Chem. Inf. Model.* **2023**, *63*, 1401–1405. DOI: [10.1021/acs.jcim.3c00218](https://doi.org/10.1021/acs.jcim.3c00218).

13. Ahmed, F.; Soellner, M. B.; Brooks, C. L. Real-World Assessment of Machine-Learned Docking Using Bioassay-Derived Benchmarks. *J. Chem. Inf. Model.* **2026**. DOI: [10.1021/acs.jcim.5c03020](https://doi.org/10.1021/acs.jcim.5c03020). (ACS publication date 2026-07-21; volume/pages not captured from Cite-this at draft time. Zenodo supporting receptors: 10.5281/zenodo.20030384.)

14. Schaller, D. A.; Christ, C. D.; Chodera, J. D.; Volkamer, A. Benchmarking Cross-Docking Strategies in Kinase Drug Discovery. *J. Chem. Inf. Model.* **2024**, *64*, 8848–8858. DOI: [10.1021/acs.jcim.4c00905](https://doi.org/10.1021/acs.jcim.4c00905). (ACS Cite this: **2024**, *64*, 23, 8848–8858.)

15. Sindt, F.; Bret, G.; Rognan, D. On the Difficulty to Rescore Hits from Ultralarge Docking Screens. *J. Chem. Inf. Model.* **2025**, *65*, 5553–5566. DOI: [10.1021/acs.jcim.5c00730](https://doi.org/10.1021/acs.jcim.5c00730). (HAL Cite this: *65*(11), 5553–5566.)

Discussion 复用 Introduction 条目时不改号：

- 5–7：DUD / DUD-E / LIT-PCBA（decoy 构建与化学偏倚）。
- 9：Zhou, Li, Hou 2013 双靶激酶对接评价（inhibitor vs noninhibitor；结构依赖；高 FP）。不是 Dual vs A-only/B-only。
- 10–11：DualDiff / FuseDiff（生成式双靶评测口径；本文未重对接其生成物）。

---

## 每条允许写什么 / 禁止写什么

| # | 正文位置 | 该文实际主张（核验） | Discussion 允许的用法 | 禁止的用法 |
|---|----------|----------------------|----------------------|------------|
| 12 | 4.1–4.2 | 数据 unbiasing 与简单方法（如 IFP）可能高估 VS；需要比 AVE / 少数 generic ML SF 更宽的比较；IFP 可过拟合回顾性基准。 | ligand-distribution / chemical-composition signal 可以产生表观 VS 性能；双靶评测必须把 ligand-level baseline 写进主文。 | 声称 DualFourClass 使用了 IFP、AVE 或 DEKOIS；把该 Viewpoint 写成双靶论文。 |
| 13 | 4.1 | DiffDock-Pocket 在 bioassay/HTS 衍生基准上的真实世界评价；配套 Zenodo 含 LIT-PCBA / DUD-E 受体准备。强调 constructed decoy 与真实 assay 数据之间的差异。 | assay-derived benchmark 可以揭示 constructed decoy 集掩盖的局限。 | 声称本文评价了 DiffDock-Pocket；声称 DualFourClass 使用了 Ahmed 等的配体/受体；编造卷页。 |
| 14 | 4.3 | 激酶 cross-docking：姿态回收依赖于所选激酶结构与对接方法；多结构对接提高低 RMSD 姿态机会。使用 OpenEye Fred / Hybrid / Posit，不是 Vina。 | receptor representation 是独立性能变量；pose QC ≠ screening 可迁移。 | 声称本文用了 OpenEye 引擎；把其姿态回收成功率写成 DualFourClass 的结果；写成双靶 docking 论文。 |
| 15 | 4.4 | 十个已成功的超大规模对接 hit list 上，八种再打分方法均不能在全部 assay 上稳健区分已知结合分子与无活性分子。 | 单靶 VS 在 hit 选择/后处理阶段仍有可靠性边界；双靶更不能只报两端有利分数。 | 声称 DualFourClass 做了超大规模对接或再打分 bake-off；把 Sindt 的失败写成 Vina 双靶失败的直接证据。 |

---

## 有意不引入 Discussion 的文献

| 候选 | 原因 |
|------|------|
| Cieplinski et al., *J. Chem. Inf. Model.* **2023**（generative models should dock well） | 生成分子对接成功率基准，不是实验 dual vs A-only/B-only；Intro 已有 DualDiff/FuseDiff。 |
| Hall 2025 超大规模对接结果库 | 不是双靶硬负评测。 |
| Vu 2025 DUDE-Z 姿态/打分 | 单靶 VS；易被读成同类双靶评测。 |
