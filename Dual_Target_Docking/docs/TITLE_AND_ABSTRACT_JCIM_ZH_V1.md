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
