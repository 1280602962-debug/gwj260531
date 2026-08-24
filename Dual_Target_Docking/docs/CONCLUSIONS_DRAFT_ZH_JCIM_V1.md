# Conclusions（中文工作稿 · JCIM Articles）

## 5. 结论

本研究建立 DualFourClass-Bench，作为针对 A-selective、B-selective 硬负样本的实验锚定计算判别环境。在四对冻结靶标上，最差方向判别摘要高度依赖靶对（`summary_min` 0.430–0.692），且把 docking 加到 ECFP4 后只产生很小的支架分组 CV AUROC 增量改善。在这四对面板中，主终点对最大 pChEMBL 与重复测定中位数聚合的选择总体不敏感。PIK3CA/mTOR 的点估计最高，并在未使用配体池稳定性检查中保持同向趋势；但其不确定性与受体敏感性排除了可推广的决策规则。

更广泛的分析表明，表观双靶对接判别取决于任务设定、配体化学组成和受体实现。EGFR/HER2 上，Dual-versus-neither comparator 为 0.756，而方向性最差一臂为 0.430；这是依赖靶对的描述性对照，不是配对显著性检验或四对定律。若干靶对上，ligand-only reference 达到或超过 docking；unused-pool holdout 还暴露出配对区间包含零的未解决 wrong-pocket point-estimate reversal。受体实现改变了表观判别的幅度和方向。这些结果支持把实验定义的选择性硬负样本、配体层混淆对照、面板外稳定性检查与受体敏感性作为双靶对接评价的互补要求。DualFourClass-Bench 的主要贡献是界定证据与可靠性边界的系统协议，而不是普适 docking winner。
