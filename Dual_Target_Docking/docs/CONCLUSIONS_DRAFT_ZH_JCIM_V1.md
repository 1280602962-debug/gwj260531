# Conclusions（中文工作稿 · JCIM Articles）

## 5. 结论

本研究建立了 DualFourClass-Bench，用于在实验定义的 A-selective 和 B-selective 硬负样本条件下评价双靶点 docking 的方向性判别。四对冻结靶标中，`summary_min` 为 0.430–0.692，且将 docking 分数加入 ECFP4 后在支架分组交叉验证中的增量有限。活性聚合规则的敏感性分析未改变主要结论，而受体结构替换则可显著改变判别表现，包括在两个 PIK3CA 相关靶对中出现方向相反的变化。

这些结果表明，双靶点 docking 的表观判别具有明显的任务、配体化学和受体实现依赖性；仅依据两个口袋中的有利 docking 分数不足以建立双靶活性的充分证据。对于双靶虚拟筛选及将 docking 用作下游评价的生成式设计流程，应同时考虑选择性硬负样本、配体层化学对照、未使用配体池和受体结构敏感性。DualFourClass-Bench 的主要用途因此是界定 docking-based dual-target recognition 的证据边界，而非提供普适的 docking 排名规则。
