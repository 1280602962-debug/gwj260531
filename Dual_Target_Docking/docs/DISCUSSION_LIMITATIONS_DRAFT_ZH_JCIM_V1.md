# Discussion · Limitations（中文工作稿片段）

> 位置：Discussion 末段或 Discussion 下 *Limitations* 小节（JCIM 常见）。  
> 仅列仓库可支撑的局限；不编造未做实验。与 Methods 2.1、2.5、Results 3.2 / 3.4–3.6 交叉引用。

---

本文有若干局限需明确。第一，评价集仅含四对靶标，其中 EGFR/HER2 为供给受限案例而非严格厚面板；结论不得外推为“全部双靶对接任务”的通论。第二，活性标签取自 ChEMBL 的**最大** pChEMBL 聚合，存在 activity inflation：同一配体–靶标若跨 assay 波动，取最大会抬高表观活性并可能扩大 dual 类。冻结缓存未保留逐条 assay 的中位数、置信度或物种字段，因而**现有数据包无法**完成 max 对 median、confidence≥8 或 Homo sapiens 过滤的敏感性重算（`T0_SKIPS.md`）；该比较需重新拉取 assay 级记录后方可补做，本稿不虚构其结果。第三，面板建造规则在靶对间不完全统一（严格规则与 θ = 6.0），存在“为凑供给而选阈值”的观感；我们因此将**统一标签重标**（Table S4）作为跨对主稳健分析，建造规则仅作 construction readout，且统一规则下 EGFR/HER2 与 PIK3CA/mTOR 的 B_only 过少时标记 underpowered。第四，尚未从严格供给池按同一配额重抽大量独立 panel 以报告 summary_min 的抽样分布——现有配体/支架 bootstrap 描述的是固定面板内的不确定度，不是 panel composition 重抽；该分析依赖对未用池分子的对接，列入后续工作。第五，RTMScore 取 Vina 九姿态最优，而 GNINA 仅对 mode 1 rescore，通道对照不对称，不能解释为三引擎公平赛马。第六，受体准备采用 meeko 默认路径，未另做 PDBFixer 补全或 Reduce 质子化枚举，组氨酸互变异构与辅因子处理边界见 Methods 2.4。第七，PIK3CA/mTOR 相对重原子数基线的优势在点估计上存在，但其 Δ 的 95% 置信区间仍包含 0，不足以支持“显著优于平凡描述符”的主张。第八，支架分组交叉验证中多数支架接近单例，对同系物泄漏的压制有限；二维指纹基线高 AUROC 表明标签与化学型相关，不能单独证明口袋物理特异性。第九，本研究为计算评测，未包含湿实验验证。
