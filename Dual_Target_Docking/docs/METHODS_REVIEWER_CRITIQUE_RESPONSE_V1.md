# Methods 审稿式批评响应（内部工作笔记）

> 对应用户转述的 8 条 Methods 质疑。原则：**有则改之，无则不编造未做实验**。

| # | 质疑 | 是否合理 | 处理 |
|---|------|----------|------|
| 1 | 标签规则像数据驱动调整 | **合理** | 提升统一重标为跨对主稳健分析；建面板规则降为 construction protocol（见 Methods 2.1 修订）。**不**把未计算过的“hybrid 6.0/5.5” invent 成主表。 |
| 2 | max pChEMBL 膨胀 | **合理** | 已有 Limitation；加强措辞。现缓存无 median 字段（`T0_SKIPS.md`），**禁止编造** max vs median 表；本地可补抓 ChEMBL assay 后再做。 |
| 3 | 缺 panel 重抽样稳定性 | **合理** | 现有配体/支架 bootstrap ≠ 从供给池重抽 1000 个 panel。需对接未用池才能做；列入待做，不伪造分布。 |
| 4 | score 方向要写死 | **合理** | Methods 2.6 明确 \(S=-\Delta G_{\mathrm{Vina}}\)，dual=positive。 |
| 5 | GNINA 只打 mode1 不公平 | **合理** | Methods/Limitations 写明不对称。2026-08-21 已核实：本环境无 gnina、主面板无 9 姿态坐标，**不能**在此补打（`GNINA_NINE_POSE_SKIP_V1.md`）。**不**假装已公平。 |
| 6 | 蛋白准备不足 | **部分合理** | 按**实际** meeko 路径补写；**不**写入未跑的 PDBFixer/Reduce。缺 histidine/辅因子专项枚举写入局限。 |
| 7 | decoy 匹配窗口无引用 | **合理** | 引用 DUD-E / property-matched decoy 惯例并写明窗口。 |
| 8 | Methods 过长 | **合理** | cognate 数字表压到 SI；正文只留门槛与结论句。 |
