# Methods 审稿式批评响应（内部工作笔记）

> 对应用户转述的 8 条 Methods 质疑。原则：**有则改之，无则不编造未做实验**。

| # | 质疑 | 是否合理 | 处理 |
|---|------|----------|------|
| 1 | 标签规则像数据驱动调整 | **合理** | 提升统一重标为跨对主稳健分析；建面板规则降为 construction protocol（见 Methods 2.1 修订）。**不**把未计算过的“hybrid 6.0/5.5” invent 成主表。 |
| 2 | max pChEMBL 膨胀 | **合理** | 已有 Limitation；加强措辞。现缓存无 median 字段（`T0_SKIPS.md`），**禁止编造** max vs median 表；本地可补抓 ChEMBL assay 后再做。 |
| 3 | 缺 panel 重抽样稳定性 | **合理** | 现有配体/支架 bootstrap ≠ 从供给池重抽 1000 个 panel。需对接未用池才能做；列入待做，不伪造分布。 |
| 4 | score 方向要写死 | **合理** | Methods 2.6 明确 \(S=-\Delta G_{\mathrm{Vina}}\)，dual=positive。 |
| 5 | GNINA 只打 mode1 不公平 | **合理，已解决** | 2026-08-21 曾核实本环境无法补打（`GNINA_NINE_POSE_SKIP_V1.md`）；2026-08-24 用户本地补做全 9 姿态公平重打并推送真实结果（`GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`）：口袋匹配 summary_min 变化 −0.04–+0.08，三对上不超过 Vina，PIK3CA/PIK3CB 上略高于 Vina（0.533 对 0.500，mode01 时已如此，二者均近随机，非 best9 新产生的现象），姿态覆盖不对称不是 GNINA 偏弱的主要原因。Methods 2.5/2.7、Limitations 已据此更新，不再列为未做项。 |
| 6 | 蛋白准备不足 | **部分合理** | 按**实际** meeko 路径补写；**不**写入未跑的 PDBFixer/Reduce。缺 histidine/辅因子专项枚举写入局限。 |
| 7 | decoy 匹配窗口无引用 | **合理** | 引用 DUD-E / property-matched decoy 惯例并写明窗口。 |
| 8 | Methods 过长 | **合理** | cognate 数字表压到 SI；正文只留门槛与结论句。 |
