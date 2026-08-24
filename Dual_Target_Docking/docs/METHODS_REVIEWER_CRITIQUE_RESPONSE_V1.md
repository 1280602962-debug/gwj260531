# Methods 审稿式批评响应（内部工作笔记）

> 对应用户转述的 8 条 Methods 质疑。原则：**有则改之，无则不编造未做实验**。

| # | 质疑 | 是否合理 | 处理 |
|---|------|----------|------|
| 1 | 标签规则像数据驱动调整 | **合理** | 提升统一重标为跨对主稳健分析；建面板规则降为 construction protocol（见 Methods 2.1 修订）。**不**把未计算过的“hybrid 6.0/5.5” invent 成主表。 |
| 2 | max pChEMBL 膨胀 | **合理；A4 已关** | 主策展仍用 max。全面板 median：标签一致率 93.6% / 98.9% / 99.0% / 100%；pair-level `summary_min` 基本不动。现写 **controlled limitation**（pChEMBL 非 assay-equivalent），不是 unresolved validity threat。禁止把 27 配体诊断当 SI。 |
| 3 | 缺 panel 重抽样稳定性 | **合理，但“对接未用池即可做 1000 独立 panel”不成立** | 现有配体/支架 bootstrap ≠ 供给池重抽。主面板+holdout 后严格硬负剩余见 `C_CLASS_EXPERIMENT_NECESSITY_VERDICT_V1.md`：最多再抽 1 个不重叠 30/30/30（PM、AChE），PIK3CB/EGFR 为 0。Holdout 已是一次 unused-pool 检验；禁止把有放回 bootstrap 或未对接的 BindingDB 合并写成 1000-panel 分布。 |
| 4 | score 方向要写死 | **合理** | Methods 2.8 明确 \(S=-E_{\mathrm{Vina}}\)，dual=positive。 |
| 5 | GNINA 只打 mode1 不公平 | **合理，已解决** | 2026-08-21 曾核实本环境无法补打（`GNINA_NINE_POSE_SKIP_V1.md`）；2026-08-24 用户本地补做全 9 姿态公平重打并推送真实结果（`GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`）：口袋匹配 summary_min 变化 −0.04–+0.08，三对上不超过 Vina，PIK3CA/PIK3CB 上略高于 Vina（0.533 对 0.500，mode01 时已如此，二者均近随机，非 best9 新产生的现象），姿态覆盖不对称不是 GNINA 偏弱的主要原因。Methods 2.7、Limitations 已据此更新，不再列为未做项。数字在 Results，不在 Methods。 |
| 6 | 蛋白准备不足 | **部分合理** | 按**实际** meeko 路径补写；**不**写入未跑的 PDBFixer/Reduce。缺 histidine/辅因子专项枚举写入局限。 |
| 7 | decoy 匹配窗口无引用 | **合理** | 引用 DUD-E / property-matched decoy 惯例并写明窗口。 |
| 8 | Methods 过长 | **合理** | cognate 数字表压到 SI；正文只留门槛与结论句。 |

---

## 第二轮（协议化 / 与 Results 分离）

原则仍是：**有则改之，不编造、不重抽已冻结面板。**

| # | 意见 | 判断 | 处理 |
|---|------|------|------|
| R2.1 | Methods 混入 Results | **合理** | 49→4、cognate 谁失败、GNINA Δ、HOAP_028、holdout/换晶点估计全部移出 Methods |
| R2.2 | 四类任务 vs 三类 primary | **合理** | 统一为 four-state benchmark + two directional pairwise primary；资源名 DualFourClass 保留 |
| R2.3 | summary_min 不是 novel metric | **合理** | 写成 worst-direction discrimination summary / 任务约束 |
| R2.4 | best descriptor 选择偏倚 | **合理** | 四描述符全报；max 只作 descriptive strong baseline |
| R2.5 | logistic 是 covariate-adjusted | **合理** | Methods 2.9.4 |
| R2.6 | 统一 Murcko 重抽全部 panel | **合理但不做** | 重抽会改冻结 AUROC；AChE/PIK3CB 现改为 class quotas + deterministic shuffle，**删除** identifier 前缀叙述；事后报告 Murcko，不重抽 |
| R2.7 | structure robustness → sensitivity | **合理** | Methods 2.12；不把 PM 写成预设 positive case |
| R2.8 | holdout ≠ external validation | **已部分做到；加强** | unused-pool, panel-external |
| R2.9 | wrong-pocket = falsification | **合理** | Methods 2.9.1 |
| R2.10 | cognate = pose-generation QC | **合理** | Methods 2.5；mode-1 偏离的数字在 Table S3 |
| R2.11 | RTM/GNINA 非三引擎竞赛 | **已有；保持** | Methods 2.7 |
| R2.12 | contact_count 写入 Methods | **合理** | Methods 2.9.6 |

