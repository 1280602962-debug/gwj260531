# PM110 vs PM48 — pocket-matched (Vina / RTM / GNINA)

> 主指标：`summary_min = min(AUROC D vs A on pocket B, AUROC D vs B on pocket A)`。
> Vina：−affinity；RTM：raw（越高越好）；GNINA：CNNscore（越高越好）。
> PM48 carryover 复用既有姿态；PM110_* 新对接后 RTM best-of-9 + GNINA mode_01（原始版本，见下方更新）。
> **2026-08-24 更新：** GNINA 已补做全 9 姿态公平重打（`data/jcim_bench_v0/analysis/GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`），下表 GNINA 行改为 best-of-9 数值；mode_01 历史值移至表后单独列出，不再作为正文引用值。

| Panel | Arm | n(D/A/B) | D↔A (B) | D↔B (A) | **summary_min** | 95% CI |
|-------|-----|----------|---------|---------|-----------------|--------|
| PM48 | vina | 18/14/12 | 0.7143 | 0.6921 | **0.6921** | [0.48, 0.80] |
| PM48 | rtm | 18/14/12 | 0.6151 | 0.6574 | **0.6151** | [0.39, 0.76] |
| PM48 | gnina (best9) | 18/14/12 | 0.6548 | 0.6852 | **0.6548** | [0.4277, 0.8087] |
| PM110 | vina | 30/30/30 | 0.7311 | 0.6483 | **0.6483** | [0.51, 0.76] |
| PM110 | rtm | 30/30/30 | 0.5756 | 0.7433 | **0.5756** | [0.43, 0.72] |
| PM110 | gnina (best9) | 30/30/30 | 0.6133 | 0.6822 | **0.6133** | [0.4613, 0.7361] |

**GNINA mode_01 历史值（供追溯，不再作正文引用）**

| Panel | D↔A (B) | D↔B (A) | summary_min | 95% CI |
|-------|--------:|--------:|-------------:|--------|
| PM48 | 0.5794 | 0.6713 | 0.5794 | [0.361, 0.746] |
| PM110 | 0.5222 | 0.7133 | 0.5222 | [0.3777, 0.6611] |

**Δ Vina (PM110−PM48)**: -0.0438
**Δ RTM**: -0.0395
**Δ GNINA (best9)**: -0.0415

## 结论

扩面后 Vina summary_min 与 PM48 接近（|Δ|<0.05）；CI 收窄。声称仍受 CLAIM_CEILING 约束。
RTM/GNINA 不改变「非通用决策臂」上限；即便 GNINA 换用 best-of-9 公平重打后数值上升（PM48 0.579→0.655；PM110 0.522→0.613），仍不高于同面板 Vina 口袋匹配（PM48 0.692；PM110 0.648），正文以 Vina 口袋匹配为主报告。
