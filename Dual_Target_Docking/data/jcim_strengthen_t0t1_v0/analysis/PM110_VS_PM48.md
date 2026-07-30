# PM110 vs PM48 — pocket-matched (Vina / RTM / GNINA)

> 主指标：`summary_min = min(AUROC D vs A on pocket B, AUROC D vs B on pocket A)`。
> Vina：−affinity；RTM：raw（越高越好）；GNINA：CNNscore（越高越好）。
> PM48 carryover 复用既有姿态；PM110_* 新对接后 RTM best-of-9 + GNINA mode_01。

| Panel | Arm | n(D/A/B) | D↔A (B) | D↔B (A) | **summary_min** | 95% CI |
|-------|-----|----------|---------|---------|-----------------|--------|
| PM48 | vina | 18/14/12 | 0.7143 | 0.6921 | **0.6921** | [0.48, 0.80] |
| PM48 | rtm | 18/14/12 | 0.6151 | 0.6574 | **0.6151** | [0.39, 0.76] |
| PM48 | gnina | 18/14/12 | 0.5794 | 0.6713 | **0.5794** | [0.36, 0.75] |
| PM110 | vina | 30/30/30 | 0.7311 | 0.6483 | **0.6483** | [0.51, 0.76] |
| PM110 | rtm | 30/30/30 | 0.5756 | 0.7433 | **0.5756** | [0.43, 0.72] |
| PM110 | gnina | 30/30/30 | 0.5222 | 0.7133 | **0.5222** | [0.38, 0.66] |

**Δ Vina (PM110−PM48)**: -0.0438
**Δ RTM**: -0.0395
**Δ GNINA**: -0.0571

## 结论

扩面后 Vina summary_min 与 PM48 接近（|Δ|<0.05）；CI 收窄。声称仍受 CLAIM_CEILING 约束。
RTM/GNINA 不改变「非通用决策臂」上限；若 RTM/GNINA summary_min 不高于 Vina，正文以 Vina 口袋匹配为主报告。
