# B 组裁决 — exhaustiveness / enrichment / PM110

> 包：`data/jcim_strengthen_t0t1_v0/`  
> 约束：CLAIM_CEILING；主指标 = 口袋匹配 `summary_min`

## 完成清单

| ID | 产出 | 状态 |
|----|------|------|
| B1 T1.1 | `EXHAUSTIVENESS_E8_VS_E16.md` + `scores_vina_E8_best.csv` | ✅ |
| B2 T1.3 | `SINGLE_TARGET_ENRICHMENT_SANITY.md` + `single_target_enrichment_v1.csv` | ✅ |
| B3 T1.4 | `pik3ca_mtor_panel110_rdkit_v0` Vina E=16 → RTM best-of-9 → GNINA mode_01 | ✅ |
| B3 对照 | `PM110_VS_PM48.md` + `pm110_vs_pm48_pocket_matched_v1.csv` | ✅ |

## 关键数字

| 对照 | 结果 |
|------|------|
| E16 vs E8 summary_min | 0.692 vs 0.660（Δ +0.032）→ exhaustiveness **不足以单独解释** PM |
| 单靶 enrichment | 4L23 AUROC **0.603**；4JT6 **0.629** → 弱辨别 |
| PM110 Vina | summary_min **0.648** [0.51, 0.76]（PM48 0.692 [0.48, 0.80]；Δ −0.044） |
| PM110 RTM | **0.576** [0.43, 0.72]（低于 Vina） |
| PM110 GNINA | **0.522** [0.38, 0.66]（近随机） |

## 对投稿的含义

1. 正文并列报告 E=8；不得把 PM 优势归因于 exhaustiveness。
2. 单靶 AUROC≈0.60–0.63：对接有弱信号，但不是强 VS 引擎。
3. 扩面后 Vina 同向、|Δ|<0.05，CI 收窄；**ceiling 不变**（非通用决策臂）。
4. 正文以 **Vina 口袋匹配** 作 PM 主报告；RTM/GNINA 作重打分对照（未增强，甚至更弱）。

## 未做（显式）

- C 组（decoy 受体 / 换构象 / holdout）未开
- EGFR 扩面、LigPrep 主面板仍禁止
