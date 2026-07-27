# Decision ablation v0 — frozen thresholds
## Frozen a priori (before case inspection)
- `shortfall_lambda` = **0.5**
- `consensus_top_frac` = **0.25** → top_n = **12** / 48
- Flags/chemotype warnings **do not** enter gated scores (`flags_enter_score=0`).
- No clash retune; no new docking.

## Arms
| arm | definition |
|-----|------------|
| A vina_mean | (−aff_A − aff_B)/2 |
| B rtm_min_z | min of per-target RTM z |
| C rtm_shortfall | rtm_min_z − λ\|zA−zB\| |
| D1 consensus_rank_mean | −mean(rank_vina, rank_rtm) |
| D2 consensus_and_top25 | pass if both ranks ≤ top_n; else deep penalty |

## Metrics (Dual vs rest)

| arm | auroc_dual_vs_rest | top10_dual | top10_A_only | top10_B_only | top10_neither | top10_hardneg | top1 | frozen_shortfall_lambda | frozen_consensus_top_frac | frozen_top_n |
|---|---|---|---|---|---|---|---|---|---|---|
| vina_mean | 0.6333333333333333 | 5 | 1 | 3 | 1 | 5 | PM48_10 | 0.5 | 0.25 | 12 |
| rtm_min_z | 0.6851851851851852 | 6 | 4 | 0 | 0 | 4 | PM48_26 | 0.5 | 0.25 | 12 |
| rtm_shortfall | 0.687037037037037 | 6 | 4 | 0 | 0 | 4 | PM48_26 | 0.5 | 0.25 | 12 |
| consensus_rank_mean | 0.6675925925925926 | 4 | 4 | 1 | 1 | 6 | PM48_01 | 0.5 | 0.25 | 12 |
| consensus_and_top25 | 0.6962962962962963 | 6 | 4 | 0 | 0 | 4 | PM48_01 | 0.5 | 0.25 | 12 |

## Key ligand ranks (1=best)

| ligand | class | pref_name | vina_mean | rtm_min_z | rtm_shortfall | consensus_rank_mean | consensus_and_top25 |
|---|---|---|---|---|---|---|---|
| PM48_01 | dual | PI-103 | 9 | 4 | 2 | 1 | 1 |
| PM48_02 | dual | OMIPALISIB | 3 | 30 | 40 | 19 | 30 |
| PM48_10 | dual | TORIN1 | 1 | 31 | 33 | 16 | 31 |
| PM48_20 | A_only | nan | 14 | 2 | 3 | 2 | 4 |
| PM48_21 | A_only | nan | 18 | 5 | 5 | 8 | 6 |
| PM48_26 | A_only | nan | 20 | 1 | 1 | 4 | 3 |
| PM48_34 | B_only | WYE-132 | 10 | 40 | 44 | 25 | 40 |

## Verdict vs success criteria
Success = hardneg Top10 ↓ vs B **and** Torin1/Omipalisib (PM48_10/02) not clearly worse than under B.

- **rtm_shortfall**: hardneg Top10 4 vs baseline 4; PM48_10 rank 33 (B=31); PM48_02 rank 40 (B=30) → FAIL
- **consensus_rank_mean**: hardneg Top10 6 vs baseline 4; PM48_10 rank 16 (B=31); PM48_02 rank 19 (B=30) → FAIL
- **consensus_and_top25**: hardneg Top10 4 vs baseline 4; PM48_10 rank 31 (B=31); PM48_02 rank 30 (B=30) → FAIL

### Bottom line

**无法同时满足**：在冻结阈值下，shortfall / consensus 未能在降低硬负 Top10 的同时保护 Torin1/Omipalisib；或硬负 Top10 未下降。主文应并列报告 `vina_mean` 与 `rtm_min_z`，并用化学型警告层标注 T2，而不是宣称决策规则已闭环。

Do **not** claim C4 extrapolation success. Do **not** retune clash to drop PM48_26.
