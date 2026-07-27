# Cross-pair bootstrap CI v0 — decision arms

Resampling: ligand bootstrap with replacement, **N=2000**, seed=20260727.
Frozen (unchanged): `shortfall_lambda=0.5`, `consensus_top_frac=0.25`.
Scores/flags: arms precomputed on full panel; flags **do not** enter scores; no re-docking / no threshold retune.
Significance for improvement vs `vina_mean`: paired-bootstrap 95% CI of ΔAUROC (arm − vina_mean) **excludes 0**.

## EGFR/HER2 panel40

| arm | AUROC [95% CI] | Top10 hardneg [95% CI] | Top10 dual [95% CI] |
|-----|----------------|------------------------|---------------------|
| vina_mean | 0.552 [0.369, 0.737] | 6 [3, 9] | 4 [1, 7] |
| rtm_min_z | 0.723 [0.536, 0.882] | 3 [1, 7] | 7 [3, 9] |
| rtm_shortfall | 0.712 [0.542, 0.876] | 4 [1, 7] | 6 [3, 9] |
| consensus_rank_mean | 0.655 [0.470, 0.829] | 5 [2, 8] | 5 [2, 8] |
| consensus_and_top25 | 0.707 [0.530, 0.875] | 3 [1, 7] | 7 [3, 9] |

### Δ vs vina_mean (paired bootstrap)

| compare | ΔAUROC [95% CI] | sig? | ΔTop10 hardneg [95% CI] | sig? |
|---------|----------------|------|-------------------------|------|
| rtm_min_z_minus_vina_mean | +0.171 [-0.015, +0.356] | no | -3 [-6, +0] | no |
| rtm_shortfall_minus_vina_mean | +0.160 [-0.035, +0.348] | no | -2 [-6, +1] | no |
| consensus_rank_mean_minus_vina_mean | +0.103 [+0.002, +0.205] | YES | -1 [-4, +0] | no |
| consensus_and_top25_minus_vina_mean | +0.155 [-0.021, +0.333] | no | -3 [-6, +0] | no |

## PIK3CA/mTOR panel48

| arm | AUROC [95% CI] | Top10 hardneg [95% CI] | Top10 dual [95% CI] |
|-----|----------------|------------------------|---------------------|
| vina_mean | 0.633 [0.455, 0.793] | 5 [1, 8] | 5 [2, 9] |
| rtm_min_z | 0.685 [0.519, 0.840] | 4 [1, 7] | 6 [3, 9] |
| rtm_shortfall | 0.687 [0.514, 0.850] | 4 [1, 7] | 6 [3, 9] |
| consensus_rank_mean | 0.668 [0.496, 0.816] | 6 [2, 8] | 4 [2, 8] |
| consensus_and_top25 | 0.696 [0.534, 0.850] | 4 [1, 7] | 6 [3, 9] |

### Δ vs vina_mean (paired bootstrap)

| compare | ΔAUROC [95% CI] | sig? | ΔTop10 hardneg [95% CI] | sig? |
|---------|----------------|------|-------------------------|------|
| rtm_min_z_minus_vina_mean | +0.052 [-0.144, +0.251] | no | -1 [-5, +4] | no |
| rtm_shortfall_minus_vina_mean | +0.054 [-0.154, +0.261] | no | -1 [-5, +4] | no |
| consensus_rank_mean_minus_vina_mean | +0.034 [-0.080, +0.139] | no | +1 [-3, +4] | no |
| consensus_and_top25_minus_vina_mean | +0.063 [-0.126, +0.255] | no | -1 [-5, +4] | no |

## EGFR panel40 — five-arm point estimates (new)

Same frozen decision arms as PIK3CA/mTOR; full table in `data/egfr_her2_panel40_v0/analysis/decision_ablation_v0/`.

| arm | AUROC | Top10 dual | Top10 hardneg |
|-----|-------|------------|---------------|
| vina_mean | 0.552 | 4 | 6 |
| rtm_min_z | 0.723 | 7 | 3 |
| rtm_shortfall | 0.712 | 6 | 4 |
| consensus_rank_mean | 0.655 | 5 | 5 |
| consensus_and_top25 | 0.707 | 7 | 3 |

## Updated conclusion (CI-aware)

1. **EGFR/HER2 — `rtm_min_z`:** point ΔAUROC = +0.171 [-0.015, +0.356] → **not significant at 95%** (CI includes 0; lower bound near 0 / borderline) vs `vina_mean`.
2. **EGFR/HER2 — other arms:** shortfall sig=False; **consensus_rank_mean sig=True** (Δ=+0.103 [+0.002, +0.205]); AND-top25 sig=False. Rank-mean consensus is the only arm with ΔAUROC CI excluding 0 on this pair; hardneg Top10 Δ still not significant.
3. **PIK3CA/mTOR — all five arms vs `vina_mean`:** `rtm_min_z` Δ=+0.052 [-0.144, +0.251] → **not significant**; shortfall=False; rank-mean=False; AND-top25=False. Point lifts exist but are unstable under ligand bootstrap.
4. **Top10 hardneg:** no arm on either pair has a paired-bootstrap Δ hardneg CI that excludes 0; apparent hardneg drops (e.g. EGFR 6→3) are **not** CI-significant.
5. **How to write it:** keep dual readout (`vina_mean` + `rtm_min_z`); do **not** claim a universally significant RTM upgrade across pairs; EGFR shows a larger point lift that is borderline for `rtm_min_z` and significant only for rank-mean consensus; PIK3CA/mTOR lifts are non-significant. Still no C4-closed / clash-retune / flags-in-score claim.
