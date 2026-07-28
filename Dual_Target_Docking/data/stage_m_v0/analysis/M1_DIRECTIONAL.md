# M1 — Directional metric (primary readout)

## Definition (frozen)

- **Primary:** AUROC(dual vs A_only) and AUROC(dual vs B_only), reported separately.
- **Summary (secondary):** `min(D/A, D/B)` and `mean(D/A, D/B)` for arm ranking only.
- **Deprecated as sole headline:** pooled Dual vs A∪B AUROC (cancels opposing directions; appendix only).
- Top10 hardneg counts reported **split by A_only / B_only**.

## EGFR/HER2 (N=110, all; prep mixed — interpret cautiously)

| arm | D/A | D/B | pooled (appendix) | Top10 A/B/dual |
|-----|-----|-----|-------------------|----------------|
| `vina_mean` | 0.689 | 0.311 | 0.516 | 5/4/1 |
| `rtm_min_z` | 0.587 | 0.317 | 0.464 | 5/2/3 |
| `heavy_atoms` | 0.700 | 0.369 | 0.549 | 2/5/3 |
| `MW` | 0.648 | 0.416 | 0.542 | 4/2/4 |
| `TPSA` | 0.703 | 0.427 | 0.577 | 3/5/2 |
| `morgan_dual_medsim` | 0.555 | 0.585 | 0.569 | 6/1/3 |

### Subsets (prep confound flagged)

- `old40`: LigPrep as-run; `new70`: RDKit+meeko. Do **not** treat old/new RTM split as a method conclusion until M4 unified prep.

- **old40** vina_mean D/A=0.628 D/B=0.453; rtm_min_z D/A=0.800 D/B=0.607
- **new70** vina_mean D/A=0.725 D/B=0.238; rtm_min_z D/A=0.444 D/B=0.101

## PIK3CA/mTOR (N=48)

| arm | D/A | D/B | pooled (appendix) | Top10 A/B/dual |
|-----|-----|-----|-------------------|----------------|
| `vina_mean` | 0.698 | 0.597 | 0.652 | 1/3/6 |
| `rtm_min_z` | 0.611 | 0.792 | 0.694 | 4/0/6 |
| `heavy_atoms` | 0.464 | 0.463 | 0.464 | 4/2/4 |
| `MW` | 0.448 | 0.472 | 0.459 | 5/2/3 |
| `TPSA` | 0.260 | 0.451 | 0.348 | 7/2/1 |
| `morgan_dual_medsim` | 0.639 | 0.694 | 0.664 | 3/1/6 |

## Gate

**M1 = Go** (definitional). Numbers rechecked vs `plan_v2_redteam_v0` (±0.005). Full table: `tables/m1_directional_auroc.csv`.
