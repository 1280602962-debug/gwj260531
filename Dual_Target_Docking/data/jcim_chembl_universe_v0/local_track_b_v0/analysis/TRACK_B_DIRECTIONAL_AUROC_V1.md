# Track B directional AUROC (local Vina pack)

Engine: AutoDock Vina 1.2.7 mode-1; `S = −E`. Seed 20260727. Exhaustiveness 8.
Estimand: pocket-matched directional AUROC; `summary_min = min(D/A, D/B)`.
Bootstrap: B=2000, seed 20260729, class-preserving. **Does not replace Table 2.**
Count **three systems** (coagulation, JAK, PPAR), not five pairs.

## Primary (panel strict 6.5/5.5 class)

| pair | system | n_scored | D/A (pocket B) | D/B (pocket A) | summary_min [95% CI] |
|------|--------|---------:|---------------:|---------------:|---------------------:|
| F2/F10 | coagulation | 107 | 0.4133 | 0.3448 | 0.3448 [0.2077, 0.4476] |
| JAK1/TYK2 | JAK | 109 | 0.5751 | 0.3649 | 0.3649 [0.2339, 0.5071] |
| JAK1/JAK2 | JAK | 110 | 0.5884 | 0.7275 | 0.5884 [0.4404, 0.7139] |
| PPARG/PPARA | PPAR | 109 | 0.6492 | 0.7061 | 0.6492 [0.5101, 0.75] |
| PPARA/PPARD | PPAR | 110 | 0.6465 | 0.4463 | 0.4463 [0.293, 0.5801] |

## Unified θ = 6.0 labels (same scores)

| pair | system | n_scored | D/A (pocket B) | D/B (pocket A) | summary_min [95% CI] |
|------|--------|---------:|---------------:|---------------:|---------------------:|
| F2/F10 | coagulation | 107 | 0.4133 | 0.3448 | 0.3448 [0.2077, 0.4466] |
| JAK1/TYK2 | JAK | 109 | 0.5751 | 0.3649 | 0.3649 [0.2339, 0.501] |
| JAK1/JAK2 | JAK | 110 | 0.5884 | 0.7275 | 0.5884 [0.4458, 0.7168] |
| PPARG/PPARA | PPAR | 109 | 0.6492 | 0.7061 | 0.6492 [0.5039, 0.7393] |
| PPARA/PPARD | PPAR | 110 | 0.6465 | 0.4463 | 0.4463 [0.2988, 0.5879] |

Companion four-descriptor reference (same scored ligands): `analysis/TRACK_B_DESCRIPTOR_REFERENCE_V1.md`. Vina does not beat the best single descriptor by more than +0.022 on any pair.

Artifacts: `tables/track_b_summary_min_v1.csv`, `tables/track_b_directional_auroc_v1.csv`, `tables/scores_vina_mode1_v1.csv`, `tables/job_status.csv`, `tables/layer3_cognate_rmsd_v1.csv`.
