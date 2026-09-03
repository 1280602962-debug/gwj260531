# Four-pair multi-seed Vina sensitivity (v2; Table 3 estimand)

Source scores: `data/jcim_multiseed_v0/tables/multiseed_scores_long_v1.csv`
Frozen seeds: 20260727 (primary, reused) + 20260811–20260814.
Protocol otherwise identical to production (receptors, boxes, exhaustiveness, modes, energy_range).

## Estimands

- Directional Dual vs A-only / Dual vs B-only / `summary_min`: pocket-matched; unchanged from v1.
- **Primary Dual-versus-neither:** per-ligand `vina_mean = (S_A+S_B)/2`, then one AUROC.
  This is the Table 3 estimand. Primary seed 20260727 recovered
  0.7560 / 0.6494 / 0.5592 / 0.5139.
- **Sensitivity only:** `mean_marginal_pocket_auroc_D_vs_neither = mean(AUC_A, AUC_B)`
  (v1 column; do not cite as Table 3 or as the Table S54 Dual-versus-neither value).
- Formulation gap = Dual-versus-neither (`vina_mean`) − `summary_min`.

## Per-seed metrics (primary Dual-versus-neither = vina_mean)

| pair | seed | dual_vs_A | dual_vs_B | summary_min | dual_vs_neither_vina_mean | gap | mean_marginal (legacy) |
|---|---:|---:|---:|---:|---:|---:|---:|
| AChE/BChE | 20260727 | 0.6504 | 0.6058 | 0.6058 | 0.6494 | 0.0436 | 0.6494 |
| AChE/BChE | 20260811 | 0.5988 | 0.6396 | 0.5988 | 0.6564 | 0.0576 | 0.659 |
| AChE/BChE | 20260812 | 0.5645 | 0.6104 | 0.5645 | 0.6427 | 0.0781 | 0.6373 |
| AChE/BChE | 20260813 | 0.6017 | 0.6104 | 0.6017 | 0.656 | 0.0543 | 0.6387 |
| AChE/BChE | 20260814 | 0.5527 | 0.6119 | 0.5527 | 0.656 | 0.1033 | 0.6507 |
| EGFR/HER2 | 20260727 | 0.6664 | 0.4297 | 0.4297 | 0.756 | 0.3263 | 0.7641 |
| EGFR/HER2 | 20260811 | 0.6598 | 0.3728 | 0.3728 | 0.7768 | 0.404 | 0.7798 |
| EGFR/HER2 | 20260812 | 0.6391 | 0.3214 | 0.3214 | 0.7619 | 0.4405 | 0.7634 |
| EGFR/HER2 | 20260813 | 0.6532 | 0.3672 | 0.3672 | 0.7589 | 0.3917 | 0.7693 |
| EGFR/HER2 | 20260814 | 0.6692 | 0.3956 | 0.3956 | 0.7679 | 0.3722 | 0.7768 |
| PIK3CA/PIK3CB | 20260727 | 0.6905 | 0.5 | 0.5 | 0.5592 | 0.0592 | 0.5798 |
| PIK3CA/PIK3CB | 20260811 | 0.6587 | 0.5019 | 0.5019 | 0.5737 | 0.0717 | 0.5926 |
| PIK3CA/PIK3CB | 20260812 | 0.6634 | 0.4783 | 0.4783 | 0.5893 | 0.111 | 0.5977 |
| PIK3CA/PIK3CB | 20260813 | 0.6541 | 0.4681 | 0.4681 | 0.5759 | 0.1078 | 0.5915 |
| PIK3CA/PIK3CB | 20260814 | 0.6772 | 0.4707 | 0.4707 | 0.5882 | 0.1175 | 0.5993 |
| PIK3CA/mTOR | 20260727 | 0.7143 | 0.6921 | 0.6921 | 0.5139 | -0.1782 | 0.5278 |
| PIK3CA/mTOR | 20260811 | 0.746 | 0.713 | 0.713 | 0.5278 | -0.1852 | 0.5486 |
| PIK3CA/mTOR | 20260812 | 0.7262 | 0.7778 | 0.7262 | 0.5556 | -0.1706 | 0.5486 |
| PIK3CA/mTOR | 20260813 | 0.7381 | 0.6759 | 0.6759 | 0.4861 | -0.1898 | 0.5347 |
| PIK3CA/mTOR | 20260814 | 0.744 | 0.7037 | 0.7037 | 0.5278 | -0.1759 | 0.5347 |

## Consistency vs primary seed (vina_mean gap)

- **EGFR/HER2**: summary_min median 0.3728 (range 0.3214–0.4297); gap-sign match 5/5; positive gap 5/5; neither>summary_min order match 5/5.
- **AChE/BChE**: summary_min median 0.5988 (range 0.5527–0.6058); gap-sign match 5/5; positive gap 5/5; neither>summary_min order match 5/5.
- **PIK3CA/PIK3CB**: summary_min median 0.4783 (range 0.4681–0.5019); gap-sign match 5/5; positive gap 5/5; neither>summary_min order match 5/5.
- **PIK3CA/mTOR**: summary_min median 0.7037 (range 0.6759–0.7262); gap-sign match 5/5; positive gap 0/5; neither>summary_min order match 5/5.

## Claim ceiling

- Allowed: report median/IQR/range across frozen seeds; state whether the primary qualitative pattern held.
- Forbidden: picking a favorable seed; replacing primary Table 2 with a multi-seed mean; claiming seed robustness beyond these four pairs; citing v1 Dual-versus-neither as Table 3.

v1 tables remain as a dated wrong-estimand snapshot and must not be copied into the article.

