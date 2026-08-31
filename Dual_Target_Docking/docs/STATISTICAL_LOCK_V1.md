# Statistical lock (submission)

Date: 2026-08-31  
Branch: `cursor/jcim-final-integration-0b1a`  
Parent freeze: `cursor/pik3ca-mtor-structure-freeze-0b1a@236fd60c`  
ChatGPT snapshot (not merged): `chatgpt/jcim-manuscript-review-20260828@d6cc6c12`

This page is the only manuscript-facing definition of primary estimands.  
`PRIMARY_METRIC_V2.md` is a dated 2026-07-29 snapshot and is not Table 2.

Any script that computes a different estimand must not be labeled primary, Table 2, or Table 3.

---

## Primary labels

θ = 6.0 four-state assignment (dual / A-only / B-only / neither), freeze date 2026-07-23.

## Primary scores

Pocket-matched directional AutoDock Vina mode-1 scores on the four frozen panels.  
\(S = -E_{\mathrm{Vina}}\) (higher is more favorable). Do not refresh the four ablation score tables.

## Directional arms (Table 2)

- Dual vs A-only: pocket B
- Dual vs B-only: pocket A
- Dual is always the positive class

## Primary descriptive summary

\(\mathrm{summary}_{\min} = \min(\mathrm{AUROC}_{D/A},\;\mathrm{AUROC}_{D/B})\)

Both arms are always reported. `summary_min` is not a scoring function.

## Table 2 uncertainty (canonical; do not swap)

Ligand-level **non-stratified** percentile bootstrap, \(B = 2000\), deterministic SHA-256-derived sub-seed.

Canonical source: `data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv`  
(`label_rule = theta_6.0`; script `build_t0_strengthen_v1.py`, offset key `stable_offset(pair, "theta_6.0")`).

Locked manuscript rounding:

| Pair | summary_min | 95% CI |
|---|---:|---|
| EGFR/HER2 | 0.430 | [0.282, 0.578] |
| AChE/BChE | 0.606 | [0.437, 0.730] |
| PIK3CA/PIK3CB | 0.500 | [0.350, 0.650] |
| PIK3CA/mTOR | 0.692 | [0.470, 0.813] |

Coeval, **not** Table 2: `pocket_matched_directional_v1.csv` (different hash-offset key).  
Deprecated: `PRIMARY_METRIC_V2.md` (seed=20260729, no offset).  
Legacy pooled control: `bootstrap_directional_ci_v1.csv`.

Future Table 2 / Figure 3 / Table S4 / forest Vina CIs must **read** this CSV. They must not each re-bootstrap.

## Dual vs neither (Table 3)

Per-ligand \(\mathrm{vina\_mean} = (S_A+S_B)/2\), then **one** AUROC.

Canonical source: `data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv`

Locked points: 0.7560 / 0.6494 / 0.5592 / 0.5139.

`mean(\mathrm{AUC}_A, \mathrm{AUC}_B)` is **not** Table 3. It may be stored only as `mean_marginal_pocket_auroc_D_vs_neither`.

## Multi-seed (Table S54)

Same estimands as Table 2 and Table 3. Report median / IQR / range across five frozen Vina seeds  
(20260727 primary, reused, plus 20260811–20260814).

Canonical sources:

- scores: `data/jcim_multiseed_v0/tables/multiseed_scores_long_v1.csv` (do not re-dock)
- analyzer: `data/jcim_multiseed_v0/scripts/analyze_multiseed_vina_v2.py`
- tables: `multiseed_auroc_by_seed_v2.csv`, `multiseed_auroc_aggregate_v2.csv`, `multiseed_consistency_v2.csv`

Primary seed Dual-versus-neither **must** recover Table 3. The v2 analyzer exits non-zero otherwise.

v1 Dual-versus-neither columns are a dated wrong-estimand snapshot. Do not copy them into the article.

## Sensitivity only

- Cluster bootstrap (document / scaffold)
- Receptor realization
- Independent GNINA pose generation
- θ grid, holdout, time split, BindingDB-native gate

## Forbidden

- Replacing Table 2 with a multi-seed mean or a later API refetch
- Packaging BindingDB as external validation
- Adding MCL1/Bcl-xL to Table 2
- Switching Table 2 CIs to pocket_matched or PRIMARY_METRIC without a new statistical freeze
