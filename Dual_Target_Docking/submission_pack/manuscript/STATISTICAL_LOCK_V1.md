# Statistical lock (submission)

Date: 2026-09-09  
Branch: `cursor/jcim-results-revision-0b1a`  
This page replaces the earlier K=4 lock. PIK3CA/PIK3CB is withdrawn and is not a Table 2 row.

This page is the only manuscript-facing definition of primary estimands.  
`PRIMARY_METRIC_V2.md` is a dated 2026-07-29 snapshot and is not Table 2.

Any script that computes a different estimand must not be labeled primary, Table 2, or Table 3.

---

## Primary labels

θ = 6.0 four-state assignment (dual / A-only / B-only / neither), freeze date 2026-07-23.

Candidate-pool rules, class quotas, and scaffold caps differ by pair (manuscript Table 1). They do not replace these analysis labels.

## Primary scores

Pocket-matched directional AutoDock Vina mode-1 scores on the eight primary pairs.  
\(S = -E_{\mathrm{Vina}}\) (higher is more favorable).

## Directional arms (Table 2)

- Dual vs A-only: pocket B
- Dual vs B-only: pocket A
- Dual is always the positive class

## Primary descriptive summary

\(\mathrm{summary}_{\min} = \min(\mathrm{AUROC}_{D/A},\;\mathrm{AUROC}_{D/B})\)

Both arms are always reported. `summary_min` is not a scoring function. Eight pairs are not averaged.

## Table 2 uncertainty (canonical; do not swap)

Ligand-level **non-stratified** percentile bootstrap, \(B = 2000\), deterministic SHA-256-derived sub-seed.

Canonical sources (do not mix):

| Pairs | CSV |
|---|---|
| EGFR/HER2, AChE/BChE, PIK3CA/mTOR | `data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv` (`label_rule = theta_6.0`) |
| F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, PPARA/PPARD | `data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv` |

Read-only loader: `scripts/primary/bootstrap_primary.py`.

Locked manuscript rounding:

| Pair | summary_min | 95% CI |
|---|---:|---|
| EGFR/HER2 | 0.430 | [0.282, 0.578] |
| AChE/BChE | 0.606 | [0.437, 0.730] |
| PIK3CA/mTOR | 0.692 | [0.470, 0.813] |
| F2/F10 | 0.345 | [0.211, 0.477] |
| JAK1/TYK2 | 0.365 | [0.231, 0.503] |
| JAK1/JAK2 | 0.588 | [0.444, 0.725] |
| PPARG/PPARA | 0.649 | [0.504, 0.751] |
| PPARA/PPARD | 0.446 | [0.296, 0.584] |

Coeval, **not** Table 2: `pocket_matched_directional_v1.csv` (different hash-offset key).  
Deprecated: `PRIMARY_METRIC_V2.md` (seed=20260729, no offset).  
Legacy pooled control: `bootstrap_directional_ci_v1.csv`.

Future Table 2 / Figure 2 Vina CIs must **read** these CSVs. They must not each re-bootstrap.

## Dual vs neither (Table 3)

Per-ligand \(\mathrm{vina\_mean} = (S_A+S_B)/2\), then **one** AUROC.

Canonical sources:

| Pairs | CSV |
|---|---|
| EGFR/HER2, AChE/BChE, PIK3CA/mTOR | `data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv` |
| F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, PPARA/PPARD | same five-pair `table2_comparable_theta6_v1.csv` (`D_vs_neither_vina_mean`) |

Locked points (three original pairs): 0.7560 / 0.6494 / 0.5139.  
Five-pair Dual-versus-neither points: 0.5188 / 0.7696 / 0.7299 / 0.6853 / 0.5647.

`mean(\mathrm{AUC}_A, \mathrm{AUC}_B)` is **not** Table 3.

## BindingDB remainder (Table S11b)

Canonical eight-pair source: `data/jcim_novelty_v0/tables/external_slice_summary_v1.csv`.  
No pair was packaged as an external evaluation set.

Historical three-pair remainder snapshot (do not use as Table S11):  
`data/jcim_novelty_v0/tables/external_slice_summary_202608_contract_v1.csv`.  
Historical REST counts: `bindingdb_independence_summary_v1.csv`.

## Multi-seed

Same estimands as Table 2 and Table 3. Report median / IQR / range across five frozen Vina seeds  
(20260727 primary, reused, plus 20260811–20260814). Typeset SI cites this as Table S9; the former working number was Table S54.

Canonical sources:

- scores: `data/jcim_multiseed_v0/tables/multiseed_scores_long_v1.csv` (do not re-dock)
- analyzer: `data/jcim_multiseed_v0/scripts/analyze_multiseed_vina_v2.py`
- tables: `multiseed_auroc_by_seed_v2.csv`, `multiseed_auroc_aggregate_v2.csv`, `multiseed_consistency_v2.csv`

Primary seed Dual-versus-neither **must** recover Table 3 on the pairs present in the multi-seed tables.

## Leave-exact-cognate-out

Remove the exact co-crystallized ligand only when it is a member of the main panel, then recompute the Table 2 directional AUROCs, descriptive `summary_min`, and Table 3 `AUC(vina_mean)` on seed 20260727.

Canonical source: `data/jcim_novelty_v0/tables/leave_cognate_out_v1.csv`

## Sensitivity only

- Cluster bootstrap (document / scaffold), including the two flagship fixed-channel Δ values
- Receptor realization (PIK3CA/mTOR)
- Independent GNINA pose generation (EGFR/HER2, PIK3CA/mTOR, JAK1/TYK2)
- θ grid, holdout, time split, BindingDB-native gate

## Forbidden

- Replacing Table 2 with a multi-seed mean or a later API refetch
- Packaging BindingDB as external validation
- Adding MCL1/Bcl-xL or PIK3CA/PIK3CB to Table 2
- Switching Table 2 CIs to pocket_matched or PRIMARY_METRIC without a new statistical freeze
- Treating the eight pairs as eight independent replicates
