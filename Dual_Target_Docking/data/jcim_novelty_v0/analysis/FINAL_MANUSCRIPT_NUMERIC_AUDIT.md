# FINAL_MANUSCRIPT_NUMERIC_AUDIT

One-number → source check for DualFourClass JCIM manuscript package.
Manuscript: `docs\MANUSCRIPT_JCIM_EN.md`

## Key manuscript numbers

| claim | value | in_ms | source hits | checksum16 | status |
|---|---|---:|---|---|---|
| EGFR/HER2 Dual-vs-neither | 0.756 | 1 | data\jcim_novelty_v0\tables\aggregation_min_mean_geometric_harmonic_v1.csv; data | 9fa170b2cf653d9a | PASS |
| EGFR/HER2 summary_min | 0.430 | 1 | data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv | 19316c05ddb5e205 | PASS |
| EGFR/HER2 weak arm Dual-vs-B-only | 0.430 | 1 | data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv | 19316c05ddb5e205 | PASS |
| AChE/BChE summary_min | 0.606 | 1 | data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv | 19316c05ddb5e205 | PASS |
| PIK3CA/PIK3CB summary_min | 0.500 | 1 | data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv | 19316c05ddb5e205 | PASS |
| PIK3CA/mTOR summary_min | 0.692 | 1 | data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv | 19316c05ddb5e205 | PASS |
| ECFP4 docking increment max | 0.020 | 1 | data/jcim_novelty_v0/tables/incremental_information_v1.csv | 97389149d3951686 | PASS |
| GNINA Dual-vs-neither EGFR | 0.783 | 1 | data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv | 37d0e1db02121c33 | PASS |
| GNINA summary_min EGFR | 0.220 | 1 | data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv | 37d0e1db02121c33 | PASS |
| MCL1 exploratory summary_min | 0.609 | 0 | data/mcl1_bclxl_panel_v0/tables/formulation_auroc_MBX_v1.csv | a6387cf3fc9fb6a1 | NOT_REPORTED_EXPECTED |
| BindingDB eligible pairs | 0 | 1 | data/jcim_novelty_v0/tables/external_slice_summary_v1.csv | fe1e1b61eeea1102 | PASS |
| K pairs | 4 | 1 | data\jcim_novelty_v0\tables\aggregation_min_mean_geometric_harmonic_v1.csv; data | 9fa170b2cf653d9a | PASS |

## Artifact presence / stale-claim checks

- PASS: CLAIM_CEILING MCL1 demotion language
- PASS: CLAIM_CEILING BindingDB zero pairs
- PASS: MCL1 formal demotion doc present
- PASS: detectable_effect_simulation_v1.csv
- PASS: scaffold_cluster_bootstrap_v1.csv
- PASS: bindingdb_external_feasibility_flow_v1.csv

## Uncertainty matrix (prespecified)

| source | EGFR/HER2 | AChE/BChE | PIK3CA/PIK3CB | PIK3CA/mTOR |
|---|---|---|---|---|
| ligand bootstrap | yes | yes | yes | yes |
| scaffold-cluster bootstrap | yes | yes | yes | yes |
| document-cluster bootstrap | yes | yes | yes / limited | neither unstable (n=4, 1 doc) |
| receptor realization | — | — | yes (alt crystals) | yes (alt crystals) |
| docking seed (five frozen seeds) | complete; median 0.3728 | complete; median 0.5988 | complete; median 0.4783 | complete; median 0.7037 |
| detectable-effect simulation | yes | yes | yes | yes |

## Rules

- Do not impute not-stably-estimable cells.
- Do not replace primary seed-20260727 Table 2 with multi-seed averages.
- BindingDB remains a supply-freeze negative result.
- MCL1/Bcl-xL remains exploratory stress-test.

