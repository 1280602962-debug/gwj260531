# FINAL_MANUSCRIPT_NUMERIC_AUDIT

One-number → source check for DualFourClass JCIM manuscript package.
Manuscript: `docs/MANUSCRIPT_JCIM_EN.md`

## Key manuscript numbers

| claim | value | in_ms | source hits | checksum16 | status |
|---|---|---:|---|---|---|
| EGFR/HER2 Dual-vs-neither | 0.756 | 1 | data/jcim_novelty_v0/tables/aggregation_min_mean_geometric_harmonic_v1.csv; data | 9fa170b2cf653d9a | PASS |
| EGFR/HER2 summary_min | 0.430 | 1 |  |  | CHECK |
| EGFR/HER2 weak arm Dual-vs-B-only | 0.430 | 1 |  |  | CHECK |
| AChE/BChE summary_min | 0.606 | 1 |  |  | CHECK |
| PIK3CA/PIK3CB summary_min | 0.500 | 1 |  |  | CHECK |
| PIK3CA/mTOR summary_min | 0.692 | 1 |  |  | CHECK |
| ECFP4 docking increment max | 0.020 | 1 |  |  | CHECK |
| GNINA Dual-vs-neither EGFR | 0.783 | 1 |  |  | CHECK |
| GNINA summary_min EGFR | 0.220 | 1 |  |  | CHECK |
| MCL1 exploratory summary_min | 0.609 | 0 |  |  | CHECK |
| BindingDB eligible pairs | 0 | 1 | data/jcim_novelty_v0/tables/aggregation_min_mean_geometric_harmonic_v1.csv; data | 9fa170b2cf653d9a | PASS |
| K pairs | 4 | 1 | data/jcim_novelty_v0/tables/aggregation_min_mean_geometric_harmonic_v1.csv; data | 9fa170b2cf653d9a | PASS |

## Artifact presence / stale-claim checks

- PASS: CLAIM_CEILING MCL1 demotion language
- PASS: CLAIM_CEILING BindingDB zero pairs
- PASS: MCL1 formal demotion doc present
- PASS: detectable_effect_simulation_v1.csv
- PASS: scaffold_cluster_bootstrap_v1.csv
- PASS: bindingdb_external_feasibility_flow_v1.csv
- PASS: leave_cognate_out_v1.csv

## Uncertainty matrix (prespecified)

| source | EGFR/HER2 | AChE/BChE | PIK3CA/PIK3CB | PIK3CA/mTOR |
|---|---|---|---|---|
| ligand bootstrap | yes | yes | yes | yes |
| scaffold-cluster bootstrap | yes | yes | yes | yes |
| document-cluster bootstrap | yes | yes | yes / limited | neither unstable (n=4, 1 doc) |
| receptor realization | — | — | yes (alt crystals) | yes (alt crystals) |
| docking seed (five frozen seeds, v2 AUC(vina_mean)) | complete; summary_min median 0.3728 | complete; median 0.5988 | complete; median 0.4783 | complete; median 0.7037 |
| detectable-effect simulation | yes | yes | yes | yes |

## Rules

- Do not impute not-stably-estimable cells.
- Do not replace primary seed-20260727 Table 2 with multi-seed averages.
- BindingDB remains a supply-freeze negative result.
- MCL1/Bcl-xL remains exploratory stress-test.
- Leave-cognate-out removes one exact co-crystallized ligand only; it is not a train/test leakage or chemotype-removal test.

