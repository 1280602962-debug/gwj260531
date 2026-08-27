# Cluster uncertainty verdict (v1)

Compares ligand-level, document-connected, and Bemis-Murcko scaffold-cluster
bootstrap 95% CIs on frozen high-confidence Vina pocket-matched directional AUROCs.
B = 2000, seed 20260729. Scaffold grouping frozen before score inspection.

## EGFR/HER2 weak arm (D_vs_B, pocket A)

| estimator | AUROC | CI lo | CI hi | groups | note |
|-----------|------:|------:|------:|-------:|------|
| ligand bootstrap | 0.4297 | 0.2913 | 0.5815 | — | class-preserving resample |
| document cluster | 0.4297 | 0.3214 | 0.6171 | 23 docs | resamples document-connected ligand groups, not individual ligands |
| scaffold cluster | 0.4297 | 0.2777 | 0.595 | 51 scaffolds | resamples Bemis-Murcko scaffold groups, not individual ligands |

The weak EGFR/HER2 arm (D_vs_B = 0.4297) stays near chance under all three
resampling schemes; scaffold and document CIs are wider than ligand bootstrap
because correlated chemotypes/documents are kept together.

## PIK3CA/mTOR issues

| contrast | ligand CI | document CI | scaffold CI | document status | scaffold note |
|----------|-----------|-------------|-------------|-----------------|---------------|
| D_vs_A | [0.5078, 0.8929] | [0.4, 0.8864] | [0.5052, 0.8941] | ok | resamples Bemis-Murcko scaffold groups, not individual ligands |
| D_vs_B | [0.4977, 0.8657] | [0.0, 0.8179] | [0.495, 0.8714] | cannot_stably_estimate (1 valid doc-blocked fold) | resamples Bemis-Murcko scaffold groups, not individual ligands |
| summary_min | [0.4702, 0.8133] | — | [0.4635, 0.8111] | ok | resamples Bemis-Murcko scaffold groups, not individual ligands |

## Cross-estimator summary (all K=4 pairs)

| pair | contrast | ligand CI width | document CI width | scaffold CI width |
|------|----------|----------------:|------------------:|------------------:|
| EGFR/HER2 | D_vs_A | 0.272 | 0.331 | 0.296 |
| EGFR/HER2 | D_vs_B | 0.290 | 0.296 | 0.317 |
| AChE/BChE | D_vs_A | 0.310 | 0.312 | 0.314 |
| AChE/BChE | D_vs_B | 0.299 | 0.324 | 0.333 |
| PIK3CA/PIK3CB | D_vs_A | 0.312 | 0.408 | 0.337 |
| PIK3CA/PIK3CB | D_vs_B | 0.304 | 0.388 | 0.347 |
| PIK3CA/mTOR | D_vs_A | 0.385 | 0.486 | 0.389 |
| PIK3CA/mTOR | D_vs_B | 0.368 | 0.818 | 0.376 |

Scaffold-cluster bootstrap generally widens CIs relative to ligand bootstrap and
is comparable to or slightly narrower than document-cluster bootstrap depending on
whether document co-reporting or scaffold reuse dominates correlation structure.
PIK3CA/mTOR remains the most uncertainty-sensitive pair because of small n and
heavy document/scaffold concentration; document-blocked CV already flags D_vs_B as
not stably estimable under leave-document-out folds.
