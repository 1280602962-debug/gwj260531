# BindingDB-native external-slice freeze

Protocol frozen in `protocol/external_slice_contract.yaml` before docking.
This session did **not** dock and does **not** package external evaluation.
BindingDB archive `202608` md5-verified. ChEMBL document lookup: `chembl_document_api_partial`.

Primary-gate pairs (directional n≥20, ≥3 sources/class, top-document ≤50%): **0**.
Thin EGFR-style replications: **0**.
Packaged as external evaluation: **no** (no docking; gates require ≥2 primary pairs plus later pose deposit).

| pair | dual/A/B/neither after ECFP<0.70 | gate | sources dual/A/B |
|------|--------------------------------:|------|------------------|
| EGFR/HER2 | 180/10/20/6 | insufficient | 16/5/4 |
| AChE/BChE | 4/8/14/59 | insufficient | 2/6/3 |
| PIK3CA/PIK3CB | 9/0/3/100 | insufficient | 4/0/1 |
| PIK3CA/mTOR | 91/4/1/2 | insufficient | 9/2/1 |
| MCL1/Bcl-xL | 1/0/2/0 | insufficient | 1/0/1 |

Stop rule: fewer than two primary BindingDB-native pairs. Keep the manuscript as a
four-pair formulation audit. Do not call the remaining ligands a database-external set.

## MCL1/Bcl-xL

ChEMBL map at θ=6.0: dual/A/B/neither 82/77/24/122.
Frozen panel: 24/24/24/24.
B-only is exhaustive on the cached map. No same-library holdout. LC6 pose-gold gate was not run (no Vina).
Do not call this pair a disparate-fold pair. It is a PPI/BH3 groove domain shift.

## Primary receptors

- MCL1 3WIY (primary): chain `A` entity 1, 2.15 Å, mutations=0.
- BCL2L1 3WIZ (primary): chain `A` entity 1, 2.45 Å, mutations=0.
- MCL1 6UDV (alternate): chain `A` entity 1, 1.35 Å, mutations=0.
- BCL2L1 3SP7 (alternate): chain `A` entity 1, 1.4 Å, mutations=0.

No BindingDB AUROC was computed. Do not inspect these counts and then change θ, boxes, or receptors.
Archive sha256 values are in `tables/bindingdb_archive_lock_v1.csv` (4 files).
