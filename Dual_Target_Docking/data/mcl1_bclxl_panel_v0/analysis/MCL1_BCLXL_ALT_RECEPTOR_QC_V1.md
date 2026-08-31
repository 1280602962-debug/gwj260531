# MCL1/Bcl-xL alternate-receptor cognate QC + failure sensitivity

Role: applicability stress-test support. Not a fifth primary pair.

## Alternate cognate QC (element-Hungarian diagnostic)

| target | PDB | cognate | top-1 | best-top3 | best9 | gate_top3<2Å |
|---|---|---|---:|---:|---:|---:|
| MCL1 | 6UDV | Q51 | nan | nan | nan | 0 |
| BCL2L1 | 3SP7 | 03B | nan | nan | nan | 0 |

## Interpretation

- 6UDV gate_top3: 0
- 3SP7 gate_top3: 0
- Matching is legacy element-Hungarian (same class as primary LC6 diagnostic); not topology-aware.
- If an alternate fails the <2 Å diagnostic, do **not** retune the box to rescue it.
- Full-panel alternate-receptor AUROC is run only if the corresponding cognate diagnostic passes; otherwise stop at QC.

## Failure sensitivity (primary 3WIY/3WIZ panel)

See `tables/mcl1_bclxl_failure_properties_v1.csv` and `tables/mcl1_bclxl_failure_class_counts_v1.csv`.
Failures were ligand-prep/embed failures, not Vina timeouts.

