# Time-split protocol freeze and result

Frozen before seeing AUROC:

- Primary cutoff year: **2018** (train first_year < 2018; test first_year ≥ 2018).
- Sensitivity cutoffs: 2015 and 2020. These were not chosen after looking at class counts that favor AUROC.
- Ligand year = minimum `document.year` among retained high-confidence records.
- Late ligands are not used to choose thresholds, receptors, or endpoints.
- Minimum AUROC gate: dual, A-only, and B-only each ≥10 in the test split; ≥15 preferred.
- Below 10: descriptive counts only. Cutoff is not moved to recover a class.
- External-validation package requires ≥2 pairs passing the AUROC gate at the primary cutoff.

Primary cutoff evaluable/underpowered pairs: **0**.
Packaged as external validation: **no**. Keep the internal formulation-audit claim.

| cutoff | pair | test dual/A/B/neither | gate | D_vs_A | D_vs_B |
|-------:|------|----------------------:|------|--------|--------|
| 2015 | EGFR/HER2 | 9/8/15/6 | descriptive_only |  |  |
| 2015 | AChE/BChE | 11/11/24/10 | underpowered_report | 0.6033 | 0.6364 |
| 2015 | PIK3CA/PIK3CB | 17/20/4/8 | descriptive_only |  |  |
| 2015 | PIK3CA/mTOR | 6/2/1/0 | descriptive_only |  |  |
| 2018 | EGFR/HER2 | 6/3/14/2 | descriptive_only |  |  |
| 2018 | AChE/BChE | 8/5/15/6 | descriptive_only |  |  |
| 2018 | PIK3CA/PIK3CB | 12/11/0/3 | unevaluable |  |  |
| 2018 | PIK3CA/mTOR | 2/0/1/0 | unevaluable |  |  |
| 2020 | EGFR/HER2 | 6/0/14/2 | unevaluable |  |  |
| 2020 | AChE/BChE | 5/3/14/6 | descriptive_only |  |  |
| 2020 | PIK3CA/PIK3CB | 6/4/0/1 | unevaluable |  |  |
| 2020 | PIK3CA/mTOR | 1/0/0/0 | unevaluable |  |  |

If this split is insufficient, BindingDB is the next option (`docs/BINDINGDB_EXTERNAL_SOP.md`).
Do not dock new ligands until the independent set and evaluation contract are frozen.
