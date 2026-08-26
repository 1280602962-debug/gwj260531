# Time-split protocol freeze

This protocol was frozen before AUROC on the dated splits was computed. Do not change the cutoff to recover a class.

## Rule

- Ligand year = minimum ChEMBL `document.year` among retained high-confidence records for that ligand.
- Train/development: first_year < cutoff.
- Time-split test: first_year ≥ cutoff.
- Late ligands are not used to choose thresholds, receptors, metrics, or which pairs to emphasize.
- Primary cutoff: **2018**.
- Sensitivity cutoffs, also frozen: **2015** and **2020**.
- AUROC is reported only if dual, A-only, and B-only each have n ≥ 10 in the test split (≥15 preferred).
- n < 10: counts only. The cutoff is not moved.
- Package as external validation only if at least two pairs pass the n ≥ 10 gate at the primary cutoff.

## Stop

If the primary cutoff fails that pair gate, keep the paper as an internal formulation audit and consider BindingDB next (`BINDINGDB_EXTERNAL_SOP.md`). Do not shop among 2015/2018/2020 for a favorable AUROC.
