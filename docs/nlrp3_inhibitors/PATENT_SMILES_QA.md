# Patent SMILES / Activity Fix & Sample QA

Date: 2026-07-13

## What was fixed

### SMILES auto-fixes
| Fix tag | Meaning |
|---------|---------|
| `n_oxide` | `N(O)=` / `N(O)` → `[N+]([O-])=` (WO2023147468 aromatic N-oxide) |
| `n_methyl_indazole` / `manual_n_methyl_indazole` | OCR `N=N(C)` → aromatic `nn(C)` (WO2022204227 #100) |

- SMILES auto-fixed rows: **76**
- Still unparseable: **0**
- Parse rate by patent: `{'WO2022204227': 1.0, 'WO2023147468': 1.0, 'WO2024064655': 1.0, 'WO2025207644': 1.0, 'WO2026054623': 1.0}`

### Activity corrections (vs PDF activity tables)
| Patent | Compound | Old CSV | Corrected |
|--------|----------|---------|-----------|
| WO2023147468 | 32 | + | ++ |
| WO2024064655 | 32 | + | +++ |
| WO2025207644 | 193 | (empty) | A |

Activity corrections applied: **3**

## Cleaned dataset

| File | Rows / content |
|------|----------------|
| `patent_bal_compounds_merged.csv` | 939 parseable compounds (main table) |
| `patent_bal_compounds_cleaned_full.csv` | full audit columns |
| `patent_smiles_fix_log.csv` | every auto-fixed SMILES |
| `patent_smiles_sample_qa.csv` | stratified sample checklist (48 compounds) |
| `qa_smiles_samples/` | rendered structure PNGs + per-patent grids |

## Automated sample checks

BAL-like substructure hit rate (parseable): `{'WO2022204227': 0.906, 'WO2023147468': 1.0, 'WO2024064655': 0.959, 'WO2025207644': 0.997, 'WO2026054623': 1.0}`

Cross-source InChIKey agreement (new vs older PatentPak dump, after same fixes):

| Patent | Same | Diff | Skipped |
|--------|------|------|---------|
| WO2022204227 | 383 | 0 | 0 |
| WO2023147468 | 75 | 0 | 0 |
| WO2024064655 | 169 | 0 | 0 |
| WO2025207644 | 286 | 0 | 0 |
| WO2026054623 | 26 | 0 | 0 |

### How to manually confirm against the patent PDF
1. Open `qa_smiles_samples/<PATENT>_grid.png`
2. Match each compound number to the patent PDF structure
3. Check ring systems, substituents, N-oxide, stereochemistry
4. Record pass/fail in `patent_smiles_sample_qa.csv`

## Notes / residual risk
- Auto-fix restores **chemical valence** for known PatentPak OCR patterns; it does **not** by itself prove identity with the patent drawing.
- WO2026054623 SMILES retain SEM protecting groups — confirm whether the activity table refers to protected intermediates or final products.
- WO2023147468 source quality was entirely `Low Confidence` before fix; fixed N-oxides are chemically plausible restorations and should be PDF spot-checked.
