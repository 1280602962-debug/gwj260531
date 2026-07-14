# Patent BAL SAR feature tables

Generated from cleaned activity-table compounds (`patent_bal_compounds_merged.csv`, n=939).

## Activity label → IC50 range (from patent text)

### + / ++ / +++ (WO2022204227, WO2023147468, WO2024064655)
Patent wording is **nested**:
- `+++` = IC50 < 1 μM
- `++` = IC50 < 10 μM
- `+` = IC50 < 50 μM

Practical exclusive bins used in summaries:
- `+++` → IC50 < 1 μM
- `++` → 1 μM ≤ IC50 < 10 μM
- `+` → 10 μM ≤ IC50 < 50 μM

### A–E (WO2025207644)
- `A` = IC50 < 100 nM
- `B` = 100–500 nM
- `C` = 500 nM–1 μM
- `D` = 1–15 μM
- `E` = IC50 > 15 μM

### A/B/C (WO2026054623)
- `A` = IC50 < 50 nM
- `B` = 50 nM < IC50 < 100 nM
- `C` = IC50 > 100 nM (CSV header truncated after C)

## Region mapping
- **west_left**: western biaryl terminus (left in usual BAL drawing)
- **central_ring / central_alkoxy**: middle aryl + ortho alkoxy
- **linker**: amide N-substituent / benzylic linker
- **east_right**: indazole / azaindazole / N-oxide core

## Files
| File | Content |
|------|---------|
| `patent_activity_scale_definitions.csv` | Patent-literal IC50 ranges per label |
| `patent_activity_scale_practical_bins.csv` | Practical exclusive bins |
| `patent_compounds_region_annotated.csv` | All 939 compounds with region tags + IC50 bins |
| `patent_sar_by_region_feature.csv` | Per-region feature × activity stats |
| `patent_sar_feature_activity_counts.csv` | Counts/% per activity label |
| `patent_sar_west_alkoxy_east_combinations.csv` | west\|alkoxy\|east triad stats (n≥3) |
| `patent_activity_label_structure_summary.csv` | Each activity label → top structures |
