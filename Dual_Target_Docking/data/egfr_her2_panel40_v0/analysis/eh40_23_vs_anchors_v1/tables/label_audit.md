# Label audit
## EH40_23 (CHEMBL3633938)
Panel: **A_only** (EGFR=6.35, HER2=5.57; threshold 6.0)

### API (EGFR CHEMBL203 / HER2 CHEMBL1824)
| target | n | n_pchembl | max | median | >=6 | types |
|--------|---|-----------|-----|--------|-----|-------|
| EGFR | 1 | 1 | 6.35 | 6.35 | 1 | IC50 |
| HER2 | 1 | 1 | 5.57 | 5.57 | 0 | IC50 |

Both ends: biochemical (B) HotSpot IC50; origin annotated unknown; no mutant keywords in descriptions.

### Conclusion: **label_ok**
- API pchembl exactly matches panel (6.35 / 5.57) → A_only at threshold 6
- Both ends measured; not mutant-driven artifact in these rows
- Caveat: only 1 public activity/end (sparse), but values are consistent
- **Not dual**

## EH40_01 (TAK-285)
- panel: dual EGFR=9.0 HER2=8.52
- API max EGFR/HER2: 9.0 / 8.52 (n_pchembl 5/4)
- dual by API max≥6 both ends: **True** → dual anchor **robust**

## EH40_02 (LAPATINIB)
- panel: dual EGFR=10.22 HER2=8.8
- API max EGFR/HER2: 10.22 / 8.8 (n_pchembl 138/68)
- dual by API max≥6 both ends: **True** → dual anchor **robust**

