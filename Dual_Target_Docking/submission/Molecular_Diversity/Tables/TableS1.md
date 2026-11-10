# Table S1 Evaluation-set construction and activity provenance

## A. Selection conditions

| Target pair | Selection threshold | Chemical selection record | Deposited quota / composition | Eligible n | Selection replay status |
| --- | --- | --- | --- | --- | --- |
| EGFR/HER2 | 6 | Murcko cap 5 within class | deposited N=110 with classes 28/38/32/12 before the later unresolved-arm exclusion | 109 | Documented; not replayed |
| JAK1/JAK2 | 6.5 / 5.5 | Small-molecule / size filters; no Murcko cap | dual:32 / A-only: 32 / B-only: 32 / neither: 14 | 110 | Documented; not replayed |
| JAK1/TYK2 | 6.5 / 5.5 | Small-molecule / size filters; no Murcko cap | dual:32 / A-only: 32 / B-only: 32 / neither: 14 | 110 | Documented; not replayed |
| PIK3CA/mTOR | 6 | Fixed 110-molecule membership | 32 / 32 / 32 / 14 | 110 | Not independently replayed in this compilation |
| AChE/BChE | 6.5 / 5.5 | ChEMBL ID-prefix diversity key in deposited script | N=100; dual / A-only / B-only / neither 28/28/28/16 | 98 | Documented; not replayed |
| F2/F10 | 6.5 / 5.5 | Small-molecule / size filters; no Murcko cap | dual:32 / A-only: 32 / B-only: 32 / neither: 14 | 110 | Documented; not replayed |
| PPARG/PPARA | 6.5 / 5.5 | Small-molecule / size filters; no Murcko cap | dual:32 / A-only: 32 / B-only: 32 / neither: 14 | 110 | Documented; not replayed |
| PPARA/PPARD | 6.5 / 5.5 | Small-molecule / size filters; no Murcko cap | dual:32 / A-only: 32 / B-only: 32 / neither: 14 | 110 | Documented; not replayed |

## B. Activity-source coverage

| Target pair | n | Both ends traced | Complete identifiers | Missing identifiers | Source / endpoint mix | Threshold crossing (n) | Value / label mismatch |
| --- | --- | --- | --- | --- | --- | --- | --- |
| EGFR/HER2 | 109 | 109 | 101 | 8 | 22 | 13 | 0 |
| JAK1/JAK2 | 110 | 110 | 109 | 1 | 9 | 0 | 0 |
| JAK1/TYK2 | 110 | 110 | 108 | 2 | 13 | 4 | 0 |
| PIK3CA/mTOR | 110 | 110 | 109 | 1 | 24 | 3 | 0 |
| AChE/BChE | 98 | 98 | 98 | 0 | 8 | 2 | 0 |
| F2/F10 | 110 | 110 | 110 | 0 | 10 | 0 | 0 |
| PPARG/PPARA | 110 | 110 | 110 | 0 | 42 | 9 | 0 |
| PPARA/PPARD | 110 | 110 | 110 | 0 | 40 | 4 | 0 |

All final labels use pChEMBL >= 6.0. Provenance flags overlap. Complete identifiers require DOI, PubMed or patent identifiers for all tied maxima at both endpoints. Across all pairs, cross-document aggregation is present in 81 observations and endpoint heterogeneity in 123 (union 168). The full source crosswalk, raw activity rows, selection evidence and exceptions are included in Supplementary_Data.zip. CAP5 is a post hoc sensitivity analysis and does not redefine initial selection
