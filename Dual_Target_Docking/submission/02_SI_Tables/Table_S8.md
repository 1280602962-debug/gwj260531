# Table S8. Random-seed and receptor-structure sensitivity

随机种子与受体结构敏感性

## A. Five-seed distributions

| Target pair | Metric | Seed 42 | Median | Minimum | Maximum | IQR | Finite seeds |
| --- | --- | --- | --- | --- | --- | --- | --- |
| EGFR/HER2 | auc_B | 0.6274 | 0.6375 | 0.6274 | 0.6477 | 0.0111 | 5 |
| EGFR/HER2 | auc_A | 0.2467 | 0.2388 | 0.2305 | 0.2467 | 0.0089 | 5 |
| EGFR/HER2 | summary_min | 0.2467 | 0.2388 | 0.2305 | 0.2467 | 0.0089 | 5 |
| JAK1/JAK2 | auc_B | 0.5908 | 0.5840 | 0.5684 | 0.5908 | 0.0093 | 5 |
| JAK1/JAK2 | auc_A | 0.7485 | 0.7485 | 0.7422 | 0.7573 | 0.0112 | 5 |
| JAK1/JAK2 | summary_min | 0.5908 | 0.5840 | 0.5684 | 0.5908 | 0.0093 | 5 |
| JAK1/TYK2 | auc_B | 0.6475 | 0.6475 | 0.6328 | 0.6631 | 0.0156 | 5 |
| JAK1/TYK2 | auc_A | 0.4219 | 0.4229 | 0.4121 | 0.4302 | 0.0015 | 5 |
| JAK1/TYK2 | summary_min | 0.4219 | 0.4229 | 0.4121 | 0.4302 | 0.0015 | 5 |
| PIK3CA/mTOR | auc_B | 0.7363 | 0.7354 | 0.6997 | 0.7402 | 0.0029 | 5 |
| PIK3CA/mTOR | auc_A | 0.5801 | 0.5752 | 0.5596 | 0.5869 | 0.0190 | 5 |
| PIK3CA/mTOR | summary_min | 0.5801 | 0.5752 | 0.5596 | 0.5869 | 0.0190 | 5 |
| AChE/BChE | auc_B | 0.6341 | 0.6422 | 0.6049 | 0.6533 | 0.0163 | 5 |
| AChE/BChE | auc_A | 0.6442 | 0.6442 | 0.6356 | 0.6448 | 0.0026 | 5 |
| AChE/BChE | summary_min | 0.6341 | 0.6356 | 0.6049 | 0.6448 | 0.0075 | 5 |
| F2/F10 | auc_B | 0.4082 | 0.4062 | 0.3955 | 0.4082 | 0.0015 | 5 |
| F2/F10 | auc_A | 0.3623 | 0.3647 | 0.3594 | 0.3682 | 0.0029 | 5 |
| F2/F10 | summary_min | 0.3623 | 0.3647 | 0.3594 | 0.3682 | 0.0029 | 5 |
| PPARG/PPARA | auc_B | 0.6553 | 0.6553 | 0.6465 | 0.6680 | 0.0127 | 5 |
| PPARG/PPARA | auc_A | 0.6021 | 0.6313 | 0.6021 | 0.6406 | 0.0225 | 5 |
| PPARG/PPARA | summary_min | 0.6021 | 0.6313 | 0.6021 | 0.6406 | 0.0225 | 5 |
| PPARA/PPARD | auc_B | 0.6787 | 0.6929 | 0.6787 | 0.6953 | 0.0039 | 5 |
| PPARA/PPARD | auc_A | 0.5107 | 0.5088 | 0.4800 | 0.5151 | 0.0146 | 5 |
| PPARA/PPARD | summary_min | 0.5107 | 0.5088 | 0.4800 | 0.5151 | 0.0146 | 5 |

## B. One-receptor substitutions

| Target pair | Alternative PDB | Replaced side | Common dual n | Common A-only n | Common B-only n | Original summary min | Alternative summary min | Summary difference | Difference CI lower | Difference CI upper | Affected AUROC difference | Unreplaced AUROC difference |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AChE/BChE | 1P0M | B | 26 | 23 | 26 | 0.6037 | 0.6003 | -0.0033 | -0.0803 | 0.0836 | -0.0033 | 0 |
| AChE/BChE | 4EY6 | A | 27 | 23 | 26 | 0.6184 | 0.5670 | -0.0514 | -0.1852 | 0 | -0.0940 | 0 |
| F2/F10 | 2Y5F | B | 30 | 32 | 31 | 0.3828 | 0.3828 | 0 | -0.0469 | 0.0229 | -0.0146 | 0 |
| F2/F10 | 3SHC | A | 30 | 32 | 31 | 0.3828 | 0.4333 | 0.0505 | 0 | 0.1441 | 0.0780 | 0 |
| PPARA/PPARD | 5U46 | B | 32 | 32 | 32 | 0.5107 | 0.5107 | 0 | -0.0225 | 0 | -0.0181 | 0 |
| PPARA/PPARD | 6KAX | A | 32 | 32 | 32 | 0.5107 | 0.4438 | -0.0669 | -0.1484 | 0.0044 | -0.0669 | 0 |


Five Vina seeds: 17, 29, 42, 71 and 101. Each seed uses its own valid population. Receptor comparisons use paired common members and replace only one target. Full direction-specific estimates, confidence intervals and members are in the source data. PIK3CA/mTOR Q110 has no alternative-receptor calculation.
