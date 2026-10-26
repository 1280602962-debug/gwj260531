# Table S9. Chemical-diversity and activity-threshold sensitivity

化学多样性与活性阈值敏感性

## A. CAP5 membership

| Target pair | Original n | CAP5 n | Excluded n | Excluded dual | Excluded A-only | Excluded B-only | Excluded neither |
| --- | --- | --- | --- | --- | --- | --- | --- |
| EGFR/HER2 | 109 | 109 | 0 | 0 | 0 | 0 | 0 |
| JAK1/JAK2 | 110 | 107 | 3 | 0 | 3 | 0 | 0 |
| JAK1/TYK2 | 110 | 110 | 0 | 0 | 0 | 0 | 0 |
| PIK3CA/mTOR | 110 | 110 | 0 | 0 | 0 | 0 | 0 |
| AChE/BChE | 98 | 98 | 0 | 0 | 0 | 0 | 0 |
| F2/F10 | 110 | 110 | 0 | 0 | 0 | 0 | 0 |
| PPARG/PPARA | 110 | 110 | 0 | 0 | 0 | 0 | 0 |
| PPARA/PPARD | 110 | 107 | 3 | 0 | 1 | 2 | 0 |

## B. Vina summary min before and after CAP5

| Target pair | Primary summary min | Primary CI lower | Primary CI upper | CAP5 summary min | CAP5 CI lower | CAP5 CI upper | Point change |
| --- | --- | --- | --- | --- | --- | --- | --- |
| EGFR/HER2 | 0.2467 | 0.1328 | 0.3783 | 0.2467 | 0.1328 | 0.3783 | 0.0000 |
| JAK1/JAK2 | 0.5908 | 0.4512 | 0.7246 | 0.6315 | 0.4871 | 0.7597 | 0.0406 |
| JAK1/TYK2 | 0.4219 | 0.2812 | 0.5635 | 0.4219 | 0.2812 | 0.5635 | 0.0000 |
| PIK3CA/mTOR | 0.5801 | 0.4346 | 0.7129 | 0.5801 | 0.4346 | 0.7129 | 0.0000 |
| AChE/BChE | 0.6341 | 0.4607 | 0.7540 | 0.6341 | 0.4607 | 0.7540 | 0.0000 |
| F2/F10 | 0.3623 | 0.2305 | 0.4805 | 0.3623 | 0.2305 | 0.4805 | 0.0000 |
| PPARG/PPARA | 0.6021 | 0.4561 | 0.7051 | 0.6021 | 0.4561 | 0.7051 | 0.0000 |
| PPARA/PPARD | 0.5107 | 0.3652 | 0.6475 | 0.4969 | 0.3427 | 0.6354 | -0.0139 |

## C. PIK3CA/mTOR strict 6.5/5.5 subset

| Target pair | Scheme | Dual n | A-only n | B-only n | AUROC B | B CI lower | B CI upper | AUROC A | A CI lower | A CI upper | Summary min | Summary CI lower | Summary CI upper |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PIK3CA/mTOR | S1 | 25 | 12 | 14 | 0.7900 | 0.6200 | 0.9300 | 0.4571 | 0.2743 | 0.6429 | 0.4571 | 0.2743 | 0.6429 |
| PIK3CA/mTOR | S2 | 25 | 12 | 14 | 0.6200 | 0.4133 | 0.8133 | 0.5343 | 0.3429 | 0.7229 | 0.5343 | 0.3343 | 0.6829 |
| PIK3CA/mTOR | S3 | 25 | 12 | 14 | 0.7500 | 0.5767 | 0.8933 | 0.6286 | 0.4457 | 0.7971 | 0.6286 | 0.4429 | 0.7800 |
| PIK3CA/mTOR | S4 | 25 | 12 | 14 | 0.6267 | 0.4267 | 0.8167 | 0.6057 | 0.4143 | 0.7829 | 0.6057 | 0.3867 | 0.7229 |
| PIK3CA/mTOR | S5 | 25 | 12 | 14 | 0.5933 | 0.4033 | 0.7833 | 0.6686 | 0.4828 | 0.8371 | 0.5933 | 0.4000 | 0.7367 |

## D. Chemistry-model sensitivity to CAP5

| Target pair | Direction | Original OOF n | CAP5 OOF n | Original ECFP4 AUROC | CAP5 ECFP4 AUROC | ECFP4 AUROC change | Original Vina increment | CAP5 Vina increment |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EGFR/HER2 | D_vs_A | 65 | 65 | 0.8002 | 0.8002 | 0.0000 | -0.0164 | -0.0164 |
| EGFR/HER2 | D_vs_B | 59 | 59 | 0.8790 | 0.8790 | 0.0000 | -0.0023 | -0.0023 |
| JAK1/JAK2 | D_vs_A | 64 | 61 | 0.9185 | 0.9079 | -0.0106 | 0.0015 | 0.0048 |
| JAK1/JAK2 | D_vs_B | 64 | 64 | 0.9678 | 0.9678 | 0.0000 | 0.0029 | 0.0029 |
| JAK1/TYK2 | D_vs_A | 63 | 63 | 0.8488 | 0.8488 | 0.0000 | -0.0050 | -0.0050 |
| JAK1/TYK2 | D_vs_B | 63 | 63 | 0.9042 | 0.9042 | 0.0000 | -0.0010 | -0.0010 |
| PIK3CA/mTOR | D_vs_A | 64 | 64 | 0.8213 | 0.8213 | 0.0000 | 0.0215 | 0.0215 |
| PIK3CA/mTOR | D_vs_B | 64 | 64 | 0.7285 | 0.7285 | 0.0000 | -0.0029 | -0.0029 |
| AChE/BChE | D_vs_A | 52 | 52 | 0.8993 | 0.8993 | 0.0000 | -0.0044 | -0.0044 |
| AChE/BChE | D_vs_B | 55 | 55 | 0.8201 | 0.8201 | 0.0000 | -0.0053 | -0.0053 |
| F2/F10 | D_vs_A | 63 | 63 | 0.9405 | 0.9405 | 0.0000 | -0.0071 | -0.0071 |
| F2/F10 | D_vs_B | 63 | 63 | 0.6633 | 0.6633 | 0.0000 | 0.0363 | 0.0363 |
| PPARG/PPARA | D_vs_A | 63 | 63 | 0.8206 | 0.8206 | 0.0000 | -0.0141 | -0.0141 |
| PPARG/PPARA | D_vs_B | 64 | 64 | 0.6670 | 0.6670 | 0.0000 | -0.0059 | -0.0059 |
| PPARA/PPARD | D_vs_A | 64 | 63 | 0.9316 | 0.9345 | 0.0028 | -0.0029 | -0.0030 |
| PPARA/PPARD | D_vs_B | 64 | 62 | 0.8740 | 0.8688 | -0.0053 | 0.0049 | 0.0042 |


CAP5 is post hoc supporting sensitivity: at most five ligands per nonempty Murcko scaffold within pair and activity class; acyclic ligands are exempt. Six ligands are excluded without replacement. Complete CAP5 method comparisons, chemical models, controls and ranking outputs are retained in Supplementary_Data.zip. The PPARA/PPARD S5 dual-versus-A-only CI changes from [0.4921875, 0.767578125] to [0.5030241935, 0.7772177419]; this boundary change is retained, without a general significance claim. The strict Q110 subset has 56 members (25 dual, 12 A-only, 14 B-only, 5 neither). Main and subset populations differ; no new paired CI for their difference is calculated.
