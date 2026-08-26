# Document-blocked CV verdict

Grouping rule: ligands sharing any retained high-confidence ChEMBL `document_id` are one group.
The same folds are used for ECFP4, physicochemical, and docking logistic models.
Grouping was not changed after seeing AUROC.

| pair | contrast | n_pos/n_neg | groups | docs | valid folds | rank full | rank blocked mean-fold | ECFP4 | physchem | dock logistic | status |
|------|----------|-------------|-------:|-----:|------------:|----------:|-----------------------:|------:|---------:|--------------:|--------|
| EGFR/HER2 | D_vs_A | 28/38 | 28 | 156 | 5 | 0.6664 | 0.7344 | 0.6269 | 0.7256 | 0.5818 | ok |
| EGFR/HER2 | D_vs_B | 28/32 | 23 | 113 | 5 | 0.4297 | 0.4971 | 0.6228 | 0.3359 | 0.4688 | ok |
| AChE/BChE | D_vs_A | 27/25 | 39 | 97 | 5 | 0.6504 | 0.6617 | 0.8474 | 0.5956 | 0.6148 | ok |
| AChE/BChE | D_vs_B | 27/28 | 41 | 95 | 5 | 0.6058 | 0.6538 | 0.7672 | 0.7262 | 0.5397 | ok |
| PIK3CA/PIK3CB | D_vs_A | 28/27 | 33 | 43 | 5 | 0.6905 | 0.7339 | 0.6521 | 0.4974 | 0.5357 | ok |
| PIK3CA/PIK3CB | D_vs_B | 28/28 | 30 | 38 | 5 | 0.5 | 0.5007 | 0.8342 | 0.7551 | 0.3367 | ok |
| PIK3CA/mTOR | D_vs_A | 18/14 | 8 | 138 | 3 | 0.7143 | 0.6044 | 0.5139 | 0.5208 | 0.5903 | ok |
| PIK3CA/mTOR | D_vs_B | 18/12 | 9 | 129 | 1 | 0.6921 |  |  |  |  | cannot_stably_estimate |

Estimable directional arms: 7/8.
If an arm is not stably estimable, that is a result, not a reason to regroup.
