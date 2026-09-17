# Document-blocked CV verdict

Grouping rule: ligands sharing any retained high-confidence ChEMBL `document_id` are one group.
The same folds are used for ECFP4, physicochemical, and docking logistic models.
Grouping was not changed after seeing AUROC.

| pair | contrast | n_pos/n_neg | groups | docs | valid folds | rank full | rank blocked mean-fold | ECFP4 | physchem | dock logistic | status |
|------|----------|-------------|-------:|-----:|------------:|----------:|-----------------------:|------:|---------:|--------------:|--------|
| EGFR/HER2 | D_vs_A | 28/38 | 28 | 156 | 5 | 0.6607 | 0.6466 | 0.6344 | 0.7359 | 0.594 | ok |
| EGFR/HER2 | D_vs_B | 28/32 | 23 | 113 | 5 | 0.3237 | 0.4694 | 0.6161 | 0.3415 | 0.5346 | ok |
| AChE/BChE | D_vs_A | 27/26 | 40 | 97 | 5 | 0.6524 | 0.6581 | 0.849 | 0.5769 | 0.6353 | ok |
| AChE/BChE | D_vs_B | 27/28 | 41 | 95 | 5 | 0.6058 | 0.6029 | 0.8095 | 0.7513 | 0.5013 | ok |
| PIK3CA/mTOR | D_vs_A | 18/14 | 8 | 138 | 3 | 0.7143 | 0.7711 | 0.5347 | 0.4722 | 0.6736 | ok |
| PIK3CA/mTOR | D_vs_B | 18/12 | 9 | 129 | 1 | 0.6921 |  |  |  |  | cannot_stably_estimate |

Estimable directional arms: 5/6.
If an arm is not stably estimable, that is a result, not a reason to regroup.
