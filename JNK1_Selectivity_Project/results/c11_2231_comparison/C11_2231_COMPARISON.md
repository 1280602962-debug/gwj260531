# C11 In-silico Comparison: 2231 (unbought) vs 690/2157 (purchased)

## Table

|   compound_id | HIT_ID       | purchased   | group   | smiles                                     |   score_JNK1 |   score_JNK2 |   score_JNK3 |   delta_sel_dock_table27 |   mmgbsa_JNK1 |   md_jnk1_bias_score |   hinge_JNK1 |   hinge_JNK2 |   hinge_JNK3 |   RMSD_JNK1 |   RMSD_JNK2 |   RMSD_JNK3 | pass_md_overall   | pose_grade   |   chemotype_sim |   AD_maxTc |   price_CNY_table27 |
|--------------:|:-------------|:------------|:--------|:-------------------------------------------|-------------:|-------------:|-------------:|-------------------------:|--------------:|---------------------:|-------------:|-------------:|-------------:|------------:|------------:|------------:|:------------------|:-------------|----------------:|-----------:|--------------------:|
|           690 | HIT103871685 | True        | G1      | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 |     -7.762   |     -6.68452 |     -4.47584 |                    1.077 |      -56.4616 |                0.97  |         1    |         0.51 |         0.77 |        0.72 |        1.98 |        0.61 | 是                 | A            |        0.232877 |   0.232877 |                 588 |
|          2157 | HIT101201113 | True        | G1      | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              |     -8.45736 |     -9.50699 |     -4.71357 |                   -1.05  |      -52.7413 |                1.03  |         0.85 |         0.46 |         0.02 |        0.49 |        1.13 |        0.35 | 是                 | A            |        0.225352 |   0.225352 |                2060 |
|          2231 | HIT100544184 | False       | G2      | COc1nc(NCc2ccccc2CN2CCCC2=O)ncc1F          |    -11.2204  |     -7.84746 |     -6.70699 |                    3.373 |      -54.4623 |                1.946 |         0.91 |         0    |         0.1  |        0.48 |        1.17 |        0.66 | 否                 | C            |        0.217949 |   0.229885 |                 832 |

## Narrative (locked with Option A)

- **Attractive for RQ-B:** Highest MD JNK1-bias score among shortlist; strongest hinge asymmetry (hinge JNK1≈0.91 vs JNK2≈0); best Glide score_JNK1 (−11.22).
- **Excluded for RQ-A gate:** pass_md_overall = NO / pose_grade C — fails project MD overall gate (requires JNK1 pass AND (JNK2 OR JNK3) pass). Purchase set prioritized pose-credible family binders (grade A), not the strongest MD-bias hypothesis.
- **Tradeoff:** Omitting 2231 weakens prospective test of MD-predicted JNK1 preference; 2157 is only a secondary bias hypothesis and has anti-JNK1 Δsel_dock (−1.05).
- **Paper use:** SI/Discussion: document opportunity cost; do not claim 2231 inactivity; optional future buy if budget allows.
