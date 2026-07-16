# C11 Purchase Panel: 690 + 2231 (updated)

## Table

|   compound_id | HIT_ID       | purchased   | role                        | group   | smiles                                     |   score_JNK1 |   score_JNK2 |   score_JNK3 |   delta_sel_dock_table27 |   mmgbsa_JNK1 |   md_jnk1_bias_score |   hinge_JNK1 |   hinge_JNK2 |   hinge_JNK3 |   RMSD_JNK1 |   RMSD_JNK2 |   RMSD_JNK3 | pass_md_overall   | pose_grade   |   chemotype_sim |   price_CNY_table27 |
|--------------:|:-------------|:------------|:----------------------------|:--------|:-------------------------------------------|-------------:|-------------:|-------------:|-------------------------:|--------------:|---------------------:|-------------:|-------------:|-------------:|------------:|------------:|------------:|:------------------|:-------------|----------------:|--------------------:|
|           690 | HIT103871685 | True        | activity/pan-leaning anchor | G1      | Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1 |     -7.762   |     -6.68452 |     -4.47584 |                    1.077 |      -56.4616 |                0.97  |         1    |         0.51 |         0.77 |        0.72 |        1.98 |        0.61 | 是                 | A            |        0.232877 |                 588 |
|          2231 | HIT100544184 | True        | MD JNK1-bias hypothesis     | G2      | COc1nc(NCc2ccccc2CN2CCCC2=O)ncc1F          |    -11.2204  |     -7.84746 |     -6.70699 |                    3.373 |      -54.4623 |                1.946 |         0.91 |         0    |         0.1  |        0.48 |        1.17 |        0.66 | 否                 | C            |        0.217949 |                 832 |
|          2157 | HIT101201113 | False       | not purchased (former alt)  | G1      | Cc1cnc(NCc2cccc3c2OCCCO3)nc1C              |     -8.45736 |     -9.50699 |     -4.71357 |                   -1.05  |      -52.7413 |                1.03  |         0.85 |         0.46 |         0.02 |        0.49 |        1.13 |        0.35 | 是                 | A            |        0.225352 |                2060 |

## Narrative

- **690:** G1 grade-A, pass_md_overall, pan-leaning hinge (high on all isoforms) — RQ-A family-activity anchor with credible pose QC.
- **2231:** Strongest MD JNK1-bias score; hinge J1≫J2; favorable Δsel_dock (+3.37) and best score_JNK1 (−11.22). Purchased to give RQ-B a real prospective test.
- **Risk:** 2231 pose_grade C and pass_md_overall = NO under archived MD gate. Extended 200 ns used ligand restraints — must not be cited as unrestrained proof. Week 2–3 priority: unrestrained MD replicas (C3) before over-interpreting bias.
- **2157:** Budget/slot tradeoff; 2157 had MD bias #2 but anti-JNK1 Δsel_dock (−1.05). Kept as in-silico comparator only.
- **Option A:** Core contribution remains selectivity-predictor failure + family pipeline. 2231 upgrades secondary RQ-B; null/mixture outcomes remain publishable.
