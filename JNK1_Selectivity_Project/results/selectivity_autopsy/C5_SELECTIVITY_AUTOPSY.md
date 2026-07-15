# C5 Selectivity-Method Autopsy (Main-Text Ready)

> **Purchase decoupling:** Δsel / Gly87 / ML selectivity labels were **not** purchase hard gates.

## Main-text summary table

| method                                                                | metric                       | value                                                                              | pass_threshold                            | verdict                                         | used_for_purchase                                 |
|:----------------------------------------------------------------------|:-----------------------------|:-----------------------------------------------------------------------------------|:------------------------------------------|:------------------------------------------------|:--------------------------------------------------|
| Glide Δsel_dock direction (archived project report)                   | direction accuracy           | 43% (3/7) VSW single-PDB; 29% (2/7) ensemble (project report)                      | ≥55%                                      | FAIL                                            | NO — decoupled; family shortlist only             |
| Glide Δsel_dock direction (recomputed from benchmark_deltas_51c1.csv) | direction accuracy           | 2/8 = 25.0% on benchmark_deltas CSV (excl. NA dirs)                                | ≥55%                                      | FAIL                                            | NO                                                |
| Glide Δsel_dock (key controls SP600125/TCS/CC-930/E1)                 | direction accuracy           | 2/4 = 50.0%                                                                        | ≥55%                                      | FAIL                                            | NO                                                |
| Gly87 (KLIFS b.l.37) occupancy heuristic                              | discriminative power         | occ_JNK1 True for 5/5 benchmarks; pred_JNK1_sel any=False; d_occ range 0.59–1.18 Å | separates JNK1-preferring vs pan/opposite | FAIL (non-discriminative)                       | NO                                                |
| ML isoform-selectivity classifier (ChEMBL paired)                     | test F1 (selective class)    | 0 (positives n_train≈8; reported in training_report / project report)              | usable precision/recall for purchase      | FAIL                                            | NO — ML used only as family pActivity recall gate |
| ML family activity gate p_family≥6.0                                  | benchmark recall / decoy FPR | recall 9/9; Taosu decoy FPR 95.3%; EF1%=9.20                                       | high recall OK for early filter           | PASS as recall filter; NOT a selectivity filter | YES (activity recall only)                        |

## Per-compound docking direction

| compound   | expected_profile    |   delta_sel_dock | exp_dir_pIC50   | pred_dir_dock   | direction_match   |   jnk1_ic50_nM |   jnk2_ic50_nM |   jnk3_ic50_nM |
|:-----------|:--------------------|-----------------:|:----------------|:----------------|:------------------|---------------:|---------------:|---------------:|
| AS602801   | pan-JNK             |          -3.0975 | JNK1            | JNK23           | False             |           80   |           90   |          230   |
| CC-401     | unknown-isoform     |          -4.902  | nan             | JNK23           | False             |          nan   |          nan   |          nan   |
| SP600125   | pan-JNK             |          -2.5735 | JNK1            | JNK23           | False             |           40   |           40   |           90   |
| TCS JNK 6O | JNK1-preferring     |          -1.181  | JNK1            | JNK23           | False             |           45   |          160   |          nan   |
| CC-930     | JNK2/JNK3-biased    |          -4.896  | JNK23           | JNK23           | True              |           61   |            7   |            6   |
| JNK-IN-8   | JNK3-preferring     |           0.4125 | JNK1            | WEAK            | False             |            4.7 |           18.7 |            1   |
| CC-90001   | pan-JNK             |          -0.925  | JNK1            | JNK23           | False             |           11   |           31   |          nan   |
| Q63        | JNK1/JNK3-over-JNK2 |          -1.0325 | JNK1            | JNK23           | False             |           33.5 |          112.9 |           33.2 |
| E1         | JNK1-preferring     |           3.045  | JNK1            | JNK1            | True              |            2.7 |           19   |            9   |

## Gly87 self-check

| ligand     |   d_occ | occ_JNK1   | occ_method           |   d_clash_JNK2 | clash_JNK2   |   d_clash_JNK3 | clash_JNK3   | pose_ok   | pose_notes                                                | pred_JNK1_sel   | exp_profile   | exp_fold                 | match   |
|:-----------|--------:|:-----------|:---------------------|---------------:|:-------------|---------------:|:-------------|:----------|:----------------------------------------------------------|:----------------|:--------------|:-------------------------|:--------|
| E1         |   0.744 | True       | superpose_JNK2_Ser87 |         10.112 | False        |          9.856 | False        | True      | hinge_HB_dist=2.91Å(OK); centroid_dist=8.39Å(limit=36.8)  | False           | JNK1偏好        | JNK2/1=7.0x; JNK3/1=3.3x | False   |
| TCS JNK 6O |   0.899 | True       | superpose_JNK2_Ser87 |          8.691 | False        |          9.14  | False        | True      | hinge_HB_dist=3.04Å(OK); centroid_dist=10.15Å(limit=36.8) | False           | JNK1偏好        | JNK2/1=3.6x              | False   |
| CC-930     |   1.18  | True       | superpose_JNK2_Ser87 |          9.362 | False        |          9.107 | False        | True      | hinge_HB_dist=2.80Å(OK); centroid_dist=7.22Å(limit=36.8)  | False           | JNK2/3偏好(反向)  | JNK2/1=0.1x; JNK3/1=0.1x | True    |
| SP600125   |   1.009 | True       | superpose_JNK2_Ser87 |          9.329 | False        |         11.204 | False        | True      | hinge_HB_dist=3.09Å(OK); centroid_dist=7.14Å(limit=36.8)  | False           | pan           | JNK2/1=1.0x; JNK3/1=2.2x | True    |
| CC-90001   |   0.59  | True       | superpose_JNK2_Ser87 |          8.868 | False        |          9.016 | False        | True      | hinge_HB_dist=3.04Å(OK); centroid_dist=11.74Å(limit=36.8) | False           | 近pan(弱2.8x)   | JNK2/1=2.8x              | True    |

## Suggested Results paragraph (English draft)
On a literature JNK benchmark panel, Glide-derived Δsel_dock failed a 55% direction-accuracy usability threshold (archived VSW single-PDB 43% [3/7]; recomputed ensemble table 2/8 = 25%). A Gly87 occupancy heuristic labeled all tested benchmarks as JNK1-occupying and did not separate JNK1-preferring controls from pan or reverse profiles. An ML selective-class model trained on sparse ChEMBL positives yielded test F1 = 0. Accordingly, these filters were retained only as negative controls and were not used to purchase the shortlist; candidate selection prioritized family-activity and pose/MD QC.
On a literature JNK benchmark panel, Glide-derived Δsel_dock reproduced experimental isoform direction for only 2/8 compounds (25% if available else 'low accuracy'). A Gly87 occupancy heuristic labeled all tested benchmarks as JNK1-occupying and did not separate JNK1-preferring controls from pan or reverse profiles. An ML selective-class model trained on sparse ChEMBL positives yielded test F1 = 0. Accordingly, these filters were retained only as negative controls and were not used to purchase the shortlist; candidate selection prioritized family-activity and pose/MD QC.
