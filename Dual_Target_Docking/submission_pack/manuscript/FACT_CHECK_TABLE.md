# 事实核查表（润色轮）

核对日期：2026-09-12。来源为当前分支有效 CSV，不是对话记忆。

| location | claim/value | source file or script | status | required action |
|---|---|---|---|---|
| Abstract; Table 2 EGFR/HER2 | 0.430 [0.282, 0.578]; D vs A 0.666; D vs B 0.430; n 28/38/32 | `unified_threshold_sensitivity_v2.csv` `theta_6.0` (0.4297 / 0.6664 / 0.4297; CI 0.2818–0.5775) | verified | none |
| Abstract; Table S4 EGFR pocket A | 0.430 vs 0.808; Δ 0.378 [0.205, 0.547] | `formulation_equal_score_negative_v1.csv` (0.4297, 0.808, 0.3783, 0.205–0.5469) | verified | none |
| Abstract; Table S4 JAK1/TYK2 pocket A | 0.444 [0.263, 0.620] | `equal_score_negative_s34_v1.csv` (0.4438, 0.2626–0.6201) | verified | none |
| Table 2 five pairs | F2 0.345 [0.211, 0.477]; JAK1/TYK2 0.365 [0.231, 0.503]; JAK1/JAK2 0.588 [0.444, 0.725]; PPARG/PPARA 0.649 [0.504, 0.751]; PPARA/PPARD 0.446 [0.296, 0.584] | `table2_comparable_theta6_v1.csv` | verified | none |
| Table 3 | EGFR 0.756 [0.562, 0.920] n_neither=12; JAK1/TYK2 0.770 [0.597, 0.906] n=14; PM 0.514 [0.222, 0.806] n=4 | `formulation_conventional_vs_directional_v1.csv`; five-pair Table 2 CSV | verified | none |
| Abstract / Table S5 | 16-direction max \|Δ\| = 0.023 | `ecfp4_incremental_s20s24_v1.csv` + original-pair increment; PPARA/PPARD D vs B −0.0234 (unscaled ECFP4 vs ECFP4+docking) | verified | none |
| Results 3.3 TPSA | 0.733 / 0.801 | `descriptor_all_four_directional_v1.csv` AChE/BChE `tpsa_D_vs_A` 0.7333, `tpsa_D_vs_B` 0.8009 | verified | none |
| Results 3.4 matched pocket | EGFR 0.170 [0.060, 0.280]; AChE 0.161 [0.037, 0.269] | `wrong_pocket_paired_delta_bootstrap_v1.csv` (0.1697, 0.06–0.2803; 0.1614, 0.037–0.269) | verified | none |
| Results 3.4 GNINA EGFR | 0.783 [0.610, 0.922]; 0.220 [0.109, 0.343] is dual–B-only | `independent_dock_formulation_v1.csv` (0.7825 [0.6104, 0.9222]; 0.2199 [0.1093, 0.3427]; `summary_min` row has empty CI) | verified | none |
| Results 3.1 census | 2,164,618 / 63,790 / 5,253 / 86 | `data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv` (`n_pairs_n_both_ge_1` / `n_pairs_n_both_ge_10` / `n_directional_n10` / `n_strict_thick`) | verified | none |
| Table 1 n_scored | JAK1/TYK2 31/32/32; F2 31/32/32; PPARG 32/31/32; AChE 27/25/28 | Table 2 CSVs | verified | none |
| Results 3.5 API max/median | EGFR 0.430 / 0.417 / 0.424, 93.6%; AChE median 0.629 | `assay_max_vs_median_agreement_v1.csv` | verified | none |
| Results 3.5 dump max/median | 1 flip CHEMBL121; summary_min unchanged to 3 d.p. | `five_pair_dump_gated_v1/max_vs_median_auroc_v1.csv` | verified | none |
| Results 3.6 BindingDB | no pair packaged as external set | `external_slice_summary_v1.csv` | verified | none |
| Results 3.6 year split | JAK1/TYK2 and JAK1/JAK2 reportable at 2018 | `five_pair_dump_gated_v1/time_split_v1.csv` | verified | none |
| Figure S5 simulation | 0.109 / 0.119 / 0.128 at true AUROC 0.50 | `detectable_effect_simulation_v1.csv`; SI Table S10 | verified in prior round | none |
| Results 3.4 holdout Δ range | −0.079 to +0.150 | five-pair `wrong_pocket_by_channel_v1.csv` `holdout_vina_20260727` (F2 −0.079; PPARA/PPARD +0.150) plus original-pair holdouts | verified | none |
| Results 3.4 PPARG rescore | RTM 0.369 [0.233, 0.475]; GNINA CNN 0.500 | `table2_comparable_by_channel_v1.csv` `rtm_best9` 0.3691 [0.233, 0.4753]; `gnina_cnn_affinity` 0.5 | verified | none |
| Results 3.4 PPARA RMSD | top-3 7.848 Å; best 1.401 Å | `layer3_cognate_rmsd_calcrrms_v1.csv` 6LXA `calcrrms_top3_A` 7.848, `calcrrms_best_A` 1.401 | verified | none |
| Results 3.2 EGFR filter | Top-10 1/5/4; precision 0.298; n=98 | `and_filter_operating_point_v1.csv` dual_percentile=50: n_dual_pass=14, n_A=9, n_B=24, precision_dual=0.2979; Top-10 from SI Table S13 / prior lock | verified (filter); Top-10 unchanged from SI | none |

未改任何科学数据文件。
