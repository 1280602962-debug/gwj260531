# Experiment scripts packed for upload

These scripts reproduce typeset tables and figures from deposited scores or saved poses.
Run them from the `Dual_Target_Docking/` tree. This folder is a flat copy for upload.
Do not upload assemble / audit / validate / checksum / freeze / pack helpers.

| Script | Repository path | Paper role |
|---|---|---|
| `build_t0_strengthen_v1.py` | `data/jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py` | Table 2 / S3 θ grid and ligand-chemistry baselines for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR from frozen Vina scores. |
| `benchmark_formulation_v1.py` | `data/jcim_novelty_v0/scripts/benchmark_formulation_v1.py` | Table 3, Table S4 fixed-channel Δ, and Table S5 incremental analyses for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR. |
| `analyze_five_pair_stack_v1.py` | `data/jcim_chembl_universe_v0/scripts/analyze_five_pair_stack_v1.py` | Table 2 / Table 3 / Table S4 / Table S5 for F2/F10, JAK1/JAK2, JAK1/TYK2, PPARG/PPARA, and PPARA/PPARD from frozen scores. Does not redock. |
| `build_p0_missing_tables_v1.py` | `data/jcim_strengthen_t0t1_v0/scripts/build_p0_missing_tables_v1.py` | Matched-versus-mismatched pocket bootstrap for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR (Table S6 / Figure 4). |
| `analyze_five_pair_local_channels_v1.py` | `data/jcim_chembl_universe_v0/scripts/analyze_five_pair_local_channels_v1.py` | Pocket channels, holdout join, JAK independent GNINA, and five-seed aggregate for F2/F10, JAK, and PPAR pairs (Tables S6–S7 / Figure 5). |
| `analyze_holdout_v1.py` | `data/jcim_holdout_v0/scripts/analyze_holdout_v1.py` | Unused-pool holdout AUROCs in Table S6 / Figure 4. Not external validation. |
| `analyze_independent_dock_v1.py` | `data/jcim_independent_dock_v0/scripts/analyze_independent_dock_v1.py` | Independent GNINA formulation AUROCs for EGFR/HER2 and PIK3CA/mTOR in Table S7 / Figure 5. |
| `analyze_multiseed_vina_v2.py` | `data/jcim_multiseed_v0/scripts/analyze_multiseed_vina_v2.py` | Five-seed AUROC aggregate for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR (Figure 5C). |
| `multiseed_fixed_membership_v1.py` | `data/jcim_chembl_universe_v0/scripts/multiseed_fixed_membership_v1.py` | Fixed-membership five-seed summary for F2/F10, JAK1/JAK2, JAK1/TYK2, PPARG/PPARA, and PPARA/PPARD (Figure 5C). Does not redock. |
| `equal_score_cluster_bootstrap_v1.py` | `data/jcim_novelty_v0/scripts/equal_score_cluster_bootstrap_v1.py` | Document/scaffold cluster resampling of the two flagship fixed-channel Δ values (Table S4 / Figure 6). |
| `claim_hardening_v1.py` | `data/jcim_novelty_v0/scripts/claim_hardening_v1.py` | Four-descriptor directional AUROCs for the Table 2 descriptor column and Figure S2. |
| `bindingdb_native_slice_eight_pairs_v1.py` | `data/jcim_novelty_v0/scripts/bindingdb_native_slice_eight_pairs_v1.py` | Independence-filtered BindingDB/PubChem eligibility counts for Table S8 / Figure 6. No external docking. |
| `eight_pair_ranking_operating_point_v1.py` | `data/jcim_novelty_v0/scripts/eight_pair_ranking_operating_point_v1.py` | Eight-pair vina_mean Top-10 and AND-filter operating points for Table S10. Does not redock. Figure 2D still uses operating_point_examples_review_v1.csv. |
| `assay_aggregation_max_vs_median_v1.py` | `data/jcim_novelty_v0/scripts/assay_aggregation_max_vs_median_v1.py` | API max-versus-median label sensitivity for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR (Table S3). |
| `high_confidence_label_rebuild_v1.py` | `data/jcim_novelty_v0/scripts/high_confidence_label_rebuild_v1.py` | High-confidence field screen for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR (Table S3). |
| `analyze_five_pair_dump_gated_v1.py` | `data/jcim_chembl_universe_v0/scripts/analyze_five_pair_dump_gated_v1.py` | Dump-gated max-versus-median footnote for F2/F10, JAK, and PPAR pairs (Table S3). |
| `analyze_eight_pair_dump_gated_v1.py` | `data/jcim_chembl_universe_v0/scripts/analyze_eight_pair_dump_gated_v1.py` | ChEMBL 37 dump join, max-versus-median, and leftover counts for all eight Table 2 pairs (Table S3). |
| `chembl_exhaustive_pair_census_v1.py` | `data/jcim_chembl_universe_v0/scripts/chembl_exhaustive_pair_census_v1.py` | ChEMBL pair-census summary plotted in Figure 1. |
| `ecfp4_docking_scaler_sensitivity_v1.py` | `data/jcim_novelty_v0/scripts/ecfp4_docking_scaler_sensitivity_v1.py` | Table S5 StandardScaler sensitivity on the same GroupKFold splits. |
| `reaudit_all14_cognate_rmsd_calcrrms_v1.py` | `data/jcim_novelty_v0/scripts/reaudit_all14_cognate_rmsd_calcrrms_v1.py` | Unified 14-receptor chemically mapped CalcRMS for Table S2 / Figure S4. Does not redock. |
| `reaudit_layer3_cognate_rmsd_v1.py` | `data/jcim_chembl_universe_v0/scripts/reaudit_layer3_cognate_rmsd_v1.py` | Atom-mapping / CalcRMS helper imported by the all-14 cognate RMSD script. |
| `replay_track_b_vina_mode1_v1.py` | `data/jcim_chembl_universe_v0/scripts/replay_track_b_vina_mode1_v1.py` | Re-reads REMARK VINA RESULT from the committed production pose tree. Does not redock. |
| `update_figures_pr32.py` | `figures/jcim_article/scripts/update_figures_pr32.py` | Official generator for Figures 1–6, S1–S4, and the TOC graphic from pinned CSVs. |
| `plot_jcim_article_figures_v3.py` | `figures/jcim_article/scripts/plot_jcim_article_figures_v3.py` | Panel functions imported by the official figure generator. |
| `jcim_figure_style.py` | `figures/jcim_article/scripts/jcim_figure_style.py` | Shared figure style used by the official figure generator. |
