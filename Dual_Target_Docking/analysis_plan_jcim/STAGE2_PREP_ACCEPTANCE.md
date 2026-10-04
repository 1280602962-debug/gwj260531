# STAGE2_PREP_ACCEPTANCE

This report is read-only prep/acceptance only. It does not prove module implementations.
Implementation evidence is `jcim_stage2_impl_check.py` on artificial data.

SECOND_STAGE_COMPUTE_EXECUTED=NO

| check | result | detail |
|---|---|---|
| layer_counts_938_934_928 | PASS | 938/934/928 |
| fold_unassigned_labeled_10_7 | PASS | 10/7:AB_001,AB_053,AB_054,EH40_31,F2F10_080,J1TYK2_092,PGPA_030 |
| fold_unassigned_after_m0_6_4 | PASS | 6/4:EH40_31,F2F10_080,J1TYK2_092,PGPA_030 |
| alias_ab040_only | PASS | ['AB_040'] |
| ab040_not_in_pop | PASS |  |
| published_layer_table_matches_independent | PASS | /tmp/v4_2_final_audited_release_20260930/Dual_Target_Docking/analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv |
| wrongpocket_ache_extra_drop_0 | PASS | lab=27 m0=25 both=25 extra=0 |
| primary_metrics_120 | PASS | 120 |
| primary_delta_32 | PASS | 32 |
| population_807 | PASS | 807 |
| seed42_master_matches_formal_m0 | PASS | mismatches=0 |
| fiveseed_seed42_1614 | PASS | 1614 |
| fiveseed_seed42_matches_m0_master | PASS | mismatches=0 |
| ab046_seed17_4bds_timeout | PASS | {'job_id': 'AChE_BChE__4BDS__ZKOHODJSCICEFD-UHFFFAOYSA-N__seed17', 'pair': 'AChE/BChE', 'pdb_id': '4BDS', 'global_ligand_entity_id': 'ZKOHODJSCICEFD-UHFFFAOYSA-N', 'canonical_ligand_id': 'AB_046', 'seed': '17', 'primary_score': '0', 'status': 'TIMEOUT', 'runtime_s': '1586.067', 'vina_rc': 'TIMEOUT', 'n_retries': '1', 'mode1_affinity': '', 'n_modes_returned': '0', 'reason': 'JOB_TIMEOUT', 'timeout_s': '1586'} |
| alt_primary_scores_match_formal_m0 | PASS | mismatches=0 |
| alt_pdbs_frozen | PASS | ['1P0M', '2Y5F', '3SHC', '4EY6', '4JT5', '4L2Y', '5U46', '6KAX'] |
| m1b_m0_membership_equals_m1_m0 | PASS |  |
| m1b_three_pairwise_membership_identical | PASS | all_pairs_identical |
| fold_keys_unique | PASS | [] |
| scaffold_not_split_across_folds | PASS | [] |
| smiles_parse_808 | PASS | 808/808 |
| no_historical_imports_jcim_stage2_lib.py | PASS |  |
| no_historical_imports_jcim_stage2_audit.py | PASS |  |
| no_historical_imports_jcim_stage2_compute.py | PASS |  |
| synthetic_unit_checks | PASS |  |
| blocked_modules_gjk | PASS | ('G', 'J', 'K') |
| executable_modules_bcdefhi | PASS | ('B', 'C', 'D', 'E', 'F', 'H', 'I') |
| authority_files_present | PASS | presence_only_not_an_edit_check |
| authority_files_unchanged_vs_d365c901 | PASS |  |
| valid_flag_matches_finite_score | PASS | 0 |
| fold_class_matches_population | PASS | 0 |
| mapping_pair_id_unique_among_independent | PASS | [] |

Independent layer totals: labeled=938 m0=934 chemistry_oof=928
Labeled fold-unassigned records/ligands: 10/7 (AB_001, AB_053, AB_054, EH40_31, F2F10_080, J1TYK2_092, PGPA_030)
After-M0 fold-unassigned records/ligands: 6/4 (EH40_31, F2F10_080, J1TYK2_092, PGPA_030)

Population members were recounted in this audit file from the raw CSVs.
The audit does not call `jcim_stage2_lib.build_layer_members`.

A second raw-CSV recount that imported neither `jcim_stage2_lib` nor the audit helpers also returned 938/934/928.
Official eight pairs currently have identical dual / A_only / B_only members across M1b−M0, M1−M0, and M1b−M1. That is a read-only observation, not a ranking result.

Prep checks are not implementation evidence. See STAGE2_IMPL_ACCEPTANCE.md.

Blocked modules remain G, J, K.
