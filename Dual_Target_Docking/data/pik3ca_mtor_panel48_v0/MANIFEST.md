# MANIFEST — pik3ca_mtor_panel48_v0

- exhaustiveness: **16** (`exhaustiveness_v0_1`)
- seed_fixed_global: 20260727
- cognate QC E=8: No-Go (archived) — `analysis/cognate_redock_v0/COGNATE_QC_VERDICT.md`
- cognate QC E=16: **Go** — `analysis/cognate_redock_v0/COGNATE_QC_VERDICT_E16.md`
- E=8 cognate poses archived at: `poses/cognate_E8_archive/`
- cognate E16 poses: `poses/cognate_E16/`
- full panel poses: `poses/4L23|4JT6/<PM48_XX>/`
- note: 4JT6 mode1 may fail redock RMSD while best_of_9 passes → keep 9 modes; plan RTM best-of-9

## Full panel status
- 96/96 jobs attempted @ E=16 seed=20260727
- 95 jobs with 9 modes; **4JT6/PM48_34 has 8 valid modes** (Vina phantom mode9 +75.84; see logs/vina/4JT6_PM48_34_NOTE.md)
- scores: tables/scores_vina.csv, tables/scores_vina_long.csv, tables/job_status.csv

## RTMScore + ablation
- pockets: `receptors/4L23_pocket_10.0.pdb`, `receptors/4JT6_pocket_10.0.pdb` (10 Å around X6K)
- model: RTMScore model1
- poses scored: 863 (4L23 432 + 4JT6 431)
- script: `scripts/run_rtm_and_ablation.py`
- tables: `scores_rtm.csv`, `scores_rtm_all_poses.csv`, `scores_rtm_4L23.csv`, `scores_rtm_4JT6.csv`
- ablation: `ablation_metrics.csv`, `ablation_ligand_scores.csv`, `ablation_ranks.csv`, `ablation_key_ligand_ranks.csv`
- primary arm (AUROC Dual vs rest): **rtm_min_z = 0.685** (vina_mean 0.633)
