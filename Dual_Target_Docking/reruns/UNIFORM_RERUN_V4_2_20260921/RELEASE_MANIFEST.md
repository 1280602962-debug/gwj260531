# RELEASE_MANIFEST

branch: cursor/v4-2-final-audited-results-20260930
kind: orphan_release
source_worktree: /tmp/pr39_fiveseed (read-only)
run_root: Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921

## Stage
PRIMARY metrics completed.
POST_PRIMARY_METRICS_QA = PASS.
Sensitivity analyses have not been run.

## Formal methods
- M0: official primary seed42 Vina mode1 (`vina_score = -minimizedAffinity`)
- M1: VERIFIED Vina mode1 pose GNINA CNNscore
- M1b: VERIFIED max CNNscore on the same saved Vina pose set
- M2: RTMScore model1 on Vina mode1 pose
- M3: TOPOLOGY_VERIFIED GNINA independent docking seed42

Old M1/M3 tables inside 09/10 are HISTORICAL_NON_FORMAL provenance only.

## Gates (copied from source)
- END_TO_END_SCIENTIFIC_INTEGRITY_GATE = PASS
- PRE_AUROC_GATE = PASS
- PRE_METRICS_ANALYSIS_GATE = PASS
- METRICS_IMPLEMENTATION_GATE = PASS
- POST_PRIMARY_METRICS_QA = PASS

## Ligands
02_ligands in this release is a real directory.
785 SDF + 785 PDBQT copied from V4.1 ligand store and byte-compared.
Source V4.2 02_ligands symlink was not modified.

## Forbidden for formal V4.2 analysis
- scripts/analysis/*
- results/canonical/current_score_master.csv
- old AUROC / old sensitivity
These forbidden references are not copied into this orphan branch.

## Formal result tables
- Dual_Target_Docking/results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv (120 rows)
- Dual_Target_Docking/results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv (32 rows)
