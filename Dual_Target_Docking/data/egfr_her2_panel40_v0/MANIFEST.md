# MANIFEST — egfr_her2_panel40_v0 freeze
- created_utc: 2026-07-27T04:44:16.501083+00:00
- PROJECT_ROOT: `/mnt/d/CADD paper exercise/dual target docking`
- OUT: `/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0`
- source vina: `/mnt/d/CADD paper exercise/dual target docking/Maestro doc/vina_docking`

## Completeness
- panel: OK
- scores_vina: OK
- scores_rtm: OK
- poses 3POZ ligands with files: 40/40
- poses 3RCD ligands with files: 40/40
- vina log files copied: 162
- vina seeds recorded: 41

## Pose issues
- none (40×2 ligands, 9 modes each)

## Missing sources
- none

## Copied files (abbrev)
- total copy records: 1080
- key tables: panel_v0_40.csv, scores_vina.csv, scores_rtm.csv, scores_rtm_all_poses.csv
- poses layout: `poses/<3POZ|3RCD>/<EH40_XX>/mode_01.pdbqt` … `mode_09.pdbqt` + `*_all_modes.pdbqt`
- protocol: see `protocol/protocol.yaml`

## Ablation / architecture / gating outputs
- tables/panel_v0_40_arch.csv
- tables/ablation_metrics.csv
- tables/ablation_ranks.csv
- tables/ablation_ranks_wide.csv
- tables/ablation_ligand_scores.csv
- tables/hardneg_cases.md
- scripts/build_ablation_table.py
- protocol/protocol.yaml

## Known gaps
- Global docking seed was NOT fixed (per-job random; see logs/vina_seeds.json — partial capture from logs that contain 'random seed').
- PoseBusters / GNINA / PLIP not available; structure gate = simple steric clash @ 2.2 Å (triggered 0 ligands).
- Architecture: only EH40_01/02/05 labeled merged; remaining 37 = unknown (no invented fused/linked cuts).
- No second target pair started (no PIK3CA/mTOR).

## Update 2026-07-27 — exhaustiveness sensitivity v1
- Added `protocol/seeds_as_run.csv` (80/80 as-run seeds) and `protocol/SEED_POLICY.md`
- Updated `protocol/protocol.yaml`: fixed_global seed=20260727; exhaustiveness_v0_1=8
- Added `analysis/exhaustiveness_sensitivity_v1/` (Exp A/B poses, tables, SENSITIVITY_VERDICT.md)
- No second target pair; no LigPrep redo; no full-panel rerun

## Decision ablation + bootstrap CI v0
- `analysis/decision_ablation_v0/` — five frozen arms (same thresholds as panel48)
- `analysis/bootstrap_ci_v0/` + `../cross_pair_bootstrap_v0/`
