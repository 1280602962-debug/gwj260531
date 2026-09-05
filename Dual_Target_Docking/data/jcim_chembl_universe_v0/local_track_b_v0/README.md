# Track B local Vina pack (`local_track_b_v0`)

Layer-3 cognate QC (8/8 PASS at E=8) and production AutoDock Vina 1.2.7 on five pairs (~1100 jobs).

- Protocol: `../analysis/DOCKING_PLAN_V1.md`, lock `../tables/track_b_local_run_v1.yaml`
- Results: `analysis/TRACK_B_DIRECTIONAL_AUROC_V1.md`
- Descriptor reference: `analysis/TRACK_B_DESCRIPTOR_REFERENCE_V1.md`
- Scores: `tables/scores_vina_mode1_v1.csv` (1094 successes; 6 timeout skips)
- Poses / full Vina logs are gitignored (regenerate with `scripts/dock_track_b_production_v1.py`)

Does **not** replace Table 2. Count three systems (coagulation, JAK, PPAR).
