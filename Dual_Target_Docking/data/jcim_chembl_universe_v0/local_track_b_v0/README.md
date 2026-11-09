# Track B local Vina pack (`local_track_b_v0`)

Frozen AutoDock Vina 1.2.7 production scores for five pairs (F2/F10, JAK1/JAK2, JAK1/TYK2, PPARG/PPARA, PPARA/PPARD). Production pose PDBQT files remain under `poses/`. Input `ligands_pdbqt/` is not in this pack. Vina logs are not required for the zero-dock statistical chain.

**Protocol (exists):** `data/jcim_chembl_universe_v0/tables/track_b_local_run_v1.yaml`

**Production scores (exists):** `tables/scores_vina_mode1_v1.csv`

**Current statistical truth is not this folder.** Use:

- `results/canonical/current_score_master.csv`
- `results/canonical/primary_directional_auroc.csv`
- `results/canonical/primary_summary_min.csv`
- `results/canonical/five_seed_summary_min.csv`
- `results/canonical/computational_robustness.csv`

This pack does not replace Table 2. Do not regenerate from deleted scripts (`dock_track_b_production_v1.py`, `DOCKING_PLAN_V1.md`, `LOCAL_RECOMPUTE_PACK_V1.md` and the old analysis markdowns are git history only).
