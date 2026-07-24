# Redock smoke: lesinurad × 9DKB

Protocol comparison for redocking (Vina / gnina / RTMScore), exhaustiveness 8 and 32.

## Key outputs
- `redock_results_lesinurad_9DKB.csv` — filled metrics table (P0–P5)
- `lesinurad_9DKB/` — prep, vina/gnina poses, RTMScore scores, logs, analysis

## Notes
- RMSD: heavy-atom pose RMSD in protein frame (symmetry-aware), ≤2 Å pass/fail
- gnina: `--cnn_scoring rescore`, seed 42, CPU (`--no_gpu`)
- RTMScore: model1, pocket 10 Å from crystal ligand
- **Analysis / conclusions:** [`docs/REDOCK_SMOKE_ANALYSIS.md`](../../docs/REDOCK_SMOKE_ANALYSIS.md)
