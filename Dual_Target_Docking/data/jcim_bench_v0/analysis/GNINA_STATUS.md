# GNINA rescore — STATUS: SKIPPED

**Date:** 2026-07-29  
**Pack scope:** JCIM docking phase (EH110, PM48 RDKit, AChE/BChE, PIK3CA/PIK3CB)

## Status

`STATUS: SKIPPED`

## Reason

No usable `gnina` binary found on this host (`which gnina` empty; no install under `/home/gwj/miniconda3` or `/home/gwj`). Agent command allows Skip when GNINA is unavailable.

## Impact

- Scoring channels for primary tables remain **Vina** + **RTMScore** (best-of-K).
- JCIM manuscript **Limitations** must state: sampling/scoring is single-engine (Vina) plus RTM rescore; CNN/GNINA channel not run — conclusions are not claimed to be score-function-invariant.
- Does **not** block Step 5 bench pack assembly.

## Resume condition

Install GNINA (GPU or CPU build), then CNN-rescore existing pose PDBQTs without re-docking. Cover packs listed above; write `tables/gnina_rescore_*.csv` and update `jcim_bench_v0` channel count.
