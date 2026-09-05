# C5 docking handoff status

## Live
- Queue waiter: `scripts/c5_queue_after_current_vina.sh` (PID via `pgrep -af c5_queue_after_current_vina`)
- Log: `data/campaigns/c5/00_verification/c5_handoff_nohup.out`
- Waiting on Track B Vina (`dock_track_b_production_v1.py`); will auto-start C5 when it clears.

## Order (locked)
1. Wait for current Track B to finish (no CPU overlap)
2. Task1 gate: benzbromarone@9DKA × 42/43/44; CNNscore Top-1 RMSD ≤ 2.0 Å vs R75
3. Gate fail → STOP (no Task2/3)
4. Gate pass → remaining W1 (32 new; reuse lesinurad@9DKB×3 + benzbromarone@9DKB seed42)
5. Task2 W4 @ 7ALV (≥146); decoys already locked
6. Task3 W2 IFP rescoring (0 new dock)
7. Timeouts: skip molecule, continue

## Prepared
- `9DKA/9DKC/9DK9_receptor.pdbqt` via `prepare_receptor_vina.py`
- `TD-3.pdbqt` carboxylate (Dimorphite-DL → Meeko)
- Crystal refs bond-ordered under `01_ligand_prep/w1_crystal_refs/`
- W4 decoys: `02_nlrp3_panel/w4_decoys_locked.csv` (n=40, seed 0xC5DEC0)

## Engine note
Worklist asks GNINA **1.3.1**; local binary is **1.3.2**. Settings otherwise locked: exh=32, modes=9, rescore, seeds 42/43/44, `--no_gpu`.
