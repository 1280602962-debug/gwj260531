# C5 docking handoff status

Updated 2026-09-05 after W1 gate audit. Authoritative summary:
`docs/C5_RANKING_AND_NEXT_DOCKS.md`.

## Live

W1 gate (benzbromarone@9DKA × 42/43/44) **finished and failed**
(CNNscore Top-1 `pose_rmsd` ≈ 3.59 Å). No gnina process is required
to wait on that gate anymore.

## Order (locked after audit)

1. ~~Wait for Track B~~ — not a scientific dependency for C5 W4/W2.
2. ~~Task1 gate~~ — **done, fail** (`search_ok_selection_fail`). Do not rerun.
3. **Task2 / W4 @ 7ALV — RUN NOW (146 new jobs).** Decoys locked.
4. **Task3 / W2 IFP rescoring — RUN (0 new docks).** Anchor on deposited
   crystal coordinates (R75/A1AIL/A1A45), not CNNscore Top-1 poses.
5. W1 remaining 29 cells: optional SI only; **not** an unlock.
6. Timeouts: skip molecule, continue.

## Must-run W4 jobs (146)

| subset | n | seeds | jobs |
|---|---:|---|---:|
| 8 non-cocrystal positives | 8 | 42/43/44 | 24 |
| locked decoys `w4_decoys_locked.csv` | 40 | 42/43/44 | 120 |
| REP_07837 @ 7ALV | 1 | 43, 44 | 2 |

Receptor: 7ALV. Box: `[16.756, 35.449, 125.714]` 20³.
Engine: exh=32, modes=9, `cnn_scoring=rescore`. Do not reuse exh=8 panel SDFs.

## Prepared

- `9DKA/9DKC/9DK9_receptor.pdbqt` via `prepare_receptor_vina.py`
- `TD-3.pdbqt` carboxylate (Dimorphite-DL → Meeko)
- Crystal refs bond-ordered under `01_ligand_prep/w1_crystal_refs/`
- W4 decoys: `02_nlrp3_panel/w4_decoys_locked.csv` (n=40, seed 0xC5DEC0)

## Engine note

Worklist asks GNINA **1.3.1**; local binary may be **1.3.2**. Settings otherwise locked: exh=32, modes=9, rescore, seeds 42/43/44.

## Do not

- Reopen URAT1 docking-score ranking
- Declare the 2.0 Å gate passed via GetBestRMS
- Start W5/W6 before shortlist freeze
- Treat remaining 29 W1 jobs as required

## W4 launch 2026-09-05 03:07 UTC

- Authorized despite W1 gate fail (`search_ok_selection_fail`).
- Positives+decoys re-prepped with `prepare_ligands_c1.py` (Dimorphite→Meeko); old panel Meeko-only PDBQT not used.
- Jobs: 146 (8×3 + 40×3 + REP_07837 seeds 43/44).
- Runner: `scripts/run_c5_w4_nlrp3_panel.py --cpu 3 --workers 2`.
- Status CSV: `data/campaigns/c5/02_nlrp3_panel/w4_job_status.csv`

## W2 complete 2026-09-05 08:18 UTC

- Crystal-anchored IFP gate on Phase I 9DKB SDFs (228 vs 64), 0 new docks.
- Key map: 11/12 (Q437 = LEU in 9DKB, unmatched).
- Primary IFP (CNNscore Top-1): **gate_pass=false** → fallback A1∩A2.
- Outputs: `data/campaigns/c5/02_urat1_ifp/`

## W2 final 2026-09-05 08:21 UTC

- IFP primary (CNNscore Top-1 + overlap/IFP/key/clash; Arg≤7.7027 from A1AIL): **gate_pass=true**
  OR=3.15, CI95=[1.72, 7.08], Fisher p=6.7e-4 (228 vs 64).
- Arg max locked to pre-registered 7.7027 (not loosened by Kabsch R75/A1A45 O–Arg).
- Q437 unmatched (LEU in 9DKB). Outputs under `data/campaigns/c5/02_urat1_ifp/`.

