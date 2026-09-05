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
