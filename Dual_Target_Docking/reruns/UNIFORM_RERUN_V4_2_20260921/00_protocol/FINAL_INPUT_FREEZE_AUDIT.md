# FINAL INPUT FREEZE AUDIT
RUN_ID: UNIFORM_RERUN_V4_1_20260921
PROTOCOL_VERSION: 4.1
Written: 2026-09-20T17:31:10Z
Authority: `scripts/analysis/build_current_score_master.py` + resolved unique source files listed in `analysis_freeze.yaml`.

## Formal freeze record

INPUT_IDENTITY_AUDIT_PASS = YES
LIGAND_PANEL_AUDIT_PASS = YES
FINAL_INPUT_FREEZE_PASS = YES

UNRESOLVED_P0 = NONE
UNRESOLVED_P1 = NONE

HOLD_REVIEW_REQUIRED: NONE

fold_id: AB_046 = 1 from unique file `results/canonical/model_fold_assignments.csv`. AB_040 remains DUPLICATE_PARENT_ALIAS and is not an independent fold row. No scaffold resplit.

Input source ambiguity: NONE

## Frozen 8 pairs (A/B order unchanged)

1. EGFR / HER2
2. JAK1 / JAK2
3. JAK1 / TYK2
4. PIK3CA / mTOR
5. AChE / BChE
6. F2 / F10
7. PPARG / PPARA
8. PPARA / PPARD

## Frozen 14 primary receptors

EGFR 3POZ; HER2 3RCD; JAK1 6N7A; JAK2 8BXH; TYK2 3LXP; PIK3CA 4L23; mTOR 4JT6; AChE 4EY7; BChE 4BDS; F2 4UDW; F10 2JKH; PPARG 9V8H; PPARA 6LXA; PPARD 5U3Q.

5U3Q primary copy = B (7UJ auth B / label N / 501, occupancy 1.0). Copy A = LEGACY_ONLY.
4JT6 auth A = label C; box ligand X6K A:2601.
4EY7 receptor context = AB dimer; box-defining E20 A:604; E20 B:605 equivalent non-primary; no A/B bake-off.

## Components

4L23 p85 REMOVE; 4JT6 mLST8 REMOVE; 9V8H peptide RETAIN; 4UDW H/I/L RETAIN; 2JKH A/L RETAIN.

## Mutations (POCKET_REMOTE, keep as deposited)

3LXP five mutations; 4BDS N→Q; 2JKH R372E.

## PTR

Input-stage experimental PTR is recorded. `replaceNonstandardResidues()` remains forbidden.
Compatibility test: all three of 6N7A, 8BXH, 3LXP lose PTR at Meeko (PDBFixer heavy-atom repair keeps PTR phosphate). Frozen V4.1 production status:

PTR_STATUS = DOCUMENTED_UNIFORM_PTR_TO_TYR_FALLBACK

The same PTR→TYR transformation is applied to 6N7A, 8BXH, and 3LXP after PDBFixer and before Reduce2. Decision did not use docking score, RMSD, or AUROC.

## Panel

Historical records: 808
Independent pair-level parent units: 807
Global unique ligand entities (desalted parent InChIKey): 785
AB_040 = DUPLICATE_PARENT_ALIAS of AB_046; no independent 3D/PDBQT/docking.
HOAB_026: activity_eligible factual = 1 (holdout panel pA=7.52 pB=5.44); holdout_eligible = 0; overlap AB_043.
EH120_059 and AB_087: activity_eligible = 0; missing ≠ inactive.

## Track-B source concentration

NOT_ASSESSABLE_FROM_FROZEN_INPUT

## P2 limitations

NON_BLOCKING_SCIENTIFIC_LIMITATION. No panel resampling.

## Not started in Phase 0

PDBFixer, Reduce2, ETKDGv3 production, Meeko production, redocking, Vina, GNINA, RTMScore.
