# PREPARATION INTEGRITY AUDIT

RUN_ID = UNIFORM_RERUN_V4_2_20260921

- RUN_ID: **UNIFORM_RERUN_V4_2_20260921**
- V4_1_Reduce2_branch: **ABANDONED_NOT_OVERWRITTEN**
- FINAL_INPUT_FREEZE_PASS: **YES**
- ligand_assets_reused: **785_from_V4_1**
- Reduce2_used: **NO**
- 5U3Q B used: **YES**
- 4JT6 auth A label C X6K A:2601: **YES**
- 4EY7 box E20 A:604: **YES**
- PTR status: **DOCUMENTED_UNIFORM_PTR_TO_TYR_FALLBACK**
- AB_046 fold_id: **1**
- AB_040 alias fold resplit: **NO**
- prepared_receptors: **14**
- held_receptors: **[]**
- scientific_check_issues: **[]**
- RECEPTOR_PREPARATION_INTEGRITY_PASS: **YES**
- LIGAND_PREPARATION_INTEGRITY_PASS: **YES**
- PHASE_5_REDOCKING_GO: **YES**
- new_PDB_searched: **NO**
- Cognate redocking jobs: **NOT STARTED**
- Production Vina: **NOT STARTED**
- STOP_AFTER_PHASE_4: **YES**

V4.1 Reduce2 branch abandoned. V4.1 receptor files are historical and were not overwritten.
Meeko overrides limited to crystal-verifiable CYX (SG–SG), blunt-end polymer breaks,
and observed-fragment templates for residues whose PDBFixer-added atoms created new covalent-range contacts or valence failure.
See 13_qa/HOLD_STAGE_COMPARISON.md. Forbidden flags were not used. No new PDB was searched.

## Receptor status
- 3POZ: PREPARED (v4_2_meeko_ok)
- 3RCD: PREPARED (v4_2_meeko_ok)
- 6N7A: PREPARED (v4_2_meeko_ok)
- 8BXH: PREPARED (v4_2_meeko_ok)
- 3LXP: PREPARED (v4_2_meeko_ok)
- 4L23: PREPARED (v4_2_meeko_ok)
- 4JT6: PREPARED (v4_2_meeko_ok)
- 4EY7: PREPARED (v4_2_meeko_ok)
- 4BDS: PREPARED (v4_2_meeko_ok)
- 4UDW: PREPARED (v4_2_meeko_ok)
- 2JKH: PREPARED (v4_2_meeko_ok)
- 9V8H: PREPARED (v4_2_meeko_ok)
- 6LXA: PREPARED (v4_2_meeko_ok)
- 5U3Q: PREPARED (v4_2_meeko_ok)

## Scientific checks
- RETAIN components present; REMOVE components absent
- 6N7A A:1041 / 8BXH A:839 / 4UDW I:565 contain only crystal-observed heavy atoms
- All 14 PDBQT files are non-empty
- No forbidden Meeko flags
- V4.1 historical receptors remain on disk
