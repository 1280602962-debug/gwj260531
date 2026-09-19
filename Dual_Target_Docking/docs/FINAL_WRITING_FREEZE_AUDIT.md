# Final writing-freeze audit

Date: 2026-09-19  
Scope: PR branch `cursor/scientific-freeze-c7cc` after CASE 2 EGFR uniform promotion.  
Zero-dock rebuild: `/tmp/dual_target_freeze_rebuild` → `results/canonical`.  
Verification: `result: PASS`. `check_current_chain.py`: PASS. Path audit unresolved = 0.

## P0

None.

| Item | Status | Evidence |
|------|--------|----------|
| Fresh zero-dock rebuild | CLOSED | `rebuild_freeze.py` exit 0; `verify_freeze_rebuild.py` `result: PASS` |
| Canonical / submission-pack sync | CLOSED | pack 31/31 canonical CSV identical |
| Broken current paths | CLOSED | `docs/CURRENT_PATH_AUDIT.md` unresolved = 0 |
| Unexplained protocol contradictions | CLOSED | PRIMARY / SENSITIVITY / INDEPENDENT GNINA / SAME-POSE documented; Track-B `do_not_now` not rewritten |
| EGFR uniform prep provenance | CLOSED | **confirmed**: `prep_uniform_ligands.py`, 110/110 PDBQT, Vina five-seed CSV, GNINA CSV |
| AChE five-seed membership | CLOSED | `same_protocol_as_primary=1`, `same_membership_as_primary=0`; fixed-membership n=88; `qualitative_change=no` |
| Receptor provenance boundary | CLOSED | 14 unique PDBs in `data/provenance/receptor_input_registry.csv`; prep not claimed as uniform Meeko/PPW |
| EGFR is only current production truth | CLOSED | `EGFR_SCORE_SOURCE` / master `score_source` = `data/egfr_her2_uniform_rdkit_v1/tables/scores_vina_mode1_fiveseed.csv`. Historical ablation is not Table 2. |

P0 count = **0**

## P1

| Item | Status | Note |
|------|--------|------|
| Ligand PDBQT for PIK3CA / AChE / Track-B production | CLOSED as documented | analysis reproducible from deposited scores; exact original PDBQT **not_recoverable** |
| Receptor prep method | CLOSED as documented | PDBQT **confirmed**; prep recipe **not_recoverable** except 4EY7 `metadata_only` |
| EGFR GNINA incomplete cases | CLOSED as documented | 33 timeout_skipped; EH120_109 fail (`CG0`); n=89 |
| EGFR Vina EH40_31 | CLOSED as documented | timeout all seeds; not in primary (B-only dropped; Table 2 n=28/37/31) |
| Historical audit memos | CLOSED as historical | `PR39_*`, `STALE_*`, `FINAL_REVIEWER_*` are pre-promotion snapshots, not current number authority |

## P2

| Item | Status |
|------|--------|
| Nested descriptor scaffold-CV (option B) | CLOSED (`descriptor_nested_scaffold_cv.csv`; ECFP pipeline unchanged) |
| Label kinds | CLOSED (`chembl_assay_adjudicated` / `chembl37_dump_panel` / `panel_pchembl_no_audit_rows`) |
| Analysis vs docking environment split | CLOSED (analysis Python 3.12.3 + pinned analysis packages; docking versions only where evidenced) |
| Figure 5C AChE membership caption | CLOSED |
| Figure 4 EGFR matched−mismatched now excludes 0 | CLOSED (new numbers win; 0.107 [0.006, 0.220]) |

## CLOSED

All writing-freeze closure items 1–10 for this round.

Current EGFR Table 2: `summary_min` **0.3341 [0.1970, 0.4713]**, n=28/37/31. Do not restore 0.3237.

## Gates

- P0 = 0
- broken current paths = 0
- unexplained protocol contradictions = 0
- EGFR uniform prep provenance = confirmed
- AChE five-seed membership documented + fixed-member QA passed
- receptor provenance boundary documented
- canonical / submission-pack synchronized
- fresh rebuild PASS

## Verdict

READY_FOR_WRITING_FREEZE
