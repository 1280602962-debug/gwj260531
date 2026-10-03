# ALT_RECEPTOR_READINESS

Stage 1 completeness only. No AUROC.

Authority: `00_protocol/ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml`

Formal subpanel: PIK3CA/mTOR, AChE/BChE, F2/F10, PPARA/PPARD.

## Existing masters (current release branch)

| File | Rows | Role |
|---|---:|---|
| `13_qa/alt_redocking_master.csv` | 40 | cognate redock QC, 8 PDB × 5 seeds |
| `13_qa/alt_production_seed42_master.csv` | 734 | production Vina seed42 jobs |
| `13_qa/alt_score_master_independent.csv` | 734 | alt score plus kept-primary scores |
| `13_qa/alt_production_completeness_audit.json` | — | expected 734, accounted 734, SUCCESS 706, TIMEOUT 28, PERMANENT_FAIL 0, MISSING_JOB 0 |

`AUROC_allowed: true` in that JSON is a historical completeness flag. This Stage 1 audit does **not** compute AUROC.

## Job trees

- `07_alt_cognate_redocking/jobs`: 40 directories; 5 per frozen alt PDB.
- `08_alt_vina_seed42/jobs`: 734 directories, matching the production master.

## Per-substitution production (seed 42)

| alt_id | pair | side | design | alt PDB | expected | SUCCESS | TIMEOUT | PERMANENT_FAIL | redock QC (5 seeds) |
|---|---|---|---|---|---:|---:|---:|---:|---|
| ALT_4L2Y | PIK3CA/mTOR | A | A2/B1 | 4L2Y | 48 | 48 | 0 | 0 | 5/5 status OK |
| ALT_4JT5 | PIK3CA/mTOR | B | A1/B2 | 4JT5 | 48 | 48 | 0 | 0 | 5/5 status OK |
| ALT_4EY6 | AChE/BChE | A | A2/B1 | 4EY6 | 99 | 91 | 8 | 0 | 5/5 status OK |
| ALT_1P0M | AChE/BChE | B | A1/B2 | 1P0M | 99 | 90 | 9 | 0 | 5/5 status OK |
| ALT_3SHC | F2/F10 | A | A2/B1 | 3SHC | 110 | 104 | 6 | 0 | 5/5 status OK |
| ALT_2Y5F | F2/F10 | B | A1/B2 | 2Y5F | 110 | 105 | 5 | 0 | 5/5 status OK |
| ALT_6KAX | PPARA/PPARD | A | A2/B1 | 6KAX | 110 | 110 | 0 | 0 | 5/5 status OK |
| ALT_5U46 | PPARA/PPARD | B | A1/B2 | 5U46 | 110 | 110 | 0 | 0 | 5/5 status OK |

Expected job counts match the freeze `ligand_reuse.pair_independent_prepare3d` numbers (48 / 99 / 110 / 110).

TIMEOUT is missing, not imputed. `alt_score_master_independent.csv` records `TIMEOUT_treated_as=missing_not_imputed`.

Redocking master parameters observed: scoring=vina, exhaustiveness=16, num_modes=9, energy_range=6, cpu=1, timeout_s=1586, `receptor_modified_from_rmsd=0` and related flags 0. RMSD columns exist for QC description later; they must not be used to change receptors.

## Readiness

`data_readiness = READY_EXISTING_RAW`

A later M0 paired AUROC can use `alt_score_master_independent.csv` with complete-case ligands (finite alt score and the required primary scores). Substitutions with TIMEOUTs remain in the design; missing scores stay missing.

No new PDB, no best-alternative selection, no AUROC in this phase.

## Not in this module

EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PPARG/PPARA have no formal alternative substitutions.
