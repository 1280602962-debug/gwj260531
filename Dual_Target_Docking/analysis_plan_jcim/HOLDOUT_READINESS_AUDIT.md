# HOLDOUT_READINESS_AUDIT

Membership is frozen. A complete labeled replication set is **not** present on this release branch.

This file does not dock, compute AUROC, drop ligands, or call the set external validation.

## Current formal file

`00_protocol/holdout_membership_freeze.csv`

Columns: `pair`, `holdout_id`, `parent_inchikey`, `holdout_eligible`, `overlap_main_ids`, `status`, `exclusion_reason`, `activity_eligible`.

Missing on this file: class, SMILES, pA/pB, Murcko, document/source, scores.

## Membership counts

| pair | n | holdout_eligible=1 | activity_eligible=1 |
|---|---:|---:|---:|
| JAK1/JAK2 | 58 | 58 | 58 |
| JAK1/TYK2 | 60 | 60 | 60 |
| PIK3CA/mTOR | 60 | 60 | 60 |
| AChE/BChE | 60 | 59 | 60 |
| F2/F10 | 60 | 60 | 60 |
| PPARG/PPARA | 60 | 60 | 60 |
| PPARA/PPARD | 60 | 60 | 60 |
| EGFR/HER2 | 0 | 0 | 0 |
| **total** | **418** | **417** | **418** |

Status: `INDEPENDENT` = 417; `RESOLVED_EXCLUDED` = 1.

The excluded row:

- pair AChE/BChE
- `HOAB_026`
- parent `ADEBPBSSDYVVLD-UHFFFAOYSA-N`
- `holdout_eligible=0`
- `overlap_main_ids=AB_043`
- `exclusion_reason=parent_overlap_with_main_panel`

`417 eligible` is a membership count, not 417 complete evaluation samples.

`pair_ligand_mapping.csv` has `analysis_set=main` for all 808 rows and empty `holdout_eligible`. Holdout members are not in the formal main mapping.

## Answers to the required questions

1. **Continuous pA/pB on this branch?** No.
2. **Frozen activity file on this branch?** `ligand_activity_aggregate_v1.csv` is absent and listed under `forbidden_activity_regeneration`. It is a main-set aggregate on `scientific-freeze-c7cc`, not a holdout table.
3. **Frozen before V4.2 PRIMARY?** Membership file is in the current formal protocol. Historical holdout panels exist on `a0ef2726` (2026-09-20), before `PRIMARY_DIRECTIONAL_METRICS.csv`.
4. **Can θ=6.0 be applied without ambiguity on this branch?** No. There are no continuous values here. Historical `holdout_panels_all_v1.csv` on `a0ef2726` has `pchembl_A`, `pchembl_B`, and `class`, but `label_rule=strict_6.5_5.5`. Those classes are not automatically the θ=6.0 classes.
5. **Parent overlap with main?** Encoded for one excluded AChE row (`AB_043`). Other `overlap_main_ids` fields are empty. This is a recorded membership rule, not a new Tanimoto filter.
6. **Murcko on this branch?** No.
7. **Can main-vs-holdout scaffold overlap be computed later?** Only after SMILES or scaffolds are restored under a later permission. Not possible from the membership file alone.
8. **Can nearest-main ECFP4 be computed later?** Same: needs holdout SMILES, which are not on this branch.
9. **Document/source on the membership file?** No source/document columns. Historical panel table on `a0ef2726` has `source=chembl37_unused_pool_post_panel_freeze` for the inspected header row. That is provenance, not a restored formal table.
10. **Per-pair membership:** table above. EGFR/HER2 has none.
11. **holdout_eligible=0:** one row, `HOAB_026`.
12. **Why excluded:** `parent_overlap_with_main_panel` / `AB_043`.

## Historical objects (git only; not formal inputs)

`git show a0ef2726:Dual_Target_Docking/data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panels_all_v1.csv` header includes `class`, `canonical_smiles`, `pchembl_A`, `pchembl_B`, `murcko_scaffold`, `label_rule`, `source`.

`data/jcim_holdout_v0/` on that commit contains SDF/PDBQT for some holdout IDs.

Those paths are **not** on this release branch and are **not** current score or label authority. Historical holdout Vina tables on that branch must not be treated as formal V4.2 M0.

## Readiness

`data_readiness = NEED_FROZEN_SOURCE_PERMISSION`

`evidence_status = MEMBERSHIP_FROZEN_LABELS_SMILES_ACTIVITY_ABSENT_ON_RELEASE`

A later phase may consider M0 seed42 reserved-sample checking only after a human permission to restore a frozen labeled holdout table, relabel at θ=6.0 from continuous values, and confirm independence. Role remains: within-source reserved-sample check. Not external validation.

Worth continuing? Only after that permission and a join audit. Not worth docking now.
