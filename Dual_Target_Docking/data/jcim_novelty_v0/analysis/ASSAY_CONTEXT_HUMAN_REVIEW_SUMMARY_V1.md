# Assay-context human review summary v1

Reviewer: `local_agent_metadata_pass` on 2026-08-26.
ChEMBL assay free-text API returned HTTP 500; `protein_construct` and
`wildtype_or_mutant` set to `unknown` for all reviewed rows.

- Priority ligands: 186
- include / uncertain / exclude: 179 / 7 / 0
- Ligands with exclude or class flip vs frozen: 0

## Label sensitivity

No frozen DualFourClass labels were changed in this pass.

Exclude/uncertain ligands (if any):

- EGFR/HER2 EH120_045 (dual): uncertain → dual; At least one activity row uncertain; frozen class retained.
- PIK3CA/mTOR PM48_04 (dual): uncertain → dual; At least one activity row uncertain; frozen class retained.
- PIK3CA/mTOR PM48_05 (dual): uncertain → dual; At least one activity row uncertain; frozen class retained.
- AChE/BChE AB_089 (neither): uncertain → neither; At least one activity row uncertain; frozen class retained.
- AChE/BChE AB_091 (neither): uncertain → neither; At least one activity row uncertain; frozen class retained.
- AChE/BChE AB_094 (neither): uncertain → neither; At least one activity row uncertain; frozen class retained.
- PIK3CA/mTOR PM48_22 (A_only): uncertain → A_only; At least one activity row uncertain; frozen class retained.

Recompute Table 2 only if any `human_reviewed_class` differs from `frozen_class` or excludes remove directional arms.
