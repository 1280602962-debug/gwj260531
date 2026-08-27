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

## Source-document reading (2026-08-27)

Public-text pass for the 7 uncertain ligands is in
`ASSAY_CONTEXT_SOURCE_READING_V1.md`. Drop PDFs in
`literature_sources/pdfs/`. Labels were **not** changed in this pass.

Highest-risk finding: PM48_05 frozen maxima (PIK3CA 10.00, mTOR 8.52)
come from reviews and conflict with Folkes 2008 (3 nM / ~580 nM,
already in ChEMBL as CHEMBL1140078). At θ=6.0 the Folkes pair is still
dual (8.52 / 6.24). A dual→A_only flip needs an extra rule, not the
frozen max-pChEMBL rule alone. EH120_045 HER2 9.12 (0.76 nM) conflicts
with BindingDB 76 nM on the same assay ID. AB_* stay neither. PM48_22
stays A_only (cellular mTOR surrogate). PM48_04 1.4 nM PI3Kα is cited
in Hong 2024, not measured there.
