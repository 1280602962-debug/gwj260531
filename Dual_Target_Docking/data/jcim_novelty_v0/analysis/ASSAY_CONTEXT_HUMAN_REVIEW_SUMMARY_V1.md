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

Pass for the 7 uncertain ligands is in
`ASSAY_CONTEXT_SOURCE_READING_V1.md`, with a row-level verdict table in
`tables/uncertain_ligand_source_verdicts_v1.csv`. Labels were **not**
changed. Table 2 was **not** recomputed.

What was actually read: US 9,181,263 B2 local PDF; Hong 2024 PMC HTML;
BindingDB article dumps for Folkes / Elmenier / Bass / Ma / Sang / Yang /
Cheng / Shi (DOI downloads were not publisher PDFs).

Confirmed in this pass:

- EH120_045: patent Table 1 EGFR 0.5 nM; Ma KINOMEscan HER2 Kd **0.75 nM**
  (the 76 nM page was the wrong ligand). Still dual.
- PM48_04: Hong 1.4 nM is cited, not measured. Yang mTOR 0.45 nM is
  truncated FLAG 1362–end. Still dual if 8.85 is dropped.
- PM48_05: Elmenier mTOR 3 nM is a **column shift** of Folkes PI3Kα.
  Bass PI3Kα 0.1 nM is a review transcription vs Folkes 3 nM. Folkes
  original 3 nM / Ki 580 nM remains dual at θ=6.0.
- AB_089/091/094: frozen BChE matches **human** serum BuChE, not horse.
  Still neither.
- PM48_22: Shi PI3Kα 1.70 nM confirmed; Cheng 5.83 remains a cellular
  p70S6K ELISA. Still A_only.
