# Assay-context human review SOP (local)

Machine extraction is already in:

- `data/jcim_novelty_v0/tables/assay_context_priority_ligands_v1.csv`
- `data/jcim_novelty_v0/tables/assay_context_audit.csv`

Do **not** start by reading all 352 ligands. Review only the priority set, in this order, including ligands that hurt the main claim.

## Order

1. EGFR/HER2: every dual, A-only, and B-only ligand.
2. PIK3CA/mTOR: the four neither ligands (single document `CHEMBL1240340`).
3. Ligands flagged `high_auroc_influence` (leave-one-out effect on directional AUROC).
4. Ligands flagged `mixed_endpoint`, `biochem_and_functional`, or `non_human_assay_organism`.
5. Ligands flagged `top_document_series`.

## Fields to fill

For each activity row in `assay_context_audit.csv`:

| column | allowed values |
|--------|----------------|
| `protein_construct` | short text or `unknown` |
| `wildtype_or_mutant` | `WT` / `mutant:<name>` / `unknown` |
| `human_include_exclude` | `include` / `exclude` / `uncertain` |
| `human_reviewed_label` | `dual` / `A_only` / `B_only` / `neither` / `drop` |
| `human_rationale` | one sentence |
| `incomparable_record` | `0` / `1` |

Also complete the ligand-level columns in `assay_context_priority_ligands_v1.csv`.

## Decision rules

- Keep biochemical vs cellular, IC50 vs Ki/Kd/EC50, and WT vs mutant as **context**, not silent relabeling.
- Exclude only when the record is clearly not the intended human protein/assay or is internally contradictory.
- Do not manufacture a single unified truth by majority vote.
- If any frozen class changes, recompute Table 2 and document the sensitivity. Do not keep the original number because it was more favorable.

## Local PDF drop (7 uncertain ligands)

Place publisher PDFs in
`data/jcim_novelty_v0/literature_sources/pdfs/` using the filenames in
that folder's README. Source-reading notes:
`data/jcim_novelty_v0/analysis/ASSAY_CONTEXT_SOURCE_READING_V1.md`.
PDFs are gitignored. After a PDF lands, the next pass fills construct /
mutation on the matching audit rows. Do not recompute Table 2 unless a
frozen class actually changes.

## What this cloud pass already did

Assay type, organism, relation, endpoint, document, assay ID, mixed-endpoint flags, and AUROC influence were extracted from the 2026-08-26 high-confidence view. Protein construct and mutation usually require the source paper or ChEMBL assay description and are left blank on purpose. A 2026-08-27 source-reading pass started from public full text (US9181263, Hong 2024 PMC) and ChEMBL/BindingDB cross-checks; labels were not changed.
