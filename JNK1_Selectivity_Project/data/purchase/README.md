# JNK1 MD Post-QC Purchase List

Wet-lab purchase shortlist after **48 Desmond MD pose QC jobs** (16 compounds × 3 PDBs: 3ELJ / 3E7O / 3TTI).

## Files

| File | Description |
|------|-------------|
| `purchase_after_md.csv` | **10 compounds** for same-batch JNK1/JNK2/JNK3 enzymatic IC50 |

## Selection summary

| Group | In list | Rationale |
|-------|---------|-----------|
| G3_control | 4 | Assay calibration rulers (always purchase) |
| G1 | 3 | MD `pass_md_overall` grade A (690, 2232, 2157) |
| G2 | 3 | 2231 (JNK1-only pass) + 1280/4795 (off-target pose backups) |

## SMILES validation

All SMILES were validated with RDKit (2026-07-03):

- **G3 controls**: canonical isomeric SMILES match `data/benchmarks/literature_benchmarks.csv`
- **G1/G2 hits**: achiral; no undefined stereocenters
- **CC-90001 / CC-930**: stereochemistry preserved (`@` markers required for correct isomer)

## Ordering notes

### G3 benchmarks (commercial)

| Compound | Suggested source |
|----------|------------------|
| SP600125 | Selleck / MCE / Cayman (CAS **129-56-6**) |
| CC-90001 | MedChemExpress or custom synthesis |
| CC-930 (Tanzisertib) | MedChemExpress / Cayman |
| E1 | Literature supplier / custom synthesis |

### G1/G2 screening hits (Enamine library)

`ligand_id` (690, 2232, …) is the **internal ID from the original Enamine VSW library**. Order by **SMILES** (column `smiles`) if catalog ID lookup fails.

## Assay plan

Run **JNK1 + JNK2 + JNK3 biochemical IC50 in the same batch** for all 10 compounds.

MD pass indicates **pose credibility in the binding site**, not confirmed activity or isoform selectivity.

## Related results

MD QC source report: user workspace `MD_QC_report_cf26.md` (not in repo).
