# Receptor identity audit — P0 finding

**Script:** `scripts/receptor_identity_audit_v1.py`  
**Table:** `tables/receptor_identity_audit_v1.csv`  
**Method:** for every receptor actually used for docking, resolve the RCSB entry title, polymer-entity SIFTS UniProt accession, source organism, and cognate ligand, then compare the accession against the UniProt the DualFourClass labels were harvested from.

## Result

| Pair | End | Expected | PDB | Resolved accession | Organism | Verdict |
|------|-----|----------|-----|--------------------|----------|---------|
| PIK3CA/mTOR | A | PIK3CA P42336 | 4L23 | P42336 | Homo sapiens | OK |
| PIK3CA/mTOR | B | MTOR P42345 | 4JT6 | P42345 | Homo sapiens | OK |
| AChE/BChE | A | ACHE P22303 | 4EY7 | P22303 | Homo sapiens | OK |
| AChE/BChE | B | BCHE P06276 | 4BDS | P06276 | Homo sapiens | OK |
| EGFR/HER2 | A | EGFR P00533 | 3POZ | P00533 | Homo sapiens | OK |
| EGFR/HER2 | B | ERBB2 P04626 | 3RCD | P04626 | Homo sapiens | OK |
| PIK3CA/PIK3CB | A | PIK3CA P42336 | 4L23 | P42336 | Homo sapiens | OK |
| **PIK3CA/PIK3CB** | **B** | **PIK3CB P42338** | **2WXF** | **O35904** | **Mus musculus** | **WRONG PROTEIN AND SPECIES** |
| PIK3CA/mTOR alt (Table S30) | A/B | PIK3CA / MTOR | 4JPS, 5DXT, 4JSX | P42336 / P42345 | Homo sapiens | OK |
| MCL1/Bcl-xL stress test | A/B | Q07820 / Q07817 | 3WIY / 3WIZ | Q07820 / Q07817 | Homo sapiens | OK |

## The 2WXF problem

- RCSB entry title: *“The crystal structure of the **murine** class IA PI 3-kinase **p110delta** in complex with PIK-39.”*
- Polymer entity 1 → SIFTS **O35904**, gene **Pik3cd**, *Mus musculus*, entity sequence coverage 0.998.
- Cognate ligand **039** = PIK-39, a p110δ-preferring inhibitor.
- Human PIK3CB is **P42338**, which has **zero** PDB cross-references (see `UNIVERSE_STRUCTURE_FEASIBILITY_V1.md`).
- Local files confirm the receptor in use is that entry: `receptors/PIK3CB_protein.pdb` and `receptors/2WXF_protein.pdb` are byte-identical (md5 `2e9002f5de98e130e59cfcb06e1bd123`), as are the two PDBQT files.

So the pair reported as **PIK3CA/PIK3CB** docked ligands whose four-state labels come from **human PIK3CB (P42338)** activity into a **mouse p110δ (PIK3CD)** pocket. The A end is correct.

Methods and SI currently list this receptor as “PIK3CA/PIK3CB, 4L23 / 2WXF (X6K / 039)” with no isoform or species qualifier. `MANIFEST.md` records it as “PIK3CB | 2WXF | 039 | RMSD 0.405 @E=8/16”, i.e. it passed cognate QC — cognate redocking cannot detect a wrong-protein receptor, because the cognate ligand belongs to that same wrong protein.

Rejected candidates 2Y3A and 4BFR are mouse **Pik3cb** (Q8BTI9) — right protein, wrong species. 3T8M was already excluded as a chimera; its entity maps to human PIK3CG (P48736). The one structure that was kept is the only one that is both the wrong isoform and the wrong species.

## What this affects

1. **Table 2, PIK3CA/PIK3CB row** (`summary_min` 0.500 [0.350, 0.650]) — the B pocket is not PIK3CB.
2. **Table S30 / B5 receptor swap** — “same PIK3CA crystals with 2WXF frozen raised worst-arm AUROC from 0.500 to 0.691 / 0.685”. The frozen B end is mouse p110δ throughout.
3. Any sentence describing this pair as an **isoform control** for PI3Kα/β. Structurally it is an α-versus-δ, human-versus-mouse comparison.
4. Panel labels remain valid as *ChEMBL* labels; only the receptor is wrong.

It does **not** affect PIK3CA/mTOR, AChE/BChE, EGFR/HER2, the alternate-crystal set, or the MCL1/Bcl-xL stress test. All of those resolve to the expected human accession.

## Why it cannot simply be re-docked

Human PIK3CB has no crystal structure in the PDB. The only options are mouse Pik3cb (2Y3A 3.30 Å, 4BFR 2.80 Å — both previously failed cognate QC or PDBQT parsing) or a homology model. Neither is a drop-in replacement, and both require a written species-surrogate exception.

## Recommended handling (not executed here)

- Withdraw PIK3CA/PIK3CB from the primary directional table; report it in the SI as a **documented receptor-identity failure**, alongside the existing 3T8M chimera lesson.
- Correct the Table S30 sentence, or drop that swap arm.
- Keep the three verified pairs as the docked set.
- Add this audit to the protocol so identity is checked **before** cognate QC on any future pair.

Table 2 and the manuscript were **not** edited by this audit.
