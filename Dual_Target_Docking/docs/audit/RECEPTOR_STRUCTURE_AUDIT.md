# RECEPTOR STRUCTURE AUDIT

Date: 2026-09-19  
Machine tables: `data/provenance/receptor_input_registry.csv`, `data/provenance/docking_box_registry.csv`  
Rule: deposited PDBQT and box JSON are evidence. Old markdown is not.

## 14 slots

| Protein | PDB | Pair pocket | Unique file | Construct (from PDBQT / yaml) | Cognate | Peptide / extra | PDBQT exists | Prep method |
|---------|-----|-------------|-------------|-------------------------------|---------|-----------------|--------------|-------------|
| EGFR | 3POZ | EGFR/HER2 A | yes | chain A kinase | 03P | none in PDBQT | yes | NOT_RECOVERABLE |
| HER2 | 3RCD | EGFR/HER2 B | yes | chain A kinase | 03P | none in PDBQT | yes | NOT_RECOVERABLE |
| PIK3CA | 4L23 | PIK3CA/mTOR A | yes | chain A; protein PDB deposited | X6K | none in PDBQT | yes | NOT_RECOVERABLE |
| mTOR | 4JT6 | PIK3CA/mTOR B | yes | chain A; protein PDB deposited | X6K | none in PDBQT | yes | NOT_RECOVERABLE |
| AChE | 4EY7 | AChE/BChE A | yes | dimer AB (`ACHE_receptor.pdbqt`) | E20 | dimer AB | yes | meeko recipe in recovered script = METADATA_ONLY |
| BChE | 4BDS | AChE/BChE B | yes | chain A | THA | none in PDBQT | yes | NOT_RECOVERABLE |
| F2 | 4UDW | F2/F10 A | yes | thrombin H/I/L | N6L | HIL construct | yes | NOT_RECOVERABLE |
| F10 | 2JKH | F2/F10 B | yes | chains A/L | BI7 | none in PDBQT | yes | NOT_RECOVERABLE |
| JAK1 | 6N7A | JAK1/JAK2 A and JAK1/TYK2 A | shared | AB; yaml JH1 ATP | KEV | none in PDBQT | yes | NOT_RECOVERABLE |
| JAK2 | 8BXH | JAK1/JAK2 B | yes | A; yaml JH1 ATP | C87 | none in PDBQT | yes | NOT_RECOVERABLE |
| TYK2 | 3LXP | JAK1/TYK2 B | yes | A; yaml JH1 ATP **not JH2** | IZA | none in PDBQT | yes | NOT_RECOVERABLE |
| PPARG | 9V8H | PPARG/PPARA A | yes | AB; yaml LBD ternary keep PG08-NL | BRL | chain B peptide **present in PDBQT** | yes | NOT_RECOVERABLE |
| PPARA | 6LXA | PPARG/PPARA B and PPARA/PPARD A | shared | A; LBD | EPA | none in PDBQT | yes | NOT_RECOVERABLE |
| PPARD | 5U3Q | PPARA/PPARD B | yes | A; yaml agonist 1 not PEG | 7UJ | none in PDBQT | yes | NOT_RECOVERABLE |

Resolution / method / species in SI Table S2 are **METADATA_ONLY** relative to this audit: they were not re-fetched from RCSB in this pass. Do not treat SI Å values as independently confirmed here.

## Points that remain correct under deposited files

- **9V8H**: peptide atoms are in the PDBQT. That is confirmed by the receptor file, not by caption prose.
- **4EY7**: dimer AB atom count (~10006) is consistent with a two-chain AChE construct.
- **4UDW**: HIL construct in PDBQT; pocket extract `4UDW_pocket_10.0.pdb` is not a full protein PDB.
- **TYK2 3LXP**: yaml/site notes say JH1 ATP, not JH2. This audit did not re-annotate the kinase domain from sequence.
- **6N7A / 6LXA**: shared files across pairs. Repeated targets are not independent replicates.

## Receptor preparation (per slot)

For every slot:

1. Actual receptor PDBQT exists: **yes**
2. Used for production: **yes** (score sources and current scripts point at these files)
3. Original preparation script: AChE has a recovered QC script (METADATA_ONLY). All others **not_recoverable**
4. Protonation: **NOT_RECOVERABLE**
5. Missing loops: **NOT_RECOVERABLE**
6. Waters: none in PDBQT (confirmed absence in file, not a recorded removal protocol)
7. Cofactors/peptides: 9V8H peptide retained; 4EY7 dimer retained; others none in PDBQT
8. Charges: **NOT_RECOVERABLE** (PDBQT charges are whatever was deposited)

Must not claim “Meeko prepared all receptors” or “Protein Preparation Wizard was used consistently.”  
Current manuscript/SI already avoid that claim. If any leftover sentence appears, it is a P1 writing error.

## Boxes

All 14 current primary box JSONs exist. Min edge ≥ 20 Å on all three axes.

EGFR corrected boxes declare `cognate_heavy_atom_AABB_plus_5A_min20`.  
Legacy `3POZ_box.json` / `3RCD_box.json` still exist; in this tree their coordinates match the corrected files. Current EGFR scripts use `*_box_corrected.json`.

PIK3CA `4L23_box.json` has `n_ligand_atoms=42` and `n_heavy_atoms=26` and no construction field. Hydrogen-inclusive AABB cannot be ruled out from the JSON alone; sizes are all 20 Å (min-edge). Status: **METADATA_ONLY**.

Sensitivity boxes (4JPS/5DXT/4JSX, AChE 6ZWI/5DYW/6QAA) are not primary.
