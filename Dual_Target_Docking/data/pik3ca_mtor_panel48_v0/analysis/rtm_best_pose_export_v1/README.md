# RTM-best pose export v1 (Maestro viewing pack)

Case poses for dual-target ranking diagnosis (PIK3CA 4L23 / mTOR 4JT6).

## Contents
- `pdbqt/` — original Vina MODELs (RTM-best mode)
- `sdf/` — Maestro-friendly SDF (heavy-atom coords via meeko SMILES IDX + LigPrep template)
- `pdb/` — coordinate PDB (HETATM) for quick open
- `receptors/` — `4L23_PIK3CA_prepared.pdb`, `4JT6_mTOR_prepared.pdb`
- `pose_inventory.csv` — index

## Ligands
| ligand | role | 4L23 RTM-best | 4JT6 RTM-best |
|--------|------|---------------|---------------|
| PM48_26 | stubborn A_only (Top1) | mode_04 | mode_01 |
| PM48_20 | stubborn A_only | mode_02 | mode_05 |
| PM48_21 | stubborn A_only | mode_01 | mode_01 |
| PM48_34 | rescued B_only (WYE-132) | mode_06 | mode_01 |
| PM48_10 | injured dual (Torin1) | mode_07 | mode_01 |
| PM48_02 | injured dual (Omipalisib) | mode_07 | mode_01 |
| PM48_01 | gold standard (PI-103) | mode_01 | **mode_03** (not mode1) |

## Maestro tip
Import receptor PDB + ligand SDF into the same Workspace; keep protein/ligand as separate entries.
