# HOLD stage comparison — UNIFORM_RERUN_V4_2_20260921

Machine-readable companion: `HOLD_STAGE_COMPARISON.json`.
Meeko covalent threshold = 1.2 × (r_C + r_N) = 1.764 Å; 1.2 × (r_C + r_O) = 1.704 Å.

No atoms were moved. No residues or RETAIN components were deleted.
`--forgive_extra_bonds`, `--allow_bad_res`, `--delete_bad_res` were not used. No new PDB was searched.

## 6N7A A:1024 LEU / A:1041 ARG

| Stage | A:1041 atoms | Closest contact | Meeko covalent? |
|---|---|---|---|
| Frozen CIF / raw | N,CA,C,O,CB only; no altLoc | LEU CD2–ARG CB = 6.443 Å | No |
| PDBFixer addMissingAtoms | +CG,CD,NE,CZ,NH1,NH2 | LEU CD1–ARG NH2 = 1.690 Å | **Yes (introduced)** |
| PTR fallback | unchanged (not PTR) | 1.690 Å | Yes |
| After observed-atom revert | N,CA,C,O,CB | 6.443 Å | No |

- Chain continuity A:1024 and A:1041: peptide C–N 1.34–1.36 Å, no internal gap.
- Origin: **PDBFixer-invented sidechain**, not present in the crystal.
- Fact-driven fix: drop fixer-added atoms on A:1041; `--set_template A:1041=OBS_BB_CB` (observed heavy atoms = backbone + CB).

## 8BXH A:839 MET / A:918 TYR

| Stage | A:839 atoms | Closest contact | Meeko covalent? |
|---|---|---|---|
| Frozen CIF / raw | N,CA,C,O,CB only; no altLoc | MET CB–TYR OH = 3.058 Å | No |
| PDBFixer addMissingAtoms | +CG,SD,CE | MET CE–TYR OH = 1.512 Å | **Yes (introduced)** |
| PTR fallback | A:918 remains crystal TYR | 1.512 Å | Yes |
| After observed-atom revert | N,CA,C,O,CB | 3.058 Å | No |

- A:918 is crystal TYR (not a PTR conversion product). Coordinates of OH unchanged.
- Origin: **PDBFixer-invented CE**, not present in the crystal.
- Fact-driven fix: drop fixer-added atoms on A:839; `--set_template A:839=OBS_BB_CB`.

## 4UDW GLN I:565 (hirudin, RETAIN)

| Stage | I:565 atoms | RDKit/Meeko SanitizeMol |
|---|---|---|
| Frozen CIF / raw | N,CA,C,O,CB,CG; no altLoc | Pass |
| PDBFixer + OXT | +CD,OE1,NE2,+OXT | Fail: OE1 valence 3 (OE1–CG 1.024 Å); OXT valence 4 |
| After OXT-not-in-raw drop | +CD,OE1,NE2 | Fail: OE1 bonded to CD+CG+NE2 |
| After observed-atom revert | N,CA,C,O,CB,CG | Pass |

- I:565 is the C-terminus of chain I (prev I:564 C–N = 1.326 Å). Crystal has no OXT.
- Component I remains RETAIN.
- Origin: **PDBFixer-invented amide atoms** with impossible intra-residue geometry.
- Fact-driven fix: drop fixer-added CD/OE1/NE2; `--set_template I:565=OBS_BB_CB_CG`.

## Uniform rule applied to all 14

PDBFixer still runs existing-residue heavy-atom repair. A fixer-added atom is kept only if it does **not** create a new Meeko-covalent inter-residue contact and the completed residue sanitizes. Otherwise the residue is reverted to its frozen-raw observed heavy atoms and parameterized with an explicit observed-fragment template. Other fixer repairs on the same receptors are kept.
