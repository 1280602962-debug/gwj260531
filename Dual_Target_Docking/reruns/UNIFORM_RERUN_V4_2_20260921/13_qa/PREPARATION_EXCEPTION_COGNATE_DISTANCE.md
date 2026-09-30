# Preparation exception vs cognate ligand

Gate: incomplete residue min heavy-atom distance to cognate ≤ 5 Å → HOLD_REVIEW_REQUIRED (do not submit that receptor). Otherwise DISTAL_PREPARATION_EXCEPTION and continue.

| PDB | Residue | Crystal-observed atoms | Min heavy distance to cognate | Closest pair | Gate |
|---|---|---|---|---|---|
| 6N7A | ARG A:1041 | N,CA,C,O,CB | 13.835 Å | CB – KEV A:1203 | DISTAL_PREPARATION_EXCEPTION |
| 8BXH | MET A:839 | N,CA,C,O,CB | 21.039 Å | CB – C87 A:1201 | DISTAL_PREPARATION_EXCEPTION |
| 4UDW | GLN I:565 | N,CA,C,O,CB,CG | 28.249 Å | CG – N6L H:1249 | DISTAL_PREPARATION_EXCEPTION |

Final PDBQT checks (all three):

- Forbidden fixer-invented sidechain atoms absent (6N7A no CG/CD/NE/CZ/NH1/NH2; 8BXH no CG/SD/CE; 4UDW no CD/OE1/NE2).
- No nonstandard inter-residue covalent contacts; only sequential peptide C–N.
- Residue SanitizeMol valence OK.

`any_hold_review_required = false`. All 14 receptors remain redocking-eligible.
