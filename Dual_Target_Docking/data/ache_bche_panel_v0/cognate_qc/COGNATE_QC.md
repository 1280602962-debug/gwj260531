# Cognate QC — AChE / BChE receptor freeze (JCIM Step 2B)

| target | PDB | ligand | E | best_of_9 RMSD (Å) | pass |
|--------|-----|--------|--:|--------------------:|:----:|
| ACHE | **4EY7** | E20 (donepezil) | 8 | **0.339** | PASS |
| BCHE | **4BDS** | THA (tacrine) | 8 | **0.386** | PASS |

### Rejected / failed candidates
| target | PDB | ligand | note |
|--------|-----|--------|------|
| BCHE | 6ZWI | QRH | RMSD≈2.3–2.5 @E8/16 |
| BCHE | 6QAA / 5DYW | HUN / 5HF | ligand PDBQT parse fail (ROOT tag) |

**Frozen for Step 3:** ACHE=`4EY7`, BCHE=`4BDS`. Panel docking E TBD (default 8; raise only if cognate requires).
