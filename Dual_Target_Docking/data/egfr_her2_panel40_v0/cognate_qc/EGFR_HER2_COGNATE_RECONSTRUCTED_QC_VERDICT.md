# EGFR/HER2 cognate reconstructed QC verdict

Original nine-mode cognate PDBQTs were not found on the local calculation
disk or in git. Poses below were **re-redocked** under the frozen protocol
(Vina, seed 20260727, exhaustiveness 8, 9 modes) and must be cited
as reconstructed QC, not as the original production gate artifact.

| target | top-1 RMSD | top-3 min | best-of-9 | pass top-1 (<2Å) | pass top-3 |
|--------|------------:|----------:|----------:|:----------------:|:----------:|
| 3POZ | 7.7795 | 4.9358 | 0.6802 | 0 | 0 |
| 3RCD | 1.7294 | 1.1126 | 1.1126 | 1 | 1 |

Crystal ligand: residue `03P` (TAK-285) in both 3POZ and 3RCD.
Artifacts: `data/egfr_her2_panel40_v0/cognate_qc/`.
