# MANIFEST — pik3ca_pik3cb_panel_v0

## Role
JCIM K=4 isoform control pair (PIK3CA / PIK3CB). Narrative: near-paralogs; expect weak directional separation.

## Panel
- N=100: dual=28, A_only=28, B_only=28, neither=16
- Rule: strict 6.5 / 5.5 (same as AChE pack)
- SMILES: ChEMBL canonical (salts desalted at prep)

## Receptors (cognate QC PASS)
| End | PDB | Cognate | Notes |
|-----|-----|---------|-------|
| PIK3CA | 4L23 | (reuse PM freeze) | same as PM48 |
| PIK3CB | 2WXF | 039 | RMSD 0.405 @E=8/16 |

## Dock protocol
- Prep: RDKit ETKDG (seed 20260727) + meeko
- Vina 1.2.7, E=8, n_modes=9, seed=20260727
- Score: Vina + RTM best-of-9
- GNINA: SKIPPED

## Failed candidates
- 2Y3A/GD9, 4BFR/J82 documented under `cognate_qc/` — not used
