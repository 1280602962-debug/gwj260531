# MANIFEST — ache_bche_panel_v0

## Role
JCIM K=4 development pair (AChE / BChE). Strict labels only; no gray zone in main panel.

## Panel
- N=100: dual=28, A_only=28, B_only=28, neither=16
- Rule: dual both ≥6.5; A_only A≥6.5 & B≤5.5; B_only symmetric; neither both ≤5.5
- SMILES: ChEMBL canonical (salts desalted to largest fragment at prep)

## Receptors (cognate QC PASS)
| End | PDB | Cognate | RMSD @E=8 |
|-----|-----|---------|-----------|
| ACHE | 4EY7 | E20 | 0.339 |
| BCHE | 4BDS | THA | 0.386 |

## Dock protocol
- Prep: RDKit ETKDG (seed 20260727) + meeko
- Vina 1.2.7, E=8, n_modes=9, seed=20260727
- Score: Vina + RTM best-of-9
- GNINA: SKIPPED (see `../jcim_bench_v0/analysis/GNINA_STATUS.md`)

## Exclusions
- Covalent/PROTAC not specially filtered beyond ChEMBL activity tables used for sampling
- Failed failed cognates (6ZWI, etc.) documented under `cognate_qc/` — not used
