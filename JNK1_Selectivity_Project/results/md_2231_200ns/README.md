# Compound 2231 — 200 ns Amber MD (JNK1/2/3)

Extended MD validation for G2 candidate **2231**, addressing §12 Q5 (extended trajectory on priority G2 molecule).

## Simulation setup

| Item | Value |
|------|-------|
| Engine | AMBER (HMR, `prod_hmr.in`) |
| Length | 200 ns × 3 systems (single replica each) |
| Receptors | JNK1/3ELJ, JNK2/3E7O, JNK3/3TTI |
| Ligand restraint | 2.0 kcal/mol/Å² on `:MOL` heavy atoms (production) |
| Analysis | cpptraj (RMSD/RMSF/RoG/SASA/H-bond); MMPBSA.py igb=8 (frames 15001–20000) |
| Raw trajectories | Local: `2231_200nsMD/` (not in Git — large binary files) |

## Reproducibility

Analysis scripts: `2231_200nsMD/analyze_md_2231.py`, `analyze_md_2231_extended.py`

```bash
cd 2231_200nsMD && python3 analyze_md_2231.py && python3 analyze_md_2231_extended.py
```

## Contents

- `figures/` — RMSD, RMSF, H-bond, MM-GBSA plots
- `tables/` — CSV summaries and `MD_2231_extended_summary.md`

See **§6.5** in `docs/JNK1_PROJECT_REPORT.md` for interpretation and limitations.
