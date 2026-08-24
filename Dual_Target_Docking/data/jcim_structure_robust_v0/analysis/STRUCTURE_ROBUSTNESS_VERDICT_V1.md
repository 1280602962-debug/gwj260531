# Receptor-structure sensitivity verdict v1

> Not a robustness certificate. Cognate QC: `STRUCTURE_ROBUSTNESS_QC_V1.md`.  
> Design: replace one pocket at a time; the other pocket keeps frozen scores.  
> Metric: pocket-matched directional AUROC; `summary_min = min(D/A, D/B)`; ligand bootstrap B=2000.  
> Combined table: `tables/receptor_realization_two_pair_v1.csv`. Audit: `data/jcim_novelty_v0/analysis/A4_B5_STATISTICAL_AUDIT_V1.md`.

## Cognate QC gate

| target | PDB | ligand | best-of-9 RMSD (Å) | verdict |
|--------|-----|--------|-------------------:|---------|
| PIK3CA | 4JPS | 1LT | 0.607 | PASS |
| PIK3CA | 5DXT | 5H5 | 0.624 | PASS |
| mTOR | 4JSX | 17G (Torin2) | 0.515 | PASS |

All three polymer entities match the intended gene (not chimeras). 3T8M remains excluded.

## Two-pair PIK3CA realization effect (B pocket frozen)

| pair | PIK3CA | kept B | n att/ok/fail | D/A | D/B | summary_min [95% CI] | Δ |
|------|--------|--------|--------------:|----:|----:|----------------------|--:|
| PIK3CA/mTOR | 4L23 | 4JT6 | 48/48/0 | 0.714 | 0.692 | **0.692 [0.464, 0.802]** | — |
| PIK3CA/mTOR | 4JPS | 4JT6 | 48/48/0 | 0.714 | 0.486 | **0.486 [0.259, 0.692]** | −0.206 |
| PIK3CA/mTOR | 5DXT | 4JT6 | 48/48/0 | 0.714 | 0.505 | **0.505 [0.292, 0.696]** | −0.187 |
| PIK3CA/mTOR | 4L23 | 4JSX | 48/48/0 | 0.639 | 0.692 | **0.639 [0.418, 0.776]** | −0.053 |
| PIK3CA/PIK3CB | 4L23 | 2WXF | 100/99/1 | 0.691 | 0.500 | **0.500 [0.347, 0.648]** | — |
| PIK3CA/PIK3CB | 4JPS | 2WXF | 100/99/1 | 0.691 | 0.707 | **0.691 [0.516, 0.779]** | +0.191 |
| PIK3CA/PIK3CB | 5DXT | 2WXF | 100/99/1 | 0.691 | 0.685 | **0.685 [0.506, 0.768]** | +0.185 |

Use deposited CSV CIs. PIK3CA/PIK3CB weak arm switches: original D/B = 0.500; 4JPS bottleneck becomes frozen D/A = 0.691.

`PAB_034` (A_only, CHEMBL5089694) timed out on original 4L23 (`timeout_900s_torsdof=23`) and on both 4JPS and 5DXT (`timeout_600s`). Failure is a docking timeout, not a label filter. No 100-ligand AUROC exists for this pair under any PIK3CA crystal here.

## Honest ceiling

The same PIK3CA crystals **lower** PIK3CA/mTOR discrimination and **raise** PIK3CA/PIK3CB discrimination. Phrase this as a **receptor-realization effect**, not robustness and not a unidirectional collapse. Two opposite-direction examples under a one-pocket-at-a-time design are stronger than a single-pair anecdote; they are not a universal law (K = 4; both pairs share PIK3CA).

Do not advertise structure-invariant dual-target docking performance.

## Files

- `receptors/{4JPS,5DXT,4JSX}_*`
- `tables/scores_vina_mode1_PM48_alt*.csv`
- `tables/scores_vina_mode1_PAB_alt*.csv`
- `tables/pocket_matched_PM48_alt*_v1.csv`
- `tables/pocket_matched_PAB_alt*_v1.csv`
- `tables/receptor_realization_two_pair_v1.csv`
