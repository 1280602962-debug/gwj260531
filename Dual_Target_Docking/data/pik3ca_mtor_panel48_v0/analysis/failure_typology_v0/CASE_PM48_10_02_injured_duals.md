# CASE PM48_10 / PM48_02 — injured duals (Torin1 / Omipalisib)

## Identity
| id | name | class | pChEMBL PIK3CA / mTOR |
|----|------|-------|------------------------|
| PM48_10 | TORIN1 | dual | 7.85 / 9.54 |
| PM48_02 | OMIPALISIB | dual | 10.72 / 9.74 |

## Ranks
| ligand | vina_mean | rtm_min_z |
|--------|-----------|-----------|
| Torin1 | **#1** | **#31** |
| Omipalisib | **#3** | **#30** |

## Pose QC (failure end = 4L23)
| ligand | RTM mode | hinge | occ | centroid→X6K | note |
|--------|----------|-------|-----|--------------|------|
| Torin1 | 7 | **no** (5.59) | 0.44 | 4.19 Å | Vina mode1 affinity better (−9.55) but RTM prefers mode7 |
| Omipalisib | 7 | **no** (3.88) | 0.64 | 2.81 Å | same pattern |

4JT6 ends remain hinge-positive with mid/high RTM — injury is **asymmetric via weak PIK3CA RTM pose**.

## Typology
**T5 rescoring_injury / pose-family mismatch**  
Classic duals that Vina ranks well; RTM+min demotes because the PIK3CA RTM-best pose is off-hinge / off-cognate.

## Protocol implication
- Do **not** claim RTM-only ranking as sole primary metric in main text.
- Report **dual readout**: vina_mean and rtm_min_z.
- Optional P1: consensus arm (both must be strong) — must show it does not resurrect T2 A_only worse than baseline.
- Optional deep-dive: full 9-mode hinge table on 4L23 for these two ligands.
