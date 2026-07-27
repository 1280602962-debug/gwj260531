# CASE PM48_01 — PI-103 (gold)

## Identity
- **PM48_01** / CHEMBL573339 / **dual** / PI-103
- pChEMBL: PIK3CA 8.89 / mTOR 8.28
- Role: dual-end cognate (X6K in 4L23 & 4JT6)

## Ranks
| arm | rank |
|-----|------|
| vina_mean | 9 |
| rtm_min_z | **4** |

## Pose QC (RTM-best)
| target | mode | hinge | clash&lt;2.2 | occ vs X6K | MCS RMSD | RTM z |
|--------|------|-------|------------|------------|----------|-------|
| 4L23 | 1 | yes (2.80) | 0 | 1.00 | 1.09 | +1.02 |
| 4JT6 | **3** | yes (2.80) | 0 | 1.00 | **0.45** | +1.21 |

## Verdict
Cognate recovered. On 4JT6, Vina mode1 is not the crystal-like pose; RTM-best mode3 is. Supports **best-of-9 + rescoring** as protocol components.

## Typology
**Gold / positive control** — not a failure case.
