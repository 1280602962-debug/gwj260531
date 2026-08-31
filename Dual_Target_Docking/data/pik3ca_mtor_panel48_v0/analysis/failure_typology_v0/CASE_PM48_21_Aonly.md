# CASE PM48_21 — stubborn A_only (EH40_23 analogue)

## Identity
- **PM48_21** / CHEMBL3814414 / **A_only**
- pChEMBL: PIK3CA **8.70** / mTOR **5.92** (hard negative)
- Scaffold: morpholine-containing fused heterocycle; Tanimoto to PI-103 **0.27**, MCS **13** atoms

## Ranks
| arm | rank |
|-----|------|
| vina_mean | 18 |
| rtm_min_z | **5** |

## Pose QC (RTM-best)
| target | mode | hinge | clash | occ vs X6K | MCS RMSD | RTM z |
|--------|------|-------|-------|------------|----------|-------|
| 4L23 | 1 | yes | 0 | **1.00** | 1.58 | +0.95 |
| 4JT6 | 1 | yes | 0 | **0.97** | **0.92** | +1.40 |

## Why gates fail
Both ends are **geometrically clean ATP-site poses** (hinge + near-cognate occupancy). Clash gate = 0. Weak-end RTM is still high → min/shortfall cannot demote.

## Typology
**T2 chemotype_homolog** — same failure class as EGFR **EH40_23**, different chemotype family (morpholine-ATP vs anilinoquinazoline).

## Protocol implication
Keep as canonical PIK3CA/mTOR T2 hard negative. Report **chemotype warning**, do **not** retune clash to force-drop.
