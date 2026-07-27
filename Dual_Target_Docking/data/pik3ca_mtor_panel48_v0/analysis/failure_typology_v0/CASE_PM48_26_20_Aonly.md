# CASE PM48_26 / PM48_20 — stubborn A_only pair (Top1–2)

## Identity
| id | ChEMBL | class | pChEMBL PIK3CA / mTOR |
|----|--------|-------|------------------------|
| **PM48_26** | CHEMBL3646543 | A_only | 8.40 / **5.32** |
| **PM48_20** | CHEMBL3646632 | A_only | 8.70 / **5.92** |

- Pairwise Tanimoto **0.455**, MCS **18** → same amino-triazine-like series
- vs PI-103: Tanimoto ~0.09 (not PI-103 clones; still ATP-competitive chemotypes)

## Ranks (rtm_min_z)
- PM48_26 **#1**, PM48_20 **#2** (above PI-103 #4)

## Pose QC
Both ligands, both ends: **hinge=yes, clash=0**, solid pocket occupancy (0.65–0.78).  
Weak end (mTOR) RTM z still **+1.53 / +1.13**.

## Why shortfall fails
Problem is **false-high inactive-end score**, not imbalance.  
`min(z) − λ|Δz|` left PM48_26 at #1 for all tested λ.

## Typology
**T2 chemotype_homolog (series)** — RTM over-scores mTOR for PIK3CA-selective ATP chemotypes.

## Protocol implication
Primary evidence that **C4 is not closed** on this pair. Main text: diagnostic boundary + warning layer.
