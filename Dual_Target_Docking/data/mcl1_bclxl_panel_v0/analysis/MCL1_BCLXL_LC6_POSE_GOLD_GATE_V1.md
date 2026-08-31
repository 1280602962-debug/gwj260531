# MCL1_BCLXL_LC6_POSE_GOLD_GATE_V1

Updated: `2026-08-27T02:37:26Z`

Primary receptors frozen before panel docking: **3WIY** (MCL1) / **3WIZ** (Bcl-xL), cognate **LC6** (Tanaka compound 10).
Protocol: Vina seed `20260727`, exhaustiveness `8`, num_modes `9`, energy_range 3.
Coordinate diagnostic: Hungarian element-matched heavy-atom absolute displacement (AD4 `A`→C). This assignment is not constrained by molecular-graph isomorphism and must not be interpreted as topology-aware symmetry-corrected RMSD.

| target | PDB | top1 Å | best-of-top3 Å | best-of-9 Å | top3 gate (<2Å) |
|--------|-----|-------:|---------------:|------------:|:---------------:|
| MCL1 | 3WIY | 1.689 | 1.689 | 1.677 | 1 |
| BCL2L1 | 3WIZ | 4.17 | 2.011 | 2.011 | 0 |

**Formal gate: UNMET / not validly completed.** The preliminary coordinate diagnostic is insufficient for a formal pose-gold decision because topology-aware atom mapping, physical-validity checks, interaction recovery, and a second seed were not completed. Independently, its 3WIZ point value is above the predeclared 2.0 Å cutoff. Do **not** package this pair as standard screening-performance evidence. Panel docking is a predeclared **applicability stress-test** only.

`panel_role=applicability_stress_test`
