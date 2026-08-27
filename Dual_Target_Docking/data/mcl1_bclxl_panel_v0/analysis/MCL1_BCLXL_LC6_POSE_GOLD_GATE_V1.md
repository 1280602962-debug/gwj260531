# MCL1_BCLXL_LC6_POSE_GOLD_GATE_V1

Updated: `2026-08-27T02:37:26Z`

Primary receptors frozen before panel docking: **3WIY** (MCL1) / **3WIZ** (Bcl-xL), cognate **LC6** (Tanaka compound 10).
Protocol: Vina seed `20260727`, exhaustiveness `8`, num_modes `9`, energy_range 3.
RMSD: Hungarian element-matched heavy-atom absolute RMSD (AD4 `A`→C).

| target | PDB | top1 Å | best-of-top3 Å | best-of-9 Å | top3 gate (<2Å) |
|--------|-----|-------:|---------------:|------------:|:---------------:|
| MCL1 | 3WIY | 1.689 | 1.689 | 1.677 | 1 |
| BCL2L1 | 3WIZ | 4.17 | 2.011 | 2.011 | 0 |

**Gate: FAIL.** At least one end best-of-top3 ≥ 2.0 Å. Per `JCIM_NO_WETLAB_DEEP_PLAN_V2`, do **not** package this pair as standard screening-performance evidence. Panel docking (if run) is a predeclared **applicability stress-test** only.

`panel_role=applicability_stress_test`
