# GNINA_POCKET_MATCHED_BEST9_VERDICT_V1

Zero new docking. Uses the already-pushed best-of-9 GNINA CNN rescore
(`scores_gnina_best.csv`) and the preserved mode_01 backups
(`scores_gnina_best_mode01_backup.csv`). Computes the **same directional
pocket-matched definition as Vina/RTM** (Methods 2.6): D vs A_only scored
on pocket B; D vs B_only scored on pocket A.

## Correction of the 2026-08-24 push's own terminology

`GNINA_BEST9_STATUS.md` / `compare_gnina_mode01_vs_best9.py` label their
`min(score_A, score_B)`-for-both-contrasts metric "pocket-matched". That is
actually **worst-pocket** (the same convention already used for
`gnina_cnn_min` / `vina_worst` / `rtm_worst` elsewhere in this repo), not the
asymmetric Methods 2.6 pocket-matched definition. Both are legitimate
diagnostics; they must not share the same name. This file reports the true
pocket-matched number; the push's own file should be read as worst-pocket.

## K=4: pocket-matched GNINA, mode_01 vs best-of-9

| pair | channel | n (D/A/B) | D vs A (pocket B) | D vs B (pocket A) | summary_min [95% CI] | Vina ref |
|------|---------|-----------|-------------------:|-------------------:|----------------------:|---------:|
| EGFR/HER2 | mode01 | 28/38/32 | 0.4897 | 0.327 | 0.327 [0.1882, 0.458] | 0.43 |
| EGFR/HER2 | best9 | 28/38/32 | 0.4709 | 0.2902 | 0.2902 [0.1605, 0.426] | 0.43 |
| AChE/BChE | mode01 | 27/25/28 | 0.4859 | 0.4418 | 0.4418 [0.2833, 0.5391] | 0.606 |
| AChE/BChE | best9 | 27/25/28 | 0.5526 | 0.4127 | 0.4127 [0.263, 0.5503] | 0.606 |
| PIK3CA/PIK3CB | mode01 | 28/27/28 | 0.6071 | 0.5536 | 0.5536 [0.3915, 0.6738] | 0.5 |
| PIK3CA/PIK3CB | best9 | 28/27/28 | 0.5701 | 0.5332 | 0.5332 [0.3672, 0.6392] | 0.5 |
| PIK3CA/mTOR | mode01 | 18/14/12 | 0.5794 | 0.6713 | 0.5794 [0.35, 0.7531] | 0.692 |
| PIK3CA/mTOR | best9 | 18/14/12 | 0.6548 | 0.6852 | 0.6548 [0.4238, 0.8088] | 0.692 |

## Does best-of-9 change the qualitative ranking vs Vina?

- **EGFR/HER2**: GNINA mode01=0.327 → best9=0.2902 (Δ=-0.0368); Vina pocket-matched reference=0.43. GNINA remains below Vina.
- **AChE/BChE**: GNINA mode01=0.4418 → best9=0.4127 (Δ=-0.0291); Vina pocket-matched reference=0.606. GNINA remains below Vina.
- **PIK3CA/PIK3CB**: GNINA mode01=0.5536 → best9=0.5332 (Δ=-0.0204); Vina pocket-matched reference=0.5. GNINA at/above Vina reference.
- **PIK3CA/mTOR**: GNINA mode01=0.5794 → best9=0.6548 (Δ=+0.0754); Vina pocket-matched reference=0.692. GNINA remains below Vina.

### One-line verdict

**Moving GNINA from mode-1-only to best-of-9 does not change which pair looks best.** PIK3CA/mTOR remains the strongest pair under GNINA (mode01 0.579 → best9 0.655), still below its own Vina pocket-matched reference (0.692). EGFR/HER2 and AChE/BChE GNINA best9 pocket-matched values are **below 0.5** (0.290 and 0.413), i.e. best-of-9 does not rescue GNINA on the pairs where Vina is also weak or descriptor-explained.

**Correction:** on PIK3CA/PIK3CB, GNINA pocket-matched is marginally *above* the Vina reference under **both** mode01 (0.554 vs 0.500) and best9 (0.533 vs 0.500) — this is not new with best-of-9, and the margin (0.03–0.05) is well within the bootstrap CI width (±0.13–0.15) on either side, so it is not a meaningful "GNINA beats Vina" result. Do not write "GNINA never exceeds Vina on any pair"; write instead: GNINA best-of-9 stays close to chance or below the same-panel Vina value on three of four pairs, and is statistically indistinguishable from Vina (both near chance, overlapping CIs) on PIK3CA/PIK3CB. `RTMScore 与 GNINA 未改变这一格局` remains supported, now with a directional (not just pooled/worst-pocket) GNINA number.

## Stability-check panels (PM48 / PM110)

| panel | channel | n (D/A/B) | D vs A (pocket B) | D vs B (pocket A) | summary_min [95% CI] |
|-------|---------|-----------|-------------------:|-------------------:|----------------------:|
| PM48 | mode01 | 18/14/12 | 0.5794 | 0.6713 | 0.5794 [0.3684, 0.75] |
| PM48 | best9 | 18/14/12 | 0.6548 | 0.6852 | 0.6548 [0.4316, 0.8117] |
| PM110 | mode01 | 30/30/30 | 0.5222 | 0.7133 | 0.5222 [0.3675, 0.6703] |
| PM110 | best9 | 30/30/30 | 0.6133 | 0.6822 | 0.6133 [0.4583, 0.7419] |

PM48 mode01 (0.5794) and PM110 mode01 (0.5222) match the pre-existing frozen
`PM110_VS_PM48.md` / `B_GROUP_VERDICT.md` numbers exactly, confirming the
mode01 backup is faithful to the original rescore. The best9 numbers here
(PM48 0.6548; PM110 0.6133) **supersede** those mode01 GNINA entries in the
Results/SI text describing the PM48↔PM110 stability check and must be used
going forward; the mode01 values are retained here only as a consistency check.

## What this is not

- Not a new docking run; not a change to the frozen K=4 ligand sets.
- Not a claim that GNINA is now a validated general-purpose score.
- Does not touch the primary Vina-based Table 2; GNINA remains a secondary channel.

## Reproduce

```bash
python3 Dual_Target_Docking/data/jcim_bench_v0/scripts/gnina_pocket_matched_best9_v1.py
```

