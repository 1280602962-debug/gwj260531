# Structure robustness verdict v1

> Cognate QC: `STRUCTURE_ROBUSTNESS_QC_V1.md`  
> Panel redock: replace one pocket at a time on frozen PM48 ligands; other pocket keeps original scores.  
> Metric: pocket-matched directional AUROC; `summary_min = min(D/A, D/B)`; ligand bootstrap B=2000.

## Cognate QC gate

| target | PDB | ligand | best-of-9 RMSD (Å) | verdict |
|--------|-----|--------|-------------------:|---------|
| PIK3CA | 4JPS | 1LT | 0.607 | PASS |
| PIK3CA | 5DXT | 5H5 | 0.624 | PASS |
| mTOR | 4JSX | 17G (Torin2) | 0.515 | PASS |

All three polymer entities match the intended gene (not chimeras). 3T8M remains excluded.

## Panel-level summary_min after pocket swap

Main panel (4L23 + 4JT6): **summary_min = 0.692** (D/A = 0.714, D/B = 0.692).

| alt receptor | replaced pocket | summary_min [95% CI] | D vs A | D vs B | Δ vs main |
|--------------|-----------------|---------------------:|-------:|-------:|----------:|
| 4JPS | A (PIK3CA) | **0.486** [0.259, 0.692] | 0.714 | 0.486 | −0.206 |
| 5DXT | A (PIK3CA) | **0.505** [0.292, 0.696] | 0.714 | 0.505 | −0.187 |
| 4JSX | B (mTOR) | **0.639** [0.418, 0.776] | 0.639 | 0.692 | −0.053 |

Notes:
- Swapping PIK3CA (4JPS/5DXT) leaves D/A unchanged (still scored on frozen 4JT6) but collapses D/B to ~0.49–0.50.
- Swapping mTOR (4JSX) leaves D/B unchanged (still scored on frozen 4L23) and mildly lowers D/A (0.714 → 0.639); summary_min remains above 0.5 at the point estimate but the CI includes values near chance.

## Honest ceiling (Plan §3.4)

Alternative **PIK3CA** crystals that pass cognate QC **do not preserve** the main-panel summary_min. This is **receptor dependence** on the PIK3CA end: the favourable PM signal on 4L23 is not automatically reproduced on 4JPS/5DXT with the same ligands and protocol.

Alternative **mTOR** crystal 4JSX is closer to the main result (Δ ≈ −0.05) but does not strengthen the claim; CI remains wide.

Claim implication: phrase PM directional AUROC as tied to the frozen cognate-QC’d pair **4L23/4JT6**, with explicit sensitivity under alternate deposited structures—do not advertise structure-invariant dual-target docking performance.

## Files

- `receptors/{4JPS,5DXT,4JSX}_*`
- `tables/scores_vina_mode1_PM48_alt*.csv`
- `tables/pocket_matched_PM48_alt*_v1.csv`
- `scripts/cognate_qc_alt_pm_v1.py`, `scripts/redock_pm48_alt_pik3ca_v1.py`
