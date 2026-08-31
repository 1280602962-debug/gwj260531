# P0 missing-content fill (frozen scores only)

No new docking. Bootstrap B=2000, seed 20260729, ligand units.
Point AUROCs are checksummed against `unified_threshold_sensitivity_v2.csv` (θ=6.0),
`pocket_matched_directional_v1.csv`, `holdout_pocket_matched_v1.csv`, and
`forest_summary_min_ci_v1.csv`.

## What was filled

1. Endpoint hierarchy (primary / pre-specified secondary / robustness / exploratory).
2. Formal frozen vs holdout table; EGFR/HER2 holdout = not eligible.
3. Paired Δ = summary_min(matched) − summary_min(wrong) on the **same** bootstrap sample.
4. Pose-fairness table (9 Vina poses; RTM and GNINA best-of-9; GNINA mode-1 historical).
5. Crystal-swap criteria remain pre-specified (true gene, non-chimera, cognate best-of-9 < 2 Å; 3T8M excluded).
6. Scaffold vs random ECFP4 merged; leakage is small — not a hunt for a leakier split.
7. Pocket-matched (not vina_mean) paired Δ vs the strongest trivial descriptor.
8. Exploratory top-10 hard-negative counts on **vina_mean** (not Table 2).

## Paired Δ (matched − wrong)

| set | pair | matched | wrong | Δ | 95% CI | excludes 0? |
|-----|------|--------:|------:|--:|--------|-------------|
| main_panel | EGFR/HER2 | 0.4297 | 0.26 | 0.1697 | [0.06, 0.2803] | True |
| main_panel | AChE/BChE | 0.6058 | 0.4444 | 0.1614 | [0.037, 0.269] | True |
| main_panel | PIK3CA/PIK3CB | 0.5 | 0.3489 | 0.1511 | [-0.0215, 0.3105] | False |
| main_panel | PIK3CA/mTOR | 0.6921 | 0.6019 | 0.0902 | [-0.1222, 0.2626] | False |
| unused_pool_holdout | AChE/BChE | 0.6175 | 0.6425 | -0.025 | [-0.1119, 0.0714] | False |
| unused_pool_holdout | PIK3CA/PIK3CB | 0.425 | 0.52 | -0.095 | [-0.2814, 0.1143] | False |
| unused_pool_holdout | PIK3CA/mTOR | 0.765 | 0.7875 | -0.0225 | [-0.1165, 0.079] | False |

Main-panel point Δ is positive on all four pairs. Holdout point Δ is negative on all three
eligible pairs (wrong ≥ matched). Whether the holdout Δ CI excludes 0 is an empirical
result of this bootstrap, not a claim that the paradox is explained.

## Pocket-matched Δ (Vina − best descriptor)

| pair | descriptor | Vina | descriptor | Δ | 95% CI | excludes 0? |
|------|------------|-----:|-----------:|--:|--------|-------------|
| EGFR/HER2 | clogp | 0.4297 | 0.4821 | -0.0524 | [-0.2, 0.1155] | False |
| AChE/BChE | tpsa | 0.6058 | 0.7333 | -0.1275 | [-0.3039, 0.0493] | False |
| PIK3CA/PIK3CB | heavy | 0.5 | 0.6217 | -0.1217 | [-0.3197, 0.0891] | False |
| PIK3CA/mTOR | heavy | 0.6921 | 0.463 | 0.2291 | [-0.0105, 0.4352] | False |

This is **not** `baseline_gate_bootstrap_v1.csv` (that file uses pooled `vina_mean`;
EGFR 0.2824 ≠ Table 2 0.4297).

## Frozen vs holdout

See `tables/frozen_vs_holdout_v1.csv`. EGFR/HER2 holdout remains not eligible.

## ML leakage check

Mean (random − scaffold) on the eight directional contrasts:
 0.0112.
Do not hunt a split that inflates the gap.

## Not done (and not invented)

- EGFR unused-pool holdout
- BindingDB/`as_is` docking panel
- 1000 independent panels
- LigPrep, PLIF, B=10000 as if it changed conclusions
- 4-class accuracy as primary
- receptor-only scorer as hard P0

Holdout rows written: 8.
