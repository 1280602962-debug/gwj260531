# Directional AUROC — docking-phase forest (v1, with baselines)

Supersedes the v0 snippet that omitted EGFR and trivial baselines.  
Primary table: [`../tables/directional_with_baselines_v1.csv`](../tables/directional_with_baselines_v1.csv)  
Verdict: [`POST_DOCKING_VERDICT.md`](POST_DOCKING_VERDICT.md)

## Docking channels — `summary_min = min(D/A, D/B)`

| pair | vina_mean | rtm_min_z | gnina_cnn_min | best baseline | docking vs baseline |
|------|----------:|----------:|-------------:|---------------|---------------------|
| PIK3CA/mTOR | **0.671** | 0.520 | 0.563 | heavy 0.463 | **PASS** (all three) |
| AChE/BChE | 0.530 | 0.409 | 0.372 | TPSA 0.753 | FAIL |
| PIK3CA/PIK3CB | 0.412 | 0.439 | 0.506 | heavy 0.599 | FAIL |
| EGFR/HER2 | 0.282 | 0.253 | 0.263 | cLogP 0.482 | FAIL |

Directional detail (vina): EGFR 0.680/0.282; PIK3CA/PIK3CB 0.703/0.412 — one-sided inversion on both kinase-like pairs.

## Note on AChE/BChE TPSA

TPSA alone reaches `summary_min` 0.753 on this strict panel (far above docking). Treat as a chemotype/polarity shortcut to discuss in Limitations; even vs heavy_atoms (0.547), vina (0.530) still fails the baseline gate.
