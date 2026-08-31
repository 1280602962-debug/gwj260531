# jcim_bench_v0 — DualFourClass-Bench pack

JCIM evaluation/benchmark aggregation after docking-phase scoring.

Authorization: [`../protocols/PAIR_ROLES_APPROVED_JCIM.yaml`](../protocols/PAIR_ROLES_APPROVED_JCIM.yaml)  
Claim ceiling: [`CLAIM_CEILING.md`](CLAIM_CEILING.md)  
GNINA: [`analysis/GNINA_STATUS.md`](analysis/GNINA_STATUS.md) (**DONE**, best-of-9 CNN rescore)  
Post-dock verdict: [`analysis/POST_DOCKING_VERDICT.md`](analysis/POST_DOCKING_VERDICT.md)  
**CI analysis pack (v1):** [`analysis/BENCHMARK_ANALYSIS_V1.md`](analysis/BENCHMARK_ANALYSIS_V1.md) · [`analysis/FIGURE_PLAN_V1.md`](analysis/FIGURE_PLAN_V1.md)

| Pair | Pack | Vina | RTM | GNINA |
|------|------|------|-----|-------|
| EGFR/HER2 | `../egfr_her2_panel120_v0/` | existing | existing | best-of-9 done |
| PIK3CA/mTOR | `../pik3ca_mtor_panel48_rdkit_v0/` | 96/96 | done + LigPrep Δ | best-of-9 done |
| AChE/BChE | `../ache_bche_panel_v0/` | 191/200 | done | best-of-9 done |
| PIK3CA/PIK3CB | `../pik3ca_pik3cb_panel_v0/` | 199/200 | done | best-of-9 done |

Primary large pose workspaces remain under  
`/mnt/d/CADD paper exercise/dual target docking/results/` (not all poses committed).

## Reproduce analysis pack (Zenodo-ready)

```bash
cd Dual_Target_Docking
python3 data/jcim_bench_v0/scripts/build_benchmark_analysis_v1.py
python3 data/jcim_bench_v0/scripts/plot_forest_ci_v1.py
```

Key outputs:

| Artifact | Path |
|----------|------|
| Assembled ligands | `tables/assembled_all_pairs_long.csv` |
| Bootstrap CIs | `tables/bootstrap_directional_ci_v1.csv`, `tables/forest_summary_min_ci_v1.csv` |
| Baseline gate Δ±CI | `tables/baseline_gate_bootstrap_v1.csv` |
| Forest / gate figures | `figures/forest_summary_min_ci_v1.{png,pdf}`, `figures/baseline_gate_delta_ci_v1.{png,pdf}` |
| Meta + failure modes | `tables/analysis_meta_v1.json` |

Older point-estimate tables (kept):  
`tables/directional_forest_v0.csv` · `tables/directional_with_baselines_v1.csv`

## Status
Docking-phase **scoring complete** (Vina + RTM + GNINA).  
**Post-dock + CI pack:** only PIK3CA/mTOR beats trivial baselines on the **point estimate**; its Δ CI still spans 0. EGFR and PIK3CA/PIK3CB are significantly below their best trivial baselines. Supports an evaluation/benchmark JCIM narrative, not a universal scorer claim.

**Next (writing, not docking):** English manuscript + Zenodo deposit; optional PM48 expand under strict quotas to tighten CIs.
