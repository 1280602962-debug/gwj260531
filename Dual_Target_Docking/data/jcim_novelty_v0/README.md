# DualFourClass-Bench novelty analyses (frozen scores)

Scripts in `scripts/` re-use pocket-matched AutoDock Vina affinities already in the four main-panel ablation tables. No new docking.

| Script | Outputs |
|---|---|
| `benchmark_formulation_v1.py` | Dual-vs-neither comparator vs directional AUROC; chemotype-constrained hard-negatives; incremental logistic models; mixed-library EF |
| `claim_hardening_v1.py` | min/arithmetic/geometric/harmonic aggregation; all-four descriptors; docking failure census |
| `plot_formulation_comparison_v1.py` | SI figure from the formulation CSV |
| `assay_aggregation_max_vs_median_v1.py` | full-panel max vs median pChEMBL (Table S29; caches under `cache/chembl_activity/`) |
| `build_master_results_table_v1.py` | `MASTER_RESULTS_TABLE.csv`, A4 agreement table, two-pair receptor table |

Run from the `Dual_Target_Docking` directory:

```bash
python3 data/jcim_novelty_v0/scripts/benchmark_formulation_v1.py
python3 data/jcim_novelty_v0/scripts/claim_hardening_v1.py
python3 data/jcim_novelty_v0/scripts/build_master_results_table_v1.py
python3 data/jcim_novelty_v0/scripts/plot_formulation_comparison_v1.py
```

Interpretation: `analysis/FORMULATION_COMPARISON_VERDICT_V1.md` and `docs/JCIM_NOVELTY_GAP_AND_PLAN_V1.md`.
