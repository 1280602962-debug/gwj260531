# Local Run Guide

This repository has two practical usage levels:

1. **Zero-dock reproducibility**: rebuild the manuscript-ready analysis tables and figures from committed CSVs.
2. **Docking / rescoring workflows**: rerun Vina, RTMScore, GNINA, holdout, or crystal-swap experiments. These require heavier local tooling and, in some cases, external local workspaces under `/mnt/d/...`.

## 1. Minimal local analysis run

From the repository root:

```bash
python3 scripts/check_local_env.py
bash scripts/run_local_repro.sh
```

Minimal Python packages are listed in `requirements-analysis.txt`.

Expected outputs include:

- `data/jcim_strengthen_t0t1_v0/analysis/PRIMARY_METRIC_V2.md`
- `data/jcim_bench_v0/tables/forest_summary_min_ci_v1.csv`
- `data/jcim_bench_v0/figures/forest_summary_min_ci_v1.png`

## 2. What the minimal run does not do

The local reproducibility entry point intentionally does **not**:

- rebuild docking poses
- call AutoDock Vina
- call GNINA
- call RTMScore
- fetch large pose workspaces from `/mnt/d/CADD paper exercise/dual target docking/results/`

Those workflows remain in the per-pack scripts under `data/*/scripts/`.

## 3. Optional heavy dependencies

These are only needed for extended workflows:

- `vina`
- `gnina`
- `obabel`
- `meeko`
- RTMScore local environment / weights
- `biopython` for some structure-context / superposition scripts

See `data/jcim_strengthen_t0t1_v0/ENV_PIN.md` for the as-run environment snapshot.

## 4. Typical commands beyond the minimal run

Main CI / benchmark pack:

```bash
python3 data/jcim_bench_v0/scripts/build_benchmark_analysis_v1.py
python3 data/jcim_bench_v0/scripts/plot_forest_ci_v1.py
```

GNINA best-of-9 comparison summary:

```bash
python3 data/jcim_bench_v0/scripts/compare_gnina_mode01_vs_best9.py
```

Supply audit:

```bash
python3 data/jcim_j0j1_v0/scripts/run_j0_supply_audit.py
```
