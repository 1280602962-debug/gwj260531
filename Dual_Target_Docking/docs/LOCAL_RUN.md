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

- `data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv` (Table 2 canonical CIs)
- `data/jcim_strengthen_t0t1_v0/analysis/PRIMARY_METRIC_V2.md` (deprecated; do not cite as Table 2)
- `data/jcim_bench_v0/tables/forest_summary_min_ci_v1.csv` (`vina_mean` forest, not pocket-matched Table 2)

Manuscript-facing numeric regression:

```bash
python3 data/jcim_novelty_v0/scripts/validate_revision_v1.py
python3 data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py --check
python3 scripts/primary/bootstrap_primary.py
```

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

Detectable-effect simulation (zero docking; ~10 min):

```bash
python3 data/jcim_novelty_v0/scripts/detectable_effect_simulation_v1.py
python3 data/jcim_novelty_v0/scripts/plot_detectable_effect_and_workflow_v1.py
```

Track B five-pair local Vina (F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, PPARA/PPARD) is specified in `data/jcim_chembl_universe_v0/analysis/DOCKING_PLAN_V1.md` and `tables/track_b_local_run_v1.yaml`. Cloud has no Vina. Do not run F2/F10 protonation sensitivity or 2Y3A E=32 cognate retest in this pack.

Independent pose-generation and optional PLIF **results are deposited**:
`data/jcim_independent_dock_v0/` (GNINA docking search, not Vina-pose rescore)
and `data/jcim_structure_robust_v0/analysis/plif_v1/` (geometric occupancy).
SOPs remain in `docs/AGENT_COMMAND_INDEPENDENT_POSE_GENERATION_V1.md` and
`docs/AGENT_COMMAND_RECEPTOR_PLIF_V1.md`. Do not rerun them in a cloud
environment without GNINA/smina and frozen pose workspaces.

## 5. Phase-1 revision analyses (still zero docking)

```bash
bash scripts/run_phase1_revision.sh
```

This rebuilds document-blocked CV, the assay-context machine extract, the frozen literature-year split, cognate inventory, assembled manuscripts, and the checksum manifest. It does not dock new ligands or mint a Zenodo DOI.

## 6. Docking scripts vs deposited scores

Many historical dock/RTM/GNINA drivers still point at the original local workspace (`/mnt/d/CADD paper exercise/...`, `/home/gwj/miniconda3/...`). They are the as-run recipes, not portable cloud commands. Publication numbers are regenerated from the deposited ligand-level CSVs under `data/*/tables/` by the zero-dock analysis scripts.

- Canonical Table 2 CIs: `unified_threshold_sensitivity_v2.csv` (`label_rule=theta_6.0`).
- Do not cite `pocket_matched_directional_v1.csv` or `PRIMARY_METRIC_V2.md` as Table 2.
- Multi-seed Table S54: `analyze_multiseed_vina_v2.py`. `analyze_multiseed_vina_v1.py` uses a different Dual-versus-neither estimand and refuses to overwrite unless `--legacy` is passed.


## Track B local Vina (five pairs)

See `data/jcim_chembl_universe_v0/local_track_b_v0/analysis/TRACK_B_DIRECTIONAL_AUROC_V1.md`.
