# Benchmark MM-GBSA Calibration Report

## Data sources
- Per-PDB Prime MM-GBSA: `jnk_docking_export/mmgbsa_results/ddg_selectivity_detail.tsv`
- Glide XP scores: `analysis_summary/benchmark_scores_by_pdb.csv`
- Experimental IC50 / ΔpIC50: `results/validation/benchmark_deltas.csv`

## Isoform aggregation
| Isoform | PDBs (ensemble) | VSW primary |
|---------|-----------------|-------------|
| JNK1 | mean(3ELJ, 4L7F) | 3ELJ |
| JNK2 | 3E7O | 3E7O |
| JNK3 | mean(3TTI, 4WHZ) | 3TTI |

## Δsel definition (same as docking pipeline)
`delta_sel = min(iso2, iso3) - iso1` — **positive => computational JNK1 preference**
(For MM-GBSA dG_bind: more negative = stronger binding; formula unchanged.)

## Potency reference (MM-GBSA dG_bind, kcal/mol)
- Non-covalent benchmark **median mmgbsa_JNK1 (ensemble)**: **-46.35**
- E1 @ 3ELJ (strong JNK1): **-63.44**
- Shortlist F2 gate used historically: **-51.60** (benchmark median)

## Selectivity magnitude (|Δsel_mmgbsa_vsw|, non-covalent)
- Median |Δsel|: **8.13** kcal/mol
- 75th percentile |Δsel|: **14.79** kcal/mol
- Max |Δsel|: **22.00** kcal/mol

## Correlation with experimental selectivity (ΔpIC50_sel)
| Metric | Spearman ρ | p | n |
|--------|------------|---|---|
| delta_sel_dock_ens | 0.679 | nan | 7 |
| delta_sel_mmgbsa_ens | 0.643 | nan | 7 |
| delta_sel_dock_vsw | 0.750 | nan | 7 |
| delta_sel_mmgbsa_vsw | 0.786 | nan | 7 |

## Direction accuracy vs experiment (non-covalent, n with IC50)
- **dock_ens**: 29% (2/7)
- **mmgbsa_ens**: 29% (2/7)
- **dock_vsw**: 43% (3/7)
- **mmgbsa_vsw**: 43% (3/7)

## Suggested calibration thresholds
- **Potency gate (JNK1 MM-GBSA @ 3ELJ)**: ≤ **-46.4** kcal/mol (benchmark median)
- **Conservative Δsel_MMGBSA gate**: ≥ **22.2** kcal/mol (vs current VSW **2.0** — benchmark noise median 8.1)
- **Do not use MM-GBSA Δsel alone for isoform ranking** — direction accuracy ~25–40% on benchmarks; align with wet IC50.

## Per-ligand calibration table (key columns)

| Ligand | Profile | mmgbsa J1/J2/J3 (VSW) | Δsel_mmgbsa | Δsel_dock | exp_dir | match_mmgbsa |
|--------|---------|------------------------|-------------|-----------|---------|--------------|
| AS602801 | pan-JNK | -47.2/-54.3/-62.6 | -15.43 | -4.97 | JNK1 | False |
| CC-90001 | pan-JNK | -50.7/-51.3/-40.0 | -0.59 | -0.25 | JNK1 | False |
| CC-930 | JNK2/JNK3-biased | -53.4/-47.9/-75.4 | -22.00 | -5.31 | JNK23 | True |
| E1 | JNK1-preferring | -63.4/-61.0/-59.5 | 2.49 | 3.46 | JNK1 | True |
| Q63 | JNK1/JNK3-over-JNK2 | -42.4/-46.1/-48.4 | -6.00 | 0.95 | JNK1 | False |
| SP600125 | pan-JNK | -47.4/-55.6/-43.5 | -8.13 | -1.53 | JNK1 | False |
| TCS JNK 6O | JNK1-preferring | -52.5/-38.3/-31.6 | 14.14 | -1.15 | JNK1 | True |

## Output files
- `benchmark_mmgbsa_by_pdb.csv` — long format (ligand × PDB)
- `benchmark_mmgbsa_calibration.csv` — wide calibration master table
- `benchmark_selectivity_correlations.csv` — Spearman vs ΔpIC50
- `suggested_thresholds.json` — machine-readable threshold suggestions

## Raw Prime outputs (per PDB)
- `benchmarks_3ELJ_prime_mmgbsa_1/benchmarks_3ELJ_prime_mmgbsa_1-out.csv`
- `benchmarks_3E7O_prime_mmgbsa_2/benchmarks_3E7O_prime_mmgbsa_2-out.csv`
- `benchmarks_3TTI_prime_mmgbsa_3/benchmarks_3TTI_prime_mmgbsa_3-out.csv`
- `benchmarks_4L7F_prime_mmgbsa_1/benchmarks_4L7F_prime_mmgbsa_1-out.csv`
- `benchmarks_4WHZ_prime_mmgbsa_3/benchmarks_4WHZ_prime_mmgbsa_3-out.csv`