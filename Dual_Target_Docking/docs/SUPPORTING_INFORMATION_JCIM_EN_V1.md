# Supporting Information (English, merged submission draft)

**Article companion:** `MANUSCRIPT_JCIM_EN.md`  
**Numeric rule:** cells below are read from frozen CSVs and rounded to three decimals, matching the main text. Per-ligand long tables, exploratory slices, and demoted pairs are archived with code and SHA-256 checksums (Note S14).  
**Legacy map:** `data/manuscript_lock/SI_TABLE_MERGE_MAP_v1.csv` (former Tables S1–S54 → Tables S1–S13).

**Main-text tables (not repeated here):** Table 1 panel composition; Table 2 eight-pair directional AUROC and four descriptors; Table 3 Dual-versus-neither versus directional contrast.

---

## Table S1. Software, seeds, and docking settings

| Item | Value |
|------|------|
| Production Python / RDKit / Meeko | Frozen local env; RDKit 2026.3.1; Meeko 0.7.1 |
| Zero-dock reanalysis | Python 3.12.13; RDKit 2026.3.5; NumPy 2.5.2; pandas 3.0.5; SciPy 1.18.1; scikit-learn 1.9.0 (`requirements-analysis.txt`) |
| AutoDock Vina | 1.2.7; default `vina`; n_modes = 9; energy_range = 3 kcal mol⁻¹ |
| GNINA | 1.3.2; CPU `--no_gpu` |
| RTMScore | `rtmscore_model1` |
| Ligand prep | Desalt (largest organic) → AddHs → ETKDGv3 (seed 20260727) → MMFF ≤200 → Meeko PDBQT |
| Panel / holdout seeds | 20260729 / 20260731 |
| Vina production seed | 20260727; five-seed adds 20260811–20260814 |
| Exhaustiveness | PIK3CA/mTOR = 16; all other primary panels and holdouts = 8 |
| Bootstrap | B = 2000; ligand-level non-stratified percentile; SHA-256 sub-seeds |
| Box | Cognate AABB + 5 Å; minimum edge 20 Å |
| Cognate gate | best-of-9 heavy-atom RMSD < 2.0 Å (search coverage, not top-1 ranking) |

Ligands need both-end scores for directional AUROC; n_scored may be below n_panel (Table 1). Failures concentrate among large or flexible ligands and are not silent missingness. Source: `ENV_PIN.md`; `docking_failure_census_v1.csv`.

---

## Table S2. Primary receptor boxes and cognate redocking RMSD

The eight pairs use 14 PDB slots (JAK1 6N7A and PPARA 6LXA are shared). PIK3CA/mTOR used E = 16 because 4JT6 failed the gate at E = 8. EGFR 3POZ is reconstructed QC (original nine-mode production files were not recovered). 9V8H is a PPARγ LBD–BRL–PG08-NL ternary complex; the peptide was retained.

**S2a. Docking boxes (Å)**

| Protein | PDB | Cognate | center (x, y, z) | size (x, y, z) |
|------|-----|----------|------------------|----------------|
| PIK3CA | 4L23 | X6K | 32.443, 45.431, 42.139 | 20.000, 20.000, 20.000 |
| mTOR | 4JT6 | X6K | 51.949, 0.065, −47.707 | 20.332, 20.000, 20.000 |
| AChE | 4EY7 | E20 | −13.988, −43.906, 27.108 | 23.341, 20.000, 20.355 |
| BChE | 4BDS | THA | 133.076, 116.113, 41.335 | 20.000, 20.000, 20.000 |
| EGFR | 3POZ | 03P | 18.680, 32.127, 11.865 | 22.189, 20.000, 22.836 |
| HER2 | 3RCD | 03P | 12.463, 3.371, 27.619 | 23.222, 23.155, 20.000 |
| F2 | 4UDW | N6L | 129.739, −14.402, 164.946 | 20.600, 20.000, 20.000 |
| F10 | 2JKH | BI7 | −3.643, 10.770, 22.771 | 20.000, 20.000, 21.118 |
| JAK1 | 6N7A | KEV | 5.911, −24.317, 21.495 | 20.000, 20.597, 20.000 |
| TYK2 | 3LXP | IZA | −4.498, 25.578, −31.064 | 20.000, 20.000, 20.000 |
| JAK2 | 8BXH | C87 | −12.589, 14.594, 21.263 | 24.069, 20.000, 20.000 |
| PPARG | 9V8H | BRL | −5.454, −5.373, −23.591 | 20.000, 20.140, 20.000 |
| PPARA | 6LXA | EPA | 11.998, 5.742, −7.443 | 20.000, 20.200, 23.432 |
| PPARD | 5U3Q | 7UJ | 40.599, 0.533, 135.353 | 20.223, 20.000, 21.346 |

**S2b. Cognate redocking (production exhaustiveness)**

| Protein | PDB | E | top-1 RMSD (Å) | top-3 (Å) | best-of-9 (Å) | Gate |
|------|-----|--:|---------------:|----------:|--------------:|------|
| PIK3CA | 4L23 | 16 | 0.624 | 0.624 | 0.624 | pass |
| mTOR | 4JT6 | 16 | 7.118 | 0.445 | 0.445 | pass (E = 8 best-of-9 = 5.003, fail) |
| AChE | 4EY7 | 8 | 0.339 | 0.339 | 0.339 | pass |
| BChE | 4BDS | 8 | 4.794 | 0.386 | 0.386 | pass |
| EGFR | 3POZ | 8 | 9.505 | 6.227 | 0.760 | search coverage pass; top-1/top-3 fail |
| HER2 | 3RCD | 8 | 1.855 | 1.394 | 1.394 | pass |
| F2 | 4UDW | 8 | 0.382 | 0.382 | 0.382 | pass |
| F10 | 2JKH | 8 | 0.658 | 0.658 | 0.658 | pass |
| JAK1 | 6N7A | 8 | 0.459 | 0.459 | 0.459 | pass |
| TYK2 | 3LXP | 8 | 0.197 | 0.197 | 0.197 | pass |
| JAK2 | 8BXH | 8 | 4.064 | 0.807 | 0.807 | pass |
| PPARG | 9V8H | 8 | 6.493 | 1.459 | 1.459 | pass |
| PPARA | 6LXA | 8 | 7.508 | 1.098 | 1.098 | pass |
| PPARD | 5U3Q | 8 | 1.452 | 1.452 | 1.452 | pass |

Source: panel `boxes/*.json`; `cognate_rank_rmsd_reaudit_v1.csv`; `layer3_cognate_rmsd_v1.csv`.

---

## Table S3. Activity-threshold and pChEMBL-aggregation sensitivity

Relabeling on frozen Vina scores. Primary analysis is θ = 6.0 (Table 2). Intervals widen when a selective class shrinks.

| Pair | Label rule | n (D / A / B) | summary_min | 95% CI |
|------|----------|--------------:|------------:|--------|
| EGFR/HER2 | θ = 5.5 | 69 / 22 / 10 | 0.425 | [0.242, 0.626] |
| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 | 0.430 | [0.282, 0.578] |
| EGFR/HER2 | θ = 6.5 | 26 / 29 / 29 | 0.460 | [0.304, 0.609] |
| EGFR/HER2 | strict 6.5/5.5 | 26 / 17 / 7 | 0.324 | [0.138, 0.525] |
| AChE/BChE | θ = 5.5 / 6.0 / 6.5 / strict | 27 / 25 / 28 | 0.606 | all CIs include 0.5 |
| PIK3CA/mTOR | θ = 5.5 | 33 / 9 / 5 | 0.502 | [0.257, 0.625] |
| PIK3CA/mTOR | θ = 6.0 | 18 / 14 / 12 | 0.692 | [0.470, 0.813] |
| PIK3CA/mTOR | θ = 6.5 | 17 / 15 / 12 | 0.674 | [0.438, 0.797] |
| PIK3CA/mTOR | strict 6.5/5.5 | 17 / 7 / 4 | 0.639 | [0.317, 0.792] |
| F2/F10 | θ = 5.5 | 33 / 31 / 31 | 0.353 | [0.216, 0.483] |
| F2/F10 | θ = 6.0 / 6.5 / strict | 31 / 32 / 32 | 0.345 | [0.211, 0.477] |
| JAK1/TYK2 | θ = 5.5 | 32 / 31 / 33 | 0.385 | [0.254, 0.526] |
| JAK1/TYK2 | θ = 6.0 / 6.5 / strict | 31 / 32 / 32 | 0.365 | [0.231, 0.503] |
| JAK1/JAK2 | θ = 5.5 | 34 / 30 / 32 | 0.572 | [0.419, 0.711] |
| JAK1/JAK2 | θ = 6.0 / 6.5 / strict | 32 / 32 / 32 | 0.588 | [0.444, 0.725] |
| PPARG/PPARA | θ = 5.5 | 33 / 30 / 32 | 0.644 | [0.498, 0.751] |
| PPARG/PPARA | θ = 6.0 / 6.5 / strict | 32 / 31 / 32 | 0.649 | [0.504, 0.751] |
| PPARA/PPARD | θ = 5.5 | 34 / 33 / 30 | 0.454 | [0.307, 0.594] |
| PPARA/PPARD | θ = 6.0 / 6.5 / strict | 32 / 32 / 32 | 0.446 | [0.296, 0.584] |

**pChEMBL max versus median (2026-08-26 API snapshot; does not replace Table 2):** EGFR/HER2 label agreement 93.6%, `summary_min` 0.430 frozen vs API-max 0.417 and median 0.424; AChE/BChE agreement 98.9%, median `summary_min` 0.629 (Δ = +0.023); PIK3CA/mTOR agreement 100%, `summary_min` unchanged. Source: `unified_threshold_sensitivity_v2.csv`; `threshold_grid_v1.csv`; `assay_max_vs_median_agreement_v1.csv`.

**High-confidence human SINGLE PROTEIN view (same-day API snapshot; does not replace Table 2):** 253/253 scored ligands with that snapshot match the four-state labels and directional AUROCs are unchanged. The view covers EGFR/HER2, AChE/BChE, and PIK3CA/mTOR. Source: `high_confidence_summary_v1.csv`.

---

## Table S4. Fixed score channel, negative class only

The pocket score is held fixed. Dual ligands are resampled once; selective and neither negatives are resampled independently. Δ = AUROC(dual vs neither) − AUROC(dual vs the matched selective class). Positive Δ means the neither negative looks easier. PIK3CA/mTOR neither n = 4 is underpowered.

| Pair | Score channel | dual vs selective | dual vs neither | Δ | 95% CI | neither underpowered |
|------|----------|---------------:|----------------:|--:|--------|:----------------:|
| EGFR/HER2 | pocket A (vs B-only) | 0.430 | 0.808 | 0.378 | [0.205, 0.547] | no |
| EGFR/HER2 | pocket B (vs A-only) | 0.666 | 0.720 | 0.054 | [−0.157, 0.246] | no |
| AChE/BChE | pocket A | 0.606 | 0.590 | −0.016 | [−0.154, 0.120] | no |
| AChE/BChE | pocket B | 0.650 | 0.709 | 0.058 | [−0.085, 0.196] | no |
| PIK3CA/mTOR | pocket A | 0.692 | 0.472 | −0.220 | [−0.514, 0.019] | yes |
| PIK3CA/mTOR | pocket B | 0.714 | 0.583 | −0.131 | [−0.433, 0.179] | yes |
| F2/F10 | pocket A | 0.345 | 0.508 | 0.163 | [−0.005, 0.339] | no |
| F2/F10 | pocket B | 0.413 | 0.528 | 0.115 | [−0.018, 0.254] | no |
| JAK1/TYK2 | pocket A | 0.365 | 0.809 | 0.444 | [0.263, 0.620] | no |
| JAK1/TYK2 | pocket B | 0.575 | 0.705 | 0.130 | [−0.060, 0.295] | no |
| JAK1/JAK2 | pocket A | 0.728 | 0.723 | −0.004 | [−0.177, 0.170] | no |
| JAK1/JAK2 | pocket B | 0.588 | 0.730 | 0.142 | [−0.058, 0.334] | no |
| PPARG/PPARA | pocket A | 0.706 | 0.759 | 0.053 | [−0.139, 0.228] | no |
| PPARG/PPARA | pocket B | 0.649 | 0.644 | −0.005 | [−0.208, 0.205] | no |
| PPARA/PPARD | pocket A | 0.446 | 0.484 | 0.038 | [−0.139, 0.216] | no |
| PPARA/PPARD | pocket B | 0.646 | 0.665 | 0.019 | [−0.172, 0.198] | no |

Source: `formulation_equal_score_negative_v1.csv`; `equal_score_negative_s34_v1.csv`. Cluster resampling of the two flagship pocket-A Δ values (shared dual draws) is in `equal_score_cluster_bootstrap_v1.csv`: EGFR/HER2 document-cluster [0.083, 0.529] and scaffold-cluster [0.168, 0.562] both exclude 0; JAK1/TYK2 scaffold-cluster [0.234, 0.633] excludes 0, document-cluster [−0.034, 0.682] includes 0. Dual-versus-neither with two-pocket mean scores is main-text Table 3 and does not isolate the negative-class effect.

---

## Table S5. Ligand chemistry: ECFP4 and docking increment

Bemis–Murcko scaffold `GroupKFold` logistic out-of-fold AUROC. The docking column is logistic AUROC from the directional docking score on the same splits; it is **not** the raw ranking AUROC in Table 2. Δ = (ECFP4+docking) − ECFP4. Across 16 arms the largest |Δ| is 0.023.

| Pair | Arm | ECFP4 | ECFP4+docking | Δ | Vina ranking AUROC (Table 2) |
|------|------|------:|--------------:|--:|--------------------------:|
| EGFR/HER2 | D vs A | 0.745 | 0.751 | +0.006 | 0.666 |
| EGFR/HER2 | D vs B | 0.890 | 0.887 | −0.002 | 0.430 |
| AChE/BChE | D vs A | 0.895 | 0.893 | −0.002 | 0.650 |
| AChE/BChE | D vs B | 0.821 | 0.808 | −0.013 | 0.606 |
| PIK3CA/mTOR | D vs A | 0.762 | 0.742 | −0.020 | 0.714 |
| PIK3CA/mTOR | D vs B | 0.889 | 0.898 | +0.009 | 0.692 |
| F2/F10 | D vs A | 0.937 | 0.929 | −0.007 | 0.413 |
| F2/F10 | D vs B | 0.665 | 0.687 | +0.021 | 0.345 |
| JAK1/TYK2 | D vs A | 0.846 | 0.840 | −0.006 | 0.575 |
| JAK1/TYK2 | D vs B | 0.902 | 0.903 | +0.001 | 0.365 |
| JAK1/JAK2 | D vs A | 0.915 | 0.916 | +0.001 | 0.588 |
| JAK1/JAK2 | D vs B | 0.968 | 0.969 | +0.001 | 0.728 |
| PPARG/PPARA | D vs A | 0.833 | 0.820 | −0.013 | 0.649 |
| PPARG/PPARA | D vs B | 0.668 | 0.676 | +0.008 | 0.706 |
| PPARA/PPARD | D vs A | 0.932 | 0.928 | −0.004 | 0.646 |
| PPARA/PPARD | D vs B | 0.858 | 0.835 | −0.023 | 0.446 |

Four single-descriptor `summary_min` values are in main-text Table 2. AChE/BChE TPSA directional AUROCs are 0.733 / 0.801. Source: `incremental_information_v1.csv`; `ecfp4_incremental_s20s24_v1.csv`; `ligand_ml_baseline_scaffold_cv_v1.csv`.

---

## Table S6. Matched minus mismatched pocket

Δ = matched `summary_min` − mismatched `summary_min`. Positive Δ means the correct pocket is higher. On the main panels only EGFR/HER2 and AChE/BChE have 95% CIs excluding 0; all seven scored holdouts include 0.

| Pair | Set | Δ | 95% CI | CI excludes 0 |
|------|------|--:|--------|:---------:|
| EGFR/HER2 | main | 0.170 | [0.060, 0.280] | yes |
| AChE/BChE | main | 0.161 | [0.037, 0.269] | yes |
| PIK3CA/mTOR | main | 0.090 | [−0.122, 0.263] | no |
| F2/F10 | main | −0.031 | [−0.117, 0.040] | no |
| JAK1/TYK2 | main | −0.065 | [−0.152, 0.038] | no |
| JAK1/JAK2 | main | −0.019 | [−0.097, 0.054] | no |
| PPARG/PPARA | main | 0.030 | [−0.081, 0.163] | no |
| PPARA/PPARD | main | 0.012 | [−0.092, 0.145] | no |
| AChE/BChE | holdout | −0.025 | [−0.112, 0.071] | no |
| PIK3CA/mTOR | holdout | −0.023 | [−0.117, 0.079] | no |
| F2/F10 | holdout | −0.079 | [−0.251, 0.075] | no |
| JAK1/TYK2 | holdout | 0.025 | [−0.087, 0.133] | no |
| JAK1/JAK2 | holdout | 0.008 | [−0.073, 0.111] | no |
| PPARG/PPARA | holdout | 0.006 | [−0.168, 0.183] | no |
| PPARA/PPARD | holdout | 0.150 | [−0.056, 0.294] | no |

Source: `wrong_pocket_paired_delta_bootstrap_v1.csv` (`set=main_panel` / `unused_pool_holdout`); `wrong_pocket_by_channel_v1.csv`.

---

## Table S7. Unused-pool holdout

Ligands come from the same ChEMBL harvest after excluding main-panel members, then frozen quotas. EGFR/HER2 leftovers were too thin for a matched holdout. JAK1/JAK2 drew 20 / 20 / 18.

| Pair | Main summary_min [95% CI] | holdout n (D / A / B) | holdout summary_min [95% CI] |
|------|----------------------------:|----------------------:|------------------------------|
| AChE/BChE | 0.606 [0.437, 0.730] | 20 / 20 / 20 | 0.618 [0.422, 0.759] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 20 / 20 / 20 | 0.765 [0.603, 0.891] |
| F2/F10 | 0.345 [0.211, 0.477] | 19 / 20 / 20 | 0.392 [0.214, 0.573] |
| JAK1/TYK2 | 0.365 [0.231, 0.503] | 20 / 20 / 20 | 0.475 [0.282, 0.660] |
| JAK1/JAK2 | 0.588 [0.444, 0.725] | 20 / 20 / 18 | 0.619 [0.420, 0.749] |
| PPARG/PPARA | 0.649 [0.504, 0.751] | 20 / 19 / 20 | 0.535 [0.350, 0.717] |
| PPARA/PPARD | 0.446 [0.296, 0.584] | 20 / 20 / 20 | 0.445 [0.241, 0.559] |

Source: `holdout_pocket_matched_v1.csv`; `table2_comparable_by_channel_v1.csv` (`holdout_vina_20260727`). This is an internal membership sensitivity, not external validation.

---

## Table S8. PIK3CA/mTOR crystal substitution

One pocket at a time: the other pocket keeps frozen main-panel scores. Only this identity-verified pair received the prespecified alternate-crystal docking.

| Replacement | Pocket replaced | D vs A | D vs B | summary_min [95% CI] |
|------|------------|-------:|-------:|----------------------|
| Main 4L23 / 4JT6 | — | 0.714 | 0.692 | 0.692 [0.470, 0.813] |
| PIK3CA → 4JPS | A | 0.714 | 0.486 | 0.486 [0.259, 0.692] |
| PIK3CA → 5DXT | A | 0.714 | 0.505 | 0.505 [0.292, 0.696] |
| mTOR → 4JSX | B | 0.639 | 0.692 | 0.639 [0.418, 0.776] |

Source: `pocket_matched_PM48_alt4JPS_v1.csv`, `..._alt5DXT_v1.csv`, `..._alt4JSX_v1.csv`. Rigid Cα superposition was exploratory and is archived with occupancy snapshots (Note S14).

---

## Table S9. Independent GNINA, alternative rescoring, and five-seed Vina

Independent GNINA searches new poses; it is not a Vina rescore. Scope is EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2. RTMScore / GNINA CNN are best-of-9 rescores of the same Vina poses. Five-seed results do not replace Table 2.

**S9a. Independent GNINA pose generation**

| Pair | Engine | summary_min [95% CI] | Dual vs neither |
|------|------|----------------------|----------------:|
| EGFR/HER2 | Vina primary | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] |
| EGFR/HER2 | GNINA independent | 0.220 [0.109, 0.343] | 0.783 [0.610, 0.922] |
| PIK3CA/mTOR | Vina primary | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] |
| PIK3CA/mTOR | GNINA independent | 0.633 | 0.569 [0.222, 0.889] |
| JAK1/TYK2 | Vina primary | 0.365 [0.231, 0.503] | 0.770 [0.597, 0.906] |
| JAK1/TYK2 | GNINA independent | 0.317 [0.183, 0.463] | 0.705 |

**S9b. PPARG/PPARA same-pose rescoring (primary advantage is unstable)**

| Channel | summary_min [95% CI] | Dual vs neither |
|------|----------------------|----------------:|
| Vina primary | 0.649 [0.504, 0.751] | 0.685 [0.493, 0.848] |
| RTMScore best-of-9 | 0.369 [0.233, 0.475] | 0.817 |
| GNINA CNN affinity | 0.500 [0.356, 0.623] | 0.884 |
| unused-pool holdout | 0.535 [0.350, 0.717] | — |

**S9c. Five-seed Vina `summary_min` range (production 20260727 + 20260811–20260814)**

| Pair | Production seed | Five-seed median | Range |
|------|-------:|------------:|------|
| EGFR/HER2 | 0.430 | 0.373 | 0.321–0.430 |
| AChE/BChE | 0.606 | 0.599 | 0.553–0.606 |
| PIK3CA/mTOR | 0.692 | 0.704 | 0.676–0.726 |
| F2/F10 | 0.345 | 0.366 | 0.345–0.385 |
| JAK1/TYK2 | 0.365 | 0.396 | 0.365–0.399 |
| JAK1/JAK2 | 0.588 | 0.588 | 0.574–0.592 |
| PPARG/PPARA | 0.649 | 0.651 | 0.649–0.691 |
| PPARA/PPARD | 0.446 | 0.454 | 0.446–0.469 |

No five-seed range on F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, or PPARA/PPARD crossed 0.5. Source: `independent_dock_formulation_v1.csv`; `table2_comparable_by_channel_v1.csv`; `multiseed_auroc_aggregate_v2.csv`; `fiveseed_summary_min_aggregate_v1.csv`.

---

## Table S10. Document-cluster uncertainty and detectable-effect simulation

Document-cluster bootstrap is reported only on pairs with complete `document_id` coverage. It does not replace Table 2 ligand-level intervals. PIK3CA/mTOR Dual versus B-only is not stably estimable under document-blocked CV.

| Pair | Contrast | Ligand-level point | Document-cluster 95% CI | n_document groups |
|------|------|------------:|---------------|---------:|
| EGFR/HER2 | D vs A | 0.666 | [0.464, 0.795] | 28 |
| EGFR/HER2 | D vs B | 0.430 | [0.321, 0.617] | 23 |
| AChE/BChE | D vs A | 0.650 | [0.497, 0.809] | 39 |
| AChE/BChE | D vs B | 0.606 | [0.429, 0.753] | 41 |
| PIK3CA/mTOR | D vs A | 0.714 | [0.400, 0.886] | 8 |
| PIK3CA/mTOR | D vs B | 0.692 | [0.000, 0.818] | 9 |

**Fixed-channel Δ cluster resampling:** see the Table S4 continuation; source `equal_score_cluster_bootstrap_v1.csv`. That analysis does not replace the Table S4 ligand-level intervals.

**Detectable effect:** at the observed class sizes, if both true directional AUROCs are 0.70, the probability that the `summary_min` 95% CI excludes 0.5 is 0.219–0.621; if both are 0.60, that probability is 0.025–0.065. Source: `document_cluster_bootstrap_v1.csv`; `detectable_effect_simulation_v1.csv`. A CI that includes 0.5 means the estimate is imprecise; it is not interpreted here as equivalence to chance (see Discussion).

---

## Table S11. BindingDB / PubChem supply counts and the external-docking gate (zero pairs pass)

BindingDB and PubChem were counted for all eight pairs under the same four-state rules used in the main evaluation. External docking further required dropping shared literature sources, duplicate structures, and ECFP4 Tanimoto ≥ 0.70 molecules, plus dual / A-only / B-only each n ≥ 20 with at least three sources per class. No pair was packaged as an external evaluation set or docked externally.

**S11a. BindingDB equal-quantity strict 6.5/5.5 supply counts (zero docking)**

| Pair | Both-end measured | Strict dual / A-only / B-only |
|------|------------------:|------------------------------:|
| PIK3CA/mTOR | 2739 | 1579 / 76 / 96 |
| AChE/BChE | 2711 | 698 / 181 / 92 |
| EGFR/HER2 | 2269 | 1336 / 34 / 31 |
| F2/F10 | 1985 | 376 / 129 / 314 |
| JAK1/TYK2 | 4184 | 2455 / 95 / 165 |
| JAK1/JAK2 | 9761 | 7134 / 131 / 54 |
| PPARG/PPARA | 2026 | 413 / 84 / 95 |
| PPARA/PPARD | 1155 | 253 / 71 / 103 |

**S11b. Remainder after independence filters (external-docking gate)**

| Pair | After filters dual / A-only / B-only | n_sources (D / A / B) | Gate |
|------|------------------------------:|-------------------:|------|
| EGFR/HER2 | 180 / 10 / 20 | 16 / 5 / 4 | fail (A-only n = 10 < 20) |
| AChE/BChE | 4 / 8 / 14 | 2 / 6 / 3 | fail |
| PIK3CA/mTOR | 91 / 4 / 1 | 9 / 2 / 1 | fail |

Source: `crossdb_strict_supply_v1.csv` (S11a; BindingDB `equal_only`); `external_slice_summary_v1.csv` (S11b). S11a is a supply count, not external validation. S11b is the remainder after dropping shared sources, duplicate structures, and high-similarity molecules. The other five pairs were counted in S11a under the same rules and likewise did not form an independence-filtered external set.

---

## Table S12. Literature-year split (primary cutoff 2018)

The test set is earliest `document.year` ≥ 2018. Directional AUROC is reported only if dual, A-only, and B-only each have n ≥ 10. Pairs with reportable year-split counts did not meet that gate on the 2018 test set, so no two-direction time-split test set was packaged.

| Pair | 2018 test n (D / A / B / neither) | Gate |
|------|------------------------------------:|------|
| EGFR/HER2 | 6 / 3 / 14 / 2 | counts only |
| AChE/BChE | 8 / 5 / 15 / 6 | counts only |
| PIK3CA/mTOR | 2 / 0 / 1 / 0 | not evaluable |

Source: `time_split_class_counts_v1.csv`. 2015 / 2020 are prespecified sensitivity cutoffs and were likewise not packaged as external validation.

---

## Table S13. EGFR/HER2 ranking operating points (exploratory)

The mixed library is all 110 evaluation ligands. Top-10 uses `vina_mean`. The AND filter is Dual+A-only+B-only (n = 98) at the dual-median `vina_worst` cutoff.

| Rule | dual | A-only | B-only | neither | precision | selective fraction |
|------|-----:|-------:|-------:|--------:|----------:|----------:|
| Top-10 (`vina_mean`) | 1 | 5 | 4 | 0 | 0.100 | 0.900 |
| AND filter (dual-median `vina_worst`) | 14 | 9 | 24 | — | 0.298 | 0.702 |

Source: `mixed_library_enrichment_v1.csv`; `and_filter_operating_point_v1.csv`. The table shows that current ranking false positives are mainly experimental selectives; it is not a new screening method.

---

## Note S14. Archived to repository / Zenodo (not typeset here)

These files answer questions that do not need a separate typeset SI table. Full CSVs and scripts ship with the GitHub Release and a later Zenodo pack:

- per-ligand docking scores, holdout membership, multi-seed long tables;
- property-caliper 1:1 matching, chemotype hard-negatives, aggregation means (arithmetic / geometric / harmonic), scaffold versus random splits;
- complete-case coverage, measurement frequency, assay-context ledger, J0 candidate-pair θ = 6.0 census;
- historical BindingDB REST counts, leave-cognate-out, PIK3CA occupancy snapshots, contact counts, and whole-chain sequence identity;
- MCL1/Bcl-xL applicability stress test (LC6 pose-gold was not established; not Table 2);
- SHA-256 manifest `REVISION_CHECKSUM_MANIFEST_v1.csv` and evaluation contract `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`.

The one-to-one map from the former 54-table numbering is `SI_TABLE_MERGE_MAP_v1.csv`.
