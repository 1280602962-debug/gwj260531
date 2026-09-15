# Supporting Information (English, merged submission draft)

**Article companion:** `MANUSCRIPT_JCIM_EN.md`  
**Numeric rule:** cells below are read from frozen CSVs. Displayed values are three-decimal rounded from the CSV string (half-up or, where already typeset, half-even / banker's). A trailing digit of exactly 5 can therefore appear as either neighbor (for example 0.5045 as 0.504 or 0.505). The CSV is authoritative. EGFR/HER2 Table 2 TPSA is 0.4275 → 0.428.  
**Scope:** typeset Tables S1–S10 answer four questions: how data and docking were done; where the core fixed-channel results come from; whether ligand chemistry and pocket correspondence support attribution; and whether the main conclusions are obviously unstable to labels, sample, or computational realization. Table S10 reports the eight-pair candidate-ranking operating points that correspond to Figure 2D / Figure S1. Historical coordinate-assignment RMSD, five-seed long tables, fixed-membership intersections, class-stratified bootstrap alternatives, sample-size simulations, raw BindingDB/PubChem supply counts, literature-year splits, and per-ligand or exploratory slices remain in the public repository with SHA-256 checksums. They are not repeated as typeset tables.  
**Legacy maps:** former 54-table numbering is `data/manuscript_lock/SI_TABLE_MERGE_MAP_v1.csv`. The immediately previous typeset set was Tables S1–S14; the remap is listed after Table S10.

**Main-text tables (not repeated here):** Table 1 panel composition; Table 2 eight-pair directional AUROC; Table 3 Dual-versus-neither versus directional contrast.

---

## Table S1. Computational settings and statistical definitions

| Item | Value |
|------|------|
| Production Python / RDKit / Meeko | Frozen local env; RDKit 2026.3.1; Meeko 0.7.1 |
| Zero-dock reanalysis | Python 3.12.13; RDKit 2026.3.5; NumPy 2.5.2; pandas 3.0.5; SciPy 1.18.1; scikit-learn 1.9.0 (`requirements-analysis.txt`) |
| AutoDock Vina | 1.2.7; default `vina`; n_modes = 9; energy_range = 3 kcal mol⁻¹ |
| GNINA | 1.3.2; CPU `--no_gpu` |
| RTMScore | `rtmscore_model1` |
| Ligand prep | Desalt (largest organic) → AddHs → ETKDGv3 (seed 20260727) → MMFF ≤200 → Meeko PDBQT |
| Panel / holdout seeds | 20260729 / 20260731 |
| Vina production seed | 20260727; five-seed archive adds 20260811–20260814 |
| Exhaustiveness | PIK3CA/mTOR = 16; all other primary panels and holdouts = 8 |
| Bootstrap | B = 2000. Table 2 reports pointwise directional intervals and a replicate-wise `summary_min` interval from the locked ligand-level non-stratified percentile bootstrap. Dual-versus-neither uses class-stratified percentile intervals; matched-minus-mismatched uses paired bootstrap. A class-stratified `summary_min` sensitivity is archived (`summary_min_stratified_sensitivity_review_v1.csv`) and does not replace Table 2 |
| GroupKFold / logistic | splits = min(5, n_scaffolds, class counts); C = 1.0; max_iter = 4000; default non-shuffled GroupKFold |
| Box | Cognate AABB + 5 Å; minimum edge 20 Å |
| Cognate gate | lowest heavy-atom RMSD among all saved poses < 2.0 Å (search coverage, not top-1 ranking) |

Ligands need both-end scores for directional AUROC; n_scored may be below n_panel (Table 1). Source: `ENV_PIN.md`; `docking_failure_census_v1.csv`.

---

## Table S2. Receptors, docking boxes, and unified cognate-redocking QC

The eight pairs use 14 primary receptor structures (JAK1 6N7A and PPARA 6LXA are shared). PIK3CA/mTOR used E = 16 because 4JT6 failed the gate at E = 8. EGFR 3POZ and HER2 3RCD were redocked under the canonical cognate heavy-atom box. 9V8H is a PPARγ LBD–BRL–PG08-NL ternary complex; the peptide was retained. Figure S4 and the near-native calls use the unified chemically mapped CalcRMS table (S2b) for all 14 slots. AChE 4EY7 and TYK2 3LXP deposited 8 poses, so their lowest RMSD is not a uniform best-of-nine. An earlier coordinate Hungarian assignment is archived in the repository and is not a second official RMSD table.

**S2a. Docking boxes (Å)**

| Protein | PDB | Cognate | Resolution (Å) | center (x, y, z) | size (x, y, z) |
|------|-----|----------|-----------:|------------------|----------------|
| PIK3CA | 4L23 | X6K | 2.50 | 32.443, 45.431, 42.139 | 20.000, 20.000, 20.000 |
| mTOR | 4JT6 | X6K | 3.60 | 51.949, 0.065, −47.707 | 20.332, 20.000, 20.000 |
| AChE | 4EY7 | E20 | 2.35 | −13.988, −43.906, 27.108 | 23.341, 20.000, 20.355 |
| BChE | 4BDS | THA | 2.10 | 133.076, 116.113, 41.335 | 20.000, 20.000, 20.000 |
| EGFR | 3POZ | 03P | 1.50 | 18.680, 32.127, 11.865 | 22.189, 20.000, 22.836 |
| HER2 | 3RCD | 03P | 3.21 | 12.463, 3.371, 27.619 | 23.222, 23.155, 20.000 |
| F2 | 4UDW | N6L | 1.16 | 129.739, −14.402, 164.946 | 20.600, 20.000, 20.000 |
| F10 | 2JKH | BI7 | 1.25 | −3.643, 10.770, 22.771 | 20.000, 20.000, 21.118 |
| JAK1 | 6N7A | KEV | 1.33 | 5.911, −24.317, 21.495 | 20.000, 20.597, 20.000 |
| TYK2 | 3LXP | IZA | 1.65 | −4.498, 25.578, −31.064 | 20.000, 20.000, 20.000 |
| JAK2 | 8BXH | C87 | 1.30 | −12.589, 14.594, 21.263 | 24.069, 20.000, 20.000 |
| PPARG | 9V8H | BRL | 1.39 | −5.454, −5.373, −23.591 | 20.000, 20.140, 20.000 |
| PPARA | 6LXA | EPA | 1.23 | 11.998, 5.742, −7.443 | 20.000, 20.200, 23.432 |
| PPARD | 5U3Q | 7UJ | 1.50 | 40.599, 0.533, 135.353 | 20.223, 20.000, 21.346 |

**S2b. Chemically mapped RMSD for all 14 primary receptors (RDKit CalcRMS)**

No redocking. Atom mapping used one recipe: Meeko topology or SDF/CCD-graph CalcRMS when the prepared ligand maps onto the crystal; PIK3CA/mTOR used graph-automorphism CalcRMS because the prepared ligand is not in the crystal frame. `2JKH/BI7` used CCD SMILES because the OpenBabel SDF had invalid nitrogen valence. EGFR 3POZ is reconstructed QC.

| Protein | PDB | E | top-1 (Å) | top-3 (Å) | lowest saved (Å) | Coverage | top-1 < 2 Å |
|------|-----|--:|----------:|----------:|-----------------:|:--------:|:-----------:|
| EGFR | 3POZ | 8 | 1.019 | 1.019 | 1.019 | pass | yes |
| HER2 | 3RCD | 8 | 1.855 | 1.394 | 1.394 | pass | yes |
| JAK1 | 6N7A | 8 | 0.459 | 0.459 | 0.459 | pass | yes |
| JAK2 | 8BXH | 8 | 10.596 | 0.807 | 0.807 | pass | no |
| TYK2 | 3LXP | 8 | 0.196 | 0.196 | 0.196 | pass | yes |
| PIK3CA | 4L23 | 16 | 0.624 | 0.624 | 0.624 | pass | yes |
| mTOR | 4JT6 | 16 | 7.118 | 0.445 | 0.445 | pass | no |
| AChE | 4EY7 | 8 | 0.339 | 0.339 | 0.339 | pass | yes |
| BChE | 4BDS | 8 | 4.794 | 0.386 | 0.386 | pass | no |
| F2 | 4UDW | 8 | 0.382 | 0.382 | 0.382 | pass | yes |
| F10 | 2JKH | 8 | 0.658 | 0.658 | 0.658 | pass | yes |
| PPARG | 9V8H | 8 | 7.085 | 1.636 | 1.636 | pass | no |
| PPARA | 6LXA | 8 | 7.857 | 7.848 | 1.401 | pass | no |
| PPARD | 5U3Q | 8 | 1.510 | 1.510 | 1.510 | pass | yes |

Source: `all14_cognate_rmsd_calcrrms_v1.csv`; per-pose `all14_cognate_rmsd_calcrrms_modes_v1.csv`. Search coverage passed for all 14 receptors. JAK2, PPARG, PPARA, BChE, mTOR, and EGFR fail the 2 Å top-1 cutoff.

---

## Table S3. Activity-label and aggregation sensitivity

Relabeling on frozen Vina scores. Primary analysis is θ = 6.0 (Table 2). The table lists every threshold that changes class composition on EGFR/HER2 and PIK3CA/mTOR, and the locked θ = 6.0 row for the other six pairs. Those six pairs keep the same main class composition at θ = 6.5 and under the strict 6.5/5.5 rule; full threshold × pair expansions remain in the repository.

| Pair | Label rule | n (D / A / B) | summary_min | 95% CI |
|------|----------|--------------:|------------:|--------|
| EGFR/HER2 | θ = 5.5 | 69 / 22 / 10 | 0.425 | [0.242, 0.626] |
| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 | 0.324 | [0.195, 0.471] |
| EGFR/HER2 | θ = 6.5 | 26 / 29 / 29 | 0.460 | [0.304, 0.609] |
| EGFR/HER2 | strict 6.5/5.5 | 26 / 17 / 7 | 0.324 | [0.138, 0.525] |
| PIK3CA/mTOR | θ = 5.5 | 33 / 9 / 5 | 0.502 | [0.257, 0.625] |
| PIK3CA/mTOR | θ = 6.0 | 18 / 14 / 12 | 0.692 | [0.470, 0.813] |
| PIK3CA/mTOR | θ = 6.5 | 17 / 15 / 12 | 0.674 | [0.438, 0.797] |
| PIK3CA/mTOR | strict 6.5/5.5 | 17 / 7 / 4 | 0.639 | [0.317, 0.792] |
| AChE/BChE | θ = 6.0 | 27 / 26 / 28 | 0.606 | [0.443, 0.735] |
| F2/F10 | θ = 6.0 | 31 / 32 / 32 | 0.345 | [0.211, 0.477] |
| JAK1/TYK2 | θ = 6.0 | 31 / 32 / 32 | 0.365 | [0.231, 0.503] |
| JAK1/JAK2 | θ = 6.0 | 32 / 32 / 32 | 0.588 | [0.444, 0.725] |
| PPARG/PPARA | θ = 6.0 | 32 / 31 / 32 | 0.649 | [0.504, 0.751] |
| PPARA/PPARD | θ = 6.0 | 32 / 32 / 32 | 0.446 | [0.296, 0.584] |

**pChEMBL maximum versus median (ChEMBL 37; frozen Vina scores; does not replace Table 2):** The same ChEMBL 37 records used for Table 2 were re-aggregated by median at θ = 6.0. Maximum pChEMBL matched every scored ligand (0 missing ends; 0 mismatches). Max-to-median class flips: EGFR/HER2 6/110 (primary `summary_min` 0.324); AChE/BChE 1/96 (CHEMBL659; 0.606 → 0.629); PPARA/PPARD 1/110 (CHEMBL121; A-only 32→31; dual-versus-A-only 0.646 → 0.636; `summary_min` remained 0.446). PIK3CA/mTOR and the other four pairs kept class composition and `summary_min` point estimates. Source: `eight_pair_dump_gated_v1/max_vs_median_auroc_v1.csv`; `eight_pair_dump_gated_v1/parity_v1.csv`.

---

## Table S4. Fixed-score control-class comparisons and cluster resampling of the two largest effects

The pocket score is held fixed. Dual ligands are resampled once; selective and neither negatives are resampled independently. Δ = AUROC(dual vs neither) − AUROC(dual vs the matched selective class). Positive Δ means the neither negative looks easier. PIK3CA/mTOR neither n = 4 is underpowered. Dual versus all non-duals, using the two-pocket mean score, is a mixed-library descriptive reference only and is archived with the score tables; it does not support the main-text conclusions.

| Pair | Score channel | dual vs selective | dual vs neither | Δ | 95% CI | neither underpowered |
|------|----------|---------------:|----------------:|--:|--------|:----------------:|
| EGFR/HER2 | pocket A (vs B-only) | 0.324 | 0.786 | 0.462 | [0.262, 0.651] | no |
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

**Cluster resampling of the two pocket-A differences (dual-versus-neither minus dual-versus-B-only).** Ligand-level intervals above do not treat ligands from the same paper or the same Bemis–Murcko scaffold as independent. Cluster intervals do not replace Table 2.

| Pair | Resampling unit | Δ point | 95% CI | CI excludes 0 |
|------|-----------------|--------:|--------|:-------------:|
| EGFR/HER2 | scaffold cluster | 0.378 | [0.168, 0.562] | yes |
| EGFR/HER2 | document cluster | 0.378 | [0.083, 0.529] | yes |
| JAK1/TYK2 | scaffold cluster | 0.444 | [0.234, 0.633] | yes |
| JAK1/TYK2 | document cluster | 0.444 | [−0.034, 0.682] | no |

Source: `formulation_equal_score_negative_v1.csv`; `equal_score_negative_s34_v1.csv`; `equal_score_cluster_bootstrap_v1.csv`. Dual-versus-neither with two-pocket mean scores is main-text Table 3.

---

## Table S5. Ligand-chemistry baselines and incremental docking information

ECFP4 and ECFP4+docking AUROCs are out-of-fold predictions under the same scaffold-grouped cross-validation. The last column is the raw Vina ranking AUROC from Table 2, shown as a descriptive reference. Δ = (ECFP4+docking) − ECFP4. Across 16 arms the largest |Δ| is 0.023. These primary values use unscaled logistic regression.

| Pair | Arm | ECFP4 | ECFP4+docking | Δ | Vina ranking AUROC (Table 2) |
|------|------|------:|--------------:|--:|--------------------------:|
| EGFR/HER2 | D vs A | 0.745 | 0.751 | +0.006 | 0.666 |
| EGFR/HER2 | D vs B | 0.890 | 0.887 | −0.002 | 0.324 |
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

**Paired Δ of Vina `summary_min` minus the best single descriptor.** The four single-descriptor matrices are archived; only the per-pair best descriptor is typeset. AChE/BChE TPSA directional AUROCs are 0.733 / 0.801. Six of eight 95% CIs include 0; F2/F10 and JAK1/TYK2 exclude 0. Figure S2 plots Vina CIs and descriptor points, not the difference CIs.

| Pair | Best descriptor | Descriptor summary_min | Δ | 95% CI | CI excludes 0 |
|------|-----------------|-----------------------:|--:|--------|:-------------:|
| EGFR/HER2 | cLogP | 0.482 | −0.052 | [−0.200, 0.116] | no |
| AChE/BChE | TPSA | 0.733 | −0.128 | [−0.304, 0.049] | no |
| PIK3CA/mTOR | heavy | 0.463 | 0.229 | [−0.011, 0.435] | no |
| F2/F10 | cLogP | 0.515 | −0.170 | [−0.324, −0.012] | yes |
| JAK1/TYK2 | cLogP | 0.580 | −0.215 | [−0.374, −0.001] | yes |
| JAK1/JAK2 | heavy | 0.578 | 0.010 | [−0.084, 0.171] | no |
| PPARG/PPARA | TPSA | 0.627 | 0.022 | [−0.166, 0.193] | no |
| PPARA/PPARD | cLogP | 0.564 | −0.117 | [−0.331, 0.107] | no |

Source: `pocket_matched_vs_best_descriptor_delta_v1.csv`; `descriptor_paired_delta_s19_v1.csv`; `incremental_information_v1.csv`; `ecfp4_incremental_s20s24_v1.csv`; `ligand_ml_baseline_scaffold_cv_v1.csv`.

**Feature-scaling sensitivity.** The same GroupKFold splits were repeated with `StandardScaler` fitted on each training fold only. Across the 16 arms the largest |Δ| was 0.008 (PIK3CA/mTOR D vs A). Scaling does not replace the unscaled 0.023 primary result. Source: `ecfp4_docking_scaler_sensitivity_v1.csv`.

---

## Table S6. Matched versus mismatched pocket and unused-pool holdout

Δ = matched `summary_min` − mismatched `summary_min`. Positive Δ means the weaker matched arm is higher. On the main panels only EGFR/HER2 and AChE/BChE have 95% CIs excluding 0; all seven scored holdouts include 0. An interval that includes 0 does not prove that no advantage exists. EGFR/HER2 has no holdout. `weaker arm switched = yes` means the weaker directional arm differs between matched and mismatched scoring, so Δ`summary_min` cannot represent both directions.

| Pair | Set | Δ | 95% CI | CI excludes 0 | Weaker arm switched |
|------|------|--:|--------|:---------:|:-------------------:|
| EGFR/HER2 | main | 0.170 | [0.060, 0.280] | yes | no |
| AChE/BChE | main | 0.161 | [0.037, 0.269] | yes | yes |
| PIK3CA/mTOR | main | 0.090 | [−0.122, 0.263] | no | no |
| F2/F10 | main | −0.031 | [−0.117, 0.040] | no | no |
| JAK1/TYK2 | main | −0.065 | [−0.152, 0.038] | no | no |
| JAK1/JAK2 | main | −0.019 | [−0.097, 0.054] | no | no |
| PPARG/PPARA | main | 0.030 | [−0.081, 0.163] | no | yes |
| PPARA/PPARD | main | 0.012 | [−0.092, 0.145] | no | yes |
| AChE/BChE | holdout | −0.025 | [−0.112, 0.071] | no | yes |
| PIK3CA/mTOR | holdout | −0.023 | [−0.117, 0.079] | no | yes |
| F2/F10 | holdout | −0.079 | [−0.251, 0.075] | no | yes |
| JAK1/TYK2 | holdout | 0.025 | [−0.087, 0.133] | no | no |
| JAK1/JAK2 | holdout | 0.008 | [−0.073, 0.111] | no | no |
| PPARG/PPARA | holdout | 0.006 | [−0.168, 0.183] | no | yes |
| PPARA/PPARD | holdout | 0.150 | [−0.056, 0.294] | no | no |

Holdout ligands come from the same ChEMBL 37 source after excluding main-panel members, then frozen quotas. JAK1/JAK2 drew 20 / 20 / 18. Holdout panels are directional (dual / A-only / B-only) and do not include a neither class. This is an internal membership sensitivity, not external validation.

| Pair | Main summary_min [95% CI] | holdout n (D / A / B) | holdout summary_min [95% CI] |
|------|----------------------------:|----------------------:|------------------------------|
| AChE/BChE | 0.606 [0.437, 0.730] | 20 / 20 / 20 | 0.618 [0.422, 0.759] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 20 / 20 / 20 | 0.765 [0.603, 0.891] |
| F2/F10 | 0.345 [0.211, 0.477] | 19 / 20 / 20 | 0.392 [0.214, 0.573] |
| JAK1/TYK2 | 0.365 [0.231, 0.503] | 20 / 20 / 20 | 0.475 [0.282, 0.660] |
| JAK1/JAK2 | 0.588 [0.444, 0.725] | 20 / 20 / 18 | 0.619 [0.420, 0.749] |
| PPARG/PPARA | 0.649 [0.504, 0.751] | 20 / 19 / 20 | 0.535 [0.350, 0.717] |
| PPARA/PPARD | 0.446 [0.296, 0.584] | 20 / 20 / 20 | 0.445 [0.241, 0.559] |

Source: `wrong_pocket_paired_delta_bootstrap_v1.csv` (`set=main_panel` / `unused_pool_holdout`); `wrong_pocket_by_channel_v1.csv`; `pocket_unidirectional_delta_v1.csv`; `holdout_pocket_matched_v1.csv`; `table2_comparable_by_channel_v1.csv` (`holdout_vina_20260727`).

---

## Table S7. Computational realization: receptor substitution, independent GNINA, and PPARG rescoring

Independent GNINA searches new poses; it is not a Vina rescore. Scope is EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2. Independent GNINA did not return both-end scores for every ligand (EGFR/HER2 EH120_109; PIK3CA/mTOR PM48_19; JAK1/TYK2 one dual and three B-only). EGFR/HER2 dual-versus-neither therefore uses n_neither = 11 versus 12 in the primary Vina Table 3. EGFR/HER2 and PIK3CA/mTOR independent-GNINA `summary_min` rows have empty CI columns in the source file; intervals below are labeled on the corresponding single arm and are not min-of-two bootstrap intervals. Five-seed Vina ranges, fixed-membership intersections, and complete-case seed tables are archived in the repository. Across those seeds the numerical values fluctuated, but the main task and pair-level patterns did not change. The EGFR/HER2 task difference was positive on all five Vina seeds.

**S7a. Independent GNINA pose generation**

| Pair | Engine | n_dual / n_A / n_B / n_neither | summary_min | Weaker-arm AUROC [95% CI] | Dual vs neither |
|------|------|------|------------:|---------------------------|----------------:|
| EGFR/HER2 | Vina primary | 28 / 38 / 32 / 12 | 0.324 [0.195, 0.471] | dual–B-only (pocket A) 0.324 [0.188, 0.464] | 0.759 [0.557, 0.923] |
| EGFR/HER2 | GNINA independent | 28 / 38 / 32 / 11 | 0.220 | dual–B-only (pocket A) 0.220 [0.109, 0.343] | 0.783 [0.610, 0.922] |
| PIK3CA/mTOR | Vina primary | 18 / 14 / 12 / 4 | 0.692 [0.470, 0.813] | dual–B-only (pocket A) 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] |
| PIK3CA/mTOR | GNINA independent | 18 / 13 / 12 / 4 | 0.633 | dual–A-only (pocket B) 0.633 [0.427, 0.825] | 0.569 [0.222, 0.889] |
| JAK1/TYK2 | Vina primary | 31 / 32 / 32 / 14 | 0.365 [0.231, 0.503] | dual–B-only (pocket A) 0.365 [0.231, 0.503] | 0.770 [0.597, 0.906] |
| JAK1/TYK2 | GNINA independent | 30 / 32 / 29 / 14 | 0.317 [0.183, 0.463] | dual–B-only (pocket A) 0.317 [0.183, 0.463] | 0.705 [0.517, 0.876] |

**S7b. PPARG/PPARA same-pose rescoring (primary advantage is unstable)**

| Channel | summary_min [95% CI] | Dual vs neither |
|------|----------------------|----------------:|
| Vina primary | 0.649 [0.504, 0.751] | 0.685 [0.493, 0.848] |
| RTMScore (all saved poses) | 0.369 [0.233, 0.475] | 0.817 |
| GNINA CNN affinity | 0.500 [0.356, 0.623] | 0.884 |
| unused-pool holdout | 0.535 [0.350, 0.717] | — |

**S7c. PIK3CA/mTOR crystal substitution**

One pocket at a time: the other pocket keeps frozen main-panel scores. Only this identity-verified pair received the prespecified alternate-crystal docking.

| Replacement | Pocket replaced | D vs A | D vs B | summary_min [95% CI] |
|------|------------|-------:|-------:|----------------------|
| Main 4L23 / 4JT6 | — | 0.714 | 0.692 | 0.692 [0.470, 0.813] |
| PIK3CA → 4JPS | A | 0.714 | 0.486 | 0.486 [0.259, 0.692] |
| PIK3CA → 5DXT | A | 0.714 | 0.505 | 0.505 [0.292, 0.696] |
| mTOR → 4JSX | B | 0.639 | 0.692 | 0.639 [0.418, 0.776] |

Source: `independent_dock_formulation_v1.csv`; `table2_comparable_by_channel_v1.csv`; `pocket_matched_PM48_alt4JPS_v1.csv`, `..._alt5DXT_v1.csv`, `..._alt4JSX_v1.csv`. Rigid Cα superposition was exploratory and is archived in the repository.

---

## Table S8. External-data eligibility after independence filters

BindingDB and PubChem were searched for all eight pairs. The typeset table and Figure 6C,D are the BindingDB independence-filtered remainder labeled at \(\theta=6.0\), not a raw supply census. PubChem served as an additional paired-data availability check and is not merged into these counts. External docking further required dropping shared literature sources, duplicate structures, and ECFP4 Tanimoto ≥ 0.70 molecules, plus dual / A-only / B-only each n ≥ 20 with at least three sources per class. The development-molecule set includes main panels, the expanded PIK3CA/mTOR PM110 panel, and internal holdouts. No pair met the independent external-evaluation eligibility criteria, so no external docking was performed. Raw BindingDB/PubChem supply counts and a publication-year subset on the already-built panels remain in the repository; neither is treated as external validation.

| Pair | After filters dual / A-only / B-only | n_sources (D / A / B) | Gate |
|------|------------------------------:|-------------------:|------|
| EGFR/HER2 | 180 / 10 / 20 | 16 / 5 / 4 | fail (A-only n = 10 < 20) |
| AChE/BChE | 4 / 8 / 14 | 2 / 6 / 3 | fail |
| PIK3CA/mTOR | 91 / 4 / 1 | 9 / 2 / 1 | fail |
| F2/F10 | 46 / 15 / 16 | 4 / 1 / 3 | fail (A/B n = 15/16; A-only has 1 source) |
| JAK1/TYK2 | 323 / 7 / 103 | 20 / 3 / 6 | fail (A-only n = 7) |
| JAK1/JAK2 | 928 / 40 / 14 | 28 / 7 / 4 | fail (B-only n = 14) |
| PPARG/PPARA | 0 / 0 / 1 | 0 / 0 / 1 | fail |
| PPARA/PPARD | 0 / 1 / 0 | 0 / 1 / 0 | fail |

Source: `external_slice_summary_v1.csv`. This table is an eligibility screen, not external validation.

---

## Table S9. Final target-pair eligibility audit

Seven pairs reached the G5 protocol-compatibility gate and were included because each has conventional noncovalent pockets representable under the common rigid-receptor Vina protocol: PIK3CA/mTOR, AChE/BChE, F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, and PPARA/PPARD. The table lists the supply-limited EGFR/HER2 exception and the pairs that failed the final structure- or protocol-compatibility gates. It is an audit of pair selection, not a docking-performance table. Source: `pair_eligibility_audit_s14_v1.csv`; `FEASIBLE_PAIR_LADDER_V1.md`; `TIER1_DOCKING_ROSTER_V1.md`; `pair_ligand_identity_qc_v1.csv`.

| Pair | Last gate reached | Included/excluded | Reason | Evidence used |
|------|-------------------|-------------------|--------|---------------|
| EGFR/HER2 | supply-limited exception | included | Did not meet the strict 6.5/5.5 selective-supply criterion (minimum selective-class count = 7) but had suitable human holo structures, a cognate-defined docking site, and sufficient dual, A-only, and B-only ligands at \(\theta=6.0\) for directional evaluation. | `TIER1_DOCKING_ROSTER_V1.md`; Table 1 |
| CTSK/CTSS | G4 ligand identity | excluded | Both holos are reversible-covalent cysteine-protease complexes and therefore require treatment outside the common noncovalent rigid-receptor Vina protocol. | `TIER1_DOCKING_ROSTER_V1.md`; `FEASIBLE_PAIR_LADDER_V1.md` |
| CREBBP/BRD4 | G4 ligand identity | excluded | CREBBP has both a HAT catalytic site and a bromodomain; the intended docking domain is not uniquely defined under the common protocol. | `TIER1_DOCKING_ROSTER_V1.md`; `FEASIBLE_PAIR_LADDER_V1.md` |
| F2/PRSS1 | G4 ligand identity | excluded | Trypsin (PRSS1) is a pharmacological antitarget rather than a designed dual-target partner under the common protocol. | `TIER1_DOCKING_ROSTER_V1.md`; `FEASIBLE_PAIR_LADDER_V1.md` |
| CNR1/CNR2 | G4 ligand identity | excluded | Membrane GPCR pair requiring construct and conformational-state treatment outside the common soluble rigid-receptor Vina protocol. | `TIER1_DOCKING_ROSTER_V1.md`; `FEASIBLE_PAIR_LADDER_V1.md` |
| HCRTR1/HCRTR2 | G4 ligand identity | excluded | Membrane GPCR pair requiring construct and conformational-state treatment outside the common soluble rigid-receptor Vina protocol. | `TIER1_DOCKING_ROSTER_V1.md`; `FEASIBLE_PAIR_LADDER_V1.md` |
| OPRM1/OPRD1 | G4 ligand identity | excluded | Membrane GPCR pair requiring construct and conformational-state treatment outside the common soluble rigid-receptor Vina protocol. | `TIER1_DOCKING_ROSTER_V1.md`; `FEASIBLE_PAIR_LADDER_V1.md` |
| OPRD1/OPRK1 | G4 ligand identity | excluded | Membrane GPCR pair requiring construct and conformational-state treatment outside the common soluble rigid-receptor Vina protocol. | `TIER1_DOCKING_ROSTER_V1.md`; `FEASIBLE_PAIR_LADDER_V1.md` |
| S1PR3/S1PR1 | G4 ligand identity | excluded | Membrane GPCR pair requiring construct and conformational-state treatment outside the common soluble rigid-receptor Vina protocol. | `TIER1_DOCKING_ROSTER_V1.md`; `FEASIBLE_PAIR_LADDER_V1.md` |
| SLC6A4/SLC6A3 | G4 ligand identity | excluded | Membrane SLC6 transporter pair requiring treatment outside the common soluble rigid-receptor Vina protocol. | `TIER1_DOCKING_ROSTER_V1.md`; `FEASIBLE_PAIR_LADDER_V1.md` |
| SLC6A2/SLC6A4 | G4 ligand identity | excluded | Membrane SLC6 transporter pair requiring treatment outside the common soluble rigid-receptor Vina protocol. | `TIER1_DOCKING_ROSTER_V1.md`; `FEASIBLE_PAIR_LADDER_V1.md` |
| OPRM1/OPRK1 | G3 human holo supply | excluded | Drug-like small-molecule filter reduced the minimum strict selective-class count from 56 to 46; the pair therefore failed the G4 ligand-identity gate. | `pair_ligand_identity_qc_v1.csv`; `FEASIBLE_PAIR_LADDER_V1.md` |
| JAK3/TYK2 | G3 human holo supply | excluded | Drug-like small-molecule filter reduced the minimum strict selective-class count from 51 to 48; the pair therefore failed the G4 ligand-identity gate. | `pair_ligand_identity_qc_v1.csv`; `FEASIBLE_PAIR_LADDER_V1.md` |

---

## Table S10. Eight-pair candidate-ranking operating points

All eight primary panels were ranked by the two-pocket mean Vina score \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\), with ties broken by ascending ligand ID. The ranking readout is the top 10% of the full four-state panel, including neither, with \(k=\lceil 0.10\,n\rceil\). That fixed screening fraction, not a fixed count, is what makes the eight pairs comparable. The top-10% dual fraction is \(\mathrm{dual}/k\). Enrichment versus the panel dual base rate is \(\mathrm{EF}_{\mathrm{dual},10\%}=(\mathrm{dual}/k)/(n_{\mathrm{dual}}/n)\); values below 1 mean dual ligands were less common in the top 10% than in the full panel. The AND filter uses Dual+A-only+B-only, excludes neither, and retains ligands with \(S_{\mathrm{worst}}\geq\) the median dual \(S_{\mathrm{worst}}\). These rows are descriptive operating points; they are not an eight-pair ranking of docking quality. JAK1/TYK2 top 10% is Figure 2D; the JAK1/TYK2 AND filter is Figure S1. Source: `eight_pair_ranking_operating_point_v1.csv`; `review_scored_membership_v1.csv`.

| Pair | n_ranked (D / A / B / N) | k | Top 10% D / A / B / N | dual / k | n_dual / n | EF_dual,10% | AND input → pass (D / A / B) | AND dual precision |
|------|-------------------------:|--:|----------------------:|---------:|-----------:|------------:|------------------------------:|-------------------:|
| EGFR/HER2 | 110 (28 / 38 / 32 / 12) | 11 | 1 / 5 / 5 / 0 | 0.091 | 0.255 | 0.357 | 98 → 47 (14 / 9 / 24) | 0.298 |
| JAK1/JAK2 | 110 (32 / 32 / 32 / 14) | 11 | 6 / 5 / 0 / 0 | 0.545 | 0.291 | 1.875 | 96 → 35 (16 / 13 / 6) | 0.457 |
| JAK1/TYK2 | 109 (31 / 32 / 32 / 14) | 11 | 1 / 3 / 7 / 0 | 0.091 | 0.284 | 0.320 | 95 → 50 (16 / 12 / 22) | 0.320 |
| PIK3CA/mTOR | 48 (18 / 14 / 12 / 4) | 5 | 4 / 1 / 0 / 0 | 0.800 | 0.375 | 2.133 | 44 → 17 (9 / 4 / 4) | 0.529 |
| AChE/BChE | 95 (27 / 25 / 28 / 15) | 10 | 5 / 3 / 1 / 1 | 0.500 | 0.284 | 1.759 | 80 → 32 (14 / 7 / 11) | 0.438 |
| F2/F10 | 107 (31 / 32 / 32 / 12) | 11 | 4 / 1 / 6 / 0 | 0.364 | 0.290 | 1.255 | 95 → 59 (16 / 20 / 23) | 0.271 |
| PPARG/PPARA | 109 (32 / 31 / 32 / 14) | 11 | 7 / 3 / 0 / 1 | 0.636 | 0.294 | 2.168 | 95 → 31 (16 / 9 / 6) | 0.516 |
| PPARA/PPARD | 110 (32 / 32 / 32 / 14) | 11 | 5 / 2 / 3 / 1 | 0.455 | 0.291 | 1.562 | 96 → 58 (16 / 22 / 20) | 0.276 |

---

## Remap from the previous typeset S1–S14 set

| Previous | Typeset now | Repository only |
|----------|-------------|-----------------|
| S1 settings | S1 (compressed) | — |
| S2a boxes + S2c CalcRMS | S2 | former S2b Hungarian RMSD |
| S3 threshold / aggregation | S3 (changed pairs + locked θ = 6.0 rows) | full threshold × pair grid |
| S4 fixed-score channel | S4 + flagship cluster resampling | dual-versus-all-nonduals side table |
| S5 ligand chemistry | S5 (16 ECFP4 arms + best descriptor) | four single-descriptor matrices |
| S6a + S7 holdout | S6 | former S6b unidirectional split |
| S8 receptor + S9a/S9b | S7 | S9c–S9e five-seed / intersection tables |
| S10 cluster (flagship Δ only) | S4 | class-stratified bootstrap; sample-size scenario; Figure S5 simulation |
| S11b independence remainder | S8 | S11a raw supply counts |
| S12 literature-year split | — | year-split CSVs |
| S13 EGFR/HER2 operating points | Table S10 (eight pairs) | Figure 2D / Figure S1 illustrate JAK1/TYK2; former three-pair mixed-library CSV remains archived |
| S14 pair audit | S9 (excluded pairs + EGFR/HER2 exception) | repeated included-pair rows |

Archived files that answer questions not typeset here include: per-ligand docking scores and multi-seed long tables; property-caliper matching; complete-case coverage; historical BindingDB supply counts; leave-cognate-out and occupancy snapshots; MCL1/Bcl-xL applicability stress test; SHA-256 manifest `REVISION_CHECKSUM_MANIFEST_v1.csv`; and evaluation contract `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`. A dated ChEMBL API snapshot and a three-pair high-confidence field screen remain in the repository and are not typeset.
