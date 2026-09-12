# Results (English working draft)

## 3. Results

### 3.1 Both-end experimental supply and four-state panel construction

Both-end ChEMBL activity that can support four-state evaluation declined rapidly as sample requirements increased. Among pairs with at least one dual-measured ligand, 2,164,618 pairs had at least 1 such ligand, and 63,790 pairs had at least 10. Requiring at least 10 dual, A-only, and B-only ligands at \(\theta=6.0\) left 5,253 pairs. A strict 6.5/5.5 bidirectional selective-supply criterion reduced that number to 86 pairs. This supply change is shown in Figure 1C.

The census describes four-state supply, not the docking menu. After the panel and structure gates in section 2.3, the primary evaluation retained EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, and PPARA/PPARD (Table 1). EGFR/HER2 and PIK3CA/mTOR were drawn from the \(\theta=6.0\) pool; the other six pairs were drawn from the strict 6.5/5.5 pool. Four-state classification requires both-end measurements of the same ligand, so both-end coverage and bidirectional selective supply jointly limited panel size.

![Figure 1](../figures/jcim_article/Fig1_four_state_and_supply.png)

**Figure 1.** Four-state dual-target evaluation and data supply. (A) Four experimental states defined by threshold \(\theta\); (B) two directional tasks: dual versus A-only uses the target B score, and dual versus B-only uses the target A score; (C) ChEMBL supply falls from 2,164,618 pairs to 86 pairs, after which the primary evaluation retains eight pairs.

### 3.2 Experimental-state definition and directional docking evaluation

With the same score channel held fixed, AUROC differences from changing the experimental-state comparison varied by pair and direction. Using the target A (EGFR) score, the EGFR/HER2 dual-versus-B-only AUROC was 0.430 and rose to 0.808 against neither (difference 0.378 [0.205, 0.547]); the other direction was smaller. Using the target A (JAK1) score, the JAK1/TYK2 difference was 0.444 [0.263, 0.620]. Most remaining differences had intervals that crossed 0, so those estimates remain imprecise (Figure 2A; Table S4). Cluster resampling is reported in section 3.5.

Under unified \(\theta=6.0\) labels, directional \(\mathrm{summary}_{\min}\) ranged from 0.345 to 0.692. The PPARG/PPARA interval lay entirely above 0.5 (0.649 [0.504, 0.751]), the F2/F10 interval lay entirely below 0.5 (0.345 [0.211, 0.477]), and the remaining six pairs crossed 0.5 (Figure 2B; Table 2). An interval that crosses 0.5 means the present estimate does not establish performance above or below chance; it does not show equivalence to a random ranking.

**Table 2.** Pocket-matched directional AUROC on the eight primary pairs (Vina; unified \(\theta=6.0\)). Class sizes are n_scored (dual / A-only / B-only). Physicochemical descriptor baselines are in Table S5.

| Pair | n_scored (dual / A-only / B-only) | dual vs A_only (pocket B) | dual vs B_only (pocket A) | summary_min [95% CI] |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 | 0.728 | 0.588 [0.444, 0.725] |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 | 0.365 | 0.365 [0.231, 0.503] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] |
| F2/F10 | 31 / 32 / 32 | 0.413 | 0.345 | 0.345 [0.211, 0.477] |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 | 0.706 | 0.649 [0.504, 0.751] |
| PPARA/PPARD | 32 / 32 / 32 | 0.646 | 0.446 | 0.446 [0.296, 0.584] |

As a descriptive two-pocket mean-score comparison, dual-versus-neither AUROCs were 0.756 [0.562, 0.920] on EGFR/HER2 and 0.770 [0.597, 0.906] on JAK1/TYK2. The directional \(\mathrm{summary}_{\min}\) values were 0.430 and 0.365 (Figure 2C; Table 3). Because that comparison changes both score aggregation and control definition, the isolated effect of control composition was taken from the fixed-score-channel analysis.

![Figure 2](../figures/jcim_article/Fig2_negative_class_formulation.png)

**Figure 2.** Docking performance depends on the experimental-state comparison. (A) \(\Delta\)AUROC after changing the control class at a fixed score channel; error bars are ligand-level bootstrap 95% confidence intervals; (B) the two directional AUROCs; (C) directional \(\mathrm{summary}_{\min}\) versus dual-versus-neither. The diamond marks the PIK3CA/mTOR neither sample (n = 4). Panel C is not a fixed-score-channel comparison.

**Table 3.** Same Vina scores under directional versus dual-versus-neither settings (unified \(\theta=6.0\)). Dual-versus-neither uses the two-pocket mean score \(S_{\mathrm{mean}}\). The PIK3CA/mTOR neither sample is small (n = 4). Dual versus all non-duals is in Table S4.

| Pair | directional summary_min [95% CI] | Dual vs neither (vina_mean) | n_neither |
|------|--------------------------------:|------------------------------:|----------:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 |
| JAK1/JAK2 | 0.588 [0.444, 0.725] | 0.730 [0.547, 0.875] | 14 |
| JAK1/TYK2 | 0.365 [0.231, 0.503] | 0.770 [0.597, 0.906] | 14 |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 |
| F2/F10 | 0.345 [0.211, 0.477] | 0.519 [0.350, 0.688] | 12 |
| PPARG/PPARA | 0.649 [0.504, 0.751] | 0.685 [0.493, 0.848] | 14 |
| PPARA/PPARD | 0.446 [0.296, 0.584] | 0.565 [0.368, 0.766] | 14 |

Panel ranking reflected the screening consequence of this setup. When all 110 EGFR/HER2 ligands were ranked by the two-pocket mean score, the Top-10 contained 1 dual, 5 A-only, and 4 B-only ligands, and no neither ligand; the denominator is 110 (Figure S1A; Table S13). The higher dual-versus-neither AUROC therefore did not correspond to fewer highly ranked single-target selectives. A two-pocket filter at the median dual \(S_{\mathrm{worst}}\) was applied to Dual+A-only+B-only (n = 98), excluding neither, and retained 14 dual ligands together with 9 A-only and 24 B-only ligands (dual precision 0.298; Figure S1B; Table S13).

### 3.3 Ligand-chemistry baselines and directional discrimination

On some pairs, ligand physicochemical features alone separated dual ligands from single-target selectives. On AChE/BChE, TPSA alone gave dual-versus-A-only and dual-versus-B-only AUROCs of 0.733 and 0.801 (Figure 3C; Table S5). Single-descriptor discrimination varied across pairs. On PIK3CA/mTOR, the best single descriptor (heavy-atom count) had \(\mathrm{summary}_{\min}\) 0.463 (Table S5). Across the eight pairs, Vina and the best single descriptor differed in direction and magnitude; six of eight difference 95% intervals included 0, whereas F2/F10 and JAK1/TYK2 excluded 0 (Table S5). Vina \(\mathrm{summary}_{\min}\) intervals and best-descriptor points are in Figure S2.

Under Bemis–Murcko scaffold-grouped cross-validation, receptor-free ECFP4 fingerprints matched or exceeded docking AUROC on several directions (Figure 3A). Adding the corresponding docking score to an ECFP4 logistic model changed AUROC by at most 0.023 across all 16 directions (Figure 3B; Table S5). No stable incremental discrimination was observed under the present panels and cross-validation setup.

![Figure 3](../figures/jcim_article/Fig3_ligand_chemistry.png)

**Figure 3.** Ligand chemistry as a competing explanation. (A) Rank AUROC from raw Vina scores versus out-of-fold ECFP4 AUROC under scaffold-grouped cross-validation; blue, Vina; orange, ECFP4; (B) AUROC change after adding the corresponding Vina score to ECFP4 on the same splits; lines join the two points for one direction and are not confidence intervals; (C) TPSA distributions on AChE/BChE.

### 3.4 Pocket correspondence and computational sensitivity

Figure 4A and Tables S6 and S7 compare matched-pocket versus mismatched-pocket \(\mathrm{summary}_{\min}\). Only EGFR/HER2 (0.170 [0.060, 0.280]) and AChE/BChE (0.161 [0.037, 0.269]) had main-panel difference intervals that excluded 0; the other six pairs crossed 0. On the seven internal holdouts, the intervals all included 0 (differences from −0.079 to +0.150). A matched-pocket advantage was not stably recovered (Figure 4A; Table S7). EGFR/HER2 had too few leftover candidates for an equivalent unused-pool holdout.

![Figure 4](../figures/jcim_article/Fig4_mismatched_pocket.png)

**Figure 4.** Matched-pocket versus mismatched-pocket scores. (A) Matched−mismatched \(\Delta\mathrm{summary}_{\min}\) on the main panels and holdouts; points and bars are estimates and ligand-level bootstrap 95% confidence intervals; (B) \(\mathrm{summary}_{\min}\) on the main panels and holdouts. \(\dagger\) marks the missing EGFR/HER2 unused-pool holdout. The bottom legend distinguishes main panels from holdouts.

After independent GNINA 1.3.2 pose generation and scoring, the EGFR/HER2 dual-versus-neither AUROC was 0.783 [0.610, 0.922] with n_neither = 11. The weaker directional arm dual-versus-B-only was 0.220 [0.109, 0.343] (n_dual = 28, n_B = 32). That interval belongs to dual-versus-B-only scored in pocket A; the independent-GNINA `summary_min` row has no bootstrap interval (Figure 5A; Table S9). JAK1/TYK2 showed the same pattern, with dual-versus-neither 0.705 and directional \(\mathrm{summary}_{\min}\) 0.317 [0.183, 0.463]. Similar differences were still observed under independent pose generation.

On PIK3CA/mTOR, replacing PIK3CA 4L23 with 4JPS lowered \(\mathrm{summary}_{\min}\) from 0.692 [0.470, 0.813] to 0.486 [0.259, 0.692]. Replacement with 5DXT gave 0.505 [0.292, 0.696]. Replacing mTOR 4JT6 with 4JSX gave 0.639 [0.418, 0.776] (Figure 5B; Table S8).

Five fixed Vina random seeds produced comparatively limited numerical fluctuation. Figure 5C shows the complete-case \(\mathrm{summary}_{\min}\) range across seeds for all eight pairs. The EGFR/HER2 task difference was positive on all five Vina seeds (Table S9). On the five added pairs, a fixed-membership intersection (both-end finite scores on every seed) stayed close to the complete-case ranges. That intersection analysis does not cover all eight pairs. PPARG/PPARA switched its weaker arm from dual-versus-A-only to dual-versus-B-only on seeds 20260811, 20260812, and 20260814 (Table S9e). Smaller seed-to-seed changes on a fixed membership do not imply robustness to experimental labels or data sources.

PPARG/PPARA was the only pair whose primary Vina \(\mathrm{summary}_{\min}\) interval lay entirely above 0.5 (0.649 [0.504, 0.751]). Same-pose RTMScore rescoring lowered it to 0.369 [0.233, 0.475], and GNINA CNN rescoring lowered it to 0.500. The unused-pool holdout was 0.535 [0.350, 0.717] (Table S7; Table S9).

PIK3CA/mTOR panel-size and exhaustiveness checks are in Figure S3. PM48 is the primary panel (quota n = 48, exhaustiveness = 16); PM110 is a larger protocol-sensitivity panel on the same pair.

For all 14 primary receptors, the lowest saved-pose RMSD under the method corresponding to each structure remained below 2.0 Å. That result shows search coverage only. AChE 4EY7 and TYK2 3LXP saved 8 poses; the other primary receptors saved 9. EGFR 3POZ is reconstructed QC, with top-1 9.505 Å and lowest saved-pose RMSD 0.760 Å. Among the eight added receptors, JAK2, PPARG, and PPARA failed the 2 Å top-1 cutoff; PPARA 6LXA top-3 was 7.848 Å (Figure S4; Table S2c).

![Figure 5](../figures/jcim_article/Fig5_computational_realization.png)

**Figure 5.** Computational-implementation sensitivity. (A) Independent GNINA versus Vina: filled, Vina; open, GNINA; blue, directional \(\mathrm{summary}_{\min}\); orange, dual-versus-neither; gray lines join the two tasks for one engine and are not confidence intervals. (B) PIK3CA/mTOR receptor substitution; error bars are ligand-level bootstrap 95% confidence intervals. (C) \(\mathrm{summary}_{\min}\) range, median, and production seed across five Vina seeds on the eight pairs; this panel does not show the task difference.

### 3.5 Label and sample-composition sensitivity

Panels drawn from the strict 6.5/5.5 candidate pool kept the same main class composition under several thresholds, so the corresponding AUROCs changed little. Class composition on EGFR/HER2 and PIK3CA/mTOR shifted more with threshold, and their directional estimates also changed (Figure 6A; Table S3). Maximum-versus-median aggregation and the high-confidence field screen are not the same check, and differences between the production cache and another API snapshot cannot all be attributed to aggregation. The same-day API snapshot covers only EGFR/HER2, AChE/BChE, and PIK3CA/mTOR. EGFR/HER2 label agreement was 93.6%; the production `summary_min` was 0.430, the same-day API-max 0.417, and the median 0.424. AChE/BChE API-max remained 0.606, and the median was 0.629. PIK3CA/mTOR label agreement was 100%, and `summary_min` was unchanged. The high-confidence human SINGLE PROTEIN field screen left class assignments and directional results unchanged on those three scored panels (253/253); it is not paper-by-paper reading. The five added pairs have a separate ChEMBL 37 dump-based max-versus-median relabel, which is not that API snapshot: only PPARA/PPARD had one class flip, and all five `summary_min` point estimates were unchanged to three decimals (Table S3).

Unused-pool holdouts built from remaining candidates, after excluding main-panel molecules, showed some dependence on sample composition. AChE/BChE, PIK3CA/mTOR, and JAK1/JAK2 stayed close to the main evaluation. JAK1/TYK2 increased. F2/F10 and PPARA/PPARD remained low. PPARG/PPARA fell from 0.649 to 0.535 [0.350, 0.717] (Figure 4B; Table S7). EGFR/HER2 has no holdout. For the fixed-score difference on target A, EGFR/HER2 scaffold-cluster and document-cluster intervals were [0.168, 0.562] and [0.083, 0.529], both excluding 0. JAK1/TYK2 scaffold-cluster was [0.234, 0.633] (excludes 0) and document-cluster was [−0.034, 0.682] (includes 0) (Figure 6B; Table S10).

### 3.6 Availability of external evaluation data

BindingDB[16] and PubChem were counted for all eight pairs as a strict 6.5/5.5 supply census (Table S11a). Independent-source remainders used \(\theta=6.0\) labels and excluded shared literature, duplicate structures, and highly similar molecules against a development set that included main panels, expanded panels, and internal holdouts (Table S11b). No pair met the independent external-evaluation eligibility criteria of at least 20 compounds and at least 3 independent sources in each of the dual, A-only, and B-only classes, so no external docking was performed (Figure 6C,D; Table S11). In an internal 2018 year split on the already-built panels, JAK1/TYK2 and JAK1/JAK2 met the two-direction sample requirement; the other six pairs did not. That analysis is a sensitivity check on literature-year distribution and is not external validation (Table S12). Unused-pool holdouts and year splits are internal sensitivity analyses, not external validation. Under the data sources, independence filters, and eligibility gates used here, no pair formed an external evaluation set. This is not a claim that no external data exist.

![Figure 6](../figures/jcim_article/Fig6_evidence_boundary.png)

**Figure 6.** Evidence boundaries. (A) Activity-threshold sensitivity; (B) ligand-, scaffold-cluster, and document-cluster resampling of the pocket A dual-versus-neither minus dual-versus-B-only difference; (C) BindingDB class counts after filtering, with color saturating at n = 20 per class; (D) independent-source counts after filtering, with color saturating at 3 sources per class. \(\dagger\) marks a class with n < 10.
