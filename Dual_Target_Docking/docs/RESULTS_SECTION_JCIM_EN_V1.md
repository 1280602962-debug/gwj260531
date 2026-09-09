# Results (English working draft)

## 3. Results

### 3.1 Both-end experimental supply and four-state panel construction

Both-end ChEMBL activity that can support four-state evaluation declined rapidly as sample requirements increased. Among pairs with at least one dual-measured ligand, 2,164,618 pairs had at least 1 such ligand, and 63,790 pairs had at least 10. Requiring at least 10 dual, A-only, and B-only ligands at \(\theta=6.0\) left 5,253 pairs. A strict 6.5/5.5 bidirectional selective-supply criterion reduced that number to 86 pairs. This supply change is shown in Figure 1C.

Among pairs that met the supply criteria, some systems were unsuitable for unified noncovalent docking. Reasons included shared high-throughput assay backgrounds, target class, site chemistry, or lack of suitable human experimental structures. The pairs retained for primary evaluation were EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, and PPARA/PPARD (Table 1). Candidate-pool rules, quotas, and scaffold caps were set separately by supply size and structure feasibility (Table 1). Four-state classification requires experimental measurements of the same ligand at both targets. Both-end coverage and bidirectional selective counts therefore jointly limited the size of strict four-state panels.

### 3.2 Experimental-state definition and directional docking evaluation

With the same score channel held fixed, AUROC differences from changing the experimental-state comparison were not uniform across pairs or directions. EGFR/HER2 and JAK1/TYK2 showed the clearest differences. Using the target A (EGFR) score, the EGFR/HER2 dual-versus-B-only AUROC was 0.430. It rose to 0.808 when the control class was replaced by both-end low-activity neither ligands (difference 0.378 [0.205, 0.547]). The other direction showed a smaller difference. Using the target A (JAK1) score, the JAK1/TYK2 difference was 0.444 [0.263, 0.620]. Scaffold-cluster resampling still excluded 0 on both flagship pairs, whereas the JAK1/TYK2 document-cluster interval included 0 (Figure 2A; Figure 6B; Tables S4 and S10). AChE/BChE, JAK1/JAK2, PPARG/PPARA, and PPARA/PPARD showed smaller same-channel differences whose ligand-level 95% intervals included 0 (Figure 2A; Table S4).

Under unified \(\theta=6.0\) labels, the weaker-arm descriptive summary \(\mathrm{summary}_{\min}\) ranged from 0.345 to 0.692 across the eight pairs. Only PPARG/PPARA had a ligand-level bootstrap 95% interval entirely above 0.5 (0.649 [0.504, 0.751]). The remaining pairs had intervals that included 0.5 or lay entirely below 0.5 (Figure 2B; Table 2).

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

Panel ranking reflected the screening consequence of this setup. When all 110 EGFR/HER2 ligands were ranked by the two-pocket mean score, the Top-10 contained 1 dual, 5 A-only, and 4 B-only ligands, and no neither ligand (Figure S1A; Table S13). The higher dual-versus-neither AUROC therefore did not correspond to fewer highly ranked single-target selectives. A two-pocket filter at the median dual \(S_{\mathrm{worst}}\) retained 14 dual ligands together with 9 A-only and 24 B-only ligands (dual precision 0.298; Figure S1B; Table S13).

### 3.3 Ligand-chemistry baselines and directional discrimination

On some pairs, ligand physicochemical features alone separated dual ligands from single-target selectives. On AChE/BChE, TPSA alone gave dual-versus-A-only and dual-versus-B-only AUROCs of 0.733 and 0.801 (Figure 3C; Table S5). Single-descriptor discrimination varied across pairs. On PIK3CA/mTOR, the best single descriptor (heavy-atom count) had \(\mathrm{summary}_{\min}\) 0.463 (Table S5). Across the eight pairs, Vina and the best single descriptor differed only modestly in \(\mathrm{summary}_{\min}\), and most difference intervals included 0 (Figure S2; Table S5).

Under Bemis–Murcko scaffold-grouped cross-validation, receptor-free ECFP4 fingerprints matched or exceeded docking AUROC on several directions (Figure 3A). Adding the corresponding docking score to an ECFP4 logistic model changed AUROC by at most 0.023 across all 16 directions (Figure 3B; Table S5). No stable incremental discrimination was observed under the present panels and cross-validation setup.

### 3.4 Pocket correspondence and computational sensitivity

To examine whether directional discrimination corresponds to the pocket of the experimental activity difference, Figure 4A and Tables S6 and S7 compare matched-pocket versus mismatched-pocket \(\mathrm{summary}_{\min}\). Only EGFR/HER2 (0.170 [0.060, 0.280]) and AChE/BChE (0.161 [0.037, 0.269]) had main-panel difference intervals that excluded 0. The other six pairs had 95% intervals that included 0. On seven unused holdout sets, matched-versus-mismatched \(\mathrm{summary}_{\min}\) intervals all included 0 (differences from −0.079 to +0.150; Figure 4A; Table S7). A matched-pocket advantage was not stably recovered on the holdouts.

On PIK3CA/mTOR, replacing PIK3CA 4L23 with 4JPS lowered \(\mathrm{summary}_{\min}\) from 0.692 [0.470, 0.813] to 0.486 [0.259, 0.692]. Replacement with 5DXT gave 0.505 [0.292, 0.696]. Replacing mTOR 4JT6 with 4JSX gave 0.639 [0.418, 0.776] (Figure 5B; Table S8).

After independent GNINA 1.3.2 pose generation and scoring, the EGFR/HER2 dual-versus-neither AUROC was 0.783 [0.610, 0.922] with n_neither = 11. The weaker directional arm dual-versus-B-only was 0.220 [0.109, 0.343] (n_dual = 28, n_B = 32). JAK1/TYK2 showed the same pattern, with dual-versus-neither 0.705 and directional \(\mathrm{summary}_{\min}\) 0.317 [0.183, 0.463] (Figure 5A; Table S9). Similar differences were still observed under independent pose generation.

PPARG/PPARA was the only pair whose primary Vina \(\mathrm{summary}_{\min}\) interval lay entirely above 0.5 (0.649 [0.504, 0.751]). Same-pose RTMScore rescoring lowered it to 0.369 [0.233, 0.475], and GNINA CNN rescoring lowered it to 0.500. The unused-pool holdout was 0.535 [0.350, 0.717] (Table S7; Table S9). Five fixed Vina random seeds produced comparatively limited numerical fluctuation (Figure 5C; Table S9). PIK3CA/mTOR panel-size and exhaustiveness checks are in Figure S3. The EGFR/HER2 task difference was positive on all five Vina seeds. Cognate redocking served as protocol quality control. All main receptors produced a saved pose with heavy-atom RMSD < 2.0 Å. For EGFR 3POZ, the top-1 RMSD was 9.505 Å and the lowest saved-pose RMSD was 0.760 Å (Figure S4; Table S2).

### 3.5 Label and sample-composition sensitivity

After activity thresholds were changed (\(\theta=5.5\), 6.0, 6.5, and strict 6.5/5.5), \(\mathrm{summary}_{\min}\) changed little on better-supplied pairs (AChE/BChE remained 0.606). Estimates fluctuated more when one-sided selective counts decreased (Figure 6A; Table S3). Replacing maximum pChEMBL with the median left the main results broadly stable. On panels where high-confidence records could be rechecked, further filtering also left the corresponding directional results unchanged (Table S3).

Unused-pool holdouts built from remaining candidates, after excluding main-panel molecules, showed some dependence on sample composition. AChE/BChE, PIK3CA/mTOR, and JAK1/JAK2 stayed close to the main evaluation. JAK1/TYK2 increased. F2/F10 and PPARA/PPARD remained low. PPARG/PPARA fell from 0.649 to 0.535 [0.350, 0.717] (Figure 4B; Table S7). Results for the seven pairs with unused-pool holdouts are in Figure S5 and Table S7. For EGFR/HER2, ligand, scaffold-cluster, and document-cluster intervals for the fixed-score difference all excluded 0; the JAK1/TYK2 document-cluster interval included 0 (Figure 6B; Table S10).

### 3.6 Availability of external evaluation data

BindingDB and PubChem were counted for all eight pairs under the same four-state rules. Independent-source remainders excluded shared literature, duplicate structures, and highly similar molecules. No pair met the independent external-evaluation eligibility criteria of at least 20 compounds and at least 3 independent sources in each of the dual, A-only, and B-only classes, so no external docking was performed (Figure 6C,D; Table S11). Literature-year splits were likewise limited by bidirectional selective supply. No year-split test set could support independent evaluation of both directions at once (Table S12). Unused-pool holdouts and year splits are internal sensitivity analyses, not external validation. Current public data were insufficient to build independent bidirectional four-state external evaluation sets for these pairs.
