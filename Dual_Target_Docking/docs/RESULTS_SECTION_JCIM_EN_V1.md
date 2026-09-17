## 3. Results

### 3.1 Paired experimental supply and docking evaluation panels

Paired experimental labels that can support a four-state docking evaluation became progressively more limited as sample requirements increased. Among unordered pairs of human single-component SINGLE PROTEIN targets, 2,164,618 pairs had at least one ligand measured at both targets, and 63,790 pairs had at least 10. Requiring at least 10 dual, A-only, and B-only ligands at \(\theta=6.0\) left 5,253 pairs. Applying the strict 6.5/5.5 activity criterion left 86 pairs. Figure 1C summarizes these counts. The later structure and docking-site gates are defined in section 2.3.

After those supply, structure, and docking-compatibility checks, the primary evaluation comprised eight target pairs: EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, and PPARA/PPARD (Table 1). EGFR/HER2 and PIK3CA/mTOR were sampled from the \(\theta=6.0\) pools; the other six pairs were sampled from the strict 6.5/5.5 pools. Four-state classification requires experimental measurements of the same ligand at both targets, so paired coverage jointly limited the docking panels.

![Figure 1](../figures/jcim_article/Fig1_four_state_and_supply.png)

**Figure 1.** Four-state dual-target evaluation and data supply. (A) Four experimental states defined by threshold \(\theta\); (B) two directional tasks: dual versus A-only uses the target B score, and dual versus B-only uses the target A score, with \(\mathrm{summary}_{\min}=\min[\mathrm{AUROC}_{D/A}(B),\;\mathrm{AUROC}_{D/B}(A)]\); (C) paired experimental coverage decreases under four-state requirements. The primary evaluation comprised eight target pairs selected using the panel-construction and structural criteria in Table 1. That eight-pair set is not a direct continuation of the census counts.

### 3.2 Directional docking performance across virtual-screening control classes

With the same target-specific docking score held fixed, AUROC differences associated with changing the experimental-state control varied across pairs and directions. EGFR/HER2 and JAK1/TYK2 showed the largest fixed-score differences. Using the EGFR-pocket score, the EGFR/HER2 dual-versus-B-only AUROC was 0.324 and rose to 0.786 against neither (difference 0.462 [0.260, 0.641]). Using the JAK1-pocket score, replacing B-only with neither increased the AUROC by 0.444 [0.261, 0.630]. For the remaining pairs, the differences were smaller or their confidence intervals included zero (Figure 2A; Table S4). Cluster resampling is reported in section 3.5.

Under unified \(\theta=6.0\) labels, the two directional AUROCs ranged from 0.345 to 0.728, whereas the descriptive weaker-arm \(\mathrm{summary}_{\min}\) ranged from 0.345 to 0.692. The PPARG/PPARA \(\mathrm{summary}_{\min}\) interval lay entirely above 0.5 (0.649 [0.511, 0.746]), the EGFR/HER2 and F2/F10 intervals lay entirely below 0.5 (0.324 [0.195, 0.474] and 0.345 [0.216, 0.482]), and the remaining 5 pairs crossed 0.5 (Figure 2B; Table 2). These values are pair-specific and are not an eight-pair ranking. Under a binormal simulation that reused these class sizes and the Table 2 bootstrap, the probability that a `summary_min` CI excludes 0.5 was at most 0.068 when the true weaker-arm AUROC was 0.60, and at least 0.817 at 0.75 except for PIK3CA/mTOR (0.453). The simulation is not observed power (Table S4).

**Table 2.** Pocket-matched directional AUROC on the eight primary pairs (Vina; unified \(\theta=6.0\)). Class sizes are n_scored (dual / A-only / B-only). Physicochemical descriptor baselines are in Table S5.

| Pair | n_scored (dual / A-only / B-only) | dual vs A-only (pocket B) [95% CI] | dual vs B-only (pocket A) [95% CI] | summary_min [95% CI] |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|
| EGFR/HER2 | 28 / 37 / 32 | 0.656 [0.520, 0.791] | 0.324 [0.195, 0.474] | 0.324 [0.195, 0.474] |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 [0.448, 0.723] | 0.728 [0.602, 0.842] | 0.588 [0.448, 0.716] |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 [0.434, 0.717] | 0.365 [0.233, 0.508] | 0.365 [0.233, 0.505] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 [0.528, 0.889] | 0.692 [0.491, 0.868] | 0.692 [0.480, 0.802] |
| AChE/BChE | 27 / 26 / 28 | 0.652 [0.492, 0.806] | 0.606 [0.443, 0.749] | 0.606 [0.439, 0.735] |
| F2/F10 | 31 / 32 / 32 | 0.413 [0.271, 0.571] | 0.345 [0.218, 0.491] | 0.345 [0.216, 0.482] |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 [0.512, 0.781] | 0.706 [0.573, 0.830] | 0.649 [0.511, 0.746] |
| PPARA/PPARD | 32 / 32 / 32 | 0.647 [0.509, 0.776] | 0.446 [0.301, 0.592] | 0.446 [0.301, 0.590] |

Under the conventional two-pocket mean ranking used in dual-target virtual screening, dual-versus-neither AUROCs ranged from 0.514 to 0.770. EGFR/HER2 reached 0.759 [0.554, 0.926] and JAK1/TYK2 reached 0.770 [0.613, 0.906], whereas their directional \(\mathrm{summary}_{\min}\) values were 0.324 and 0.365 (Figure 2C; Table 3). The PIK3CA/mTOR neither sample was \(n=4\). Because Table 3 changes both score aggregation and the control class, the influence of the control class was taken from the fixed-pocket analysis in Figure 2A.

![Figure 2](../figures/jcim_article/Fig2_negative_class_formulation.png)

**Figure 2.** Docking performance depends on the experimental-state comparison. (A) \(\Delta\)AUROC after changing the control class at a fixed score channel (dual–neither minus dual–single-target-active); error bars are ligand-level bootstrap 95% confidence intervals; (B) the two directional AUROCs; (C) directional \(\mathrm{summary}_{\min}\) versus dual-versus-neither. The diamond marks the PIK3CA/mTOR neither sample (n = 4). Panel C is not a fixed-score-channel comparison. (D) JAK1/TYK2 top-10% class composition under two-pocket mean Vina ranking (\(k=11\) of 109).

**Table 3.** Same Vina scores under directional versus dual-versus-neither settings (unified \(\theta=6.0\)). Dual-versus-neither uses the two-pocket mean score \(S_{\mathrm{mean}}\). The PIK3CA/mTOR neither sample is small (n = 4).

| Pair | directional summary_min [95% CI] | Dual vs neither (vina_mean) | n_neither |
|------|--------------------------------:|------------------------------:|----------:|
| EGFR/HER2 | 0.324 [0.195, 0.471] | 0.759 [0.557, 0.923] | 12 |
| JAK1/JAK2 | 0.588 [0.444, 0.725] | 0.730 [0.547, 0.875] | 14 |
| JAK1/TYK2 | 0.365 [0.231, 0.503] | 0.770 [0.597, 0.906] | 14 |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 |
| F2/F10 | 0.345 [0.211, 0.477] | 0.519 [0.350, 0.688] | 12 |
| PPARG/PPARA | 0.649 [0.504, 0.751] | 0.685 [0.493, 0.848] | 14 |
| PPARA/PPARD | 0.446 [0.296, 0.584] | 0.565 [0.368, 0.766] | 14 |

That conventional ranking was inspected at a fixed-fraction operating point on the JAK1/TYK2 panel. When all 109 ligands were ranked by the two-pocket mean score, the top 10% (\(k=\lceil 0.10\,n\rceil=11\)) contained 1 dual, 3 A-only, and 7 B-only ligands, and no neither ligand (Figure 2D). That top-10% dual fraction (0.091) is below the panel dual base rate (31/109 = 0.284), so \(\mathrm{EF}_{\mathrm{dual},10\%}=0.320\). The higher dual-versus-neither AUROC therefore did not correspond to fewer highly ranked single-target-active ligands. A two-pocket filter at the median dual \(S_{\mathrm{worst}}\) was applied to Dual+A-only+B-only (\(n=95\)), excluding neither, and retained 16 dual ligands together with 12 A-only and 22 B-only ligands (dual precision 0.32; Figure S1).

The same two-pocket mean ranking was applied at the top 10% of each panel, with \(k=\lceil 0.10\,n\rceil\), so that pairs of different size are compared at the same screening fraction (Table S10). Dual ligands in the top 10% ranged from 1/11 (EGFR/HER2 and JAK1/TYK2; dual fraction 0.091) to 4/5 (PIK3CA/mTOR; 0.800) and 7/11 (PPARG/PPARA; 0.636). Relative to each panel's dual base rate, \(\mathrm{EF}_{\mathrm{dual},10\%}\) ranged from 0.320 (JAK1/TYK2) and 0.354 (EGFR/HER2) to 2.133 (PIK3CA/mTOR) and 2.168 (PPARG/PPARA). Values below 1 mean dual ligands were less common in the top 10% than in the full panel. Neither ligands were absent from the top 10% on five pairs and occupied one slot on AChE/BChE, PPARG/PPARA, and PPARA/PPARD. At least one single-target-active ligand remained in the top 10% on every pair. These counts are descriptive operating points on the retrospective panels; they are not an eight-pair ranking of docking quality. The AND filter on Dual+A-only+B-only is reported in the same table.

### 3.3 Ligand-chemistry baselines and incremental docking discrimination

If dual-target docking ranks experimental classes that already differ in ligand chemistry, a receptor-free baseline can separate the same screening states. On AChE/BChE, TPSA alone gave dual-versus-A-only and dual-versus-B-only AUROCs of 0.742 and 0.801 (Figure 3C; Table S5). On PIK3CA/mTOR, the best single descriptor, heavy-atom count, had \(\mathrm{summary}_{\min}\) 0.463 (Table S5). The difference between Vina and the best single descriptor varied by pair: 5 of eight 95% intervals included 0, whereas JAK1/TYK2, PIK3CA/mTOR, F2/F10 excluded 0 (Table S5). Vina \(\mathrm{summary}_{\min}\) intervals and best-descriptor points are in Figure S2.

Under Bemis–Murcko scaffold-grouped cross-validation, ligand-only ECFP4 models separated several directional comparisons (Figure 3A). Adding the corresponding docking score changed AUROC by at most 0.011 across the 16 directional comparisons, with no consistent direction of change (Figure 3B; Table S5). In this retrospective setting, the docking score did not add a stable ranking increment beyond a ligand-based model of the same experimental classes.

![Figure 3](../figures/jcim_article/Fig3_ligand_chemistry.png)

**Figure 3.** Ligand chemistry as a competing explanation, not a head-to-head predictive benchmark. (A) Rank AUROC from raw Vina scores versus out-of-fold ECFP4 AUROC under scaffold-grouped cross-validation; blue, Vina; orange, ECFP4; (B) AUROC change after adding the corresponding Vina score to ECFP4 on the same splits; lines join the two points for one direction and are not confidence intervals; (C) TPSA distributions on AChE/BChE.

### 3.4 Pocket correspondence and docking-implementation sensitivity

If the directional AUROC is a structure-based docking result, exchanging the two precomputed pocket scores without redocking should reduce the weaker-arm summary. In the main panels, only AChE/BChE had a matched-minus-mismatched \(\mathrm{summary}_{\min}\) 95% interval that excluded zero (0.177 [0.053, 0.291]). EGFR/HER2 was 0.056 [−0.037, 0.155], which includes 0, so this pair is not interpreted as having a matched-pocket advantage. The other six main-panel intervals included zero. All seven available internal-holdout intervals also included zero, with differences from \(-0.079\) to \(+0.150\) (Figure 4A; Table S6). EGFR/HER2 did not have an unused-pool holdout.

![Figure 4](../figures/jcim_article/Fig4_mismatched_pocket.png)

**Figure 4.** Matched-pocket versus mismatched-pocket scores. (A) Matched−mismatched \(\Delta\mathrm{summary}_{\min}\) on the main panels and holdouts; points and bars are estimates and ligand-level bootstrap 95% confidence intervals; (B) \(\mathrm{summary}_{\min}\) on the main panels and holdouts. \(\dagger\) marks the missing EGFR/HER2 unused-pool holdout. The bottom legend distinguishes main panels from holdouts.

Independent GNINA 1.3.2 pose generation used the same receptors, ligands, and docking boxes on three pairs (Figure 5A; Table S7). For EGFR/HER2, the dual-versus-neither AUROC was 0.737 [0.529, 0.919] ($n_{\mathrm{neither}}=11$), whereas dual-versus-B-only was 0.265 [0.144, 0.401]. For JAK1/TYK2, dual-versus-neither was 0.705 [0.524, 0.872] and directional \(\mathrm{summary}_{\min}\) was 0.317 [0.187, 0.455]. For PIK3CA/mTOR, \(\mathrm{summary}_{\min}\) was 0.633, the weaker arm was dual-versus-A-only 0.633 [0.410, 0.769], and dual-versus-neither was 0.569 [0.236, 0.889] (\(n=18/13/12/4\)). These runs are a separate pose-generation pipeline, not a confirmation of the Vina matched-minus-mismatched result.

On PIK3CA/mTOR, replacing PIK3CA 4L23 with 4JPS lowered \(\mathrm{summary}_{\min}\) from 0.692 [0.480, 0.802] to 0.486 [0.264, 0.694]. Replacement with 5DXT gave 0.505 [0.296, 0.713]. Replacing mTOR 4JT6 with 4JSX gave 0.639 [0.435, 0.783] (Figure 5B; Table S7).

Five fixed Vina random seeds produced comparatively limited numerical fluctuation relative to the larger task and receptor effects (Figure 5C). The EGFR/HER2 fixed-score task difference remained positive on all five seeds.

PPARG/PPARA was the only pair whose primary Vina \(\mathrm{summary}_{\min}\) interval lay entirely above 0.5 (0.649 [0.511, 0.746]). Same-pose RTMScore rescoring of the Vina poses lowered it to 0.369 [0.233, 0.475], and GNINA CNN rescoring lowered it to 0.500 [0.356, 0.623] (Table S7).

PIK3CA/mTOR panel-size and exhaustiveness checks are in Figure S3. PM48 is the primary panel (quota n = 48, exhaustiveness = 16); PM110 is a larger protocol-sensitivity panel on the same pair.

Each of the 14 primary receptors produced at least one saved pose below 2.0 Å under the unified chemically mapped CalcRMS table. That result is a search-coverage check. AChE 4EY7 and TYK2 3LXP saved 8 poses; the other primary receptors saved 9. EGFR 3POZ top-1 was 1.019 Å and HER2 3RCD top-1 was 1.947 Å, both below 2 Å. BChE, mTOR, JAK2, PPARG, and PPARA failed the 2 Å top-1 cutoff; PPARA 6LXA top-3 was 7.848 Å (Figure S4; Table S2).

![Figure 5](../figures/jcim_article/Fig5_computational_realization.png)

**Figure 5.** Computational-implementation sensitivity. (A) Independent GNINA versus Vina: filled, Vina; open, GNINA; blue, directional \(\mathrm{summary}_{\min}\); orange, dual-versus-neither; gray lines join the two tasks for one engine and are not confidence intervals. (B) PIK3CA/mTOR receptor substitution; error bars are ligand-level bootstrap 95% confidence intervals. (C) \(\mathrm{summary}_{\min}\) range, median, and production seed across five Vina seeds on the eight pairs; this panel does not show the task difference.

### 3.5 Label and sample-composition sensitivity

Panels drawn from the strict 6.5/5.5 candidate pool kept the same main class composition under several thresholds, so the corresponding AUROCs changed little. Class composition on EGFR/HER2 and PIK3CA/mTOR shifted more with threshold, and their directional estimates also changed (Figure 6A; Table S3). Maximum and median aggregation used the same qualified source records, the same ligand intersection, and the current docking scores; dump-missing ligands were not mixed into that denominator. Class flips: EGFR/HER2 5/109 (primary \(\mathrm{summary}_{\min}\) 0.324); AChE/BChE 1/94 (CHEMBL659; max 0.606 to 0.629); PPARA/PPARD 1/110 (CHEMBL121; \(\mathrm{summary}_{\min}\) remained 0.446). The other pairs kept class composition and \(\mathrm{summary}_{\min}\) point estimates (Table S3).

Unused-pool holdouts built from remaining candidates, after excluding main-panel molecules, showed some dependence on sample composition. AChE/BChE, PIK3CA/mTOR, and JAK1/JAK2 stayed close to the main evaluation. JAK1/TYK2 increased. F2/F10 and PPARA/PPARD remained low. PPARG/PPARA fell from 0.649 to 0.535 [0.360, 0.705] (Figure 4B; Table S6). EGFR/HER2 has no holdout. The ligand-level fixed-score difference on EGFR/HER2 pocket A is 0.462 [0.260, 0.641]. After substituting the corrected-box scores into the frozen scaffold and document groupings, the EGFR/HER2 cluster intervals were [0.235, 0.665] (scaffold) and [0.125, 0.644] (document); both exclude 0. JAK1/TYK2 scaffold-cluster was [0.220, 0.631] (excludes 0) and document-cluster could not be recomputed because the ligand–document map and ChEMBL 37 sqlite are unavailable; a previously deposited interval included 0 and is not treated as a current calculation (Table S4).

### 3.6 Availability of an independent external docking set

BindingDB[16] and PubChem were searched for paired experimental data across all eight target pairs. The \(\theta=6.0\) independence filters removed shared literature sources, duplicate structures, and molecules with ECFP4 Tanimoto similarity \(\geq 0.70\) to the development set. No pair met the independent external-evaluation eligibility criteria of at least 20 dual, A-only, and B-only ligands from at least three independent sources per class. No external docking set was therefore formed; this search is not external validation (Figure 6C,D; Table S8).

![Figure 6](../figures/jcim_article/Fig6_evidence_boundary.png)

**Figure 6.** Evidence boundaries. (A) Activity-threshold sensitivity; (B) ligand-, scaffold-cluster, and document-cluster resampling of the pocket A dual-versus-neither minus dual-versus-B-only difference; (C) BindingDB class counts after filtering, with color saturating at n = 20 per class; (D) independent-source counts after filtering, with color saturating at 3 sources per class. \(\dagger\) marks a class with n < 10.
