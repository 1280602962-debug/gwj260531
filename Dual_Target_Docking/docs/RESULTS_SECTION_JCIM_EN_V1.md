# Results (JCIM Articles draft, English)

## 3. Results

### 3.1 Data supply and four-state panels

Public bioactivity resources constrain how strictly a dual-target docking evaluation can be built. In a frozen ChEMBL supply audit of 49 candidate pairs, ligands that are potent on one target and experimentally weak on the other—directional selective hard negatives—were scarce under a strict 6.5/5.5 rule: only four pairs met a thick-panel gate of ≥50 strict hard negatives on both ends. After excluding metal-dependent HDAC1/HDAC6, PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB remained relatively well supplied, whereas EGFR/HER2 retained only seven strict B-only ligands and was kept as a supply-limited case (Table 1). A zero-docking BindingDB/PubChem recount supported the same scarcity conclusion (Table S12).

Primary docking labels used a unified θ = 6.0 four-state assignment (dual, A-only, B-only, neither), frozen before scores were inspected; the stricter rule served only supply qualification and panel construction (Methods 2.1–2.2). A post-hoc θ = 6.0 label census found 17 unique pairs with directional classes each n ≥ 10 (Table S44), but those counts are supply diagnostics only: docking evaluation remains the original four pairs.

### 3.2 Four-pair primary directional AUROCs

On the frozen K = 4 panels, AutoDock Vina mode-1 scores were evaluated with the two primary pocket-matched directional AUROCs—dual versus A-only in pocket B and dual versus B-only in pocket A—with `summary_min` reported only as a conservative descriptive summary (Figure 1B; Methods 2.4). EGFR/HER2, AChE/BChE, PIK3CA/PIK3CB, and PIK3CA/mTOR gave `summary_min` values of 0.430, 0.606, 0.500, and 0.692, respectively (Table 2; Figure 4A). All four ligand-bootstrap 95% intervals for `summary_min` included 0.5. Pair ranking was unchanged under arithmetic, geometric, and harmonic aggregation (Table S26).

As a prespecified formulation contrast on the same scores, Dual versus neither used experimental inactives (`vina_mean`; Table 3; Figure 3). On EGFR/HER2, Dual versus neither yielded AUROC 0.756 [0.562, 0.920] (n_neg = 12), whereas directional `summary_min` remained 0.430 [0.282, 0.578]; Dual versus all non-duals fell to 0.551. In a mixed-library ranking of all 110 EGFR/HER2 ligands by `vina_mean`, the Top-10 contained one dual and nine experimental selectives (hard-negative fraction 0.90; Table S25). Holding the pocket score fixed, Dual versus B-only versus Dual versus neither on EGFR/HER2 pocket A differed by 0.378 [0.205, 0.547], indicating that the weak arm reflects negative-class composition rather than mean aggregation alone (Table S34). AChE/BChE and PIK3CA/PIK3CB showed only small Dual-versus-neither increments whose intervals overlap the directional arms. PIK3CA/mTOR Dual versus neither is underpowered (neither n = 4).

Independent GNINA 1.3.2 pose generation on the same frozen EGFR/HER2 ligands, receptors, and boxes left the formulation gap intact: Dual versus neither 0.783 [0.610, 0.922] versus directional `summary_min` 0.220 [0.109, 0.343], with a Top-10 again dominated by selectives (Table S32). This check asks whether the formulation effect persists when pose search is changed; it is not an engine bake-off.

The five-seed Vina sensitivity preserved the sign of the Dual-versus-neither minus `summary_min` gap for every pair (Table S54). Median `summary_min` (range) was 0.373 (0.321–0.430) for EGFR/HER2, 0.599 (0.553–0.606) for AChE/BChE, 0.478 (0.468–0.502) for PIK3CA/PIK3CB, and 0.704 (0.676–0.726) for PIK3CA/mTOR. EGFR/HER2 retained the largest positive gap at every seed (range 0.334–0.442). These repeats support the qualitative seed stability of the formulation contrast within the frozen panels; they do not replace the production-seed estimates in Tables 2–3.

**Table 2.** Pocket-matched directional AUROC on the frozen K = 4 set (Vina; unified θ = 6.0), with all four prespecified descriptor `summary_min` values. Class sizes are n_scored (dual / A-only / B-only). The highest descriptor is a best single-descriptor reference.

| Pair | n_scored (dual / A-only / B-only) | dual vs A_only (pocket B) | dual vs B_only (pocket A) | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.350, 0.650] | 0.622 | 0.620 | 0.595 | 0.418 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] | 0.463 | 0.448 | 0.310 | 0.260 |

**Table 3.** Same Vina scores under Dual-versus-neither versus directional formulations (unified θ = 6.0). Dual-versus-neither uses experimental inactives (`vina_mean`). PIK3CA/mTOR Dual versus neither is underpowered (n_neg = 4).

| Pair | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|---|---:|---:|---:|---:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/PIK3CB | 0.500 [0.350, 0.650] | 0.559 [0.373, 0.746] | 16 | 0.556 [0.437, 0.672] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |

A detectable-effect simulation at the observed class sizes shows that these panels resolve large directional effects more readily than moderate ones (Table S31). Failure of a `summary_min` CI to exclude 0.5 therefore does not establish equivalence to chance.

### 3.3 Chemistry and source confounding

Docking was compared with four prespecified physicochemical descriptors and with ECFP4 under Bemis–Murcko scaffold GroupKFold (Figure 4B–C; Tables 2, S19–S20, S24). Relative to the best single-descriptor reference, paired `summary_min` differences included zero on all four pairs. On AChE/BChE, mean TPSA differed sharply between dual and selective classes, and TPSA alone exceeded Vina on the corresponding contrast; adding heavy-atom count and TPSA raised dual-versus-B-only AUROC from 0.606 to 0.807 while the docking odds ratio remained near one. ECFP4 scaffold-grouped fold AUROCs on several arms were well above the docking contrasts (for example ≈0.89 versus 0.43 for EGFR/HER2 dual-versus-B-only). Under the present scaffold-grouped task, adding the pocket-matched docking score to ECFP4 changed AUROC by at most 0.020 in absolute value (Table S24). That incremental result is limited to these labels, series, and receptors: it does not show that docking lacks structural information in general, only that it did not add a stable large increment beyond 2D chemistry here.

Document and scaffold correlation matter for uncertainty. Document-blocked CV left the EGFR/HER2 weak arm at 0.430 (document-cluster bootstrap 95% CI [0.321, 0.617]; Table S39). Scaffold-cluster bootstrap on the same arm gave [0.278, 0.595] (Table S39 companion; `scaffold_cluster_bootstrap_v1.csv`). All three estimators—ligand, document-cluster, and scaffold-cluster—keep the weak arm near chance with intervals that span 0.5. On PIK3CA/mTOR, all four neither ligands and their retained records came from one ChEMBL document (Table S37), and document-blocked Dual versus B-only was not stably estimable; document-cluster CI lower bounds for that pair’s B arm reach 0.0, whereas scaffold-cluster intervals remain wide but finite (Table S40; cluster uncertainty verdict). Complete-case dual-measured fractions on the maps were only 14.5%–34.0% across pairs (Table S37), so the panels enrich jointly profiled chemistry.

### 3.4 Receptor realization and docking-failure sensitivity

Holding one pocket frozen and replacing the other with alternate crystals that passed cognate QC changed apparent discrimination in opposite directions (Figure 5; Table S30). On PIK3CA/mTOR, replacing PIK3CA 4L23 with 4JPS or 5DXT while holding mTOR at 4JT6 dropped `summary_min` from 0.692 to 0.486 [0.259, 0.692] and 0.505 [0.292, 0.696]. On PIK3CA/PIK3CB, the same PIK3CA crystals with 2WXF frozen raised `summary_min` from 0.500 to 0.691 and 0.685. Receptor realization is therefore a sensitivity, not evidence of structural robustness.

Docking failures were concentrated among large or flexible ligands (Table S27). On AChE/BChE, rank-extreme lower bounds remained directionally consistent with complete-case estimates; on PIK3CA/PIK3CB, using the available pocket score for one failed A-only ligand left `summary_min` at 0.500. Unused-pool holdout and exhaustiveness/PM110 checks are internal sensitivities (Tables S8, S16); the holdout also exposed a wrong-pocket reversal relative to the main panels (Table S17), an unresolved out-of-panel failure mode rather than a robustness claim.

### 3.5 External supply failure and evidence boundary

A BindingDB-native 202608 archive rebuild under a contract frozen before docking applied literature, structure, and ECFP4 < 0.70 filters and yielded **zero pairs** meeting the pre-frozen primary external gate; remaining counts are upper bounds because ChEMBL document lookup was incomplete, and the slice was not docked (Tables S48–S49; Figure S8).[16] The pre-frozen 2018 literature-year split likewise failed the sample gate on the primary cutoff and is not packaged as external validation (Table S41). The manuscript therefore remains a data-constrained four-pair formulation audit with an explicit failed external-supply audit.

MCL1/Bcl-xL was formally demoted: LC6 topology-aware pose-gold was not established, and any panel docking is retained only as an exploratory repository archive, not as a fifth main pair or domain-extension claim. Post-hoc AND-like dual filter operating points and full-map ligand-only models are Supporting Information diagnostics only (Tables S46–S47); they reinforce that Dual versus neither is chemically easier than Dual versus selectives and do not replace Table 2.
