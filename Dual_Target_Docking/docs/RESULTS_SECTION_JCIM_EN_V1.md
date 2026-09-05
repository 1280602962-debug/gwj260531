# Results (JCIM Articles draft, English)

## 3. Results

### 3.1 Data supply and four-state panels

Public bioactivity resources constrain how strictly a dual-target docking evaluation can be built. In a frozen ChEMBL supply audit of 49 candidate pairs, ligands that are potent on one target and experimentally weak on the other—directional selective hard negatives—were scarce under a strict 6.5/5.5 rule: only four pairs met a thick-panel gate of ≥50 strict hard negatives on both ends. After excluding metal-dependent HDAC1/HDAC6, PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB remained relatively well supplied on ChEMBL labels alone, whereas EGFR/HER2 retained only seven strict B-only ligands and was kept as a supply-limited case (Table 1; Figure 1C). PIK3CA/PIK3CB was docked under this supply screen but was withdrawn from the primary evaluation set after a post-hoc receptor-identity audit found that its "PIK3CB" receptor (PDB 2WXF) is murine PIK3CD, not human PIK3CB (Methods 2.2; Supporting Information); the primary docked set therefore comprises the two supply-qualified ordinary pairs, PIK3CA/mTOR and AChE/BChE, plus the supply-limited EGFR/HER2 case. A zero-docking BindingDB/PubChem recount supported the same scarcity conclusion (Table S12; Figure S2).

Primary docking labels used a unified θ = 6.0 four-state assignment (dual, A-only, B-only, neither), frozen before scores were inspected; the stricter rule served only supply qualification and panel construction (Methods 2.1–2.2). A post-hoc θ = 6.0 label census found 17 unique pairs with directional classes each n ≥ 10 (Table S44), but those counts are supply diagnostics only: docking evaluation remains the three pairs retained after withdrawing PIK3CA/PIK3CB (Methods 2.2). Complete-case dual-measured fractions on the maps were 14.5%–34.0% (Table S37). A metadata review of 186 priority ligands gave 179 include / 7 uncertain / 0 exclude and did not change any frozen class (Table S42). Document-blocked Dual versus B-only on PIK3CA/mTOR was not stably estimable (Table S40). Within the 2026-08-26 API snapshot, max versus median aggregation gave EGFR/HER2 worst-arm AUROC 0.417 and 0.424; that comparison is independent of the frozen Table 2 estimate of 0.430 (Table S29).

### 3.2 Formulation-dependent recognition

On the frozen K = 3 panels, AutoDock Vina mode-1 scores were evaluated with the two primary pocket-matched directional AUROCs—dual versus A-only in pocket B and dual versus B-only in pocket A—with worst-arm AUROC reported as a conservative descriptive summary (Figure 1B; Methods 2.4). EGFR/HER2, AChE/BChE, and PIK3CA/mTOR gave worst-arm AUROC values of 0.430, 0.606, and 0.692, respectively (Table 2; Figure 2A). All three ligand-bootstrap 95% intervals included 0.5. Pair ranking was unchanged under arithmetic, geometric, and harmonic aggregation (Table S26).

As a prespecified formulation contrast on the same scores, Dual versus neither used experimental inactives (`vina_mean`; Table 3; Figure 2B). On EGFR/HER2, Dual versus neither yielded AUROC 0.756 [0.562, 0.920] (n_neg = 12), whereas directional worst-arm AUROC remained 0.430 [0.282, 0.578]; Dual versus all non-duals fell to 0.551. In a mixed-library ranking of all 110 EGFR/HER2 ligands by `vina_mean`, the Top-10 contained one dual and nine experimental selectives (hard-negative fraction 0.90; Table S25). Holding the pocket A score fixed, Dual versus B-only versus Dual versus neither on EGFR/HER2 differed by 0.378 [0.205, 0.547], isolating negative-class composition from mean-score aggregation (Table S34; Figure 2C). AChE/BChE showed only a small Dual-versus-neither increment whose interval overlaps the directional arm. PIK3CA/mTOR Dual versus neither is underpowered (neither n = 4).

Independent GNINA 1.3.2 pose generation on the same frozen EGFR/HER2 ligands, receptors, and boxes left the formulation contrast intact: Dual versus neither 0.783 [0.610, 0.922] versus directional worst-arm AUROC 0.220 [0.109, 0.343], with a Top-10 again dominated by selectives (Table S32; Figure 4A). PIK3CA/mTOR worst-arm AUROC remained 0.633. Across five prespecified Vina seeds, the directional worst-arm estimates showed similar pair-specific patterns (Table S54; Figure 4C).

**Table 2.** Pocket-matched directional AUROC on the frozen K = 3 set (Vina; unified θ = 6.0), with all four prespecified descriptor `summary_min` values. Class sizes are n_scored (dual / A-only / B-only). The highest descriptor is a best single-descriptor reference. A fourth pair, PIK3CA/PIK3CB, is withdrawn from this table after a receptor-identity failure (Methods 2.2; Supporting Information) and is not counted among K.

| Pair | n_scored (dual / A-only / B-only) | dual vs A_only (pocket B) | dual vs B_only (pocket A) | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] | 0.463 | 0.448 | 0.310 | 0.260 |

**Table 3.** Same Vina scores under Dual-versus-neither versus directional formulations (unified θ = 6.0). Dual-versus-neither uses experimental inactives (`vina_mean`). PIK3CA/mTOR Dual versus neither is underpowered (n_neg = 4).

| Pair | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|---|---:|---:|---:|---:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |

A detectable-effect simulation at the observed class sizes shows that these panels resolve large directional effects more readily than moderate ones (Table S31). Failure of a worst-arm AUROC CI to exclude 0.5 therefore does not establish equivalence to chance.

### 3.3 Chemistry as a competing explanation

Docking was compared with four prespecified physicochemical descriptors and with ECFP4 under Bemis–Murcko scaffold GroupKFold (Figure 3; Tables 2, S19–S20, S24). Relative to the best single-descriptor reference, paired worst-arm AUROC differences included zero on all three pairs (Figure S4). On AChE/BChE, mean TPSA differed sharply between dual and selective classes, and TPSA alone exceeded Vina on the corresponding contrast (Figure 3C); adding heavy-atom count and TPSA raised dual-versus-B-only logistic AUROC from 0.606 to 0.807 while the docking odds ratio remained near one (Figure S9). ECFP4 scaffold-grouped fold AUROCs on several arms were well above the docking contrasts (for example ≈0.89 versus 0.43 for EGFR/HER2 dual-versus-B-only; Figure 3A). Under the present scaffold-grouped task, adding the pocket-matched docking score to ECFP4 changed AUROC by at most 0.020 in absolute value (Table S24; Figure 3B). Property-caliper matching of experimental selectives on Dual-versus-B-only arms remained near chance where sample size allowed (Table S45; Figure S9).

Document and scaffold correlation matter for uncertainty. Document-blocked CV left the EGFR/HER2 weak arm at 0.430 (document-cluster bootstrap 95% CI [0.321, 0.617]; Table S39). Scaffold-cluster bootstrap on the same arm gave [0.278, 0.595]. All three estimators keep the weak arm near chance with intervals that span 0.5.

### 3.4 Receptor realization

Holding mTOR frozen and replacing PIK3CA with alternate crystals that passed cognate QC lowered apparent discrimination on PIK3CA/mTOR (Figure 4B; Table S30): replacing PIK3CA 4L23 with 4JPS or 5DXT while holding mTOR at 4JT6 dropped worst-arm AUROC from 0.692 to 0.486 [0.259, 0.692] and 0.505 [0.292, 0.696], i.e. to near chance. A parallel swap was also run on the since-withdrawn PIK3CA/PIK3CB pair and moved in the opposite direction (0.500 to 0.691 and 0.685), but because that pair's B pocket was later found to be a receptor-identity failure (murine PIK3CD, not human PIK3CB; Supporting Information), this comparison is reported only as a caution that cognate-ligand RMSD QC cannot detect a wrong-protein receptor, not as evidence that receptor substitution is direction-dependent across pairs. On the one receptor-verified pair available for this test, substitution reduced rather than reversed apparent discrimination.

Docking failures were concentrated among large or flexible ligands (Table S27). On AChE/BChE, rank-extreme lower bounds remained directionally consistent with complete-case estimates. Unused-pool holdout and mismatched-pocket scoring controls did not show a stable matched-pocket advantage on the holdout panels (Figure 5; both holdout paired 95% CIs include 0). Exhaustiveness, PM110, threshold-grid, and BindingDB-native gate checks bound how far the three-pair evidence can be pushed (Figure 6; Tables S8, S16).

### 3.5 External BindingDB slice

A BindingDB-native 202608 archive rebuild under a contract frozen before docking applied literature, structure, and ECFP4 < 0.70 filters and yielded zero pairs meeting the pre-frozen primary external gate; remaining counts are upper bounds because ChEMBL document lookup was incomplete, and the slice was not docked (Tables S48–S49; Figure 6D; Figure S8).[16] The pre-frozen 2018 literature-year split likewise failed the sample gate on the primary cutoff and is not claimed as external validation (Table S41).

### 3.6 Practical consequences for dual-target screening

An AND-like dual filter at the Dual-median `vina_worst` cutoff on EGFR/HER2 retained 14/28 dual ligands but also 33 selectives (precision 0.298; hard-negative fraction 0.702; Table S46; Figure S7). Full-map ligand-only ECFP4 models on the complete ChEMBL graphs of the three primary pairs recovered Dual versus neither more readily than Dual versus selectives (EGFR/HER2: 0.921 versus 0.864 on Dual versus B-only; Table S47; Figure S7). These diagnostics describe what a two-pocket filter does on experimentally labeled chemistry; they do not replace Table 2.
