# Results (JCIM Articles draft, English)

## 3. Results

### 3.1 Data supply and four-state panels

Public bioactivity resources constrain how strictly a dual-target docking evaluation can be built. In a frozen ChEMBL supply audit of 49 candidate pairs, ligands that are potent on one target and experimentally weak on the other—directional selective hard negatives—were scarce under a strict 6.5/5.5 rule: only four pairs met a thick-panel gate of ≥50 strict hard negatives on both ends. After excluding metal-dependent HDAC1/HDAC6, PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB remained relatively well supplied on ChEMBL labels alone, whereas EGFR/HER2 retained only seven strict B-only ligands and was kept as a supply-limited case (Table 1; Figure 1C). PIK3CA/PIK3CB was docked under this supply screen but was withdrawn from the primary evaluation set after a post-hoc receptor-identity audit found that its "PIK3CB" receptor (PDB 2WXF) is murine PIK3CD, not human PIK3CB (Methods 2.2; Note S0). A later ChEMBL 37 dump census then added five ordinary noncovalent pairs that passed pre-declared G1–G5 gates (F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, PPARA/PPARD). The primary table therefore has eight rows. A zero-docking BindingDB/PubChem recount supported the same scarcity conclusion on the original freeze and did not flip the ≥50 hard-neg gate on the five added pairs (Table S11; Figure S2).

Primary docking labels used a unified θ = 6.0 four-state assignment (dual, A-only, B-only, neither), frozen before scores were inspected; the stricter rule served only supply qualification and panel construction (Methods 2.1–2.2). A post-hoc θ = 6.0 label census found 17 unique pairs with directional classes each n ≥ 10 (Note S14), but those J0-list counts are supply diagnostics only and are not the docking menu. Complete-case dual-measured fractions on the original maps were 14.5%–34.0% (Note S14). A metadata review of 186 priority ligands gave 179 include / 7 uncertain / 0 exclude and did not change any frozen class (Note S14). Document-blocked Dual versus B-only on PIK3CA/mTOR was not stably estimable (Table S10). Within the 2026-08-26 API snapshot, max versus median aggregation gave EGFR/HER2 worst-arm AUROC 0.417 and 0.424; that comparison is independent of the frozen Table 2 estimate of 0.430 (Table S3). One post-census ligand (PPARA/PPARD CHEMBL121) flips class under median aggregation and does not replace the max-pChEMBL labels.

### 3.2 Formulation-dependent recognition

On the eight-row primary set, AutoDock Vina mode-1 scores were evaluated with the two primary pocket-matched directional AUROCs—dual versus A-only in pocket B and dual versus B-only in pocket A—with worst-arm AUROC reported as a conservative descriptive summary (Figure 1B; Methods 2.4). The 2026-07-23 pairs EGFR/HER2, AChE/BChE, and PIK3CA/mTOR gave worst-arm AUROC values of 0.430, 0.606, and 0.692; the five post-census pairs gave 0.345 (F2/F10), 0.365 (JAK1/TYK2), 0.588 (JAK1/JAK2), 0.649 (PPARG/PPARA), and 0.446 (PPARA/PPARD) (Table 2; Figure 2A). All ligand-bootstrap 95% intervals except PPARG/PPARA included 0.5. Eight pairs are not averaged into one master AUROC. Pair ranking on the original three was unchanged under arithmetic, geometric, and harmonic aggregation (Note S14).

As a prespecified formulation contrast on the same scores, Dual versus neither used experimental inactives (`vina_mean`; Table 3; Figure 2B). On EGFR/HER2, Dual versus neither yielded AUROC 0.756 [0.562, 0.920] (n_neg = 12), whereas directional worst-arm AUROC remained 0.430 [0.282, 0.578]; Dual versus all non-duals fell to 0.551. In a mixed-library ranking of all 110 EGFR/HER2 ligands by `vina_mean`, the Top-10 contained one dual and nine experimental selectives (hard-negative fraction 0.90; Table S13). Holding the pocket A score fixed, Dual versus B-only versus Dual versus neither on EGFR/HER2 differed by 0.378 [0.205, 0.547], isolating negative-class composition from mean-score aggregation (Table S4; Figure 2C). JAK1/TYK2 reproduced the same gap (Dual versus neither 0.770 [0.597, 0.906] versus directional 0.365 [0.231, 0.503]). AChE/BChE showed only a small Dual-versus-neither increment whose interval overlaps the directional arm. PIK3CA/mTOR Dual versus neither is underpowered (neither n = 4). PPARG/PPARA Dual versus neither (0.685) is not a substitute for its directional worst arm; unused-pool holdout drops that arm to 0.535 [0.350, 0.717], and same-pose RTMScore collapses it to 0.369 [0.233, 0.475] while Dual versus neither rises to 0.817.

Independent GNINA 1.3.2 pose generation on the same frozen EGFR/HER2 ligands, receptors, and boxes left the formulation contrast intact: Dual versus neither 0.783 [0.610, 0.922] versus directional worst-arm AUROC 0.220 [0.109, 0.343], with a Top-10 again dominated by selectives (Table S9; Figure 4A). PIK3CA/mTOR worst-arm AUROC remained 0.633. The same independent search on JAK1/TYK2 kept the gap (0.705 versus 0.317). Across five prespecified Vina seeds, the directional worst-arm estimates showed similar pair-specific patterns (Table S9; Figure 4C). No post-census five-seed range crossed 0.5.

**Table 2.** Pocket-matched directional AUROC on the eight-row primary set (Vina; unified θ = 6.0), with all four prespecified descriptor `summary_min` values. Class sizes are n_scored (dual / A-only / B-only). The highest descriptor is a best single-descriptor reference. The five post-census rows were added after the 2026-07-23 freeze. PIK3CA/PIK3CB is withdrawn after a receptor-identity failure (Methods 2.2; Note S0) and is not a ninth main-table row.

| Pair | n_scored (dual / A-only / B-only) | dual vs A_only (pocket B) | dual vs B_only (pocket A) | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] | 0.463 | 0.448 | 0.310 | 0.260 |
| F2/F10 | 31 / 32 / 32 | 0.413 | 0.345 | 0.345 [0.211, 0.477] | 0.432 | 0.482 | 0.515 | 0.345 |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 | 0.365 | 0.365 [0.231, 0.503] | 0.369 | 0.389 | 0.580 | 0.425 |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 | 0.728 | 0.588 [0.444, 0.725] | 0.578 | 0.565 | 0.480 | 0.570 |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 | 0.706 | 0.649 [0.504, 0.751] | 0.507 | 0.478 | 0.485 | 0.627 |
| PPARA/PPARD | 32 / 32 / 32 | 0.646 | 0.446 | 0.446 [0.296, 0.584] | 0.490 | 0.436 | 0.564 | 0.351 |

**Table 3.** Same Vina scores under Dual-versus-neither versus directional formulations (unified θ = 6.0). Dual-versus-neither uses experimental inactives (`vina_mean`). PIK3CA/mTOR Dual versus neither is underpowered (n_neg = 4). Dual versus neither is not the primary claim.

| Pair | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|---|---:|---:|---:|---:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |
| F2/F10 | 0.345 [0.211, 0.477] | 0.519 [0.350, 0.688] | 12 | 0.405 [0.273, 0.534] |
| JAK1/TYK2 | 0.365 [0.231, 0.503] | 0.770 [0.597, 0.906] | 14 | 0.527 [0.411, 0.638] |
| JAK1/JAK2 | 0.588 [0.444, 0.725] | 0.730 [0.547, 0.875] | 14 | 0.668 [0.561, 0.770] |
| PPARG/PPARA | 0.649 [0.504, 0.751] | 0.685 [0.493, 0.848] | 14 | 0.675 [0.571, 0.780] |
| PPARA/PPARD | 0.446 [0.296, 0.584] | 0.565 [0.368, 0.766] | 14 | 0.522 [0.406, 0.640] |

A detectable-effect simulation at the observed class sizes shows that these panels resolve large directional effects more readily than moderate ones (Table S10). Failure of a worst-arm AUROC CI to exclude 0.5 therefore does not establish equivalence to chance.

### 3.3 Chemistry as a competing explanation

Docking was compared with four prespecified physicochemical descriptors and with ECFP4 under Bemis–Murcko scaffold GroupKFold (Figure 3; Tables 2 and S5). Relative to the best single-descriptor reference, paired worst-arm AUROC differences included zero on the original three pairs (Figure S4). On the post-census rows, Vina minus the best descriptor is small and positive only on JAK1/JAK2 and PPARG/PPARA and negative on F2/F10, JAK1/TYK2, and PPARA/PPARD; PPARG/PPARA’s directional CI above 0.5 is matched by TPSA (0.627) and is not a docking success beyond chemistry. On AChE/BChE, mean TPSA differed sharply between dual and selective classes, and TPSA alone exceeded Vina on the corresponding contrast (Figure 3C); adding heavy-atom count and TPSA raised dual-versus-B-only logistic AUROC from 0.606 to 0.807 while the docking odds ratio remained near one (Figure S9). ECFP4 scaffold-grouped fold AUROCs on several arms were well above the docking contrasts (for example ≈0.89 versus 0.43 for EGFR/HER2 dual-versus-B-only; Figure 3A). Under the present scaffold-grouped task, adding the pocket-matched docking score to ECFP4 changed AUROC by at most 0.020 in absolute value (Table S5; Figure 3B). Property-caliper matching of experimental selectives on Dual-versus-B-only arms remained near chance where sample size allowed (Note S14; Figure S9).

Document and scaffold correlation matter for uncertainty. Document-blocked CV left the EGFR/HER2 weak arm at 0.430 (document-cluster bootstrap 95% CI [0.321, 0.617]; Table S10). Scaffold-cluster bootstrap on the same arm gave [0.278, 0.595]. All three estimators keep the weak arm near chance with intervals that span 0.5.

### 3.4 Receptor realization

Holding mTOR frozen and replacing PIK3CA with alternate crystals that passed cognate QC lowered apparent discrimination on PIK3CA/mTOR (Figure 4B; Table S8): replacing PIK3CA 4L23 with 4JPS or 5DXT while holding mTOR at 4JT6 dropped worst-arm AUROC from 0.692 to 0.486 [0.259, 0.692] and 0.505 [0.292, 0.696], i.e. to near chance. A parallel swap was also run on the since-withdrawn PIK3CA/PIK3CB pair and moved in the opposite direction (0.500 to 0.691 and 0.685), but because that pair's B pocket was later found to be a receptor-identity failure (murine PIK3CD, not human PIK3CB; Note S0), this comparison is reported only as a caution that cognate-ligand RMSD QC cannot detect a wrong-protein receptor, not as evidence that receptor substitution is direction-dependent across pairs. On the one receptor-verified pair available for this test, substitution reduced rather than reversed apparent discrimination.

Docking failures were concentrated among large or flexible ligands (Table 1). On AChE/BChE, rank-extreme lower bounds remained directionally consistent with complete-case estimates. Unused-pool holdout and mismatched-pocket scoring controls did not show a stable matched-pocket advantage on the original holdout panels (Figure 5; both holdout paired 95% CIs include 0). On the five post-census holdouts, JAK1/JAK2 stayed same-direction (0.619 [0.420, 0.749]; drawn 20/20/18), PPARG/PPARA did not replicate (0.535 [0.350, 0.717]), and F2/F10, JAK1/TYK2, and PPARA/PPARD stayed weak. Exhaustiveness, PM110, threshold-grid, and BindingDB-native gate checks bound how far the eight-row evidence can be pushed (Figure 6; Tables S3 and S11).

### 3.5 External BindingDB slice

A BindingDB-native 202608 archive rebuild under a contract frozen before docking applied literature, structure, and ECFP4 < 0.70 filters and yielded zero pairs meeting the pre-frozen primary external gate; remaining counts are upper bounds because ChEMBL document lookup was incomplete, and the slice was not docked (Table S11; Figure 6D; Figure S8).[16] The pre-frozen 2018 literature-year split likewise failed the sample gate on the original primary cutoff. On the five post-census pairs the same earliest-document-year rule is reportable only for JAK1/TYK2 and JAK1/JAK2, which share JAK1 and are not an independent external set. Unused-pool holdout and the year split are internal confirmations; neither is claimed as external validation (Table S12).

### 3.6 Practical consequences for dual-target screening

An AND-like dual filter at the Dual-median `vina_worst` cutoff on EGFR/HER2 retained 14/28 dual ligands but also 33 selectives (precision 0.298; hard-negative fraction 0.702; Table S13; Figure S7). Full-map ligand-only ECFP4 models on the complete ChEMBL graphs of the original three retained pairs recovered Dual versus neither more readily than Dual versus selectives (EGFR/HER2: 0.921 versus 0.864 on Dual versus B-only; Note S14; Figure S7). These diagnostics describe what a two-pocket filter does on experimentally labeled chemistry; they do not replace Table 2.
