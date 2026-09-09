# Docking-Based Dual-Target Recognition: Multi-Pair Evaluation Design and Sources of Discrimination

## Abstract

Molecular docking is often used to rank and prioritize candidate ligands in dual-target virtual screening. In retrospective evaluation, the control set may be compounds below the activity threshold on both targets, or experimentally selective single-target ligands. We built a four-state evaluation on eight human target pairs, comprising dual, A-only, B-only, and neither ligands. Holding the same target score channel fixed, changing the control class altered apparent docking discrimination, and the change was not uniform across pairs or evaluation directions. On EGFR/HER2, the EGFR-pocket AUROC for dual-versus-B-only was 0.430 and rose to 0.808 against neither (difference 0.378 [0.205, 0.547]). JAK1/TYK2 showed a similar difference (0.444 [0.263, 0.620]), and independent GNINA pose generation reproduced both task differences. A scaffold-grouped ligand-only ECFP4 baseline matched docking on several directions. Adding the matched-pocket docking score changed AUROC by at most 0.023. Only the EGFR/HER2 and AChE/BChE main panels showed a matched-pocket advantage on the weaker arm. A matched-pocket advantage was not stably recovered on unused holdout ligands. Among compounds with both-end experimental activities, strong dual-versus-neither discrimination does not necessarily imply strong discrimination against single-target selectives. Retrospective dual-target evaluation should report both directional selectivity comparisons, together with ligand-chemistry and pocket-correspondence controls.

**Keywords:** dual-target docking; experimental state; selectivity; AutoDock Vina; virtual screening evaluation

## 1. Introduction

Multitarget drug design uses one molecule to modulate two or more disease-related targets. Relative to combination therapy, a single multitarget ligand can reduce pharmacokinetic mismatch and formulation complexity.[1,2] Generative design and ultra-large library docking have already produced experimentally confirmed dual-target candidates. POLYGON tested 32 MEK1/mTOR compounds.[20] Ultra-large library docking has also found multitarget ligands among hundreds of millions of molecules, and experimental structures support some predicted binding modes.[19] Such work shows that computational methods can find dual-active molecules in a candidate library. What remains unanswered is how docking scores change meaning once ligands already carry both-end experimental labels.

Molecular docking searches ligand poses in a binding site and scores them. It remains a common ranking method in virtual screening.[3,4] Many benchmarks treat known actives versus artificial decoys as a binary task, as in DUD and DUD-E. LIT-PCBA instead uses high-throughput assay actives and inactives.[5–7] The source and definition of control compounds can change virtual-screening results and even reorder methods.[12,13,21] Control-set composition is therefore part of virtual-screening evaluation itself.

When both-end experimental activities are available, dual-target screening can define four experimental states: dual, A-only, B-only, and neither. Dual differs from A-only or B-only at only one target, and from neither at both targets. Dual-versus-selective and dual-versus-neither discrimination therefore do not mean the same thing. Structure-based dual-target screening usually docks each candidate into both pockets and prefers molecules that score favorably on both. Zhou and coworkers evaluated docking-based dual-kinase virtual screening and found that single-target inhibitors were an important source of predicted dual-inhibitor false positives.[9] Their analysis covered several false-positive types, but it did not treat dual-versus-A-only and dual-versus-B-only as two directional tasks tied to the corresponding targets. Kinase-Bench further showed that experimental selectives can test structural discrimination between related kinases, with the aim of identifying inhibitors selective for one kinase. Distinguishing dual from A-only and B-only ligands among compounds already measured at both ends is a different problem.[22]

Prior work has shown that single-target selectives affect dual-target screening. Less systematically examined is how evaluation-task definition changes the interpretation of docking performance when both-end experimental activities are already available. Does changing the experimental-state comparison, at a fixed score channel, change the performance judgment? How far can ligand chemistry without receptor structure explain the apparent discrimination? Is remaining discrimination consistent with pocket-level structural information, and is that correspondence stable? Answers to these questions clarify what dual-target docking performance evaluation means and where the observed discrimination comes from.

## 2. Methods

### 2.1 Study design

Ligands were assigned to dual, A-only, B-only, or neither classes from experimental activity at both targets. Directional evaluation compared dual ligands with each single-target selective class, treating dual as the positive class. Dual and A-only ligands both meet the activity threshold at target A. Their experimental difference lies at target B, so pocket B scores were used. Dual and B-only ligands differ at target A, so pocket A scores were used. The four experimental states and two directional tasks are shown in Figure 1A,B.

The same target score channel was held fixed while dual ligands were compared with single-target selectives and with neither ligands. A ligand-chemistry baseline that does not use receptor structure was then used to test whether that discrimination can be explained by ligand features alone. Pocket correspondence, receptor substitution, and an independent docking implementation were used to test whether remaining discrimination is consistent with the corresponding target structure. Analyses of data processing, sample composition, and external-data availability were used to assess stability and scope.

### 2.2 Bioactivity processing and experimental-state definition

Experimental activities were taken from ChEMBL. Standardized quantitative records, including IC50, Ki, Kd, EC50, and Potency, used pChEMBL on a negative-log molar scale. When several valid records existed for the same ligand–target pair, the primary analysis used the maximum pChEMBL as the representative value. Salt and mixture records were split by connected fragment, retaining the organic fragment with the most heavy atoms. Only ligands with usable experimental activity at both targets entered four-state classification. Absence of a measurement at either end was not treated as low activity at that target.

For any pair A/B, an activity threshold \(\theta\) defined four classes: dual, \(p_{\mathrm{A}}\geq\theta\) and \(p_{\mathrm{B}}\geq\theta\); A-only, \(p_{\mathrm{A}}\geq\theta\) and \(p_{\mathrm{B}}<\theta\); B-only, \(p_{\mathrm{A}}<\theta\) and \(p_{\mathrm{B}}\geq\theta\); neither, both \(<\theta\). Primary analysis used \(\theta=6.0\) throughout.

Bidirectional selective supply for candidate pairs was also counted under a strict 6.5/5.5 separation rule. The active end had pChEMBL \(\geq 6.5\) and the low-activity end \(\leq 5.5\). Ligands in the 5.5–6.5 gray zone were excluded from that strict class. The strict rule was used for candidate-pair supply assessment; primary analysis used \(\theta=6.0\).

To examine activity-processing choices, the maximum pChEMBL was replaced by the median. Experimental states were also reassigned from higher-confidence human SINGLE PROTEIN records on already scored ligands from a same-day API snapshot. Both analyses kept panel membership and docking scores unchanged (Table S3).

### 2.3 Pair screening and panel construction

**Pair-supply screen.** Candidate pairs came from one ChEMBL extract. Targets were restricted to human, single-component SINGLE PROTEIN records with both-end quantitative activity. Maximum pChEMBL was the representative value, and ligands were assigned to four states as in section 2.2. Bidirectional selective supply was counted under the strict 6.5/5.5 rule (active end \(\geq 6.5\), low-activity end \(\leq 5.5\); gray-zone ligands were excluded from that candidate pool). Inclusion also depended on protein class, site chemistry, and whether human experimental structures were suitable for unified noncovalent docking. A pair was retained only if a four-state panel could be drawn under those rules. The included pairs were EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, and PPARA/PPARD.

**Panel construction.** Sample-supply size and structure feasibility differed across pairs, so some pairs used different candidate-pool rules, class quotas, or scaffold caps (Table 1). EGFR/HER2 and PIK3CA/mTOR were drawn from the \(\theta=6.0\) pool. The other six pairs were drawn from the strict 6.5/5.5 pool under quota sampling. EGFR/HER2 and PIK3CA/mTOR limited a single Bemis–Murcko scaffold to at most 5 and 2 repeats within a class, respectively; the remaining pairs had no additional scaffold cap.

**Unified primary analysis.** All primary AUROCs used experimental states reassigned at \(\theta=6.0\). Only ligands with valid docking scores in the required direction were included, so n_scored can fall below n_panel.

**Table 1.** Composition and main docking settings of the dual-target evaluation panels. Quotas are construction targets (dual / A-only / B-only / neither). n_scored is the dual / A-only / B-only count with valid both-end Vina scores that entered the primary directional AUROCs. Receptor resolutions are in Table S2. The PG08-NL peptide was retained in the PPARG pocket 9V8H.

| Pair | Candidate pool | Quota (D / A / B / N) | PDB (A / B) | n_scored (dual / A-only / B-only) | Vina exhaustiveness |
|------|----------------|----------------------:|-------------|------------------------------------:|--------------------:|
| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 38 / 32 | 8 |
| JAK1/JAK2 | strict 6.5/5.5 | 32 / 32 / 32 / 14 | 6N7A / 8BXH | 32 / 32 / 32 | 8 |
| JAK1/TYK2 | strict 6.5/5.5 | 32 / 32 / 32 / 14 | 6N7A / 3LXP | 31 / 32 / 32 | 8 |
| PIK3CA/mTOR | θ = 6.0 | 18 / 14 / 12 / 4 | 4L23 / 4JT6 | 18 / 14 / 12 | 16 |
| AChE/BChE | strict 6.5/5.5 | 28 / 28 / 28 / 16 | 4EY7 / 4BDS | 27 / 25 / 28 | 8 |
| F2/F10 | strict 6.5/5.5 | 32 / 32 / 32 / 14 | 4UDW / 2JKH | 31 / 32 / 32 | 8 |
| PPARG/PPARA | strict 6.5/5.5 | 32 / 32 / 32 / 14 | 9V8H / 6LXA | 32 / 31 / 32 | 8 |
| PPARA/PPARD | strict 6.5/5.5 | 32 / 32 / 32 / 14 | 6LXA / 5U3Q | 32 / 32 / 32 | 8 |

### 2.4 Receptor and ligand preparation and docking

#### 2.4.1 Receptor preparation and docking-site definition

Human experimental structures and cognate ligands for each target are listed in Table 1. Each cognate ligand defined the docking site of its receptor. JAK1 6N7A was used for both JAK1/TYK2 and JAK1/JAK2. PPARA 6LXA was used for both PPARG/PPARA and PPARA/PPARD. TYK2 3LXP used the JH1 ATP site. In PPARG 9V8H, the bound PG08-NL peptide was retained, and only the cognate small molecule and crystal waters were removed.

Protein identity, species, binding site, and cognate-ligand position were checked before use. The initial docking box was the Cartesian range of cognate heavy atoms, expanded by 5 Å on each side along \(x\), \(y\), and \(z\). Any edge that remained shorter than 20 Å was extended to 20 Å. Crystal waters and the cognate ligand used to define the site were then removed. Meeko prepared each receptor as PDBQT. Alternate locations without an altLoc flag, or marked A, were retained. Box coordinates are in Table S2.

#### 2.4.2 Ligand preparation

Processed ChEMBL SMILES from section 2.2 were the ligand input. RDKit added hydrogens, and ETKDGv3 generated one three-dimensional starting conformer. Local geometry was then optimized with the MMFF force field for at most 200 iterations. Optimized ligands were converted to PDBQT with Meeko. Three-dimensional embedding used a fixed random seed (Table S1). Protonation states and tautomers were not enumerated systematically. Unspecified stereochemistry followed RDKit’s default handling of the given SMILES.

#### 2.4.3 AutoDock Vina docking and scoring

Primary docking used AutoDock Vina 1.2.7 with the default Vina scoring function. At most nine poses were retained for each ligand–receptor pair, with `energy_range` set to 3 kcal mol\(^{-1}\). PIK3CA/mTOR used exhaustiveness = 16 to obtain a search setting that met the near-native redocking coverage requirement. The other pairs used exhaustiveness = 8. All primary Vina analyses used the rank-1 pose affinity as the ligand–receptor score.

#### 2.4.4 Cognate-ligand redocking

Before panel docking, each main receptor was redocked with its cognate ligand to check whether the defined site and search settings could generate a near-native pose. At most nine poses were retained, and heavy-atom RMSD to the experimental cognate conformation was computed for each saved pose. The lowest RMSD among saved poses was used to assess search coverage of near-native geometry. A value below 2.0 Å was taken as evidence that the search could produce a near-native pose. Rank-1 RMSD and the lowest RMSD among the top three poses were also reported, to separate pose generation from score ranking. Boxes and redocking results are listed together in Table S2.

#### 2.4.5 Alternative scoring and independent docking

RTMScore and GNINA 1.3.2 CNN were used to rescore all available Vina poses. RTMScore took the highest score among saved poses as the ligand-level result. GNINA used CNN affinity as the main CNN rescoring readout, with CNNscore as a supplementary analysis. These results re-score poses already generated by Vina and are not independent pose generation.

GNINA 1.3.2 was also used to generate poses independently on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2. This analysis used the same receptors, ligands, and docking boxes as the primary analysis. Exhaustiveness values were 8, 16, and 8, respectively, with at most nine poses per ligand. The ligand-level score was the rank-1 minimizedAffinity, inverted so that higher values indicate more favorable predicted binding. Independent GNINA did not return both-end scores for every ligand. EGFR/HER2 dual-versus-neither used n_neither = 11 (EH120_109 missing) rather than the primary Vina dual-versus-neither count of 12; PIK3CA/mTOR A-only was n = 13 rather than 14; JAK1/TYK2 Dual/A-only/B-only/neither counts were 30/32/29/14 rather than 31/32/32/14. The analysis asks whether the main observation remains under another pose-generation and scoring pipeline. It is not a comparison of overall GNINA versus Vina performance (Table S9).

### 2.5 Evaluation metrics and statistical analysis

#### 2.5.1 Primary directional discrimination metrics

When dual was compared with A-only, both classes meet the activity threshold at target A, so pocket B scores were used: \(\mathrm{AUC}_{D/A}(B)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{A\text{-}only};\;S_{B})\). When dual was compared with B-only, pocket A scores were used: \(\mathrm{AUC}_{D/B}(A)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{B\text{-}only};\;S_{A})\). Dual was the positive class in all directional analyses.

More negative AutoDock Vina affinity indicates more favorable predicted binding. For a common AUROC direction, Vina affinity was converted to \(S_{\mathrm{Vina}}=-E_{\mathrm{Vina}}\), so higher \(S_{\mathrm{Vina}}\) is more favorable. Both directional AUROCs were reported. Their lower value was defined as \(\mathrm{summary}_{\min}\), a descriptive summary of the weaker arm, not a new ligand scoring function: \(\mathrm{summary}_{\min}=\min(\mathrm{AUC}_{D/A}(B),\;\mathrm{AUC}_{D/B}(A))\).

#### 2.5.2 Control-class and two-pocket filter comparisons

Neither ligands, which are below threshold at both targets, were also used as a control class to evaluate dual-versus-neither discrimination. For ligands scored at both targets, the mean score was \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\). AUROC was then computed with dual as the positive class and neither as the control. Merging A-only, B-only, and neither into one control set gave a dual-versus-all-nonduals AUROC that is reported only as a mixed-library descriptive reference (Table S4).

To isolate the effect of changing the control class, the target score was held fixed and only the control composition was changed. Pocket A scores were used for dual-versus-B-only and dual-versus-neither. Pocket B scores were used for dual-versus-A-only and dual-versus-neither (Table S4). The lower of the two pocket scores, \(S_{\mathrm{worst}}=\min(S_{A},S_{B})\), was also used as a two-pocket filter. The median of dual \(S_{\mathrm{worst}}\) values defined the filter threshold. Dual counts, dual recall, dual precision, and retained single-target selectives at that threshold were tabulated (Table S13). Because the two-pocket mean score changes both control composition and score form, the isolated effect of control-class change was taken from the fixed-pocket comparisons.

#### 2.5.3 Confidence intervals and resampling

The 95% confidence intervals for \(\mathrm{summary}_{\min}\) in Table 2 were estimated from 2000 ligand-level non-stratified percentile bootstrap replicates. Each replicate drew, with replacement, the original number of ligands from the pooled dual, A-only, and B-only set used in that pair’s directional analysis. Dual-versus-A-only and dual-versus-B-only AUROCs were then recomputed, and the smaller value was that replicate’s \(\mathrm{summary}_{\min}\). The 95% interval was the 2.5th and 97.5th percentiles of valid \(\mathrm{summary}_{\min}\) values. Point estimates were computed from the full analysis sample, not from bootstrap means. Dual-versus-neither AUROC intervals used class-stratified percentile bootstrap.

When different scoring methods or matched versus mismatched pocket scores were compared, each scheme reused the same ligand resample to preserve pairing. Cluster bootstrap with Bemis–Murcko scaffold groups and literature-connected groups was used to assess scaffold and source correlation. Algorithm details are in the Supporting Information (Table S4; Table S10). Table S10 also reports an independent sample-size scenario that preserves class sizes by design; that simulation is not the Table 2 ligand-level interval.

### 2.6 Baselines, controls, and sensitivity analyses

#### 2.6.1 Ligand-chemistry controls

To test discrimination from ligand chemistry without receptor structure, ECFP4 fingerprints (radius = 2, 2048 bits) were computed, together with molecular weight, heavy-atom count, cLogP, and TPSA. Directional AUROCs and \(\mathrm{summary}_{\min}\) were computed for each of the four single descriptors. The single descriptor with the highest \(\mathrm{summary}_{\min}\) on each pair served as the physicochemical baseline. Because that best descriptor was selected on the same panel, its difference from Vina is descriptive only (Table S5).

ECFP4, docking-only, and ECFP4+docking models used logistic regression without feature scaling. Bemis–Murcko scaffolds were the grouping variable, and out-of-fold predictions under the same GroupKFold splits were used to compute AUROC. The docking-only model used only the corresponding directional docking score. That AUROC is computed from out-of-fold predictions, whereas the primary analysis ranks the raw docking scores, so the two values are not identical. The AUROC difference between ECFP4 and ECFP4+docking describes incremental discrimination after adding the docking score. A sensitivity analysis applied StandardScaler on each training fold of the same splits (Table S5).

#### 2.6.2 Structural attribution and score correspondence

The main pocket-correspondence analysis compared \(\mathrm{summary}_{\min}\) under matched versus mismatched conditions and estimated the difference with a paired bootstrap. Same-direction AUROC differences were supplementary (Table S6). A positive \(\Delta=\mathrm{summary}_{\min}^{\mathrm{matched}}-\mathrm{summary}_{\min}^{\mathrm{mismatched}}\) means the weaker matched arm is higher.

Receptor-structure sensitivity was evaluated mainly on PIK3CA/mTOR. With mTOR 4JT6 held fixed, PIK3CA 4L23 was replaced by 4JPS or 5DXT. mTOR 4JT6 was also replaced by 4JSX. Ligand sets, experimental states, score definitions, and statistics were otherwise unchanged (Table S8).

#### 2.6.3 Robustness and external-data availability

Experimental states were reassigned under alternative activity thresholds and activity-aggregation rules, with panel membership and docking scores held fixed (Table S3). For pairs with enough remaining candidates, unused-pool holdout sets were built after excluding main-panel members. Holdout ligands were drawn under a fixed sampling rule that limited over-representation of a single Bemis–Murcko scaffold. Directional AUROCs were recomputed on these different ligands to assess sensitivity to panel membership. Because the ligands still come from the same source, this analysis is an internal robustness check (Table S7). Primary Vina docking was also repeated with several fixed random seeds, and different exhaustiveness settings were compared on PIK3CA/mTOR (Table S9). Alternate seeds use complete-case ligands; AChE/BChE n_complete is 95 on the production seed and 89–90 on the other four.

A literature-year split used the earliest source-document year of each ligand (Table S12). BindingDB and PubChem were counted for all eight pairs under the same four-state rules (Table S11a). External docking further required an independent-source remainder on the same eight pairs. Shared literature, duplicate structures, and molecules with ECFP4 Tanimoto \(\geq 0.70\) were removed. Dual, A-only, and B-only each needed \(n\geq 20\) with at least three sources per class (Table S11b). Only pairs meeting those conditions would have entered external evaluation and docking. Counts, eligibility criteria, and year cutoffs are in Tables S11 and S12.

### 2.7 Software and reproducibility

Molecular processing, fingerprints, descriptors, and statistics were implemented in Python, mainly with RDKit, NumPy, pandas, SciPy, and scikit-learn. Receptor and ligand PDBQT files were prepared with Meeko. Primary docking used AutoDock Vina 1.2.7. GNINA 1.3.2 was used for supplementary scoring and independent pose generation on selected pairs. RTMScore provided an additional rescoring of Vina poses.

Software versions, random seeds, and main computational parameters are summarized in Table S1. Docking boxes and cognate redocking results are in Table S2.

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

## 4. Discussion

### 4.1 Experimental-state definition and evaluation tasks

What dual-target docking evaluation answers depends on the control ligands and their experimental states. Dual-versus-neither tests whether dual-active ligands can be separated from compounds that are low-activity at both targets. Dual-versus-A-only and dual-versus-B-only test whether docking scores can separate dual ligands from selectives that remain active at one target. These comparisons answer different questions and should not be read as two difficulty levels of the same task.

A fixed score channel removes score aggregation as a confounder. Replacing B-only with neither also replaces the specific molecules, so scaffolds, size, charge, and experimental provenance may change with the class. The analysis therefore reflects a joint change in evaluation class and compound set rather than experimental activity state alone. EGFR/HER2 and JAK1/TYK2 still showed clear differences when the same target score was held fixed and only the control class and its compounds were changed. Aggregation of two-pocket scores therefore cannot explain those cases. The magnitude was not uniform across other pairs and directions. Selective ligands are not a single, uniformly harder control class; the effect of control-class choice still depends on the pair and direction. Building on Zhou and coworkers,[9] the present study treats the two single-target selective states as corresponding directional tasks. It then examines evaluation-design change at a fixed score channel, and uses ligand baselines and pocket correspondence to bound how far the discrimination can be explained.

Classic virtual-screening benchmarks have used different control types. DUD and DUD-E rely on decoys, whereas LIT-PCBA uses inactives from high-throughput assays.[5–7] Recent comparisons further show that changing the source and definition of control compounds can change virtual-screening performance and even reorder methods.[12,13,21] Dual-target evaluation also needs both selectivity directions. A higher dual-versus-neither AUROC shows better separation from both-end low-activity ligands, not automatic evidence of dual-target selectivity recognition. Retrospective evaluations that aim at dual-target selectivity still need both directional selective controls. Kinase-Bench tests selective enrichment against kinase decoys.[22] Distinguishing dual from selective ligands among compounds already measured at both ends is a different problem.

### 4.2 Ligand chemistry and apparent sources of discrimination

Retrospective activity classes can themselves carry chemical structure differences. Simple physicochemical descriptors already separated some experimental states on AChE/BChE. ECFP4 also gave high AUROC on several directions under scaffold-grouped cross-validation. Adding the corresponding docking score changed AUROC by at most 0.023 across 16 directions. Docking scores may still contain receptor-related information. A high retrospective AUROC nevertheless cannot, by its numerical value alone, be attributed to three-dimensional complementarity. Whether docking adds discrimination still needs comparison with a receptor-free baseline under the same split.

Such class-wise chemical differences may reflect how public bioactivity data are generated. Database compounds often come from medicinal-chemistry series around a limited set of lead scaffolds. Different activity classes can therefore differ in scaffold, substitution pattern, and physicochemical distribution at the same time. Caba and coworkers examined protein–ligand representations together with ligand-only Morgan fingerprints in structure-based virtual screening. Adding Morgan fingerprints improved several models but did not further improve the best PLEC model.[23] Ligand-only information can likewise carry a strong classification signal in retrospective structure-based screening. Simple physicochemical matching can reduce some low-dimensional property differences, but it cannot fully control high-dimensional differences among scaffolds and series, and strict matching further reduces available samples.

### 4.3 Structural attribution and pocket correspondence

Kinase cross-docking studies have shown that receptor-conformer choice can affect pose recovery and virtual-screening performance.[14] On PIK3CA/mTOR, directional AUROC also changed after PIK3CA or mTOR receptors were replaced. Because this analysis systematically covered only one primary pair, the result shows receptor-structure sensitivity in that system rather than a general dual-target receptor dependence.

Cognate redocking also showed that generating a near-native pose does not mean that pose receives the highest score. Lowest RMSD therefore mainly indicates whether search can produce a near-native pose. It does not evaluate whether those poses are ranked correctly, and it cannot by itself establish virtual-screening discrimination. Changing pose generation and scoring changed specific AUROC values. Clear dual-versus-selective versus dual-versus-neither differences on EGFR/HER2 and JAK1/TYK2 were still observed under independent GNINA. Evaluation-design differences are therefore not wholly dependent on the primary Vina implementation. RTMScore, GNINA CNN rescoring, and receptor substitution nevertheless showed that specific directional AUROCs are not stable to scoring and structural realization. Independent GNINA tests whether an evaluation-design difference can reappear under another pose-generation pipeline rather than whether a matched-pocket advantage is present.

Score-channel exchange provides a more direct test of structural correspondence. Figure 4 and Table S6 report matched-versus-mismatched \(\mathrm{summary}_{\min}\) differences. On the main panels, EGFR/HER2 and AChE/BChE showed a matched-pocket scoring advantage, but the other six pairs had difference intervals that included 0. All seven usable unused-ligand holdouts also had intervals that included 0. The holdouts did not provide a matched-pocket advantage consistent with the main panels. The result indicates insufficient correspondence evidence, not that structural information is absent. On homologous targets, a mismatched pocket can still produce related and biologically meaningful scores. A high directional AUROC on one panel is therefore not enough to attribute that discrimination to the corresponding pocket.

### 4.4 Practical role and evidence boundaries of dual-target docking

Under the retrospective conditions examined here, docking scores should be treated as ranking or screening signals. They are not, by themselves, evidence of true dual-target activity or selectivity. Favorable scores in both pockets mean only that the molecule received favorable predicted scores at both targets under the present computational settings. They do not confirm experimental dual activity or selectivity.

Prospective studies have shown that multitarget docking and generative design can produce experimentally confirmed multitarget candidates.[19,20] Those studies evaluate the ability to find active molecules from a candidate library, whereas the present analysis evaluates retrospective discrimination among already labeled experimental states. Retrospective AUROC therefore cannot replace prospective experimental hit rates, and it should not be used to dismiss docking for library compression. Likewise, a high retrospective AUROC cannot by itself show that a method will stably find dual-target molecules with a desired activity balance.

Retrospective dual-target docking evaluation should at least report dual-versus-A-only and dual-versus-B-only separately. Ligand-only baselines, matched versus mismatched pocket comparisons, and receptor or scoring sensitivity analyses are then needed to judge whether an apparent AUROC has a stable structural explanation. Direct experimental measurement is still required for dual-target activity, activity balance, and mechanism.

### 4.5 Scope and data limitations

The scope of the present conclusions is limited at the data, panel-construction, and computational-implementation layers.

At the data layer, scarce and heterogeneous both-end measurements are the fundamental constraint. Four-state classification requires experimental measurements of the same ligand at both targets. Public records differ in activity metric, assay conditions, and protein construct, and they often concentrate on a few well-studied medicinal-chemistry series. Activity thresholds, aggregation of replicate records, and high-confidence filtering did not change the main observations. They cannot remove heterogeneity that is inherent to mixed experimental sources and chemical series. After overlapping sources, duplicate structures, and highly similar molecules were excluded, BindingDB did not yield an independent external evaluation set that met the present eligibility criteria. Year splits likewise did not provide a test set that could evaluate both selectivity directions at once.

At the panel-construction layer, bidirectional selective supply differed substantially across pairs, so candidate pools and sampling constraints were not fully uniform. Primary statistics used unified \(\theta=6.0\) labels, with holdout and sensitivity analyses. Construction differences can still affect cross-pair comparability. The eight pairs were used to examine consistency across systems and should not be treated as independent statistical replicates. JAK1 and PPARA were each used twice, and some pairs belong to related kinase or nuclear-receptor families. Except for PPARG/PPARA, the bootstrap 95% intervals for primary Vina \(\mathrm{summary}_{\min}\) did not lie entirely above 0.5, indicating wide statistical uncertainty (Figure S6; Table S10).

At the computational layer, the study relied mainly on one representative crystal structure per target and default AutoDock Vina settings. Cognate redocking showed coverage of near-native poses. The check mainly reflects search coverage of near-native geometry and cannot replace evaluation of score ranking or screening performance. Independent GNINA pose generation, alternative receptor conformers, and other scoring functions covered only some systems. Ligand preparation did not systematically enumerate protonation states, tautomers, or conformational ensembles. The effect of those approximations on prediction accuracy still needs systematic assessment. The specific numerical values observed here should therefore not be extrapolated as a general dual-target docking performance across other protein families and broader chemical space.

Stronger evidence for this type of evaluation will require paired activity measurements of the same compound set at both targets, under more uniform experimental conditions. Paired measurements would reduce heterogeneity from mixed literature and assay conditions. Sampling across chemical space would further reduce dependence on a few medicinal-chemistry series. Such data would be better suited to testing whether four-state evaluation can extend to more protein families, and to supporting year splits, cross-database tests, and prospective experimental validation.

## 5. Conclusions

This study examined dual-target docking evaluation design and sources of discrimination on eight human target pairs with both-end experimental activities.

Holding the same target score channel fixed, changing the control class and its corresponding compound set can change the performance judgment in a pair- and direction-dependent way. Discriminating dual-active from both-end low-activity ligands does not guarantee discriminating dual-active from single-target selectives, because the two comparisons are different discrimination problems.

Some retrospective discrimination can be obtained from ligand chemistry that does not use receptor structure, and the increment after adding a docking score is limited. A matched-pocket scoring advantage was not stably recovered on unused holdout ligands. Retrospective discrimination therefore cannot be attributed simply to pocket-specific three-dimensional complementarity.

Retrospective dual-target docking evaluation should report both directional selectivity tasks and combine them with ligand-chemistry baselines and pocket-correspondence controls. The resulting metrics describe discrimination among already labeled experimental states. They are not equivalent to prospective virtual-screening hit finding.

## Data and Software Availability

Benchmark membership, experimental-state labels, receptor and docking-box definitions, per-ligand docking scores, analysis tables, and all scripts used to regenerate the reported statistics and figures are available in the `Dual_Target_Docking` directory of the public repository at https://github.com/1280602962-debug/gwj260531.

Typeset Supporting Information contains **Tables S1–S13**. The tables cover software and seeds; docking boxes and cognate RMSD; threshold and pChEMBL-aggregation sensitivity; fixed-score control-class contrast; ECFP4 increment; matched versus mismatched pockets; unused-pool holdout; PIK3CA/mTOR crystal substitution; independent GNINA and five-seed Vina; document-cluster uncertainty and sample-size scenario; BindingDB external eligibility; 2018 year-split counts; and EGFR/HER2 operating points. The map from the former working Tables S1–S54, including the old multi-seed Table S54, is `Dual_Target_Docking/data/manuscript_lock/SI_TABLE_MERGE_MAP_v1.csv`.

Per-ligand scores, holdout membership, multi-seed long tables, property-caliper matches, chemotype-matched selectives, aggregation means, complete-case and assay-context ledgers, the J0 candidate-pair census, historical BindingDB REST counts, leave-cognate-out, occupancy snapshots, contact counts, whole-chain sequence identity, and the MCL1/Bcl-xL applicability archive (not Table 2) are listed in **Note S14**. These files will ship with code and SHA-256 checksums in a GitHub Release. A Zenodo DOI will be issued from a tagged snapshot, not from this moving branch.

`data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` indexes the principal numerical results and their source CSVs. SHA-256 checksums of manuscript-facing tables are in `REVISION_CHECKSUM_MANIFEST_v1.csv`. The native-slice contract is `protocol/external_slice_contract.yaml`; the evaluation contract is `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`. The ChEMBL supply audit was frozen on 2026-07-23; the high-confidence activity view was fetched on 2026-08-26; the BindingDB-native archive lock is release 202609. The analysis environment and zero-docking reproduction commands are documented in the repository README. The BindingDB TSV archives themselves are not redistributed; CI checks committed CSVs only.

## References

(1) Anighoro, A.; Bajorath, J.; Rastelli, G. Polypharmacology: Challenges and Opportunities in Drug Discovery. *J. Med. Chem.* **2014**, *57*, 7874–7887. DOI: 10.1021/jm5006463.

(2) Proschak, E.; Stark, H.; Merk, D. Polypharmacology by Design: A Medicinal Chemist’s Perspective on Multitargeting Compounds. *J. Med. Chem.* **2019**, *62*, 420–444. DOI: 10.1021/acs.jmedchem.8b00760.

(3) Kitchen, D. B.; Decornez, H.; Furr, J. R.; Bajorath, J. Docking and Scoring in Virtual Screening for Drug Discovery: Methods and Applications. *Nat. Rev. Drug Discov.* **2004**, *3*, 935–949. DOI: 10.1038/nrd1549.

(4) Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. Model.* **2021**, *61*, 3891–3898. DOI: 10.1021/acs.jcim.1c00203.

(5) Huang, N.; Shoichet, B. K.; Irwin, J. J. Benchmarking Sets for Molecular Docking. *J. Med. Chem.* **2006**, *49*, 6789–6801. DOI: 10.1021/jm0608356.

(6) Mysinger, M. M.; Carchia, M.; Irwin, J. J.; Shoichet, B. K. Directory of Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking. *J. Med. Chem.* **2012**, *55*, 6582–6594. DOI: 10.1021/jm300687e.

(7) Tran-Nguyen, V.-K.; Jacquemard, C.; Rognan, D. LIT-PCBA: An Unbiased Data Set for Machine Learning and Virtual Screening. *J. Chem. Inf. Model.* **2020**, *60*, 4263–4273. DOI: 10.1021/acs.jcim.0c00155.

(8) Su, M.; Yang, Q.; Du, Y.; Feng, G.; Liu, Z.; Li, Y.; Wang, R. Comparative Assessment of Scoring Functions: The CASF-2016 Update. *J. Chem. Inf. Model.* **2019**, *59*, 895–913. DOI: 10.1021/acs.jcim.8b00545.

(9) Zhou, S.; Li, Y.; Hou, T. Feasibility of Using Molecular Docking-Based Virtual Screening for Searching Dual Target Kinase Inhibitors. *J. Chem. Inf. Model.* **2013**, *53*, 982–996. DOI: 10.1021/ci400065e.

(10) Zhou, X.; Guan, J.; Zhang, Y.; Peng, X.; Wang, L.; Ma, J. Reprogramming Pretrained Target-Specific Diffusion Models for Dual-Target Drug Design. In *The Thirty-eighth Annual Conference on Neural Information Processing Systems (NeurIPS 2024)*; 2024. arXiv:2410.20688.

(11) Wu, J.; Qiao, A.; Wang, Z.; Wei, Z.; Chen, S. FuseDiff: Symmetry-Preserving Joint Diffusion for Dual-Target Structure-Based Drug Design. In *Proceedings of the 32nd ACM SIGKDD Conference on Knowledge Discovery and Data Mining, Vol. 2*; ACM: New York, 2026; pp 12432–12443. DOI: 10.1145/3770855.3819050.

(12) Tran-Nguyen, V.-K.; Ballester, P. J. Beware of Simple Methods for Structure-Based Virtual Screening: The Critical Importance of Broader Comparisons. *J. Chem. Inf. Model.* **2023**, *63*, 1401–1405. DOI: 10.1021/acs.jcim.3c00218.

(13) Ahmed, F.; Soellner, M. B.; Brooks, C. L., III. Real-World Assessment of Machine-Learned Docking Using Bioassay-Derived Benchmarks. *J. Chem. Inf. Model.* **2026**, *66*, 8752–8759. DOI: 10.1021/acs.jcim.5c03020.

(14) Schaller, D. A.; Christ, C. D.; Chodera, J. D.; Volkamer, A. Benchmarking Cross-Docking Strategies in Kinase Drug Discovery. *J. Chem. Inf. Model.* **2024**, *64*, 8848–8858. DOI: 10.1021/acs.jcim.4c00905.

(15) Sindt, F.; Bret, G.; Rognan, D. On the Difficulty to Rescore Hits from Ultralarge Docking Screens. *J. Chem. Inf. Model.* **2025**, *65*, 5553–5566. DOI: 10.1021/acs.jcim.5c00730.

(16) Liu, T.; Hwang, L.; Burley, S. K.; Nitsche, C. I.; Southan, C.; Walters, W. P.; Gilson, M. K. BindingDB in 2024: a FAIR Knowledgebase of Protein-Small Molecule Binding Data. *Nucleic Acids Res.* **2025**, *53*, D1633–D1644. DOI: 10.1093/nar/gkae1075.

(17) Tanaka, Y.; Aikawa, K.; Nishida, G.; Homma, M.; Sogabe, S.; Igaki, S.; Hayano, Y.; Sameshima, T.; Miyahisa, I.; Kawamoto, T.; Tawada, M.; Imai, Y.; Inazuka, M.; Cho, N.; Imaeda, Y.; Ishikawa, T. Discovery of Potent Mcl-1/Bcl-xL Dual Inhibitors by Using a Hybridization Strategy Based on Structural Analysis of Target Proteins. *J. Med. Chem.* **2013**, *56*, 9635–9645. DOI: 10.1021/jm401170c.

(18) García-Ortegón, M.; Simm, G. N. C.; Tripp, A. J.; Hernández-Lobato, J. M.; Bender, A.; Bacallado, S. DOCKSTRING: Easy Molecular Docking Yields Better Benchmarks for Ligand Design. *J. Chem. Inf. Model.* **2022**, *62*, 3486–3502. DOI: 10.1021/acs.jcim.1c01334.

(19) Wu, Y.; Vigneron, S.; Braz, J.; Srinivasan, K.; Fink, E. A.; Huang, X.-P.; Xu, X.; Huebner, H.; Kim, J. Y.; Wang, J.; Pfeiffer, T.; Sakamoto, K.; Moroz, Y. S.; Radchenko, D. S.; Rodriguiz, R. M.; Irwin, J. J.; Gmeiner, P.; Billesboelle, C.; Roth, B. L.; Basbaum, A. I.; Manglik, A.; Wetsel, W. C.; Shoichet, B. K. Large Library Docking for Polypharmacology. *J. Med. Chem.* **2026**, *69*, 6210–6229. DOI: 10.1021/acs.jmedchem.5c03810.

(20) Munson, B. P.; Chen, M.; Bogosian, A.; Kreisberg, J. F.; Licon, K.; Kuenzi, B. M.; Ideker, T. De novo generation of multi-target compounds using deep generative chemistry. *Nat. Commun.* **2024**, *15*, 3636. DOI: 10.1038/s41467-024-47120-y.

(21) Gu, S.; Shen, C.; Zhang, X.; Sun, H.; Cai, H.; Luo, H.; Zhao, H.; Liu, B.; Du, H.; Zhao, Y.; Fu, C.; Zhai, S.; Deng, Y.; Liu, H.; Hou, T.; Kang, Y. Benchmarking AI-powered docking methods from the perspective of virtual screening. *Nat. Mach. Intell.* **2025**, *7*, 509–520. DOI: 10.1038/s42256-025-00993-0.

(22) Wei, T.-H.; Zhou, S.-S.; Jing, X.-L.; Liu, J.-C.; Sun, M.; Zhao, Z.-H.; Li, Q.-Q.; Wang, Z.-X.; Yang, J.; Zhou, Y.; Wang, X.; Ling, C.-X.; Ding, N.; Xue, X.; Yu, Y.-C.; Wang, X.-L.; Yin, X.-Y.; Sun, S.-L.; Cao, P.; Li, N.-G.; Shi, Z.-H. Kinase-Bench: Comprehensive Benchmarking Tools and Guidance for Achieving Selectivity in Kinase Drug Discovery. *J. Chem. Inf. Model.* **2024**, *64*, 9528–9550. DOI: 10.1021/acs.jcim.4c01830.

(23) Caba, K.; Tran-Nguyen, V.-K.; Rahman, T.; Ballester, P. J. Comprehensive machine learning boosts structure-based virtual screening for PARP1 inhibitors. *J. Cheminform.* **2024**, *16*, 40. DOI: 10.1186/s13321-024-00832-1.
