# Methods (English working draft)

## 2. Methods

### 2.1 Study design

Ligands were assigned to dual, A-only, B-only, or neither classes from experimental activity at both targets. Directional evaluation compared dual ligands with each single-target selective class, treating dual as the positive class. Dual and A-only ligands both meet the activity threshold at target A. Their experimental difference lies at target B, so pocket B scores were used. Dual and B-only ligands differ at target A, so pocket A scores were used. The four experimental states and two directional tasks are shown in Figure 1A,B.

The same target score channel was held fixed while dual ligands were compared with single-target selectives and with neither ligands. A ligand-chemistry baseline that does not use receptor structure was then used to test whether that discrimination can be explained by ligand features alone. Pocket correspondence, receptor substitution, and an independent docking implementation were used to test whether remaining discrimination is consistent with the corresponding target structure. Analyses of data processing, sample composition, and external-data availability were used to assess stability and scope.

### 2.2 Bioactivity processing and experimental-state definition

Experimental activities were taken from ChEMBL. Standardized quantitative records, including IC50, Ki, Kd, EC50, and Potency, used pChEMBL on a negative-log molar scale. When several valid records existed for the same ligand–target pair, the primary analysis used the maximum pChEMBL as the representative value. Salt and mixture records were split by connected fragment, retaining the organic fragment with the most heavy atoms. Only ligands with usable experimental activity at both targets entered four-state classification. Absence of a measurement at either end was not treated as low activity at that target.

For any pair A/B, an activity threshold \(\theta\) defined four classes: dual, \(p_{\mathrm{A}}\geq\theta\) and \(p_{\mathrm{B}}\geq\theta\); A-only, \(p_{\mathrm{A}}\geq\theta\) and \(p_{\mathrm{B}}<\theta\); B-only, \(p_{\mathrm{A}}<\theta\) and \(p_{\mathrm{B}}\geq\theta\); neither, both \(<\theta\). Primary analysis used \(\theta=6.0\) throughout.

Bidirectional selective supply for candidate pairs was also counted under a strict 6.5/5.5 separation rule. The active end had pChEMBL \(\geq 6.5\) and the low-activity end \(\leq 5.5\). Ligands in the 5.5–6.5 gray zone were excluded from that strict class. The strict rule was used for candidate-pair supply assessment; primary analysis used \(\theta=6.0\).

Activity-processing checks were kept separate and all held panel membership and docking scores unchanged (Table S3). On a 2026-08-26 same-day API snapshot of the three original pairs, maximum pChEMBL was replaced by the median; that snapshot is not the frozen production cache. A high-confidence human SINGLE PROTEIN field screen was applied to the same snapshot; it is an automatic filter, not paper-by-paper reading. A separate ChEMBL 37 dump-based max-versus-median relabel covers the five added pairs and is not the API snapshot.

### 2.3 Pair screening and panel construction

**Pair-supply screen.** Candidate pairs came from one ChEMBL extract. Targets were restricted to human, single-component SINGLE PROTEIN records with both-end quantitative activity. Maximum pChEMBL was the representative value, and ligands were assigned to four states as in section 2.2. Bidirectional selective supply was counted under the strict 6.5/5.5 rule (active end \(\geq 6.5\), low-activity end \(\leq 5.5\); gray-zone ligands were excluded from that candidate pool). Inclusion also depended on protein class, site chemistry, and whether human experimental structures were suitable for unified noncovalent docking. A pair was retained only if a four-state panel could be drawn under those rules. The included pairs were EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, and PPARA/PPARD.

**Panel construction.** Sample-supply size and structure feasibility differed across pairs, so some pairs used different candidate-pool rules, class quotas, or scaffold caps (Table 1). EGFR/HER2 and PIK3CA/mTOR were drawn from the \(\theta=6.0\) pool. The other six pairs were drawn from the strict 6.5/5.5 pool under quota sampling. EGFR/HER2 and PIK3CA/mTOR limited a single Bemis–Murcko scaffold to at most 5 and 2 repeats within a class, respectively; the remaining pairs had no additional scaffold cap.

**Unified primary analysis.** All primary AUROCs used experimental states reassigned at \(\theta=6.0\). Only ligands with valid both-end Vina scores entered the primary directional table, so n_scored can fall below n_panel. Each directional AUROC uses the corresponding pocket score.

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

Human experimental structures (PDB IDs) for each target are listed in Table 1. Cognate ligands and docking boxes are in Table S2. Each cognate ligand defined the docking site of its receptor. JAK1 6N7A was used for both JAK1/TYK2 and JAK1/JAK2. PPARA 6LXA was used for both PPARG/PPARA and PPARA/PPARD. TYK2 3LXP used the JH1 ATP site. In PPARG 9V8H, the bound PG08-NL peptide was retained, and only the cognate small molecule and crystal waters were removed.

Protein identity, species, binding site, and cognate-ligand position were checked before use. The initial docking box was the Cartesian range of cognate heavy atoms, expanded by 5 Å on each side along \(x\), \(y\), and \(z\). Any edge that remained shorter than 20 Å was extended to 20 Å. Crystal waters and the cognate ligand used to define the site were then removed. Meeko prepared each receptor as PDBQT. Alternate locations without an altLoc flag, or marked A, were retained. Box coordinates are in Table S2.

#### 2.4.2 Ligand preparation

Processed ChEMBL SMILES from section 2.2 were the ligand input. RDKit added hydrogens, and ETKDGv3 generated one three-dimensional starting conformer. Local geometry was then optimized with the MMFF force field for at most 200 iterations. Optimized ligands were converted to PDBQT with Meeko. Three-dimensional embedding used a fixed random seed (Table S1). Protonation states and tautomers were not enumerated systematically. Unspecified stereochemistry followed RDKit’s default handling of the given SMILES.

#### 2.4.3 AutoDock Vina docking and scoring

Primary docking used AutoDock Vina 1.2.7 with the default Vina scoring function. At most nine poses were retained for each ligand–receptor pair, with `energy_range` set to 3 kcal mol\(^{-1}\). PIK3CA/mTOR used exhaustiveness = 16 to obtain a search setting that met the near-native redocking coverage requirement. The other pairs used exhaustiveness = 8. All primary Vina analyses used the rank-1 pose affinity as the ligand–receptor score.

#### 2.4.4 Cognate-ligand redocking

Before panel docking, each main receptor was redocked with its cognate ligand to assess whether the search settings could cover the experimental bound conformation. At most nine poses were retained, and heavy-atom RMSD to the experimental cognate conformation was computed for each saved pose. Search coverage was assessed from the lowest heavy-atom RMSD among all saved poses. AChE 4EY7 and TYK2 3LXP deposited 8 poses; the other primary receptors deposited 9, so the coverage metric is not a uniform best-of-nine. A value below 2.0 Å was taken as evidence that the search could produce a near-native pose, not that rank-1 placement or selectivity prediction is reliable. Rank-1 RMSD was also reported, to separate pose generation from score ranking. The eight added receptors were rechecked with chemically mapped RMSD (Table S2c); near-native calls for those receptors use that table. The original six receptors keep their existing CalcRMS or production-QC tables and were not part of that added recheck. EGFR 3POZ is reconstructed QC, not a recovered original production output. The earlier coordinate-assignment numbers remain in Table S2b as a historical comparison. Boxes and redocking results are listed together in Table S2.[8]

#### 2.4.5 Alternative scoring and independent docking

RTMScore and GNINA 1.3.2 CNN were used to rescore all available Vina poses. RTMScore took the highest score among saved poses as the ligand-level result. GNINA used CNN affinity as the main CNN rescoring readout, with CNNscore as a supplementary analysis. These results re-score poses already generated by Vina and are not independent pose generation.

GNINA 1.3.2 was also used to generate poses independently on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2. This analysis used the same receptors, ligands, and docking boxes as the primary analysis. Exhaustiveness values were 8, 16, and 8, respectively, with at most nine poses per ligand. The ligand-level score was the rank-1 minimizedAffinity, inverted so that higher values indicate more favorable predicted binding. Independent GNINA did not return both-end scores for every ligand. EGFR/HER2 dual-versus-neither used n_neither = 11 (EH120_109 missing) rather than the primary Vina dual-versus-neither count of 12; PIK3CA/mTOR A-only was n = 13 rather than 14; JAK1/TYK2 Dual/A-only/B-only/neither counts were 30/32/29/14 rather than 31/32/32/14. The analysis asks whether the main observation remains under another pose-generation and scoring pipeline. It is not a comparison of overall GNINA versus Vina performance (Table S9).

### 2.5 Evaluation metrics and statistical analysis

#### 2.5.1 Primary directional discrimination metrics

When dual was compared with A-only, both classes meet the activity threshold at target A, so pocket B scores were used: \(\mathrm{AUC}_{D/A}(B)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{A\text{-}only};\;S_{B})\). When dual was compared with B-only, pocket A scores were used: \(\mathrm{AUC}_{D/B}(A)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{B\text{-}only};\;S_{A})\). Dual was the positive class in all directional analyses.

More negative AutoDock Vina affinity indicates more favorable predicted binding. For a common AUROC direction, Vina affinity was converted to \(S_{\mathrm{Vina}}=-E_{\mathrm{Vina}}\), so higher \(S_{\mathrm{Vina}}\) is more favorable. Both directional AUROCs were reported. Their lower value was defined as \(\mathrm{summary}_{\min}\), a descriptive summary of the weaker arm, not a new ligand scoring function: \(\mathrm{summary}_{\min}=\min(\mathrm{AUC}_{D/A}(B),\;\mathrm{AUC}_{D/B}(A))\).

#### 2.5.2 Control-class and two-pocket filter comparisons

Neither ligands, which are below threshold at both targets, were also used as a control class to evaluate dual-versus-neither discrimination. For ligands scored at both targets, the mean score was \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\). AUROC was then computed with dual as the positive class and neither as the control. Merging A-only, B-only, and neither into one control set gave a dual-versus-all-nonduals AUROC that is reported only as a mixed-library descriptive reference (Table S4).

To avoid confounding the comparison by score-channel or aggregation differences, target scores were held fixed and AUROCs were compared under different control classes. Pocket A scores were used for dual-versus-B-only and dual-versus-neither. Pocket B scores were used for dual-versus-A-only and dual-versus-neither (Table S4). Changing the control class also changes the specific compounds. The lower of the two pocket scores, \(S_{\mathrm{worst}}=\min(S_{A},S_{B})\), was also used as a two-pocket filter. The median of dual \(S_{\mathrm{worst}}\) values defined the filter threshold. Dual counts, dual recall, dual precision, and retained single-target selectives at that threshold were tabulated (Table S13). Because the two-pocket mean score changes both control composition and score form, the isolated effect of control-class change was taken from the fixed-pocket comparisons.

#### 2.5.3 Confidence intervals and resampling

The 95% confidence intervals for \(\mathrm{summary}_{\min}\) in Table 2 were estimated from 2000 ligand-level non-stratified percentile bootstrap replicates. Each replicate drew, with replacement, the original number of ligands from the pooled dual, A-only, and B-only set used in that pair’s directional analysis. Dual-versus-A-only and dual-versus-B-only AUROCs were then recomputed, and the smaller value was that replicate’s \(\mathrm{summary}_{\min}\). The 95% interval was the 2.5th and 97.5th percentiles of valid \(\mathrm{summary}_{\min}\) values. Point estimates were computed from the full analysis sample, not from bootstrap means. Dual-versus-neither AUROC intervals used class-stratified percentile bootstrap.

When different scoring methods or matched versus mismatched pocket scores were compared, each scheme reused the same ligand resample to preserve pairing. Cluster bootstrap with Bemis–Murcko scaffold groups and literature-connected groups was used to assess scaffold and source correlation. Algorithm details are in the Supporting Information (Table S10). Table S10 also reports an independent sample-size scenario that preserves class sizes by design; that simulation is not the Table 2 ligand-level interval.

### 2.6 Baselines, controls, and sensitivity analyses

#### 2.6.1 Ligand-chemistry controls

To test discrimination from ligand chemistry without receptor structure, ECFP4 fingerprints (radius = 2, 2048 bits) were computed, together with molecular weight, heavy-atom count, cLogP, and TPSA. Directional AUROCs and \(\mathrm{summary}_{\min}\) were computed for each of the four single descriptors. The single descriptor with the highest \(\mathrm{summary}_{\min}\) on each pair served as the physicochemical baseline. Because that best descriptor was selected on the same panel, its difference from Vina is descriptive only (Table S5).

ECFP4, docking-only, and ECFP4+docking models used logistic regression without feature scaling. Bemis–Murcko scaffolds were the grouping variable, and out-of-fold predictions under the same GroupKFold splits were used to compute AUROC. The docking-only model used only the corresponding directional docking score. That AUROC is computed from out-of-fold predictions, whereas the primary analysis ranks the raw docking scores, so the two values are not identical. The AUROC difference between ECFP4 and ECFP4+docking describes incremental discrimination after adding the docking score. A sensitivity analysis applied StandardScaler on each training fold of the same splits (Table S5).

#### 2.6.2 Structural attribution and score correspondence

The main pocket-correspondence analysis compared \(\mathrm{summary}_{\min}\) under matched versus mismatched conditions and estimated the difference with a paired bootstrap. Same-direction AUROC differences were supplementary (Table S6). A positive \(\Delta=\mathrm{summary}_{\min}^{\mathrm{matched}}-\mathrm{summary}_{\min}^{\mathrm{mismatched}}\) means the weaker matched arm is higher.

Receptor-structure sensitivity was evaluated mainly on PIK3CA/mTOR. With mTOR 4JT6 held fixed, PIK3CA 4L23 was replaced by 4JPS or 5DXT. mTOR 4JT6 was also replaced by 4JSX. Ligand sets, experimental states, score definitions, and statistics were otherwise unchanged (Table S8).

#### 2.6.3 Robustness and external-data availability

Experimental states were reassigned under alternative activity thresholds and activity-aggregation rules, with panel membership and docking scores held fixed (Table S3). For pairs with enough remaining candidates, unused-pool holdout sets were built after excluding main-panel members. Holdout ligands were drawn under a fixed sampling rule that limited over-representation of a single Bemis–Murcko scaffold. Directional AUROCs were recomputed on these different ligands to assess sensitivity to panel membership. Because the ligands still come from the same source, this analysis is an internal robustness check (Table S7). Primary Vina docking was also repeated with several fixed random seeds, and different exhaustiveness settings were compared on PIK3CA/mTOR (Table S9). Alternate seeds use complete-case ligands; AChE/BChE n_complete is 95 on the production seed and 89–90 on the other four. The five added pairs also have a fixed-membership intersection of ligands with both-end scores on every seed; that check does not cover all eight pairs (Table S9e). Missing scores were not treated as low activity.

A literature-year split used the earliest source-document year of each ligand (Table S12). BindingDB[16] and PubChem were counted for all eight pairs under the same four-state rules as a strict 6.5/5.5 supply census (Table S11a). External docking used a separate independence filter at \(\theta=6.0\) (Table S11b), not a row-wise filter of S11a. The development-molecule set included main panels, the expanded PIK3CA/mTOR PM110 panel, and internal holdouts. Shared literature, duplicate structures, and molecules with ECFP4 Tanimoto \(\geq 0.70\) were removed. Dual, A-only, and B-only each needed \(n\geq 20\) with at least three sources per class. Only pairs meeting those conditions would have entered external evaluation and docking. Counts, eligibility criteria, and year cutoffs are in Tables S11 and S12.

### 2.7 Software and reproducibility

Molecular processing, fingerprints, descriptors, and statistics were implemented in Python, mainly with RDKit, NumPy, pandas, SciPy, and scikit-learn. Receptor and ligand PDBQT files were prepared with Meeko. Primary docking used AutoDock Vina 1.2.7. GNINA 1.3.2 was used for supplementary scoring and independent pose generation on selected pairs. RTMScore provided an additional rescoring of Vina poses.

Software versions, random seeds, and main computational parameters are summarized in Table S1. Docking boxes and cognate redocking results are in Table S2.
