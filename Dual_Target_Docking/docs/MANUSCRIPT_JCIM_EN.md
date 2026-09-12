# Evaluating Dual-Target Molecular Docking: Experimental-State Comparisons across Multiple Target Pairs and Sources of Discrimination

## Abstract

Molecular docking is often used to rank candidates in dual-target virtual screening. In retrospective evaluation, the control set may be compounds below the activity threshold at both targets, or single-target selectives that remain active at one target. We defined a four-state evaluation on eight human target pairs—dual, A-only, B-only, and neither—and compared those states after holding the same target score channel fixed. Changing the control class and its compounds altered apparent docking discrimination, and the change varied by pair and direction. On EGFR/HER2, the EGFR-pocket AUROC for dual-versus-B-only was 0.430 and rose to 0.808 against neither (difference 0.378 [0.205, 0.547]); JAK1/TYK2 showed a similar difference (0.444 [0.263, 0.620]). Effects on the other pairs were smaller or less certain. Receptor-free ECFP4 carried substantial class information on several directions. Under the same scaffold-grouped logistic-regression setting, adding the matched-pocket docking score changed AUROC by at most 0.023. A matched-pocket advantage was supported on only two main panels and was not consistently reproduced in holdouts. These results support a narrower interpretation: among compounds with experimental measurements at both targets, strong dual-versus-neither discrimination does not establish strong discrimination against single-target selectives. Retrospective dual-target evaluation should report both directional selectivity comparisons, together with ligand-chemistry and pocket-correspondence controls.

**Keywords:** dual-target docking; experimental state; selectivity; AutoDock Vina; virtual screening evaluation

## 1. Introduction

Multitarget drug design uses one molecule to modulate two or more disease-related targets. Relative to combination therapy, a single multitarget ligand can reduce pharmacokinetic mismatch and formulation complexity.[1,2] Generative design and ultra-large library docking have already produced experimentally confirmed dual-target candidates. POLYGON tested 32 MEK1/mTOR compounds.[19] Ultra-large library docking has also found multitarget ligands among hundreds of millions of molecules, and experimental structures support some predicted binding modes.[18] Such work shows that computational methods can find dual-active molecules in a candidate library. Structure-based dual-target generative models can likewise design against two pockets at once.[10,11] How the interpretation of docking scores changes when benchmark ligands have experimental measurements at both targets requires a separate evaluation.

Molecular docking searches ligand poses in a binding site and scores them. It remains a common ranking method in virtual screening.[3,4] Many benchmarks treat known actives versus artificial decoys as a binary task, as in DUD and DUD-E. LIT-PCBA instead uses high-throughput assay actives and inactives. DOCKSTRING provides a reproducible docking workflow for comparing ligand-design methods.[5–7,17] The source and definition of control compounds can change virtual-screening results and even reorder methods.[12,13,20] Rescoring hits from ultra-large docking screens is also not always stable.[15] Control-set composition is therefore part of virtual-screening evaluation itself.

When the same ligand has usable measurements at both targets, an activity threshold assigns it to one of four states: dual-active at both ends; selective for target A or target B; or below threshold at both ends. These states are labeled dual, A-only, B-only, and neither. Dual differs from A-only or B-only at only one target, and from neither at both, so dual-versus-A-only, dual-versus-B-only, and dual-versus-neither answer different questions. Structure-based dual-target screening usually docks each candidate into both pockets and prefers molecules that score favorably on both. Zhou and coworkers evaluated docking-based dual-kinase virtual screening and found that single-target inhibitors were an important source of predicted dual-inhibitor false positives.[9] Their analysis covered several false-positive types, but it did not explicitly formulate the two selective states as separate target-aligned directional tasks. Kinase-Bench uses experimental selectives to test structural discrimination between related kinases, with the aim of identifying inhibitors selective for one kinase. Distinguishing dual from A-only and B-only ligands among compounds already measured at both ends is a different problem.[21]

Prior work has shown that single-target selectives affect dual-target screening. What remains less clear is how the evaluation task itself changes the interpretation of docking performance when both-end experimental measurements are already available. We asked three questions. First, after the score channel is held fixed, does changing the experimental-state comparison change the performance judgment? Second, what discrimination can ligand chemistry obtain without receptor structure? Third, is the observed docking discrimination consistent with the corresponding pocket, and is that correspondence stable under the available controls?

## 2. Methods

### 2.1 Study design

Ligands were assigned to dual, A-only, B-only, or neither classes from experimental activity at both targets. Directional evaluation compared dual ligands with each single-target selective class, treating dual as the positive class. Dual and A-only ligands both meet the activity threshold at target A. Their experimental difference lies at target B, so pocket B scores were used. Dual and B-only ligands differ at target A, so pocket A scores were used. The four experimental states and two directional tasks are shown in Figure 1A,B.

The primary analysis comprised the two directional AUROCs at a fixed target score channel. Their lower value, \(\mathrm{summary}_{\min}\), is a descriptive weaker-arm summary, not a separate overall-performance endpoint. Fixed-pocket dual-versus-neither comparisons were used to assess the influence of control class. Two-pocket mean-score dual-versus-neither comparisons were reported only as a conventional descriptive reference (Table 3). Ligand-chemistry baselines, matched versus mismatched pocket scores, receptor substitution, independent pose generation, label sensitivity, unused-pool holdouts, and publication-year subset analyses were used to examine possible sources of discrimination and the sensitivity of the main findings. BindingDB and PubChem counts asked whether an external evaluation set could be formed.

### 2.2 Bioactivity processing and experimental-state definition

Experimental activities were taken from ChEMBL. Table 2 primary labels used the maximum pChEMBL on the production activity tables for each pair. For EGFR/HER2, AChE/BChE, and PIK3CA/mTOR those tables were a 2026-07-23 ChEMBL REST harvest. For JAK1/JAK2, JAK1/TYK2, F2/F10, PPARG/PPARA, and PPARA/PPARD they were extracted from the ChEMBL 37 SQLite dump (release 2026-05-01). The pair-universe census used the same ChEMBL 37 dump. Later API-snapshot and dump-based checks were sensitivity analyses only (Table S3). Standardized quantitative records, including IC50, Ki, Kd, EC50, and Potency, used pChEMBL on a negative-log molar scale. When several valid records existed for the same ligand–target pair, the primary analysis used the maximum pChEMBL as the representative value. That choice favors the strongest recorded value and can be affected by a single high reading, assay condition, or endpoint type. Salt and mixture records were split by connected fragment, retaining the organic fragment with the most heavy atoms. Only ligands with usable experimental measurements at both targets entered four-state classification. Absence of a measurement at either end was not treated as low activity at that target.

For any pair A/B, an activity threshold \(\theta\) defined four classes: dual, \(p_{\mathrm{A}}\geq\theta\) and \(p_{\mathrm{B}}\geq\theta\); A-only, \(p_{\mathrm{A}}\geq\theta\) and \(p_{\mathrm{B}}<\theta\); B-only, \(p_{\mathrm{A}}<\theta\) and \(p_{\mathrm{B}}\geq\theta\); neither, both \(<\theta\). Primary analysis used \(\theta=6.0\) throughout. A-only and B-only denote experimental state relative to the current pair and threshold, not proteome-wide selectivity. Dual means only that both study targets meet the threshold; it does not imply activity balance, cellular effect, or therapeutic value. Later references to single-target selectives use this relative definition.

Bidirectional selective supply for candidate pairs was also counted under a strict 6.5/5.5 separation rule. The active end had pChEMBL \(\geq 6.5\) and the low-activity end \(\leq 5.5\). Ligands in the 5.5–6.5 gray zone were excluded from that strict class. The strict rule was used for candidate-pair supply assessment; primary analysis used \(\theta=6.0\).

Activity-processing checks were kept separate and all held panel membership and docking scores unchanged (Table S3). Maximum-versus-median aggregation was examined on two archived sources that are not interchangeable with the production tables or with each other. The 2026-08-26 ChEMBL API snapshot covered EGFR/HER2, AChE/BChE, and PIK3CA/mTOR. JAK1/JAK2, JAK1/TYK2, F2/F10, PPARG/PPARA, and PPARA/PPARD were reassessed using the ChEMBL 37 dump. Those checks only partly mitigate the bias of taking the strongest record, and differences between sources cannot be treated as equivalent sensitivity tests. A high-confidence human SINGLE PROTEIN field screen was applied to the API snapshot; it is an automatic filter, not paper-by-paper reading.

### 2.3 Pair screening and panel construction

**Pair-supply screen.** The pair universe was all unordered pairs among human, single-component SINGLE PROTEIN targets (5,869 targets; 4,672 with at least one qualifying molecule; 10,911,456 possible pairs). Self-pairs were excluded. A pair entered the census if at least one ligand had quantitative activity at both ends (2,164,618 pairs). Subsequent gates used the ChEMBL 37 dump. Requiring at least 10 ligands measured at both ends left 63,790 pairs. Requiring at least 10 dual, A-only, and B-only ligands at \(\theta=6.0\) left 5,253 pairs. A strict 6.5/5.5 bidirectional-selective rule requiring at least 50 A-only ligands and 50 B-only ligands left 86 pairs. Removing qHTS hub proteins, CYP ADME panels, and zinc-dependent metal enzymes left 26 pairs. Requiring at least five human holo structures per end (\(\leq 3.5\) Å, at least one non-polymer ligand) left 19 pairs. A drug-like small-molecule filter that retained at least 50 ligands in each strict selective class left 17 pairs spanning 12 target systems. The remaining candidates were further restricted to targets with conventional noncovalent small-molecule pockets that could be represented under the common rigid-receptor Vina protocol. Membrane GPCRs, transporters, reversible-covalent sites, and other systems requiring treatment outside this common protocol were excluded; pair-level retention and exclusion reasons are reported in Table S14.

After those protocol-compatibility restrictions, seven pairs remained: PIK3CA/mTOR, AChE/BChE, F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, and PPARA/PPARD. CTSK/CTSS reached the earlier supply and structure gates but was excluded because both holos are reversible-covalent. EGFR/HER2 was retained as a supply-limited case: it did not meet the strict 6.5/5.5 selective-supply criterion (minimum selective-class count = 7) but had suitable human holo structures, a cognate-defined docking site, and sufficient dual, A-only, and B-only ligands at \(\theta=6.0\) for directional evaluation. These eight pairs are the primary evaluation set. Pair-level retention and exclusion at the final structure- and protocol-compatibility gates is listed in Table S14. Pairs that did not satisfy the final protocol-compatibility criteria were not included in the docking evaluation; the corresponding reasons are listed in Table S14.

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

Protein identity, species, binding site, and cognate-ligand placement were verified before receptor inclusion. The initial docking box was the Cartesian range of cognate heavy atoms, expanded by 5 Å on each side along \(x\), \(y\), and \(z\). Any edge that remained shorter than 20 Å was extended to 20 Å. Crystal waters and the cognate ligand used to define the site were then removed. Meeko prepared each receptor as PDBQT. Alternate locations without an altLoc flag, or marked A, were retained. Box coordinates are in Table S2.

#### 2.4.2 Ligand preparation

Processed ChEMBL SMILES from section 2.2 were the ligand input. RDKit added hydrogens, and ETKDGv3 generated one three-dimensional starting conformer. Local geometry was then optimized with the MMFF force field for at most 200 iterations. Optimized ligands were converted to PDBQT with Meeko. Three-dimensional embedding used a fixed random seed (Table S1). Protonation states and tautomers were not enumerated systematically. Unspecified stereochemistry followed RDKit’s default handling of the given SMILES.

#### 2.4.3 AutoDock Vina docking and scoring

Primary docking used AutoDock Vina 1.2.7 with the default Vina scoring function. At most nine poses were retained for each ligand–receptor pair, with `energy_range` set to 3 kcal mol\(^{-1}\). PIK3CA/mTOR used exhaustiveness = 16 to obtain a search setting that met the near-native redocking coverage requirement. The other pairs used exhaustiveness = 8. All primary Vina analyses used the rank-1 pose affinity reported as `REMARK VINA RESULT` in the first saved model. That field is the production `*_affinity` column used for Table 2.

#### 2.4.4 Cognate-ligand redocking

Before panel docking, each main receptor was redocked with its cognate ligand to assess whether the search settings could cover the experimental bound conformation. At most nine poses were retained. Search coverage used the lowest heavy-atom RMSD among saved poses; rank-1 RMSD was reported separately. AChE 4EY7 and TYK2 3LXP deposited 8 poses, and the other primary receptors deposited 9, so the coverage metric is not a uniform best-of-nine. A value below 2.0 Å indicates that the search could produce a near-native pose, not that rank-1 placement or selectivity prediction is reliable. All 14 primary receptors were then rescored with one chemically mapped CalcRMS recipe (Table S2c); no redocking was performed. For PIK3CA/mTOR the prepared ligand is not in the crystal frame, so RMSD used graph-automorphism CalcRMS to the crystal coordinates. EGFR 3POZ is reconstructed QC, not a recovered production output. Redocking remains a search-coverage check. Coordinate-assignment values remain in Table S2b as a historical comparison.[8]

#### 2.4.5 Alternative scoring and independent docking

RTMScore and GNINA 1.3.2 CNN were used to rescore all available Vina poses. RTMScore took the highest score among saved poses as the ligand-level result. GNINA used CNN affinity as the main CNN rescoring readout, with CNNscore as a supplementary analysis. These results re-score poses already generated by Vina and are not independent pose generation.

GNINA 1.3.2 was also used to generate poses independently on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2. This analysis used the same receptors, ligands, and docking boxes as the primary analysis. Exhaustiveness values were 8, 16, and 8, respectively, with at most nine poses per ligand. The ligand-level score was the rank-1 minimizedAffinity, inverted so that higher values indicate more favorable predicted binding. Independent GNINA did not return both-end scores for every ligand. EGFR/HER2 dual-versus-neither used n_neither = 11 (EH120_109 missing) rather than the primary Vina dual-versus-neither count of 12; PIK3CA/mTOR A-only was n = 13 rather than 14; JAK1/TYK2 Dual/A-only/B-only/neither counts were 30/32/29/14 rather than 31/32/32/14. The analysis asks whether the main observation remains under another pose-generation and scoring pipeline. It is not a comparison of overall GNINA versus Vina performance (Table S9).

### 2.5 Evaluation metrics and statistical analysis

#### 2.5.1 Primary directional discrimination metrics

When dual was compared with A-only, both classes meet the activity threshold at target A, so pocket B scores were used: \(\mathrm{AUC}_{D/A}(B)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{A\text{-}only};\;S_{B})\). When dual was compared with B-only, pocket A scores were used: \(\mathrm{AUC}_{D/B}(A)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{B\text{-}only};\;S_{A})\). Dual was the positive class in all directional analyses.

More negative AutoDock Vina affinity indicates more favorable predicted binding. For a common AUROC direction, Vina affinity was converted to \(S_{\mathrm{Vina}}=-E_{\mathrm{Vina}}\), so higher \(S_{\mathrm{Vina}}\) is more favorable. Both directional AUROCs were the primary endpoints. Their lower value was defined as \(\mathrm{summary}_{\min}\), a descriptive weaker-arm summary rather than an overall-performance metric or a new ligand scoring function: \(\mathrm{summary}_{\min}=\min(\mathrm{AUC}_{D/A}(B),\;\mathrm{AUC}_{D/B}(A))\). Because the minimum of two AUROCs is selectively downward-shifted, it is not treated as an ordinary independent endpoint.

#### 2.5.2 Control-class and two-pocket filter comparisons

Neither ligands, which are below threshold at both targets, were also used as a control class to evaluate dual-versus-neither discrimination. For ligands scored at both targets, the mean score was \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\). AUROC was then computed with dual as the positive class and neither as the control. Merging A-only, B-only, and neither into one control set gave a dual-versus-all-nonduals AUROC that is reported only as a mixed-library descriptive reference (Table S4).

To avoid confounding the comparison by score-channel or aggregation differences, target scores were held fixed and AUROCs were compared under different control classes. Pocket A scores were used for dual-versus-B-only and dual-versus-neither. Pocket B scores were used for dual-versus-A-only and dual-versus-neither (Table S4). Changing the control class also changes the specific compounds. The lower of the two pocket scores, \(S_{\mathrm{worst}}=\min(S_{A},S_{B})\), was also used as a two-pocket filter. The median of dual \(S_{\mathrm{worst}}\) values defined the filter threshold. Dual counts, dual recall, dual precision, and retained single-target selectives at that threshold were tabulated (Table S13). The fixed-pocket comparison was therefore used to assess the influence of control definition while avoiding simultaneous changes in score channel or score aggregation.

#### 2.5.3 Confidence intervals and resampling

The 95% confidence intervals for \(\mathrm{summary}_{\min}\) in Table 2 were estimated from 2000 ligand-level non-stratified percentile bootstrap replicates. Each replicate drew, with replacement, the original number of ligands from the pooled dual, A-only, and B-only set used in that pair’s directional analysis, so class counts can vary across replicates. Dual-versus-A-only and dual-versus-B-only AUROCs were then recomputed, and the smaller value was that replicate’s \(\mathrm{summary}_{\min}\). If a replicate lacked a required class, it was omitted from the interval. The 95% interval was the 2.5th and 97.5th percentiles of valid \(\mathrm{summary}_{\min}\) values. Point estimates were computed from the full analysis sample, not from bootstrap means. Dual-versus-neither AUROC intervals, and the dual-versus-all-nonduals intervals, used class-stratified percentile bootstrap. The pooled interval is retained because \(\mathrm{summary}_{\min}\) is a joint function of both arms. A class-stratified alternative for \(\mathrm{summary}_{\min}\) was not recomputed; Table S10 records that choice.

When different scoring methods or matched versus mismatched pocket scores were compared, each scheme reused the same ligand resample to preserve pairing. Cluster bootstrap with Bemis–Murcko scaffold groups and literature-connected groups was used to assess scaffold and source correlation. Algorithm details are in the Supporting Information (Table S10). Table S10 also reports an independent sample-size scenario that preserves class sizes by design; that simulation is not the Table 2 ligand-level interval.

### 2.6 Baselines, controls, and sensitivity analyses

#### 2.6.1 Ligand-chemistry controls

To test discrimination from ligand chemistry without receptor structure, ECFP4 fingerprints (radius = 2, 2048 bits) were computed, together with molecular weight, heavy-atom count, cLogP, and TPSA. Directional AUROCs and \(\mathrm{summary}_{\min}\) were computed for each of the four single descriptors. The single descriptor with the highest \(\mathrm{summary}_{\min}\) on each pair served as the physicochemical baseline. Because that best descriptor was selected on the same panel, its difference from Vina is descriptive only (Table S5).

ECFP4, docking-only, and ECFP4+docking models used logistic regression without feature scaling. Bemis–Murcko scaffolds were the grouping variable, and out-of-fold predictions under the same GroupKFold splits were used to compute AUROC. The docking-only model used only the corresponding directional docking score. That AUROC is computed from out-of-fold predictions, whereas the primary analysis ranks the raw docking scores, so the two values are not identical. The AUROC difference between ECFP4 and ECFP4+docking describes incremental discrimination after adding the docking score under that shared split. The Figure 3A comparison of raw Vina ranking with scaffold-grouped ECFP4 out-of-fold predictions is a competing-explanation contrast, not a head-to-head predictive benchmark. A sensitivity analysis applied StandardScaler on each training fold of the same splits (Table S5).

#### 2.6.2 Structural attribution and score correspondence

Matched scoring used pocket B for dual-versus-A-only and pocket A for dual-versus-B-only. The mismatched control exchanged these two precomputed score channels without redocking. The analysis compared \(\mathrm{summary}_{\min}\) under the two assignments and estimated \(\Delta=\mathrm{summary}_{\min}^{\mathrm{matched}}-\mathrm{summary}_{\min}^{\mathrm{mismatched}}\) with a paired bootstrap. A positive \(\Delta\) means the weaker matched arm is higher. Same-direction AUROC differences were supplementary (Table S6).

Receptor-structure sensitivity was evaluated mainly on PIK3CA/mTOR. With mTOR 4JT6 held fixed, PIK3CA 4L23 was replaced by 4JPS or 5DXT. mTOR 4JT6 was also replaced by 4JSX. Ligand sets, experimental states, score definitions, and statistics were otherwise unchanged (Table S8).

#### 2.6.3 Robustness and external-data availability

Experimental states were reassigned under alternative activity thresholds and activity-aggregation rules, with panel membership and docking scores held fixed (Table S3). For pairs with enough remaining candidates, unused-pool holdout sets were built after excluding main-panel members. Holdout ligands were drawn under a fixed sampling rule that limited over-representation of a single Bemis–Murcko scaffold. Directional AUROCs were recomputed on these different ligands to assess sensitivity to panel membership. Because the ligands still come from the same source, this analysis is an internal robustness check (Table S7). Primary Vina docking was also repeated with several fixed random seeds, and different exhaustiveness settings were compared on PIK3CA/mTOR (Table S9). Missing scores were not treated as low activity. Complete-case counts and fixed-membership intersections are in Table S9e.

A publication-year subset analysis used the earliest source-document year of each ligand on the already-built panels (Table S12). It is an internal temporal subset, not a trained-then-tested temporal validation. BindingDB[16] and PubChem were counted for all eight pairs under the same four-state rules as a strict 6.5/5.5 supply census (Table S11a). External docking used a separate independence filter at \(\theta=6.0\) (Table S11b), not a row-wise filter of S11a. The development-molecule set included main panels, the expanded PIK3CA/mTOR PM110 panel, and internal holdouts. Shared literature, duplicate structures, and molecules with ECFP4 Tanimoto \(\geq 0.70\) were removed. Dual, A-only, and B-only each needed \(n\geq 20\) with at least three sources per class. Only pairs meeting those conditions would have entered external evaluation and docking. Counts, eligibility criteria, and year cutoffs are in Tables S11 and S12.

### 2.7 Software and reproducibility

Molecular processing, fingerprints, descriptors, and statistics were implemented in Python, mainly with RDKit, NumPy, pandas, SciPy, and scikit-learn. Receptor and ligand PDBQT files were prepared with Meeko. Primary docking used AutoDock Vina 1.2.7. GNINA 1.3.2 was used for supplementary scoring and independent pose generation on selected pairs. RTMScore provided an additional rescoring of Vina poses.

Software versions, random seeds, and main computational parameters are summarized in Table S1. Docking boxes and cognate redocking results are in Table S2.

## 3. Results

### 3.1 Both-end experimental supply and four-state panel construction

Both-end ChEMBL activity that can support four-state evaluation declined rapidly as sample requirements increased. Among unordered pairs of human single-component SINGLE PROTEIN targets, 2,164,618 pairs had at least one ligand measured at both targets, and 63,790 pairs had at least 10. Requiring at least 10 dual, A-only, and B-only ligands at \(\theta=6.0\) left 5,253 pairs. A strict 6.5/5.5 bidirectional selective-supply criterion reduced that number to 86 pairs. This supply change is shown in Figure 1C. The census counts and the later structure and pocket gates are defined in section 2.3.

The census describes four-state supply, not the docking menu. After those gates, the primary evaluation comprised EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, and PPARA/PPARD (Table 1). EGFR/HER2 and PIK3CA/mTOR were drawn from the \(\theta=6.0\) pool; the other six pairs were drawn from the strict 6.5/5.5 pool. Four-state classification requires both-end measurements of the same ligand, so both-end coverage and bidirectional selective supply jointly limited panel size.

![Figure 1](../figures/jcim_article/Fig1_four_state_and_supply.png)

**Figure 1.** Four-state dual-target evaluation and data supply. (A) Four experimental states defined by threshold \(\theta\); (B) two directional tasks: dual versus A-only uses the target B score, and dual versus B-only uses the target A score; (C) the census summarizes the availability of paired experimental data under increasing supply requirements. The primary evaluation comprised eight target pairs retained according to the panel-construction and structural criteria in Table 1. The eight-pair evaluation is shown separately from the supply census.

### 3.2 Experimental-state definition and directional docking evaluation

With the same score channel held fixed, AUROC differences from changing the experimental-state comparison varied by pair and direction. Using the EGFR-pocket score, the EGFR/HER2 dual-versus-B-only AUROC was 0.430 and rose to 0.808 against neither (difference 0.378 [0.205, 0.547]); the other direction was smaller. Using the JAK1-pocket score, the JAK1/TYK2 difference was 0.444 [0.263, 0.620]. Effects on the remaining pairs were smaller or uncertain (Figure 2A; Table S4). Cluster resampling is reported in section 3.5.

Under unified \(\theta=6.0\) labels, the two directional AUROCs, and their descriptive weaker-arm value \(\mathrm{summary}_{\min}\), ranged from 0.345 to 0.692. The PPARG/PPARA \(\mathrm{summary}_{\min}\) interval lay entirely above 0.5 (0.649 [0.504, 0.751]), the F2/F10 interval lay entirely below 0.5 (0.345 [0.211, 0.477]), and the remaining six pairs crossed 0.5 (Figure 2B; Table 2). These values are pair-specific and are not an eight-pair ranking.

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

As a descriptive two-pocket mean-score comparison, dual-versus-neither AUROCs were 0.756 [0.562, 0.920] on EGFR/HER2 and 0.770 [0.597, 0.906] on JAK1/TYK2. The directional \(\mathrm{summary}_{\min}\) values were 0.430 and 0.365 (Figure 2C; Table 3). Because that comparison changes both score aggregation and control definition, the influence of control definition was taken from the fixed-score-channel analysis.

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

Panel ranking was a retrospective operating-point diagnostic on this constructed set, not a commercial-library screen. When all 110 EGFR/HER2 ligands were ranked by the two-pocket mean score, the Top-10 contained 1 dual, 5 A-only, and 4 B-only ligands, and no neither ligand; the denominator is 110 (Figure S1A; Table S13). The higher dual-versus-neither AUROC therefore did not correspond to fewer highly ranked single-target selectives. A two-pocket filter at the median dual \(S_{\mathrm{worst}}\) was applied to Dual+A-only+B-only (n = 98), excluding neither, and retained 14 dual ligands together with 9 A-only and 24 B-only ligands (dual precision 0.298; Figure S1B; Table S13).

### 3.3 Ligand-chemistry baselines and directional discrimination

On some pairs, ligand physicochemical features alone separated dual ligands from single-target selectives. On AChE/BChE, TPSA alone gave dual-versus-A-only and dual-versus-B-only AUROCs of 0.733 and 0.801 (Figure 3C; Table S5). Single-descriptor discrimination varied across pairs. On PIK3CA/mTOR, the best single descriptor (heavy-atom count) had \(\mathrm{summary}_{\min}\) 0.463 (Table S5). Across the eight pairs, Vina and the best single descriptor differed in direction and magnitude; six of eight difference 95% intervals included 0, whereas F2/F10 and JAK1/TYK2 excluded 0 (Table S5). Vina \(\mathrm{summary}_{\min}\) intervals and best-descriptor points are in Figure S2.

Under Bemis–Murcko scaffold-grouped cross-validation, receptor-free ECFP4 fingerprints matched or exceeded docking AUROC on several directions (Figure 3A). That comparison is a competing-explanation contrast: ECFP4 uses scaffold-grouped out-of-fold predictions, whereas primary Vina uses full-panel raw ranking. Adding the corresponding docking score to an ECFP4 logistic model changed AUROC by at most 0.023 across all 16 directions (Figure 3B; Table S5). No consistent point-estimate improvement was observed across the 16 directions.

![Figure 3](../figures/jcim_article/Fig3_ligand_chemistry.png)

**Figure 3.** Ligand chemistry as a competing explanation, not a head-to-head predictive benchmark. (A) Rank AUROC from raw Vina scores versus out-of-fold ECFP4 AUROC under scaffold-grouped cross-validation; blue, Vina; orange, ECFP4; (B) AUROC change after adding the corresponding Vina score to ECFP4 on the same splits; lines join the two points for one direction and are not confidence intervals; (C) TPSA distributions on AChE/BChE.

### 3.4 Pocket correspondence and computational sensitivity

Figure 4A and Tables S6 and S7 compare matched-pocket versus mismatched-pocket \(\mathrm{summary}_{\min}\). Only EGFR/HER2 (0.170 [0.060, 0.280]) and AChE/BChE (0.161 [0.037, 0.269]) had main-panel difference intervals that excluded 0; the other six pairs crossed 0. On the seven internal holdouts, the intervals all included 0 (differences from −0.079 to +0.150). A matched-pocket advantage was not consistently reproduced across the available holdouts (Figure 4A; Table S7). EGFR/HER2 had too few leftover candidates for an equivalent unused-pool holdout.

![Figure 4](../figures/jcim_article/Fig4_mismatched_pocket.png)

**Figure 4.** Matched-pocket versus mismatched-pocket scores. (A) Matched−mismatched \(\Delta\mathrm{summary}_{\min}\) on the main panels and holdouts; points and bars are estimates and ligand-level bootstrap 95% confidence intervals; (B) \(\mathrm{summary}_{\min}\) on the main panels and holdouts. \(\dagger\) marks the missing EGFR/HER2 unused-pool holdout. The bottom legend distinguishes main panels from holdouts.

After independent GNINA 1.3.2 pose generation and scoring, the EGFR/HER2 dual-versus-neither AUROC was 0.783 [0.610, 0.922] with n_neither = 11. The weaker directional arm dual-versus-B-only was 0.220 [0.109, 0.343] (n_dual = 28, n_B = 32). That interval belongs to dual-versus-B-only scored in pocket A; the independent-GNINA `summary_min` row has no bootstrap interval (Figure 5A; Table S9). JAK1/TYK2 showed the same pattern, with dual-versus-neither 0.705 and directional \(\mathrm{summary}_{\min}\) 0.317 [0.183, 0.463]. Similar task differences persisted on EGFR/HER2 and JAK1/TYK2 under independent pose generation.

On PIK3CA/mTOR, replacing PIK3CA 4L23 with 4JPS lowered \(\mathrm{summary}_{\min}\) from 0.692 [0.470, 0.813] to 0.486 [0.259, 0.692]. Replacement with 5DXT gave 0.505 [0.292, 0.696]. Replacing mTOR 4JT6 with 4JSX gave 0.639 [0.418, 0.776] (Figure 5B; Table S8).

Five fixed Vina random seeds produced comparatively limited numerical fluctuation relative to the larger task and receptor effects (Figure 5C; Table S9). The EGFR/HER2 task difference was positive on all five Vina seeds. Complete-case counts, fixed-membership intersections, and weaker-arm switches are in Table S9e.

PPARG/PPARA was the only pair whose primary Vina \(\mathrm{summary}_{\min}\) interval lay entirely above 0.5 (0.649 [0.504, 0.751]). Same-pose RTMScore rescoring lowered it to 0.369 [0.233, 0.475], and GNINA CNN rescoring lowered it to 0.500. The unused-pool holdout was 0.535 [0.350, 0.717] (Table S7; Table S9).

PIK3CA/mTOR panel-size and exhaustiveness checks are in Figure S3. PM48 is the primary panel (quota n = 48, exhaustiveness = 16); PM110 is a larger protocol-sensitivity panel on the same pair.

Each of the 14 primary receptors produced at least one saved pose below 2.0 Å under the unified chemically mapped CalcRMS table. That result is a search-coverage check. AChE 4EY7 and TYK2 3LXP saved 8 poses; the other primary receptors saved 9. EGFR 3POZ top-1 was 9.505 Å, with lowest saved-pose RMSD 0.760 Å. BChE, mTOR, JAK2, PPARG, and PPARA also failed the 2 Å top-1 cutoff; PPARA 6LXA top-3 was 7.848 Å (Figure S4; Table S2c).

![Figure 5](../figures/jcim_article/Fig5_computational_realization.png)

**Figure 5.** Computational-implementation sensitivity. (A) Independent GNINA versus Vina: filled, Vina; open, GNINA; blue, directional \(\mathrm{summary}_{\min}\); orange, dual-versus-neither; gray lines join the two tasks for one engine and are not confidence intervals. (B) PIK3CA/mTOR receptor substitution; error bars are ligand-level bootstrap 95% confidence intervals. (C) \(\mathrm{summary}_{\min}\) range, median, and production seed across five Vina seeds on the eight pairs; this panel does not show the task difference.

### 3.5 Label and sample-composition sensitivity

Panels drawn from the strict 6.5/5.5 candidate pool kept the same main class composition under several thresholds, so the corresponding AUROCs changed little. Class composition on EGFR/HER2 and PIK3CA/mTOR shifted more with threshold, and their directional estimates also changed (Figure 6A; Table S3). The 2026-08-26 API snapshot on EGFR/HER2, AChE/BChE, and PIK3CA/mTOR, and the separate ChEMBL 37 dump check on the other five pairs, did not change the main directional judgments; agreement counts and max-versus-median values are in Table S3. Those two sources are not interchangeable.

Unused-pool holdouts built from remaining candidates, after excluding main-panel molecules, showed some dependence on sample composition. AChE/BChE, PIK3CA/mTOR, and JAK1/JAK2 stayed close to the main evaluation. JAK1/TYK2 increased. F2/F10 and PPARA/PPARD remained low. PPARG/PPARA fell from 0.649 to 0.535 [0.350, 0.717] (Figure 4B; Table S7). EGFR/HER2 has no holdout. For the fixed-score difference on target A, EGFR/HER2 scaffold-cluster and document-cluster intervals were [0.168, 0.562] and [0.083, 0.529], both excluding 0. JAK1/TYK2 scaffold-cluster was [0.234, 0.633] (excludes 0) and document-cluster was [−0.034, 0.682] (includes 0) (Figure 6B; Table S10).

### 3.6 Availability of external evaluation data

BindingDB[16] and PubChem were surveyed for all eight pairs under the strict 6.5/5.5 supply definition (Table S11a). After the separate \(\theta=6.0\) independence filter removed shared literature, duplicate structures, and molecules with ECFP4 Tanimoto similarity \(\geq 0.70\) to the development set, no pair retained at least 20 dual, A-only, and B-only ligands from at least three independent sources per class. No pair met the independent external-evaluation eligibility criteria; therefore, no external docking set was formed (Figure 6C,D; Table S11b). In the internal 2018 publication-year subset, only JAK1/TYK2 and JAK1/JAK2 met the two-direction reporting gate (Table S12). This publication-year analysis is an internal sensitivity analysis, not external validation.

![Figure 6](../figures/jcim_article/Fig6_evidence_boundary.png)

**Figure 6.** Evidence boundaries. (A) Activity-threshold sensitivity; (B) ligand-, scaffold-cluster, and document-cluster resampling of the pocket A dual-versus-neither minus dual-versus-B-only difference; (C) BindingDB class counts after filtering, with color saturating at n = 20 per class; (D) independent-source counts after filtering, with color saturating at 3 sources per class. \(\dagger\) marks a class with n < 10.

## 4. Discussion

### 4.1 Experimental-state definition and evaluation tasks

What dual-target docking evaluation answers depends on the control ligands and their experimental states. Dual-versus-neither tests whether dual-active ligands can be separated from compounds that are low-activity at both targets. Dual-versus-A-only and dual-versus-B-only test whether docking scores can separate dual ligands from selectives that remain active at one target. These comparisons answer different questions and should not be read as two difficulty levels of the same task.

A fixed score channel removes score aggregation as a confounder. Replacing B-only with neither also replaces the specific molecules, so scaffolds, size, charge, and experimental provenance may change with the class. The analysis therefore reflects a joint change in evaluation class and compound set rather than experimental activity state alone. EGFR/HER2 and JAK1/TYK2 still showed clear differences when the same target score was held fixed and only the control class and its compounds were changed. Aggregation of two-pocket scores therefore cannot explain those cases. The magnitude was not uniform across other pairs and directions. Large fixed-score differences occurred both in the supply-limited EGFR/HER2 panel and in JAK1/TYK2, which satisfied the standard supply screen. F2/F10 ranked below chance on both directional arms, further illustrating that even the direction of the docking signal can be pair dependent. Selective ligands are not a single, uniformly harder control class; the effect of control-class choice still depends on the pair and direction. Building on Zhou and coworkers,[9] the present study treats the two single-target selective states as corresponding directional tasks. It then examines evaluation-design change at a fixed score channel, and uses ligand baselines and pocket correspondence to bound how far the discrimination can be explained.

Comparisons across DUD-E-, LIT-PCBA-, and related benchmark settings have shown that control composition can change virtual-screening evaluation,[5–7,12,13,20] and Kinase-Bench addresses kinase-selective enrichment.[21] The present study asks a different question: among compounds already measured at both ends, how changing the experimental-state comparison changes the interpretation of docking discrimination. The main inferences are within-pair fixed-score comparisons, not an eight-pair leaderboard.

### 4.2 Ligand chemistry and apparent sources of discrimination

Retrospective activity classes can themselves carry chemical structure differences. Simple physicochemical descriptors already separated some experimental states on AChE/BChE. Receptor-free ECFP4 carried substantial class information on several directions under scaffold-grouped cross-validation. Adding the corresponding docking score changed AUROC by at most 0.023 across 16 directions. Docking scores may still contain receptor-related information. A high retrospective AUROC nevertheless cannot, by its numerical value alone, be attributed to three-dimensional complementarity. Whether docking adds discrimination still needs comparison with a receptor-free baseline under the same split.

Such class-wise chemical differences may reflect how public bioactivity data are generated. Database compounds often come from medicinal-chemistry series around a limited set of lead scaffolds. Different activity classes can therefore differ in scaffold, substitution pattern, and physicochemical distribution at the same time. Some dual ligands may also arise from linked, fused, or otherwise combined pharmacophore designs,[1,2] which can systematically shift molecular size, polarity, flexibility, and local structural motifs relative to single-target series. Ligand architecture was not explicitly annotated here, so the ligand-only signal cannot be assigned to any one design mechanism. Figure 3 therefore supports the presence of receptor-free class information; it does not establish that ECFP4 is a formally better predictor than docking. Caba and coworkers examined protein–ligand representations together with ligand-only Morgan fingerprints in structure-based virtual screening. Adding Morgan fingerprints improved several models but did not further improve the best PLEC model.[22] Ligand-only information can likewise carry a strong classification signal in retrospective structure-based screening. Simple physicochemical matching can reduce some low-dimensional property differences, but it cannot fully control high-dimensional differences among scaffolds and series.

### 4.3 Structural attribution and pocket correspondence

Score-channel exchange, without redocking, provides a more direct test of structural correspondence. On the main panels, only EGFR/HER2 and AChE/BChE had positive matched-minus-mismatched intervals that excluded 0. The seven internal holdouts did not consistently reproduce a matched-pocket advantage. On homologous targets, a mismatched pocket can still produce related scores, so the control is an imperfect specificity control, particularly for homologous targets. A high directional AUROC on one panel is therefore not enough to attribute that discrimination to the corresponding pocket.

Kinase cross-docking studies have shown that receptor-conformer choice can affect pose recovery and virtual-screening performance.[14] On PIK3CA/mTOR, directional AUROC also changed after PIK3CA or mTOR receptors were replaced. Because this analysis systematically covered only one primary pair, the result shows receptor-structure sensitivity in that system rather than a general dual-target receptor dependence.

Cognate redocking showed that generating a near-native pose does not mean that pose receives the highest score. Lowest RMSD therefore mainly indicates search coverage. It does not evaluate whether those poses are ranked correctly. Changing pose generation and scoring changed specific AUROC values. For EGFR/HER2 and JAK1/TYK2, similar task differences were observed under independent GNINA pose generation, indicating that these two observations were not specific to the primary Vina pose-generation workflow. RTMScore, GNINA CNN rescoring, and receptor substitution nevertheless showed that specific directional AUROCs are sensitive to scoring and structural realization. Independent GNINA tests whether an evaluation-design difference can reappear under another pose-generation pipeline, not whether a matched-pocket advantage is present.

### 4.4 Practical role and evidence boundaries of dual-target docking

Under the retrospective conditions examined here, docking scores should be treated as ranking or screening signals. They are not, by themselves, evidence of true dual-target activity or selectivity. Favorable scores in both pockets mean only that the molecule received favorable predicted scores at both targets under the present computational settings. They do not confirm experimental dual activity or selectivity.

Prospective studies have shown that multitarget docking and generative design can produce experimentally confirmed multitarget candidates.[18,19] Those studies evaluate the ability to find active molecules from a candidate library, whereas the present analysis evaluates retrospective discrimination among already labeled experimental states. Retrospective AUROC therefore cannot replace prospective experimental hit rates, and it should not be used to dismiss docking for library compression. Likewise, a high retrospective AUROC cannot by itself show that a method will stably find dual-target molecules with a desired activity balance.

For retrospective dual-target benchmarks of this type, we recommend reporting dual-versus-A-only and dual-versus-B-only separately. Ligand-only baselines, matched versus mismatched pocket comparisons, and receptor or scoring sensitivity analyses are then needed to judge whether an apparent AUROC has a stable structural explanation. Direct experimental measurement is still required for dual-target activity, activity balance, and mechanism.

### 4.5 Scope and data limitations

The main limitations arise from heterogeneous activity data, nonuniform panel construction, and partial coverage of alternative docking settings.

The main data limitation is scarce and heterogeneous experimental measurements at both targets. Four-state classification requires both-end assays of the same ligand. Endpoint types remain heterogeneous despite pChEMBL scaling: IC50, Ki, Kd, EC50, and Potency are mapped onto one negative-log axis and are not biologically equivalent. Taking the maximum pChEMBL further favors the strongest recorded value. Maximum-versus-median and high-confidence field screens are separate checks with different coverage (Results 3.5; Table S3) and only partly mitigate that bias. The resulting classes are database-derived retrospective states, not measurements from a single assay system. Under the sources, filters, and gates used here, BindingDB did not form an external evaluation set. In the 2018 publication-year subset, only JAK1/TYK2 and JAK1/JAK2 met the reporting gate.

Panel construction was not fully uniform across target pairs. EGFR/HER2 and PIK3CA/mTOR were sampled from \(\theta=6.0\) pools, whereas the other six pairs were sampled from strict 6.5/5.5 pools; scaffold caps also differed across panels. Absolute AUROCs are therefore not intended for quantitative ranking across target pairs. The main inferences instead come from within-pair comparisons that hold the score channel fixed. The eight pairs were used to examine consistency across systems and should not be treated as independent statistical replicates: JAK1 and PPARA were each used twice, and some pairs belong to related kinase or nuclear-receptor families. A separate sample-size scenario analysis is provided in Figure S5; because its percentile interval showed imperfect nominal coverage, it was treated only as a diagnostic and not as a power estimate (Table S10).

At the computational layer, the study relied mainly on one representative crystal structure per target and AutoDock Vina as the primary docking and scoring method. Cognate redocking showed coverage of near-native poses. That check mainly reflects search coverage and cannot replace evaluation of score ranking or screening performance. Independent GNINA pose generation, alternative receptor conformers, and other scoring functions covered only some systems. Ligand preparation did not systematically enumerate protonation states, tautomers, or conformational ensembles. The specific numerical values observed here should therefore not be extrapolated as general dual-target docking performance across other protein families and broader chemical space.

Stronger evidence for this type of evaluation will require paired activity measurements of the same compounds at both targets under more uniform experimental conditions, together with independent or prospective tests.

## 5. Conclusions

This study examined how experimental-state definition affects retrospective dual-target docking evaluation across eight human target pairs with experimental measurements at both targets. Holding the same target score channel fixed, changing the control class and its compounds can change the performance judgment in a pair- and direction-dependent way. Discriminating dual-active from both-end low-activity ligands does not establish discrimination against single-target selectives.

Some retrospective discrimination can be obtained from ligand chemistry that does not use receptor structure. Under the present ECFP4 model and scaffold-grouped splits, adding a docking score did not produce a consistent AUROC gain across directions. Across the seven internal holdouts, matched-pocket advantages were not consistently reproduced. A high retrospective AUROC alone therefore does not establish that the observed discrimination arises from three-dimensional information specific to the corresponding pocket. For retrospective dual-target benchmarks of this type, both directional selectivity tasks should be reported and interpreted with ligand-chemistry and pocket-correspondence controls. The resulting metrics describe discrimination among already labeled experimental states. They are not equivalent to prospective virtual-screening hit finding. Direct experiment remains required for dual-target activity and selectivity.

## Data and Software Availability

All panel definitions, experimental-state labels, receptor and docking-box specifications, per-ligand score tables, analysis scripts, and figure-generation code are available in the `Dual_Target_Docking` directory of the public repository at https://github.com/1280602962-debug/gwj260531. The reported statistical analyses can be reproduced from the deposited score and metadata tables without rerunning the full docking campaigns. Typeset Supporting Information contains Tables S1–S14.

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

(17) García-Ortegón, M.; Simm, G. N. C.; Tripp, A. J.; Hernández-Lobato, J. M.; Bender, A.; Bacallado, S. DOCKSTRING: Easy Molecular Docking Yields Better Benchmarks for Ligand Design. *J. Chem. Inf. Model.* **2022**, *62*, 3486–3502. DOI: 10.1021/acs.jcim.1c01334.

(18) Wu, Y.; Vigneron, S.; Braz, J.; Srinivasan, K.; Fink, E. A.; Huang, X.-P.; Xu, X.; Huebner, H.; Kim, J. Y.; Wang, J.; Pfeiffer, T.; Sakamoto, K.; Moroz, Y. S.; Radchenko, D. S.; Rodriguiz, R. M.; Irwin, J. J.; Gmeiner, P.; Billesboelle, C.; Roth, B. L.; Basbaum, A. I.; Manglik, A.; Wetsel, W. C.; Shoichet, B. K. Large Library Docking for Polypharmacology. *J. Med. Chem.* **2026**, *69*, 6210–6229. DOI: 10.1021/acs.jmedchem.5c03810.

(19) Munson, B. P.; Chen, M.; Bogosian, A.; Kreisberg, J. F.; Licon, K.; Kuenzi, B. M.; Ideker, T. De novo generation of multi-target compounds using deep generative chemistry. *Nat. Commun.* **2024**, *15*, 3636. DOI: 10.1038/s41467-024-47120-y.

(20) Gu, S.; Shen, C.; Zhang, X.; Sun, H.; Cai, H.; Luo, H.; Zhao, H.; Liu, B.; Du, H.; Zhao, Y.; Fu, C.; Zhai, S.; Deng, Y.; Liu, H.; Hou, T.; Kang, Y. Benchmarking AI-powered docking methods from the perspective of virtual screening. *Nat. Mach. Intell.* **2025**, *7*, 509–520. DOI: 10.1038/s42256-025-00993-0.

(21) Wei, T.-H.; Zhou, S.-S.; Jing, X.-L.; Liu, J.-C.; Sun, M.; Zhao, Z.-H.; Li, Q.-Q.; Wang, Z.-X.; Yang, J.; Zhou, Y.; Wang, X.; Ling, C.-X.; Ding, N.; Xue, X.; Yu, Y.-C.; Wang, X.-L.; Yin, X.-Y.; Sun, S.-L.; Cao, P.; Li, N.-G.; Shi, Z.-H. Kinase-Bench: Comprehensive Benchmarking Tools and Guidance for Achieving Selectivity in Kinase Drug Discovery. *J. Chem. Inf. Model.* **2024**, *64*, 9528–9550. DOI: 10.1021/acs.jcim.4c01830.

(22) Caba, K.; Tran-Nguyen, V.-K.; Rahman, T.; Ballester, P. J. Comprehensive machine learning boosts structure-based virtual screening for PARP1 inhibitors. *J. Cheminform.* **2024**, *16*, 40. DOI: 10.1186/s13321-024-00832-1.
