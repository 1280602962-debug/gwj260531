# Methods (English working draft)

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
