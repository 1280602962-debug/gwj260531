# Evaluating Dual-Target Molecular Docking: Experimental-State Comparisons across Multiple Target Pairs and Sources of Discrimination

## Abstract

Molecular docking is often used to rank candidates in dual-target virtual screening. In retrospective evaluation, the control set may be compounds below the activity threshold at both targets, or single-target selectives that remain active at one target. We defined a four-state evaluation on eight human target pairs (dual, A-only, B-only, and neither) and compared those states after holding the same target score channel fixed. Changing the control class and its compounds altered apparent docking discrimination, and the change varied by pair and direction. On EGFR/HER2, the EGFR-pocket AUROC for dual-versus-B-only was 0.324 and rose to 0.786 against neither (difference 0.462 [0.260, 0.641]); JAK1/TYK2 showed a similar difference (0.444 [0.261, 0.630]). Effects on the other pairs were smaller or less certain. Receptor-free ECFP4 carried substantial class information on several directions. Under the same scaffold-grouped logistic-regression setting, adding the matched-pocket docking score changed AUROC by at most 0.023. Only one of eight main-panel matched-minus-mismatched 95% intervals excluded zero (AChE/BChE; EGFR/HER2 included zero), and none of the seven available holdout intervals did; these intervals were not multiplicity-adjusted. These results support a narrower interpretation: among compounds with experimental measurements at both targets, strong dual-versus-neither discrimination does not establish strong discrimination against single-target selectives. Retrospective dual-target evaluation should report both directional selectivity comparisons, together with ligand-chemistry and pocket-correspondence controls.

**Keywords:** dual-target docking; experimental state; selectivity; AutoDock Vina; virtual screening evaluation

## 1. Introduction

Multitarget drug design uses one molecule to modulate two or more disease-related targets. Relative to combination therapy, a single multitarget ligand can reduce pharmacokinetic mismatch and formulation complexity.[1,2] Generative design and ultra-large library docking have already produced experimentally confirmed dual-target candidates. POLYGON tested 32 MEK1/mTOR compounds.[19] Ultra-large library docking has also found multitarget ligands among hundreds of millions of molecules, and experimental structures support some predicted binding modes.[18] Such work shows that computational methods can find dual-active molecules in a candidate library. Structure-based dual-target generative models can likewise design against two pockets at once.[10,11] How the interpretation of docking scores changes when benchmark ligands have experimental measurements at both targets requires a separate evaluation.

Molecular docking searches ligand poses in a binding site and scores them. It remains a common ranking method in virtual screening.[3,4] Many benchmarks treat known actives versus artificial decoys as a binary task, as in DUD and DUD-E. LIT-PCBA instead uses high-throughput assay actives and inactives. DOCKSTRING provides a reproducible docking workflow for comparing ligand-design methods.[5–7,17] The source and definition of control compounds can change virtual-screening results and even reorder methods.[12,13,20] Rescoring hits from ultra-large docking screens is also not always stable.[15] Control-set composition is therefore part of virtual-screening evaluation itself.

When the same ligand has usable measurements at both targets, an activity threshold assigns it to one of four states: dual-active at both ends; selective for target A or target B; or below threshold at both ends. These states are labeled dual, A-only, B-only, and neither. Dual differs from A-only or B-only at only one target, and from neither at both, so dual-versus-A-only, dual-versus-B-only, and dual-versus-neither answer different questions. Structure-based dual-target screening usually docks each candidate into both pockets and prefers molecules that score favorably on both. Zhou and coworkers evaluated docking-based dual-kinase virtual screening and found that single-target inhibitors were an important source of predicted dual-inhibitor false positives.[9] Their analysis covered several false-positive types, but it did not explicitly formulate the two selective states as separate target-aligned directional tasks. Kinase-Bench uses experimental selectives to test structural discrimination between related kinases, with the aim of identifying inhibitors selective for one kinase. Distinguishing dual from A-only and B-only ligands among compounds with experimental measurements at both targets is a different problem.[21]

Prior work has shown that single-target selectives affect dual-target screening. What remains less clear is how the evaluation task itself changes the interpretation of docking performance when experimental measurements at both targets are already available. We asked three questions. First, after the score channel is held fixed, does changing the experimental-state comparison change the performance judgment? Second, what discrimination can ligand chemistry obtain without receptor structure? Third, is the observed docking discrimination consistent with the corresponding pocket, and is that correspondence stable under the available controls?

## 2. Methods

### 2.1 Study design

Dual-target virtual screening typically docks each candidate into two binding sites and ranks molecules that score favorably in both pockets. Retrospective evaluation of that ranking depends on the experimental control set. When ligands have usable measurements at both targets, an activity threshold assigns four experimental states: dual, A-only, B-only, and neither. Dual-versus-neither is the conventional virtual-screening contrast against compounds below threshold at both ends. Dual-versus-A-only and dual-versus-B-only test whether the same docking scores can reject single-target-active ligands, a documented source of dual-target false positives. Dual was the positive class. Dual and A-only ligands both meet the activity threshold at target A; their experimental difference lies at target B, so pocket B scores were used. Dual and B-only ligands differ at target A, so pocket A scores were used. The four experimental states and the two directional docking comparisons are shown in Figure 1A,B.

The two directional AUROCs at a fixed target score were the primary docking endpoints. Their lower value, \(\mathrm{summary}_{\min}=\min\{\mathrm{AUROC}_{D/A}(B),\mathrm{AUROC}_{D/B}(A)\}\), is a descriptive weaker-arm summary, not a separate overall-performance endpoint. A fixed-pocket dual-versus-neither comparison asked how AUROC changed when the control class no longer contained single-target-active ligands (Table S4). Two-pocket mean scores, \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\), provided the conventional dual-target ranking readout against neither (Table 3). Ligand-chemistry baselines asked whether receptor-free properties already separated the same experimental classes. Matched versus mismatched pocket scores asked whether discrimination followed the assigned docking site. Receptor substitution, independent pose generation, label sensitivity, and unused-pool holdouts assessed the docking protocol and the evaluation panels. BindingDB and PubChem were searched to test whether an independent external docking set could be formed.

### 2.2 Bioactivity processing and experimental-state definition

Experimental activities were taken from ChEMBL 37. Table 2 primary labels used the maximum pChEMBL for each ligand–target pair. Quantitative records with IC50, Ki, Kd, EC50, Potency, IC50app, or Ki app were retained when a pChEMBL value was available and were used on the standard negative-log molar scale. When several valid records existed for the same ligand–target pair, the primary analysis used the maximum pChEMBL as the representative value. That choice favors the strongest recorded value and can be affected by a single high reading, assay condition, or endpoint type. Salt and mixture records were split by connected fragment, retaining the organic fragment with the most heavy atoms. Only ligands with usable experimental measurements at both targets entered four-state classification. Absence of a measurement at either end was not treated as low activity at that target. The pair-universe census used the same ChEMBL 37 records.

For any pair A/B, an activity threshold \(\theta\) defined four classes: dual, \(p_{\mathrm{A}}\geq\theta\) and \(p_{\mathrm{B}}\geq\theta\); A-only, \(p_{\mathrm{A}}\geq\theta\) and \(p_{\mathrm{B}}<\theta\); B-only, \(p_{\mathrm{A}}<\theta\) and \(p_{\mathrm{B}}\geq\theta\); neither, both \(<\theta\). Primary analysis used \(\theta=6.0\) throughout. A-only and B-only denote experimental state relative to the current pair and threshold, not proteome-wide selectivity. Dual means only that both study targets meet the threshold; it does not imply activity balance, cellular effect, or therapeutic value. Later references to single-target selectives use this relative definition.

Bidirectional selective supply for candidate pairs was also counted under a strict 6.5/5.5 separation rule. The active end had pChEMBL \(\geq 6.5\) and the low-activity end \(\leq 5.5\). Ligands in the 5.5–6.5 gray zone were excluded from that strict class. The strict rule was used for candidate-pair supply assessment; primary analysis used \(\theta=6.0\).

Activity-processing checks held panel membership and docking scores unchanged (Table S3). Maximum-versus-median aggregation was repeated on the same ChEMBL 37 records for all eight pairs, so that only the experimental labels changed. That check does not replace Table 2.

### 2.3 Pair screening and panel construction

**Pair-supply screen.** Candidate pairs were required to support both four-state experimental evaluation and noncovalent small-molecule docking. ChEMBL 37 contained 5,869 human, single-component SINGLE PROTEIN targets, of which 4,672 had at least one qualifying molecule. The pair universe comprised the 10,911,456 unordered pairs among those 4,672 targets; self-pairs were excluded. A pair entered the census if at least one ligand had quantitative activity at both ends (2,164,618 pairs). Subsequent gates used the same ChEMBL 37 records. Requiring at least 10 ligands measured at both ends left 63,790 pairs. Requiring at least 10 dual, A-only, and B-only ligands at \(\theta=6.0\) left 5,253 pairs. A strict 6.5/5.5 bidirectional-selective rule requiring at least 50 A-only ligands and 50 B-only ligands left 86 pairs. Removing qHTS hub proteins, CYP ADME panels, and zinc-dependent metal enzymes left 26 pairs. Requiring at least five human holo structures per end (\(\leq 3.5\) Å, at least one non-polymer ligand) left 19 pairs. A molecular-identity and size filter required a valid molecular structure, free-base molecular weight of 150–750 Da, 10–60 heavy atoms, and no listed metal element. Pairs retaining at least 50 ligands in each strict selective class after this filter numbered 17, spanning 12 target systems. The final manual screen required a defined noncovalent small-molecule site compatible with the docking protocol. Membrane GPCRs, transporters, selected reversible-covalent complexes, and an ambiguously assigned domain were outside this scope; F2/PRSS1 was excluded because PRSS1 served as an antitarget in the study's pharmacological framing.

After these scope and docking-protocol checks, seven pairs remained: JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, and PPARA/PPARD. CTSK/CTSS reached the earlier supply and structure gates but was excluded because both holos are reversible-covalent and would have required a different docking treatment. EGFR/HER2 was retained as a supply-limited docking case. It did not meet the strict 6.5/5.5 selective-supply criterion (minimum selective-class count = 7). It had suitable human holo structures, a cognate-defined docking site, and sufficient dual, A-only, and B-only ligands at \(\theta=6.0\) for directional evaluation. These eight pairs are the primary docking evaluation set (Table 1). The executable filters, manual decisions, and pair-level retention or exclusion reasons are listed in Table S9.

**Panel construction.** Sample-supply size and structure feasibility differed across pairs, so some pairs used different candidate-pool rules, class quotas, or scaffold caps (Table 1). EGFR/HER2 and PIK3CA/mTOR were drawn from the \(\theta=6.0\) pool. The other six pairs were drawn from the strict 6.5/5.5 pool under quota sampling. EGFR/HER2 and PIK3CA/mTOR limited a single Bemis–Murcko scaffold to at most 5 and 2 repeats within a class, respectively; the remaining pairs had no additional scaffold cap.

**Unified primary analysis.** All primary AUROCs used experimental states reassigned at \(\theta=6.0\). Only ligands with valid both-end Vina scores entered the primary directional table, so n_scored can fall below n_panel. Each directional AUROC uses the corresponding pocket score.

**Table 1.** Composition and main docking settings of the dual-target evaluation panels. Quotas are construction targets (dual / A-only / B-only / neither). n_scored is the dual / A-only / B-only count with valid both-end Vina scores that entered the primary directional AUROCs. Receptor resolutions are in Table S2. The PG08-NL peptide was retained in the PPARG pocket 9V8H.

| Pair | Candidate pool | Quota (D / A / B / N) | PDB (A / B) | n_scored (dual / A-only / B-only) | Vina exhaustiveness |
|------|----------------|----------------------:|-------------|------------------------------------:|--------------------:|
| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 38 / 32 | 8 |
| JAK1/JAK2 | strict 6.5/5.5 | 32 / 32 / 32 / 14 | 6N7A / 8BXH | 32 / 32 / 32 | 8 |
| JAK1/TYK2 | strict 6.5/5.5 | 32 / 32 / 32 / 14 | 6N7A / 3LXP | 31 / 32 / 32 | 8 |
| PIK3CA/mTOR | θ = 6.0 | 18 / 14 / 12 / 4 | 4L23 / 4JT6 | 18 / 14 / 12 | 16 |
| AChE/BChE | strict 6.5/5.5 | 28 / 28 / 28 / 16 | 4EY7 / 4BDS | 27 / 26 / 28 | 8 |
| F2/F10 | strict 6.5/5.5 | 32 / 32 / 32 / 14 | 4UDW / 2JKH | 31 / 32 / 32 | 8 |
| PPARG/PPARA | strict 6.5/5.5 | 32 / 32 / 32 / 14 | 9V8H / 6LXA | 32 / 31 / 32 | 8 |
| PPARA/PPARD | strict 6.5/5.5 | 32 / 32 / 32 / 14 | 6LXA / 5U3Q | 32 / 32 / 32 | 8 |

### 2.4 Receptor and ligand preparation and docking

#### 2.4.1 Receptor preparation and docking-site definition

The preparation scripts record receptor-specific chain selection, alternate-location handling, and residue-template choices in Table S2. The deposited preparation records do not specify a common pH-dependent protonation or missing-loop reconstruction procedure.

Human experimental structures (PDB IDs) for each target are listed in Table 1. Cognate ligands and docking boxes are in Table S2. Each cognate ligand defined the docking site of its receptor. JAK1 6N7A was used for both JAK1/TYK2 and JAK1/JAK2. PPARA 6LXA was used for both PPARG/PPARA and PPARA/PPARD. TYK2 3LXP used the JH1 ATP site. In PPARG 9V8H, the bound PG08-NL peptide was retained, and only the cognate small molecule and crystal waters were removed.

Protein identity, species, binding site, and cognate-ligand placement were verified before receptor inclusion. The initial docking box was the Cartesian range of cognate heavy atoms, expanded by 5 Å on each side along \(x\), \(y\), and \(z\). Any edge that remained shorter than 20 Å was extended to 20 Å. Crystal waters and the cognate ligand used to define the site were then removed. Meeko prepared each receptor as PDBQT. Alternate locations without an altLoc flag, or marked A, were retained. Box coordinates are in Table S2.

#### 2.4.2 Ligand preparation

Processed ChEMBL SMILES from section 2.2 were the ligand input. RDKit added hydrogens, and ETKDGv3 generated one three-dimensional starting conformer. Local geometry was then optimized with the MMFF force field for at most 200 iterations. Optimized ligands were converted to PDBQT with Meeko. Three-dimensional embedding used a fixed random seed (Table S1). Protonation states and tautomers were not enumerated systematically. Unspecified stereochemistry followed RDKit’s default handling of the given SMILES.

#### 2.4.3 AutoDock Vina docking and scoring

Primary docking used AutoDock Vina 1.2.7 with the default Vina scoring function, the engine used for the dual-target ranking evaluation. At most nine poses were retained for each ligand–receptor pair, with `energy_range` set to 3 kcal mol\(^{-1}\). AChE 4EY7 and TYK2 3LXP saved eight poses. PIK3CA/mTOR used exhaustiveness = 16 to obtain a search setting that met the near-native redocking coverage requirement. The other pairs used exhaustiveness = 8. All primary Vina analyses used the rank-1 pose affinity reported as `REMARK VINA RESULT` in the first saved model. That field is the production `*_affinity` column used for Table 2. Primary Vina settings are in Table S1.

#### 2.4.4 Cognate-ligand redocking

Cognate-ligand redocking assessed recovery of the experimental bound conformation. The original outputs report the top-ranked RMSD, the lowest RMSD among the first three ranked poses (top-3), and the lowest RMSD among all saved poses. The last metric assessed search coverage at a 2.0 Å cutoff and is not a claim that the top-ranked pose was near-native; top-1 assessed pose ranking. All 14 primary receptors were rescored from saved poses with one chemically mapped CalcRMS recipe (Table S2), without redocking. AChE 4EY7 and TYK2 3LXP deposited 8 poses; the other primary receptors deposited 9. For PIK3CA/mTOR, graph-automorphism CalcRMS was used because the prepared ligand was not in the crystal frame. EGFR 3POZ and HER2 3RCD were redocked under the canonical cognate heavy-atom box.[8]

#### 2.4.5 Alternative scoring and independent docking

The independent GNINA analyses were targeted sensitivity checks on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2, selected to contrast prominent control-class effects and the PIK3CA/mTOR protocol-sensitivity case. They were not a prespecified eight-pair engine comparison.

RTMScore and GNINA 1.3.2 CNN were used to rescore all available Vina poses. RTMScore took the highest score among saved poses as the ligand-level result. GNINA used CNN affinity as the main CNN rescoring readout, with CNNscore as a supplementary analysis. These results re-score poses already generated by Vina and are not independent pose generation.

GNINA 1.3.2 was also used to generate poses independently on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2. This analysis used the same receptors, ligands, and docking boxes as the primary analysis. Exhaustiveness values were 8, 16, and 8, respectively, with at most nine poses per ligand. The ligand-level score was the rank-1 minimizedAffinity, inverted so that higher values indicate more favorable predicted binding. Independent GNINA did not return both-end scores for every ligand. EGFR/HER2 dual-versus-neither used n_neither = 11 (EH120_109 missing) rather than the primary Vina dual-versus-neither count of 12; PIK3CA/mTOR A-only was n = 13 rather than 14; JAK1/TYK2 Dual/A-only/B-only/neither counts were 30/32/29/14 rather than 31/32/32/14. The analysis asks whether the main observation remains under another pose-generation and scoring pipeline. It is not a comparison of overall GNINA versus Vina performance (Table S7).

### 2.5 Evaluation metrics and statistical analysis

#### 2.5.1 Primary directional discrimination metrics

When dual was compared with A-only, both classes meet the activity threshold at target A, so pocket B scores were used: \(\mathrm{AUROC}_{D/A}(B)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{A\text{-}only};\;S_{B})\). When dual was compared with B-only, pocket A scores were used: \(\mathrm{AUROC}_{D/B}(A)=\mathrm{AUROC}(\mathrm{dual},\;\mathrm{B\text{-}only};\;S_{A})\). Dual was the positive class in all directional analyses.

More negative AutoDock Vina affinity indicates more favorable predicted binding. For a common AUROC direction, Vina affinity was converted to \(S_{\mathrm{Vina}}=-E_{\mathrm{Vina}}\), so higher \(S_{\mathrm{Vina}}\) is more favorable. Both directional AUROCs were the primary endpoints. Their lower value was defined as \(\mathrm{summary}_{\min}\), a descriptive weaker-arm summary rather than an overall-performance metric or a new ligand scoring function: \(\mathrm{summary}_{\min}=\min[\mathrm{AUROC}_{D/A}(B),\;\mathrm{AUROC}_{D/B}(A)]\). Because the minimum of two AUROCs is selectively downward-shifted, it is not treated as an ordinary independent endpoint.

#### 2.5.2 Conventional two-pocket ranking and control-class comparisons

Dual-target virtual screening commonly ranks candidates by combining the two pocket scores. For ligands scored at both targets, the conventional readout was the mean score \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\), with dual as the positive class and neither as the control (Table 3). Merging A-only, B-only, and neither into one control set gave a dual-versus-all-nonduals AUROC that is archived with the score tables.

To isolate the effect of the virtual-screening control class from score aggregation, the same pocket score was held fixed and AUROCs were compared under different experimental-state controls. Pocket A scores were used for dual-versus-B-only and that dual-versus-neither comparison; pocket B scores were used for dual-versus-A-only and the corresponding dual-versus-neither comparison (Table S4). Changing the control class also changes the specific compounds. The lower of the two pocket scores, \(S_{\mathrm{worst}}=\min(S_{A},S_{B})\), was used as an AND-type two-pocket filter. The median of dual \(S_{\mathrm{worst}}\) values defined the filter threshold. Dual counts, dual recall, dual precision, and retained single-target-active ligands at that threshold are reported in Table S10. The corresponding eight-pair top-10% class composition under two-pocket mean ranking is shown in Figure 2D. The cutoff is \(k=\lceil 0.10\,n\rceil\) on the full four-state panel, so pairs of different size are compared at the same screening fraction. The top-10% dual fraction is \(\mathrm{dual}/k\). Enrichment versus the panel dual base rate is \(\mathrm{EF}_{\mathrm{dual},10\%}=(\mathrm{dual}/k)/(n_{\mathrm{dual}}/n)\). The same ranking and filter rules on all eight primary pairs are in Table S10.

#### 2.5.3 Confidence intervals and resampling

Table 2 reports pointwise 95% confidence intervals for both directional AUROCs and their descriptive minimum. Primary intervals used a class-stratified nonparametric percentile bootstrap, with B = 2000 replicates and seed 20260729. For each directional AUROC, ligands were resampled with replacement within the two experimental-state classes, preserving the original class sizes. For \(\mathrm{summary}_{\min}\), dual ligands were resampled once per replicate and that same dual draw was applied to both pocket scores; A-only and B-only were resampled independently within class; the two directional AUROCs were recomputed and their minimum was taken inside the replicate. Fixed-score \(\Delta\)AUROC used the same shared dual resample for dual-versus-selective and dual-versus-neither, with the two negative classes resampled independently. Matched-versus-mismatched pocket comparisons used paired ligand resampling so that a ligand's matched and mismatched scores were drawn together. Point estimates used the full analysis sample. A ligand-level non-stratified bootstrap was retained as a sensitivity analysis and did not replace the class-stratified intervals. Cluster bootstrap by Bemis–Murcko scaffold and literature-connected document groups was a source-dependence sensitivity for the two largest fixed-score differences (Table S4) and did not replace the ligand-level primary intervals. A binormal detectable-effect simulation reused the same class-stratified shared-dual bootstrap, the current eight-pair class sizes, B = 2000, and seed 20260729; it estimates the probability that a CI excludes 0.5 under a specified true AUROC and is not observed power. Intervals were not adjusted for multiple comparisons.

When different scoring methods or matched versus mismatched pocket scores were compared, each scheme reused the same ligand resample to preserve pairing. Cluster bootstrap with Bemis–Murcko scaffold groups and literature-connected groups was used to assess scaffold and source correlation. Scaffold- and document-cluster intervals for the two largest fixed-score differences are in Table S4.

### 2.6 Baselines, controls, and sensitivity analyses

#### 2.6.1 Ligand-chemistry controls

For each directional model, the number of GroupKFold splits was the minimum of five, the number of Bemis–Murcko scaffolds, and the two class counts. Logistic regression used C = 1.0 and a maximum of 4000 iterations; GroupKFold used its default non-shuffled partitioning. Full settings are listed in Table S1.

To test whether a ligand-based ranking already separated the same experimental classes without docking, ECFP4 fingerprints (radius = 2, 2048 bits) were computed, together with molecular weight, heavy-atom count, cLogP, and TPSA. Directional AUROCs and \(\mathrm{summary}_{\min}\) were computed for each of the four single descriptors. The single descriptor with the highest \(\mathrm{summary}_{\min}\) on each pair served as the physicochemical baseline. Because that best descriptor was selected on the same panel, its difference from Vina is descriptive only (Table S5).

The descriptor was selected on the full panel and held fixed during bootstrap. Its difference from Vina and the associated interval are therefore descriptive and conditional on that selection, not a selection-adjusted test.

ECFP4, docking-only, and ECFP4+docking models used logistic regression without feature scaling. Bemis–Murcko scaffolds were the grouping variable, and out-of-fold predictions under the same GroupKFold splits were used to compute AUROC. The docking-only model used only the corresponding directional docking score. That AUROC is computed from out-of-fold predictions, whereas the primary analysis ranks the raw docking scores, so the two values are not identical. The AUROC difference between ECFP4 and ECFP4+docking measures the incremental contribution of the docking score to a ligand-based model under that shared split. The Figure 3A comparison of raw Vina ranking with scaffold-grouped ECFP4 out-of-fold predictions contrasts structure-based and ligand-based rankings; it is not a head-to-head predictive benchmark. A sensitivity analysis applied StandardScaler on each training fold of the same splits (Table S5).

#### 2.6.2 Structural attribution and score correspondence

Because \(\Delta\) compares the two lower-arm summaries, a positive value does not establish an advantage in both directions; the identity of the weaker arm can change after score-channel exchange.

Dual-target docking assigns each score to a pocket. Matched scoring used pocket B for dual-versus-A-only and pocket A for dual-versus-B-only. The mismatched control exchanged these two precomputed score channels without redocking. The analysis compared \(\mathrm{summary}_{\min}\) under the two assignments and estimated \(\Delta=\mathrm{summary}_{\min}^{\mathrm{matched}}-\mathrm{summary}_{\min}^{\mathrm{mismatched}}\) with a paired bootstrap. A positive \(\Delta\) means the weaker matched arm is higher. Whether the weaker arm switched is recorded in Table S6.

Receptor-structure sensitivity was evaluated mainly on PIK3CA/mTOR. With mTOR 4JT6 held fixed, PIK3CA 4L23 was replaced by 4JPS or 5DXT. mTOR 4JT6 was also replaced by 4JSX. Ligand sets, experimental states, score definitions, and statistics were otherwise unchanged (Table S7).

#### 2.6.3 Robustness and external-data availability

Experimental states were reassigned under alternative activity thresholds and activity-aggregation rules, with panel membership and docking scores held fixed (Table S3). For pairs with enough remaining candidates, unused-pool holdout sets were built after excluding main-panel members. Holdout ligands were drawn under a fixed sampling rule that limited over-representation of a single Bemis–Murcko scaffold. Directional AUROCs were recomputed on these different ligands to assess sensitivity to panel membership. Because the ligands still come from the same source, this analysis is an internal robustness check (Table S6). Primary Vina docking was also repeated with five fixed random seeds, and different exhaustiveness settings were compared on PIK3CA/mTOR. Complete-case counts and fixed-membership intersections are archived in the repository. A binormal detectable-effect simulation used the current eight-pair class sizes and the same class-stratified shared-dual bootstrap as Table 2.

BindingDB[16] and PubChem were searched for all eight pairs as candidate sources of an independent external docking set. An independence filter at \(\theta=6.0\) removed shared literature, duplicate structures, and molecules with ECFP4 Tanimoto \(\geq 0.70\) to the development set (Table S8). The development-molecule set included main panels, the expanded PIK3CA/mTOR PM110 panel, and internal holdouts. Dual, A-only, and B-only each needed \(n\geq 20\) with at least three sources per class. Only pairs meeting those conditions would have entered external docking evaluation.

### 2.7 Software and reproducibility

Molecular processing, fingerprints, descriptors, and statistics were implemented in Python, mainly with RDKit, NumPy, pandas, SciPy, and scikit-learn. Receptor and ligand PDBQT files were prepared with Meeko. Primary docking used AutoDock Vina 1.2.7. GNINA 1.3.2 was used for supplementary scoring and independent pose generation on selected pairs. RTMScore provided an additional rescoring of Vina poses.

Software versions, random seeds, and main computational parameters are summarized in Table S1. Docking boxes and cognate redocking results are in Table S2.

## 3. Results

### 3.1 Paired experimental supply and docking evaluation panels

Paired experimental labels that can support a four-state docking evaluation became progressively more limited as sample requirements increased. Among unordered pairs of human single-component SINGLE PROTEIN targets, 2,164,618 pairs had at least one ligand measured at both targets, and 63,790 pairs had at least 10. Requiring at least 10 dual, A-only, and B-only ligands at \(\theta=6.0\) left 5,253 pairs. Applying the strict 6.5/5.5 activity criterion left 86 pairs. Figure 1C summarizes these counts. The later structure and docking-site gates are defined in section 2.3.

After those supply, structure, and docking-compatibility checks, the primary evaluation comprised eight target pairs: EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, and PPARA/PPARD (Table 1). EGFR/HER2 and PIK3CA/mTOR were sampled from the \(\theta=6.0\) pools; the other six pairs were sampled from the strict 6.5/5.5 pools. Four-state classification requires experimental measurements of the same ligand at both targets, so paired coverage jointly limited the docking panels.

![Figure 1](../figures/jcim_article/Fig1_four_state_and_supply.png)

**Figure 1.** Four-state dual-target evaluation and data supply. (A) Four experimental states defined by threshold \(\theta\); (B) two directional tasks: dual versus A-only uses the target B score, and dual versus B-only uses the target A score, with \(\mathrm{summary}_{\min}=\min[\mathrm{AUROC}_{D/A}(B),\;\mathrm{AUROC}_{D/B}(A)]\); (C) paired experimental coverage decreases under four-state requirements. The primary evaluation comprised eight target pairs selected using the panel-construction and structural criteria in Table 1. That eight-pair set is not a direct continuation of the census counts.

### 3.2 Directional docking performance and candidate-ranking consequence

With the same target-specific docking score held fixed, AUROC differences associated with changing the experimental-state control varied across pairs and directions. EGFR/HER2 and JAK1/TYK2 showed the largest fixed-score differences. Using the EGFR-pocket score, the EGFR/HER2 dual-versus-B-only AUROC was 0.324 and rose to 0.786 against neither (difference 0.462 [0.260, 0.641]). Using the JAK1-pocket score, replacing B-only with neither increased the AUROC by 0.444 [0.261, 0.630]. For the remaining pairs, the differences were smaller or their confidence intervals included zero (Figure 2A; Table S4). Cluster resampling is reported in section 3.6.

Under unified \(\theta=6.0\) labels, the two directional AUROCs ranged from 0.324 to 0.728, whereas the descriptive weaker-arm \(\mathrm{summary}_{\min}\) ranged from 0.324 to 0.692. The PPARG/PPARA \(\mathrm{summary}_{\min}\) interval lay entirely above 0.5 (0.649 [0.511, 0.746]), the EGFR/HER2 and F2/F10 intervals lay entirely below 0.5 (0.324 [0.192, 0.463] and 0.345 [0.216, 0.482]), and the remaining 5 pairs crossed 0.5 (Figure 2B; Table 2). These values are pair-specific and are not an eight-pair ranking. Under a binormal simulation that reused these class sizes and the Table 2 bootstrap, the probability that a `summary_min` CI excludes 0.5 was at most 0.068 when the true weaker-arm AUROC was 0.60, and at least 0.817 at 0.75 except for PIK3CA/mTOR (0.453). The simulation is not observed power (Table S4).

**Table 2.** Pocket-matched directional AUROC on the eight primary pairs (Vina; unified \(\theta=6.0\)). Class sizes are n_scored (dual / A-only / B-only). Physicochemical descriptor baselines are in Table S5.

| Pair | n_scored (dual / A-only / B-only) | dual vs A-only (pocket B) [95% CI] | dual vs B-only (pocket A) [95% CI] | summary_min [95% CI] |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|
| EGFR/HER2 | 28 / 38 / 32 | 0.661 [0.524, 0.790] | 0.324 [0.192, 0.463] | 0.324 [0.192, 0.463] |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 [0.448, 0.723] | 0.728 [0.602, 0.842] | 0.588 [0.448, 0.716] |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 [0.434, 0.717] | 0.365 [0.233, 0.508] | 0.365 [0.233, 0.505] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 [0.528, 0.889] | 0.692 [0.491, 0.868] | 0.692 [0.480, 0.802] |
| AChE/BChE | 27 / 26 / 28 | 0.652 [0.492, 0.806] | 0.606 [0.443, 0.749] | 0.606 [0.439, 0.735] |
| F2/F10 | 31 / 32 / 32 | 0.413 [0.271, 0.571] | 0.345 [0.218, 0.491] | 0.345 [0.216, 0.482] |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 [0.512, 0.781] | 0.706 [0.573, 0.830] | 0.649 [0.511, 0.746] |
| PPARA/PPARD | 32 / 32 / 32 | 0.647 [0.509, 0.776] | 0.446 [0.301, 0.592] | 0.446 [0.301, 0.590] |

Under the conventional two-pocket mean ranking used in dual-target virtual screening, dual-versus-neither AUROCs ranged from 0.514 to 0.770. EGFR/HER2 reached 0.759 [0.554, 0.926] and JAK1/TYK2 reached 0.770 [0.613, 0.906], whereas their directional \(\mathrm{summary}_{\min}\) values were 0.324 and 0.365 (Figure 2C; Table 3). The PIK3CA/mTOR neither sample was \(n=4\). Because Table 3 changes both score aggregation and the control class, the influence of the control class was taken from the fixed-pocket analysis in Figure 2A.

![Figure 2](../figures/jcim_article/Fig2_negative_class_formulation.png)

**Figure 2.** Core dual-target candidate-ranking result. (A) Fixed-score \(\Delta\)AUROC (dual–neither minus dual–matched single-target-active) for eight pairs and two score directions; error bars are ligand-level bootstrap 95% confidence intervals; (B) the two directional AUROCs: dual versus A-only using pocket B, and dual versus B-only using pocket A; (C) directional \(\mathrm{summary}_{\min}\) versus two-pocket mean dual-versus-neither AUROC. Panel C is not a fixed-score comparison. The n = 4 annotation marks the PIK3CA/mTOR neither sample. (D) Eight-pair top-10% class composition under two-pocket mean Vina ranking, \(k=\lceil 0.10\,n\rceil\). Bars are 100% stacked dual / A-only / B-only / neither; \(k\) is annotated at right. \(\mathrm{EF}_{\mathrm{dual},10\%}\) is in Table 3.

**Table 3.** Candidate-ranking consequence of two-pocket mean Vina ranking (unified \(\theta=6.0\)). Dual-versus-neither uses \(S_{\mathrm{mean}}=(S_{A}+S_{B})/2\). Top 10% composition uses \(k=\lceil 0.10\,n\rceil\) on the full four-state panel. \(\mathrm{EF}_{\mathrm{dual},10\%}=(\mathrm{dual}/k)/(n_{\mathrm{dual}}/n)\). The PIK3CA/mTOR neither sample is small (n = 4).

| Pair | two-pocket mean D-vs-neither AUROC [95% CI] | n_neither | Top 10% D/A/B/N | \(\mathrm{EF}_{\mathrm{dual},10\%}\) |
|------|--------------------------------------------:|----------:|----------------:|-------------------------------------:|
| EGFR/HER2 | 0.759 [0.554, 0.926] | 12 | 1 / 5 / 5 / 0 | 0.357 |
| JAK1/JAK2 | 0.730 [0.558, 0.882] | 14 | 6 / 5 / 0 / 0 | 1.875 |
| JAK1/TYK2 | 0.770 [0.613, 0.906] | 14 | 1 / 3 / 7 / 0 | 0.320 |
| PIK3CA/mTOR | 0.514 [0.222, 0.806] | 4 | 4 / 1 / 0 / 0 | 2.133 |
| AChE/BChE | 0.649 [0.484, 0.810] | 15 | 5 / 3 / 1 / 1 | 1.778 |
| F2/F10 | 0.519 [0.341, 0.691] | 12 | 4 / 1 / 6 / 0 | 1.255 |
| PPARG/PPARA | 0.685 [0.498, 0.859] | 14 | 7 / 3 / 0 / 1 | 2.168 |
| PPARA/PPARD | 0.565 [0.373, 0.748] | 14 | 5 / 2 / 3 / 1 | 1.563 |

The same two-pocket mean ranking was applied at the top 10% of each panel, with \(k=\lceil 0.10\,n\rceil\), so that pairs of different size are compared at the same screening fraction (Figure 2D; Table 3; Table S10). Dual ligands in the top 10% ranged from 1/11 (EGFR/HER2 and JAK1/TYK2; dual fraction 0.091) to 4/5 (PIK3CA/mTOR; 0.800) and 7/11 (PPARG/PPARA; 0.636). Relative to each panel's dual base rate, \(\mathrm{EF}_{\mathrm{dual},10\%}\) ranged from 0.320 (JAK1/TYK2) and 0.357 (EGFR/HER2) to 2.133 (PIK3CA/mTOR) and 2.168 (PPARG/PPARA). EGFR/HER2 occupied 1 dual, 5 A-only, 5 B-only, and 0 neither slots (\(\mathrm{EF}_{\mathrm{dual},10\%}=0.357\)); AChE/BChE occupied 5 / 3 / 1 / 1 (\(\mathrm{EF}_{\mathrm{dual},10\%}=1.778\)). Values below 1 mean dual ligands were less common in the top 10% than in the full panel. Neither ligands were absent from the top 10% on five pairs and occupied one slot on AChE/BChE, PPARG/PPARA, and PPARA/PPARD. At least one single-target-active ligand remained in the top 10% on every pair. These counts are descriptive operating points on the retrospective panels; they are not an eight-pair ranking of docking quality and are not prospective hit rates. The AND filter on Dual+A-only+B-only is reported in Table S10.

### 3.3 Ligand-chemistry baselines and incremental docking discrimination

If dual-target docking ranks experimental classes that already differ in ligand chemistry, a receptor-free baseline can separate the same screening states. On AChE/BChE, TPSA alone gave dual-versus-A-only and dual-versus-B-only AUROCs of 0.742 and 0.801 (Figure S1; Table S5). On PIK3CA/mTOR, the best single descriptor, heavy-atom count, had \(\mathrm{summary}_{\min}\) 0.463 (Table S5). The difference between Vina and the best single descriptor varied by pair: 5 of eight 95% intervals included 0, whereas JAK1/TYK2, PIK3CA/mTOR, F2/F10 excluded 0 (Figure S1; Table S5).

Under Bemis–Murcko scaffold-grouped cross-validation, ligand-only ECFP4 models separated several directional comparisons (Figure 3A). Adding the corresponding docking score changed AUROC by at most 0.023 across the 16 directional comparisons, with no consistent direction of change (Figure 3B; Table S5). These panels compare receptor-free and structure-based rankings of the same experimental classes; they are not a formal algorithm leaderboard. In this retrospective setting, the docking score did not add a stable ranking increment beyond a ligand-based model of the same experimental classes.

![Figure 3](../figures/jcim_article/Fig3_ligand_chemistry.png)

**Figure 3.** Ligand chemistry and incremental docking information, not a head-to-head predictive benchmark. (A) Rank AUROC from raw Vina scores versus out-of-fold ECFP4 AUROC under scaffold-grouped cross-validation for eight pairs and 16 directions; (B) \(\Delta\)AUROC after adding the corresponding Vina score to ECFP4 on the same splits and ligands; lines join the two points for one direction and are not confidence intervals. AChE/BChE TPSA distributions and the single-descriptor forest are in Figure S1.

### 3.4 Pocket correspondence

If the directional AUROC is a structure-based docking result, exchanging the two precomputed pocket scores without redocking should reduce the weaker-arm summary. In the main panels, only AChE/BChE had a matched-minus-mismatched \(\mathrm{summary}_{\min}\) 95% interval that excluded zero (0.177 [0.053, 0.291]). EGFR/HER2 was 0.056 [−0.029, 0.157], which includes 0, so this pair is not interpreted as having a matched-pocket advantage. The other six main-panel intervals included zero. All seven available unused-pool internal-holdout intervals also included zero, with differences from \(-0.079\) to \(+0.150\) (Figure 4; Table S6). EGFR/HER2 did not have an unused-pool holdout. An interval that includes zero does not establish that the receptor contributes no information; it means that this control does not support a matched-pocket advantage for that pair.

![Figure 4](../figures/jcim_article/Fig4_mismatched_pocket.png)

**Figure 4.** Pocket correspondence. (A) Main-panel matched-minus-mismatched \(\Delta\mathrm{summary}_{\min}\) for eight pairs with ligand-level bootstrap 95% confidence intervals; (B) unused-pool internal holdout matched-minus-mismatched \(\Delta\mathrm{summary}_{\min}\) for seven pairs. Color distinguishes main versus holdout ligand sets, not whether an interval excludes zero. Only the AChE/BChE main-panel interval excludes zero; the EGFR/HER2 main-panel interval includes zero; all seven holdout intervals include zero. \(\dagger\) marks the missing EGFR/HER2 unused-pool holdout. Main-versus-holdout \(\mathrm{summary}_{\min}\) values are in Table S6.

### 3.5 Computational robustness

Independent GNINA 1.3.2 pose generation used the same receptors, ligands, and docking boxes on three pairs (Figure 5A; Table S7). For EGFR/HER2, the dual-versus-neither AUROC was 0.737 [0.529, 0.919] ($n_{\mathrm{neither}}=11$), whereas dual-versus-B-only was 0.265 [0.144, 0.401]. For JAK1/TYK2, dual-versus-neither was 0.705 [0.524, 0.872] and directional \(\mathrm{summary}_{\min}\) was 0.317 [0.187, 0.455]. For PIK3CA/mTOR, \(\mathrm{summary}_{\min}\) was 0.633, the weaker arm was dual-versus-A-only 0.633 [0.410, 0.769], and dual-versus-neither was 0.569 [0.236, 0.889] (\(n=18/13/12/4\)). These runs ask whether the main observation remains under another pose-generation pipeline. They are not a comparison of overall GNINA versus Vina performance.

On PIK3CA/mTOR, replacing PIK3CA 4L23 with 4JPS lowered \(\mathrm{summary}_{\min}\) from 0.692 [0.480, 0.802] to 0.486 [0.264, 0.694]. Replacement with 5DXT gave 0.505 [0.296, 0.713]. Replacing mTOR 4JT6 with 4JSX gave 0.639 [0.435, 0.783] (Figure 5B; Table S7).

Five fixed Vina random seeds produced comparatively limited numerical fluctuation relative to the larger task and receptor effects (Figure 5C). The EGFR/HER2 fixed-score task difference remained positive on all five seeds. Same-pose RTMScore and GNINA CNN rescoring of PPARG/PPARA are in Table S7. PPARG/PPARA was the only pair whose primary Vina \(\mathrm{summary}_{\min}\) interval lay entirely above 0.5 (0.649 [0.511, 0.746]); RTMScore lowered it to 0.369 [0.233, 0.475], and GNINA CNN rescoring lowered it to 0.500 [0.356, 0.623].

![Figure 5](../figures/jcim_article/Fig5_computational_realization.png)

**Figure 5.** Computational robustness of a single docking realization. (A) Independent GNINA versus Vina on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2: filled, Vina; open, GNINA; blue, directional \(\mathrm{summary}_{\min}\); orange, dual-versus-neither. Gray lines join the two tasks for one engine and are not confidence intervals. (B) PIK3CA/mTOR receptor substitution of primary 4L23 / 4JT6 by 4JPS, 5DXT, or 4JSX; error bars are ligand-level bootstrap 95% confidence intervals. (C) \(\mathrm{summary}_{\min}\) range, median, and production seed across five Vina seeds on the eight pairs. The figure asks whether the reported results depend on a single docking realization.

### 3.6 Sensitivity and external-data scope

Panels drawn from the strict 6.5/5.5 candidate pool kept the same main class composition under several thresholds, so the corresponding AUROCs changed little. Class composition on EGFR/HER2 and PIK3CA/mTOR shifted more with threshold, and their directional estimates also changed (Figure S4; Table S3). Relabeling the same ChEMBL 37 records by median rather than maximum pChEMBL flipped 6/110 EGFR/HER2 classes (primary \(\mathrm{summary}_{\min}\) 0.324), 1/96 AChE/BChE classes (CHEMBL659; 0.606 to 0.629), and 1/110 PPARA/PPARD classes (CHEMBL121; \(\mathrm{summary}_{\min}\) remained 0.446). The other five pairs kept class composition and \(\mathrm{summary}_{\min}\) point estimates (Table S3).

Unused-pool holdouts built from remaining candidates, after excluding main-panel molecules, showed some dependence on sample composition. AChE/BChE, PIK3CA/mTOR, and JAK1/JAK2 stayed close to the main evaluation. JAK1/TYK2 increased. F2/F10 and PPARA/PPARD remained low. PPARG/PPARA fell from 0.649 to 0.535 [0.360, 0.705] (Table S6). EGFR/HER2 has no holdout. The ligand-level fixed-score difference on EGFR/HER2 pocket A is 0.462 [0.260, 0.641]. After substituting the corrected-box scores into the frozen scaffold and document groupings, the EGFR/HER2 cluster intervals were [0.235, 0.665] (scaffold) and [0.125, 0.644] (document); both exclude 0. JAK1/TYK2 scaffold-cluster was [0.220, 0.631] (excludes 0) and document-cluster was [−0.034, 0.682] (includes 0) (Figure S4; Table S4).

PIK3CA/mTOR panel-size and exhaustiveness checks are in Figure S2. PM48 is the primary panel (quota n = 48, exhaustiveness = 16); PM110 is a larger protocol-sensitivity panel on the same pair. Each of the 14 primary receptors produced at least one saved pose below 2.0 Å under the unified chemically mapped CalcRMS table (Figure S3; Table S2). That result is a search-coverage check. AChE 4EY7 and TYK2 3LXP saved 8 poses; the other primary receptors saved 9. EGFR 3POZ top-1 was 1.019 Å and HER2 3RCD top-1 was 1.947 Å, both below 2 Å. BChE, mTOR, JAK2, PPARG, and PPARA failed the 2 Å top-1 cutoff; PPARA 6LXA top-3 was 7.848 Å.

BindingDB[16] and PubChem were searched for paired experimental data across all eight target pairs. The \(\theta=6.0\) independence filters removed shared literature sources, duplicate structures, and molecules with ECFP4 Tanimoto similarity \(\geq 0.70\) to the development set. No pair met the independent external-evaluation eligibility criteria of at least 20 dual, A-only, and B-only ligands from at least three independent sources per class. No external docking set was therefore formed; this search is an eligibility screen, not external validation (Figure S5; Table S8). 0/8 pairs passed the full external docking gate.


## 4. Discussion

### 4.1 The candidate-ranking problem

What dual-target docking evaluation answers depends on the control ligands and their experimental states. Dual-versus-neither tests whether dual-active ligands can be separated from compounds that are low-activity at both targets. Dual-versus-A-only and dual-versus-B-only test whether docking scores can separate dual ligands from selectives that remain active at one target. These comparisons answer different questions and should not be read as two difficulty levels of the same task.

A fixed score channel removes score aggregation as a confounder. Replacing B-only with neither also replaces the specific molecules, so scaffolds, size, charge, and experimental provenance may change with the class. The analysis therefore reflects a joint change in evaluation class and compound set rather than experimental activity state alone. EGFR/HER2 and JAK1/TYK2 still showed clear differences when the same target score was held fixed and only the control class and its compounds were changed. Aggregation of two-pocket scores therefore cannot explain those cases.

The magnitude was not uniform across other pairs and directions. Large fixed-score differences occurred both in the supply-limited EGFR/HER2 panel and in JAK1/TYK2, which satisfied the standard supply screen. F2/F10 ranked below chance on both directional arms, further illustrating that even the direction of the docking signal can be pair dependent. Selective ligands are not a single, uniformly harder control class; the effect of control-class choice still depends on the pair and direction.

Building on Zhou and coworkers,[9] this study treats the two single-target selective states as corresponding directional tasks. It then examines evaluation-design change at a fixed score channel, and uses ligand baselines and pocket correspondence to bound how far the discrimination can be explained.

Comparisons across DUD-E, LIT-PCBA, and related benchmark settings have shown that control composition can change virtual-screening evaluation.[5–7,12,13,20] Kinase-Bench addresses kinase-selective enrichment.[21] This study asks a different question: among compounds with experimental measurements at both targets, how changing the experimental-state comparison changes the interpretation of docking discrimination. The main inferences are within-pair fixed-score comparisons, not an eight-pair leaderboard.

### 4.2 Ligand chemistry and apparent sources of discrimination

Retrospective activity classes can themselves carry chemical structure differences. Simple physicochemical descriptors already separated some experimental states on AChE/BChE. Receptor-free ECFP4 carried substantial class information on several directions under scaffold-grouped cross-validation. Adding the corresponding docking score changed AUROC by at most 0.023 across 16 directions. Docking scores may still contain receptor-related information. A high retrospective AUROC nevertheless cannot, by its numerical value alone, be attributed to three-dimensional complementarity. Whether docking adds discrimination still needs comparison with a receptor-free baseline under the same split.

Such class-wise chemical differences may reflect how public bioactivity data are generated. Database compounds often come from medicinal-chemistry series around a limited set of lead scaffolds. Different activity classes can therefore differ in scaffold, substitution pattern, and physicochemical distribution at the same time. Some dual ligands may also arise from linked, fused, or otherwise combined pharmacophore designs,[1,2] which can systematically shift molecular size, polarity, flexibility, and local structural motifs relative to single-target series. Ligand architecture was not explicitly annotated here, so the ligand-only signal cannot be assigned to any one design mechanism. Figure 3 therefore supports the presence of receptor-free class information; it does not establish that ECFP4 is a formally better predictor than docking. Caba and coworkers examined protein–ligand representations together with ligand-only Morgan fingerprints in structure-based virtual screening. Adding Morgan fingerprints improved several models but did not further improve the best PLEC model.[22] Ligand-only information can likewise carry a strong classification signal in retrospective structure-based screening. Simple physicochemical matching can reduce some low-dimensional property differences, but it cannot fully control high-dimensional differences among scaffolds and series.

### 4.3 Pocket correspondence

Score-channel exchange, without redocking, provides a more direct test of structural correspondence. On the main panels, only AChE/BChE had a positive matched-minus-mismatched interval that excluded 0; the EGFR/HER2 interval included 0. The seven internal holdouts did not consistently reproduce a matched-pocket advantage. On homologous targets, a mismatched pocket can still produce related scores, so the control is an imperfect specificity control. A high directional AUROC on one panel is therefore not enough to attribute that discrimination to the corresponding pocket. An interval that includes zero does not mean that the receptor has no information.

### 4.4 Computational robustness

Kinase cross-docking studies have shown that receptor-conformer choice can affect pose recovery and virtual-screening performance.[14] On PIK3CA/mTOR, directional AUROC also changed after PIK3CA or mTOR receptors were replaced. Because this analysis systematically covered only one primary pair, the result shows receptor-structure sensitivity in that system rather than a general dual-target receptor dependence.

Cognate redocking showed that generating a near-native pose does not mean that pose receives the highest score. Lowest RMSD therefore mainly indicates search coverage. It does not evaluate whether those poses are ranked correctly. Changing pose generation and scoring changed specific AUROC values. For EGFR/HER2 and JAK1/TYK2, similar task differences were observed under independent GNINA pose generation, indicating that these two observations were not specific to the primary Vina pose-generation workflow. RTMScore, GNINA CNN rescoring, and receptor substitution nevertheless showed that specific directional AUROCs are sensitive to scoring and structural realization. Independent GNINA tests whether an evaluation-design difference can reappear under another pose-generation pipeline, not whether a matched-pocket advantage is present.

### 4.5 Scope, evidence boundaries, and prospective validation

Under the retrospective conditions examined here, docking scores should be treated as ranking or screening signals. They are not, by themselves, evidence of true dual-target activity or selectivity. Favorable scores in both pockets mean only that the molecule received favorable predicted scores at both targets under the present computational settings. They do not confirm experimental dual activity or selectivity.

Prospective studies have shown that multitarget docking and generative design can produce experimentally confirmed multitarget candidates.[18,19] Those studies evaluate the ability to find active molecules from a candidate library, whereas this study evaluates retrospective discrimination among already labeled experimental states. Retrospective AUROC therefore cannot replace prospective experimental hit rates, and it should not be used to dismiss docking for library compression. Likewise, a high retrospective AUROC cannot by itself show that a method will stably find dual-target molecules with a desired activity balance. The JAK1/TYK2 top-10% class composition comes from a fixed retrospective evaluation panel in which 109 ligands included 31 dual ligands; it is not a prospective hit rate in a commercial library. Eight-pair top-10% composition is in Figure 2D and Table 3.

For retrospective dual-target benchmarks of this type, we recommend reporting dual-versus-A-only and dual-versus-B-only separately. Ligand-only baselines, matched versus mismatched pocket comparisons, and receptor or scoring sensitivity analyses are then needed to judge whether an apparent AUROC has a stable structural explanation. Direct experimental measurement is still required for dual-target activity, activity balance, and mechanism.

The main limitations are heterogeneous activity data, nonuniform panel construction, and partial coverage of alternative docking settings. Experimental measurements at both targets remain scarce and heterogeneous. Four-state classification requires assays of the same ligand at both targets. Endpoint types remain heterogeneous despite pChEMBL scaling: IC50, Ki, Kd, EC50, and Potency are mapped onto one negative-log axis and are not biologically equivalent. Taking the maximum pChEMBL further favors the strongest recorded value. Repeating the assignment on the same ChEMBL 37 records with median aggregation produced the class flips in Table S3 and did not change the main directional judgments. That check only partly mitigates the bias of taking the strongest record. The resulting classes are database-derived retrospective states, not measurements from a single assay system. External-evaluation availability is reported in Results 3.6.

Panel construction was not fully uniform across target pairs. EGFR/HER2 and PIK3CA/mTOR were sampled from \(\theta=6.0\) pools, whereas the other six pairs were sampled from strict 6.5/5.5 pools; scaffold caps also differed across panels. Absolute AUROCs are therefore not intended for quantitative ranking across target pairs. The main inferences instead come from within-pair comparisons that hold the score channel fixed. The eight pairs were used to examine consistency across systems and should not be treated as independent statistical replicates: JAK1 and PPARA were each used twice, and some pairs belong to related kinase or nuclear-receptor families.

At the computational layer, the study relied mainly on one representative crystal structure per target and AutoDock Vina as the primary docking and scoring method. Cognate redocking showed coverage of near-native poses. That check mainly reflects search coverage and cannot replace evaluation of score ranking or screening performance. Independent GNINA pose generation, alternative receptor conformers, and other scoring functions covered only some systems. Ligand preparation did not systematically enumerate protonation states, tautomers, or conformational ensembles. The specific numerical values observed here should therefore not be extrapolated as general dual-target docking performance across other protein families and broader chemical space.

Stronger evidence for this type of evaluation will require paired activity measurements of the same compounds at both targets under more uniform experimental conditions, together with independent or prospective tests.

## 5. Conclusions

This study examined how experimental-state definition affects retrospective dual-target docking evaluation across eight human target pairs with experimental measurements at both targets. Holding the same target score channel fixed, changing the control class and its compounds can change the performance judgment in a pair- and direction-dependent way. Discriminating dual-active from both-end low-activity ligands does not establish discrimination against single-target selectives.

Some retrospective discrimination can be obtained from ligand chemistry that does not use receptor structure. Under the present ECFP4 model and scaffold-grouped splits, adding a docking score did not produce a consistent AUROC gain across directions. Across the seven internal holdouts, matched-pocket advantages were not consistently reproduced. A high retrospective AUROC alone therefore does not establish that the observed discrimination arises from three-dimensional information specific to the corresponding pocket. For retrospective dual-target benchmarks of this type, both directional selectivity tasks should be reported and interpreted with ligand-chemistry and pocket-correspondence controls. The resulting metrics describe discrimination among already labeled experimental states. They are not equivalent to prospective virtual-screening hit finding. Direct experiment remains required for dual-target activity and selectivity.

## Data and Software Availability

All panel definitions, experimental-state labels, receptor and docking-box specifications, per-ligand score tables, analysis scripts, and figure-generation code are available in the `Dual_Target_Docking` directory of the public repository at https://github.com/1280602962-debug/gwj260531. The reported statistical analyses can be reproduced from the deposited score and metadata tables without rerunning the full docking campaigns. Typeset Supporting Information contains Tables S1–S10. Additional quality-control tables remain in the public repository.

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
