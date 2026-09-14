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
| AChE/BChE | strict 6.5/5.5 | 28 / 28 / 28 / 16 | 4EY7 / 4BDS | 27 / 25 / 28 | 8 |
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

Cognate-ligand redocking assessed recovery of the experimental bound conformation. The original outputs report the top-ranked RMSD, the lowest RMSD among the first three ranked poses (top-3), and the lowest RMSD among all saved poses. The last metric assessed search coverage at a 2.0 Å cutoff and is not a claim that the top-ranked pose was near-native; top-1 assessed pose ranking. All 14 primary receptors were rescored from saved poses with one chemically mapped CalcRMS recipe (Table S2), without redocking. AChE 4EY7 and TYK2 3LXP deposited 8 poses; the other primary receptors deposited 9. For PIK3CA/mTOR, graph-automorphism CalcRMS was used because the prepared ligand was not in the crystal frame. EGFR 3POZ is reconstructed QC, not a recovered production output.[8]

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

To isolate the effect of the virtual-screening control class from score aggregation, the same pocket score was held fixed and AUROCs were compared under different experimental-state controls. Pocket A scores were used for dual-versus-B-only and that dual-versus-neither comparison; pocket B scores were used for dual-versus-A-only and the corresponding dual-versus-neither comparison (Table S4). Changing the control class also changes the specific compounds. The lower of the two pocket scores, \(S_{\mathrm{worst}}=\min(S_{A},S_{B})\), was used as an AND-type two-pocket filter. The median of dual \(S_{\mathrm{worst}}\) values defined the filter threshold. Dual counts, dual recall, dual precision, and retained single-target-active ligands at that threshold are shown in Figure S1. The corresponding Top-10 class composition under two-pocket mean ranking is shown in Figure 2D.

#### 2.5.3 Confidence intervals and resampling

Table 2 reports pointwise 95% confidence intervals for both directional AUROCs and their descriptive minimum. The original ligand-level non-stratified percentile bootstrap was retained, with 2000 replicates and the recorded deterministic seeds. Each replicate sampled the original number of ligands with replacement from the pooled dual, A-only, and B-only set, recomputed both directional AUROCs using the same resampled dual ligands, and took their minimum. Replicates missing a required class were omitted. The 2.5th and 97.5th percentiles were calculated separately for each directional AUROC and for the replicate-wise minimum. Point estimates used the full analysis sample. A supplementary class-stratified bootstrap preserved the three class sizes and reused the same dual draw in both directions. This post hoc sensitivity analysis is archived with the score tables and did not replace the original intervals. Dual-versus-neither and dual-versus-all-nonduals intervals used class-stratified percentile bootstrap. Intervals were not adjusted for multiple comparisons.

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

Experimental states were reassigned under alternative activity thresholds and activity-aggregation rules, with panel membership and docking scores held fixed (Table S3). For pairs with enough remaining candidates, unused-pool holdout sets were built after excluding main-panel members. Holdout ligands were drawn under a fixed sampling rule that limited over-representation of a single Bemis–Murcko scaffold. Directional AUROCs were recomputed on these different ligands to assess sensitivity to panel membership. Because the ligands still come from the same source, this analysis is an internal robustness check (Table S6). Primary Vina docking was also repeated with five fixed random seeds, and different exhaustiveness settings were compared on PIK3CA/mTOR. Complete-case counts and fixed-membership intersections are archived in the repository.

BindingDB[16] and PubChem were searched for all eight pairs as candidate sources of an independent external docking set. An independence filter at \(\theta=6.0\) removed shared literature, duplicate structures, and molecules with ECFP4 Tanimoto \(\geq 0.70\) to the development set (Table S8). The development-molecule set included main panels, the expanded PIK3CA/mTOR PM110 panel, and internal holdouts. Dual, A-only, and B-only each needed \(n\geq 20\) with at least three sources per class. Only pairs meeting those conditions would have entered external docking evaluation.

### 2.7 Software and reproducibility

Molecular processing, fingerprints, descriptors, and statistics were implemented in Python, mainly with RDKit, NumPy, pandas, SciPy, and scikit-learn. Receptor and ligand PDBQT files were prepared with Meeko. Primary docking used AutoDock Vina 1.2.7. GNINA 1.3.2 was used for supplementary scoring and independent pose generation on selected pairs. RTMScore provided an additional rescoring of Vina poses.

Software versions, random seeds, and main computational parameters are summarized in Table S1. Docking boxes and cognate redocking results are in Table S2.
