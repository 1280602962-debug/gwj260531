# A Four-Pair Formulation Audit of Docking-Based Dual-Target Recognition

## Abstract

Docking scores at two targets are often interpreted as evidence of dual-target potential, but the apparent performance of that strategy depends on the negative class used for evaluation. We assembled four ChEMBL-derived ligand states—dual, A-only, B-only, and neither—for four target pairs and evaluated AutoDock Vina with two pocket-matched directional AUROCs that compare dual ligands with the corresponding single-target selectives. Their minimum (`summary_min`) was used as a conservative descriptive summary. The clearest formulation effect occurred for EGFR/HER2: Dual versus neither gave an AUROC of 0.756, whereas directional `summary_min` was 0.430; independent GNINA pose generation gave 0.783 and 0.220, respectively. The other pairs did not show the same separation. Across all four pairs, the primary `summary_min` 95% bootstrap intervals included 0.5, adding docking to scaffold-grouped ECFP4 changed cross-validated AUROC by at most 0.020, and receptor replacement shifted discrimination in a target-pair-dependent manner. Five frozen Vina seeds preserved the sign of the formulation contrast for every pair. A prespecified BindingDB-native external-set construction yielded no pair meeting the eligibility gate, so external docking was not performed. These results support selectivity hard negatives and ligand- and receptor-aware controls as essential components of dual-target docking evaluation.

**Keywords:** dual-target docking; benchmark formulation; selectivity hard negatives; chemical confounding; receptor realization; virtual screening

## 1. Introduction

Multitarget drug design seeks single molecules that modulate two or more biological targets and is increasingly pursued with structure-based computational methods.[1,2] Molecular docking remains a standard component of structure-based virtual screening (SBVS), where candidate ligands are placed in a binding site and ranked by a receptor–ligand scoring function.[3,4] In a dual-target setting, an intuitive strategy is to dock each candidate against both proteins and prioritize molecules with favorable scores at both sites. The validity of that inference, however, depends on how the evaluation problem is constructed.

Benchmark design is known to influence virtual-screening conclusions. DUD and DUD-E introduced property-matched decoys to reduce trivial ligand-property separation,[5,6] whereas LIT-PCBA used experimental assay labels to reduce several biases associated with artificial decoys.[7] CASF-2016 separated scoring, ranking, docking, and screening power on protein–ligand complexes.[8] These benchmarks address important single-target evaluation problems, but none explicitly represents the four experimental states required for dual-target selectivity: active at both targets, selective for A, selective for B, or inactive at both.

For dual-target recognition, A-only and B-only ligands are particularly informative negatives because they are already active at one member of the target pair. Zhou, Li, and Hou evaluated docking-based dual-kinase screening against noninhibitors,[9] and recent dual-target generative studies also use docking-derived success criteria when ranking designed molecules.[10,11] A favorable comparison against inactive compounds does not establish that docking can separate a true dual-active ligand from a compound that is potent at only one target. That distinction motivates a directional evaluation: dual versus A-only should be judged in pocket B, and dual versus B-only in pocket A.

Here we introduce DualFourClass-Bench as a four-state formulation audit over four frozen target pairs. The primary analysis consists of the two pocket-matched directional AUROCs, with their minimum (`summary_min`) used only as a compact worst-direction summary. We compare this formulation with Dual versus neither on the same docking scores and examine whether the resulting interpretation changes with ligand chemistry, activity aggregation, document and scaffold structure, receptor realization, docking random seed, and an independent pose-generation engine. The objective is not to rank docking engines, but to determine how the evidentiary standard for dual-target docking changes when experimentally selective hard negatives are included.

## 2. Methods

### 2.1 Experimental-state definition and analysis hierarchy

Ligand activities were obtained from the public ChEMBL Web API activity endpoint. pChEMBL values place several standardized potency and affinity measurements (including IC50, EC50, Ki, Kd, and Potency) on an approximate −log10 scale. Because these endpoints and assay contexts are not equivalent, pChEMBL was used as an operational curation scale rather than as an assay-harmonized measurement.

For a ligand–target pair with multiple usable records, the maximum pChEMBL value was used for the primary curation. Ligands without a usable measurement at both targets were excluded from analyses requiring a four-state label; missing measurements were not treated as inactivity. Multicomponent records were split by connected component and the organic fragment with the largest heavy-atom count was retained.

For each target pair A/B, four states were defined: **dual** (active at both targets), **A-only** (active at A and weaker at B), **B-only** (the converse), and **neither** (weaker at both). A strict 6.5/5.5 rule was used only to assess whether a candidate target pair had sufficient label supply: dual required pChEMBL ≥ 6.5 at both targets; A-only required A ≥ 6.5 and B ≤ 5.5; B-only was symmetric; neither required both values ≤ 5.5. The 5.5–6.5 gray zone was excluded from this supply audit. All primary benchmark labels were subsequently defined with a single threshold, θ = 6.0: dual, both values ≥ θ; A-only, A ≥ θ and B < θ; B-only, B ≥ θ and A < θ; neither, both values < θ. The two thresholds therefore served different, fixed purposes: supply qualification and primary evaluation, respectively. Relabeling at θ = 5.5 and 6.5 and under the strict 6.5/5.5 rule was analyzed as sensitivity (Table S4).

The analysis hierarchy distinguishes the four-pair θ = 6.0 Vina directional analysis from sensitivity and exploratory analyses (`docs/ANALYSIS_HIERARCHY_V1.md`). Sensitivity analyses include label aggregation, grouped uncertainty, receptor swaps, docking failures, alternate panel membership, random seeds, formulation contrasts, and external-set feasibility. Exploratory diagnostics include operating-point filters, full-map ligand-only models, and geometric summaries. This hierarchy was used to determine which analyses support the main claims and which remain Supporting Information diagnostics.

### 2.2 Benchmark construction

The initial ChEMBL supply audit considered 49 candidate target pairs. Metal-dependent systems such as HDAC1/HDAC6 were excluded before the final benchmark was assembled. The frozen evaluation set contains PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB, and EGFR/HER2. EGFR/HER2 was retained as a supply-limited case because its strict B-only pool was small.

Ligands were sampled under fixed class quotas with random seed 20260729. A per-class Bemis–Murcko scaffold cap was applied where the corresponding sampling files supported it: at most two molecules per scaffold for PIK3CA/mTOR (PM48) and at most five for EGFR/HER2. AChE/BChE and PIK3CA/PIK3CB used fixed class quotas and a deterministic shuffle. Panels were not resampled after docking results were inspected. PIK3CA/mTOR PM48 contained 18/14/12/4 dual/A-only/B-only/neither ligands under θ = 6.0; AChE/BChE and PIK3CA/PIK3CB targeted 28/28/28/16 under the strict construction gate, and EGFR/HER2 contained 110 ligands (Table 1). All primary AUROCs were computed after relabeling the scored panels uniformly at θ = 6.0.

**Table 1.** DualFourClass-Bench composition and docking settings. `n_scored` gives dual/A-only/B-only ligands with valid scores in both pockets for the primary directional analysis.

| Pair | Construction rule | PDB (A / B) | Resolution (Å) | n_panel | n_scored (dual / A-only / B-only) | Vina exhaustiveness |
|---|---|---|---|---:|---:|---:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | strict 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | strict 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

### 2.3 Receptor preparation and docking

Primary receptors were PIK3CA/mTOR 4L23/4JT6 (cognate ligands X6K/PI-103), AChE/BChE 4EY7/4BDS (E20/THA), PIK3CA/PIK3CB 4L23/2WXF (X6K/039), and EGFR/HER2 3POZ/3RCD (03P/TAK-285). The docking site was defined from the cognate ligand. The axis-aligned bounding box around cognate heavy atoms was expanded by 5 Å in each direction, with each box edge set to at least 20 Å (Table S2). Water and cognate ligands were removed before receptor conversion with Meeko.

Each primary receptor was subjected to cognate redocking before production docking. Nine poses were retained, and the fixed search-coverage criterion was a best-of-nine heavy-atom RMSD < 2.0 Å to the crystallographic ligand. If the default search effort did not meet this criterion, the prespecified fallback exhaustiveness was used; this yielded production exhaustiveness 16 for PIK3CA/mTOR and 8 for the other pairs (Table S3). A later ranked-pose audit reported top-1, top-3, and best deposited-pose RMSD where topology-checked reconstruction was possible. The original EGFR/HER2 nine-pose production PDBQTs were unavailable, so their ranked-pose QC was reconstructed under the frozen protocol and is labeled accordingly in the Supporting Information.

Ligands were prepared from frozen ChEMBL SMILES by retaining the largest organic fragment, adding explicit hydrogens, generating a three-dimensional structure with ETKDGv3 (seed 20260727), performing up to 200 MMFF optimization steps, and converting to PDBQT with Meeko. Protonation, tautomer, and conformer ensembles were not systematically enumerated. Production docking used AutoDock Vina 1.2.7 with the `vina` scoring function, nine poses, `energy_range = 3` kcal mol−1, and seed 20260727. The primary readout was the mode-1 Vina energy. RTMScore and GNINA CNN rescoring of the same Vina poses were secondary scoring-function checks (Tables S14–S15).

GNINA 1.3.2 was also run in docking-search mode on EGFR/HER2 and PIK3CA/mTOR to generate poses independently from Vina using the same frozen ligand inputs, receptors, boxes, exhaustiveness values, nine retained poses, and seed 20260727. Mode-1 `minimizedAffinity` was used for this engine-level sensitivity (Table S32).

Random-seed sensitivity used five fixed Vina seeds: 20260727 (the primary run) and 20260811–20260814. Receptors, boxes, exhaustiveness, retained modes, and `energy_range` were otherwise unchanged. Per-seed scores and AUROCs are deposited in `data/jcim_multiseed_v0/`.

### 2.4 Directional endpoints and statistical analysis

For each target pair A/B, dual was the positive class. Dual versus A-only was evaluated with the pocket B score,
\[
\mathrm{AUC}_{D/A}=\mathrm{AUROC}(\mathrm{dual},\mathrm{A\!-\!only};S_B),
\]
and dual versus B-only with the pocket A score,
\[
\mathrm{AUC}_{D/B}=\mathrm{AUROC}(\mathrm{dual},\mathrm{B\!-\!only};S_A).
\]
Because a more negative Vina energy is favorable, the analyzed score was \(S_{\mathrm{Vina}}=-E_{\mathrm{Vina}}\).

The worst-direction summary was
\[
\mathrm{summary}_{\min}=\min(\mathrm{AUC}_{D/A},\mathrm{AUC}_{D/B}).
\]
`summary_min` is a descriptive compression of the two directional AUROCs; the component AUROCs remain the underlying discrimination estimates. Arithmetic, geometric, and harmonic aggregation were examined as sensitivity (Table S26).

Dual versus neither used experimental inactives as the negative class and `vina_mean` as the two-pocket score. Dual versus all non-duals was also reported. These comparisons change the target estimand and are therefore interpreted as formulation contrasts rather than as paired tests of one common performance quantity. A same-pocket comparison was additionally performed to separate negative-class composition from score aggregation (Table S34).

Ligand-level uncertainty was estimated by stratified bootstrap resampling with replacement (B = 2000, seed 20260729), reporting percentile 95% confidence intervals. Document-cluster and Bemis–Murcko scaffold-cluster bootstraps resampled correlated groups rather than individual ligands. A simulation-based detectable-effect analysis used the observed class sizes and a binormal score model to characterize how readily the panels could resolve moderate versus large AUROCs (Table S31; Figure S6).

### 2.5 Sensitivity, confounding, and external-set analyses

Four prespecified ligand descriptors—heavy-atom count, molecular weight, cLogP, and TPSA—were evaluated with the same directional workflow. Morgan/ECFP4 fingerprints (radius 2, 2048 bits) with logistic regression were evaluated under Bemis–Murcko scaffold GroupKFold. Models using ECFP4 alone and ECFP4 plus the pocket-matched docking score were compared on identical folds (Tables S20, S24). Additional diagnostics included ligand-efficiency normalization, potency- and size-constrained subsets, property-caliper matching, and a wrong-pocket falsification control.

Activity aggregation was re-evaluated from API-refetched assay-level records using the median rather than the maximum pChEMBL value while retaining the frozen docking scores (Table S29). A dated high-confidence ChEMBL view retained Homo sapiens `SINGLE PROTEIN` targets, assay confidence ≥8, exact relations, accepted quantitative endpoints, and records without validity or potential-duplicate flags (Table S36). Complete-case coverage and source-document concentration were quantified separately (Table S37). A document-blocked analysis grouped ligands connected by retained ChEMBL `document_id` values, and a literature-year split used a 2018 primary cutoff with 2015 and 2020 sensitivities (Tables S39–S41).

Panel-membership sensitivity used an unused ChEMBL pool after removing the main panels and PM110; it was available for PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB (20 dual/20 A-only/20 B-only per pair; seed 20260731). Missing docking scores were analyzed by class and ligand properties, and directional rank-extreme bounds were computed rather than silently treating failed jobs as ordinary missing observations (Table S27).

Receptor sensitivity replaced one crystal structure at a time while holding the other pocket fixed. PIK3CA alternatives 4JPS and 5DXT were evaluated in PIK3CA/mTOR and PIK3CA/PIK3CB; mTOR 4JSX was evaluated in PIK3CA/mTOR. All replacement structures satisfied the same target, cognate-site, resolution, and cognate-redocking criteria used for the primary receptors (Tables S10, S30).

For an external-set feasibility test, a BindingDB-native slice was rebuilt from the versioned 202608 article and patent archives using a contract fixed before external docking (`external_slice_contract.yaml`). The construction required human wild-type single-chain target mapping, exact IC50/Ki/Kd measurements on both targets, median aggregation within ligand–target–endpoint, θ = 6.0 four-state labels, removal of overlapping documents and structures, and maximum ECFP4 Tanimoto < 0.70 to development ligands. The eligibility gate required dual, A-only, and B-only each n ≥ 20, at least three sources per class, and no single source contributing >50% of a class. External docking was to proceed only for pairs passing this gate.

## 3. Results

### 3.1 Selectivity hard negatives constrain benchmark supply

The strict 6.5/5.5 supply audit showed that bidirectional selective hard negatives were scarce. Among 49 candidate pairs, only four met a thick-panel criterion of at least 50 strict hard negatives in each direction. After excluding metal-dependent HDAC1/HDAC6, PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB remained the best-supplied cases, whereas EGFR/HER2 had only seven strict B-only ligands and was retained as a supply-limited case. A count-only BindingDB/PubChem audit supported the same qualitative supply constraint (Table S12).

The primary docking analysis used the unified θ = 6.0 four-state labels, not the strict supply threshold. A θ = 6.0 census of the original candidate list identified 17 unique pairs with at least 10 ligands in each directional class (Table S44), but this census was not used to expand the four-pair docking benchmark after results were observed.

### 3.2 The effect of negative-class formulation is target-pair dependent

The four primary Vina panels produced different directional discrimination profiles (Table 2; Figure 4A). EGFR/HER2 yielded directional AUROCs of 0.666 and 0.430, AChE/BChE 0.650 and 0.606, PIK3CA/PIK3CB 0.691 and 0.500, and PIK3CA/mTOR 0.714 and 0.692. The corresponding `summary_min` values were 0.430, 0.606, 0.500, and 0.692. All four ligand-bootstrap 95% intervals for `summary_min` included 0.5. Aggregating the two directional AUROCs by arithmetic, geometric, or harmonic mean did not change the pair ordering (Table S26).

**Table 2.** Pocket-matched directional AUROC on the frozen four-pair set (Vina; θ = 6.0). The descriptor columns report `summary_min` for the four ligand-property references.

| Pair | n_scored (dual / A-only / B-only) | dual vs A-only (pocket B) | dual vs B-only (pocket A) | summary_min [95% CI] | heavy atoms | MW | cLogP | TPSA |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.350, 0.650] | 0.622 | 0.620 | 0.595 | 0.418 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] | 0.463 | 0.448 | 0.310 | 0.260 |

Reformulating the same scores as Dual versus neither produced the largest change for EGFR/HER2 (Table 3; Figure 3). Its Dual-versus-neither AUROC was 0.756 [0.562, 0.920], compared with directional `summary_min` 0.430 [0.282, 0.578]. When the pocket A score was held fixed, replacing B-only negatives with neither negatives changed AUROC by 0.378 [0.205, 0.547] (Table S34), showing that the result cannot be attributed only to averaging scores across two pockets. AChE/BChE and PIK3CA/PIK3CB showed much smaller differences, while the PIK3CA/mTOR neither class contained only four ligands and was underpowered. Thus, the favorable Dual-versus-neither result was not a general four-pair phenomenon.

**Table 3.** Directional and Dual-versus-neither formulations on the same frozen Vina scores (θ = 6.0).

| Pair | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|---|---:|---:|---:|---:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/PIK3CB | 0.500 [0.350, 0.650] | 0.559 [0.373, 0.746] | 16 | 0.556 [0.437, 0.672] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |

Independent GNINA pose generation preserved the EGFR/HER2 formulation gap: Dual versus neither was 0.783 [0.610, 0.922], whereas directional `summary_min` was 0.220 [0.109, 0.343] (Table S32). For PIK3CA/mTOR, independent GNINA gave a directional `summary_min` of 0.633 versus 0.692 for Vina. These results show that the EGFR/HER2 formulation effect was not specific to one Vina pose-generation run, while also reinforcing that the magnitude and direction of performance are pair dependent.

### 3.3 Ligand chemistry and data provenance contribute to the observed separation

Simple ligand properties were competitive with docking on several directional tasks (Table 2). AChE/BChE was the clearest example: TPSA alone exceeded the corresponding Vina discrimination, and adding heavy-atom count and TPSA to a logistic model increased dual-versus-B-only AUROC from 0.606 to 0.807 (Tables S19, S28). ECFP4 provided an even stronger ligand-only reference on several arms; under Bemis–Murcko scaffold GroupKFold, the EGFR/HER2 dual-versus-B-only ECFP4 AUROC was approximately 0.89 compared with 0.43 for docking. Across all directional tasks, adding the pocket-matched docking score to ECFP4 changed scaffold-grouped cross-validated AUROC by at most 0.020 in absolute value (Table S24).

This ligand-only signal was partly associated with literature and chemical-series structure. Document-blocked analysis left the EGFR/HER2 weak directional arm at 0.430, with a document-cluster bootstrap interval of [0.321, 0.617]; scaffold-cluster bootstrap gave [0.278, 0.595]. Complete-case coverage was limited: among structures with a usable value on at least one target, only 14.5%–34.0% across the four pairs had usable measurements at both targets (Table S37). PIK3CA/mTOR was especially concentrated in its neither class, where all four ligands came from one ChEMBL document.

Activity-label sensitivity did not materially alter the primary point estimates. In the API-refetched max-versus-median analysis, `summary_min` changed from 0.417 to 0.424 for EGFR/HER2, 0.606 to 0.629 for AChE/BChE, and remained 0.500 and 0.692 for PIK3CA/PIK3CB and PIK3CA/mTOR, respectively (Table S29). A high-confidence current-ChEMBL filter retained the θ = 6.0 operational state for all 352 scored ligands (Table S36). These checks support label stability under the tested record filters, while the underlying assays remain heterogeneous in endpoint and experimental context.

### 3.4 Receptor choice and random seed are distinct sources of sensitivity

Replacing one receptor while holding the other fixed produced target-pair-dependent changes (Figure 5; Table S30). For PIK3CA/mTOR, replacing PIK3CA 4L23 with 4JPS or 5DXT reduced `summary_min` from 0.692 to 0.486 [0.259, 0.692] and 0.505 [0.292, 0.696], respectively. In contrast, applying the same PIK3CA replacements to PIK3CA/PIK3CB increased `summary_min` from 0.500 to 0.691 and 0.685. Replacing mTOR 4JT6 with 4JSX gave 0.639 for PIK3CA/mTOR. The direction of receptor-induced change therefore depended on the target pair rather than following a common shift.

Random-seed sensitivity was smaller than the receptor shifts and preserved each pair's qualitative formulation pattern. Across the five fixed Vina seeds, `summary_min` ranged from 0.321–0.430 for EGFR/HER2, 0.553–0.606 for AChE/BChE, 0.468–0.502 for PIK3CA/PIK3CB, and 0.676–0.726 for PIK3CA/mTOR. For every pair, the sign of the difference between Dual versus neither and directional `summary_min` matched the primary seed in all five runs (`data/jcim_multiseed_v0/analysis/MULTISEED_VINA_VERDICT_V1.md`). Thus, the central pair-specific formulation pattern was not created by selecting a favorable Vina random seed.

Docking failures were also examined rather than treated as exchangeable missing data. Failure rates and ligand-property shifts are reported in Table S27, together with arm-available and rank-extreme directional bounds. The unused-pool panels produced mixed changes across pairs, and wrong-pocket point estimates could exceed matched-pocket estimates in that internal holdout (Table S17). These controls argue against interpreting a single receptor or panel realization as a general property of docking.

### 3.5 A prespecified external-set construction failed before docking

The BindingDB-native 202608 reconstruction was evaluated against its fixed eligibility gate before any external docking. After literature-overlap, structure-overlap, and ECFP4 < 0.70 filters, no target pair met the requirement that dual, A-only, and B-only each contain at least 20 ligands with at least three sources per class and no source contributing more than half of a class (Tables S48–S49; Figure S8). The remaining counts are upper bounds because the ChEMBL literature-overlap mapping was incomplete. In accordance with the stop rule, no BindingDB external docking AUROC was computed.

The prespecified 2018 literature-year split also did not yield enough bidirectional hard negatives for at least two evaluable pairs (Table S41). These outcomes identify an important data-supply problem: independently profiled compounds with measurements at both targets and adequate selective negatives are difficult to assemble at the scale required for an external four-state benchmark.

## 4. Discussion

### 4.1 Selectivity hard negatives change the evidentiary standard

The main result is a formulation effect, not a four-pair docking leaderboard. On EGFR/HER2, the same docking scores looked substantially more favorable when dual ligands were compared with experimentally inactive ligands than when they were compared directionally with single-target selectives. The same pattern remained under independent GNINA pose generation and across five Vina random seeds. The other three pairs did not reproduce the same magnitude or direction of separation, so the evidence supports a target-pair-dependent formulation effect rather than a universal claim that Dual-versus-neither evaluation systematically overestimates dual-target recognition.

This distinction extends a familiar lesson from single-target benchmarking. Property-matched decoys and assay-derived benchmarks were developed because the composition of the negative class can dominate apparent virtual-screening performance.[5–7,12,13] In a dual-target task, single-target selectives are the natural hard negatives: they test whether a method can identify what is missing at the opposite target, rather than simply recognizing chemistry associated with activity somewhere in the pair. A dual-pocket score is therefore more informative when it is challenged against both A-only and B-only ligands.

### 4.2 Docking should be interpreted alongside ligand and receptor controls

The ligand-only baselines show that experimental four-state labels contain substantial chemical-series information. Under scaffold-grouped evaluation, adding docking to ECFP4 changed AUROC by no more than 0.020 in this data set. This does not imply that docking contains no structural information; it means that the additional separation supplied by the present docking scores was small relative to the 2D chemical signal under these splits. For benchmark studies, a ligand-only baseline is therefore necessary to determine whether an apparent docking result exceeds what can already be inferred from chemotype and simple molecular properties.

Receptor realization provides a second, independent source of variation. Replacing PIK3CA altered the two PIK3CA-containing target pairs in opposite directions despite using the same replacement structures. This observation is consistent with cross-docking studies that treat receptor representation as part of the evaluation problem rather than as a fixed nuisance variable.[14] The multi-seed analysis helps separate this effect from stochastic search variation: the qualitative formulation pattern was stable across the five Vina seeds, whereas receptor replacement produced much larger pair-specific shifts.

### 4.3 Implications for dual-target virtual screening

A practical dual-target evaluation can therefore be organized around four checks. First, report both pocket-matched directional comparisons against A-only and B-only ligands rather than only Dual versus inactive compounds. Second, include ligand-only chemical baselines under scaffold- or source-aware splits. Third, test sensitivity to panel membership or literature grouping when the data allow it. Fourth, examine at least one alternative receptor realization for structurally flexible targets. These checks address different failure modes—formulation, chemical confounding, data provenance, and receptor dependence—and together provide a more informative interpretation of dual-pocket scores.

The same logic applies when docking is used as a filter in generative or large-scale virtual-screening workflows. A molecule can satisfy two numerical docking thresholds without being distinguishable from experimentally selective chemistry. Evaluation should therefore be tied to the intended biological contrast, not only to favorable absolute scores or comparisons with inactive compounds.

### 4.4 Limitations

The primary benchmark contains only four target pairs, two of which share PIK3CA, and three pairs involve kinase ATP-site recognition. It is consequently a case-based formulation audit rather than a representative survey of dual-target pharmacology. All four ligand-bootstrap `summary_min` intervals include 0.5, and the detectable-effect analysis shows limited sensitivity to moderate directional effects at the available class sizes. Pairwise point estimates should therefore not be interpreted as a ranking of target tractability.

The four-state labels are experimentally derived but not assay harmonized. They combine different quantitative endpoints and require measurements at both targets, which enriches the benchmark for jointly profiled chemistry and well-studied series. The high-confidence ChEMBL and max-versus-median analyses support robustness to the tested record filters, but construct, mutation, and detailed assay-condition metadata remain incompletely harmonized. The unused ChEMBL pool and document/year analyses are internal sensitivities, not a substitute for an independently docked external set.

External transfer remains untested because the fixed BindingDB eligibility gate failed before docking. Receptor and ligand preparation also define the chemical domain of the present results: cognate best-of-nine redocking measures search coverage rather than guaranteed top-ranked pose correctness, and protonation/tautomer/conformer ensembles were not systematically explored. An exploratory MCL1/Bcl-xL extension was excluded from the main benchmark after formal topology-aware pose-gold validation could not be established; its archived calculations are not used in the primary evidence.

## 5. Conclusions

Dual-target docking should be evaluated against experimentally selective hard negatives, not only against compounds inactive at both targets. In this four-pair audit, the resulting formulation effect was strongest for EGFR/HER2 and persisted under independent pose generation and multiple Vina seeds, whereas the remaining pairs showed different behavior. Ligand chemistry and receptor realization also materially influenced the apparent discrimination. These observations support a benchmark design in which directional selectivity contrasts, ligand-only baselines, provenance-aware sensitivity analyses, and receptor alternatives are reported together before favorable scores in two pockets are interpreted as evidence of dual-target activity. Broader target coverage and genuinely independent external data will be required to establish how widely these findings transfer.

## Data and Software Availability

Benchmark memberships, four-state labels, receptor and box definitions, per-ligand docking scores, analysis scripts, and manuscript-facing result tables are available in the public repository `1280602962-debug/gwj260531`, under `Dual_Target_Docking`. `data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` indexes the principal numerical results and their source files. The five-seed Vina analysis is deposited in `data/jcim_multiseed_v0/`, the BindingDB-native construction contract in `data/jcim_novelty_v0/protocol/external_slice_contract.yaml`, and the benchmark evaluation contract in `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`. Manuscript-facing table checksums are recorded in `REVISION_CHECKSUM_MANIFEST_v1.csv`, and the pinned analysis environment and reproduction commands are documented in the repository. The ChEMBL supply audit was frozen on 2026-07-23, the high-confidence ChEMBL view was fetched on 2026-08-26, and the BindingDB-native archive lock is release 202608. BindingDB source archives are not redistributed.

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