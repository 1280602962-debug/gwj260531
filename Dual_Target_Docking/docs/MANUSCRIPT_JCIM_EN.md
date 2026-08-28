# A Four-Pair Formulation Audit of Docking-Based Dual-Target Recognition

## Abstract

Benchmark construction determines what favorable scores at two targets can establish about dual activity. We assigned ChEMBL ligands for four frozen target pairs to dual, A-only, B-only, or neither states and evaluated AutoDock Vina in two pocket-matched directions: dual versus A-only in pocket B and dual versus B-only in pocket A. The smaller directional AUROC (`summary_min`) served as a conservative descriptive summary. On EGFR/HER2, Dual versus neither produced an AUROC of 0.756, whereas `summary_min` was 0.430; independent GNINA pose generation preserved the contrast (0.783 versus 0.220). The other pairs did not show the same gap. Across five frozen Vina seeds, the sign of the Dual-versus-neither minus `summary_min` gap was unchanged for every pair, with median `summary_min` values of 0.373, 0.599, 0.478, and 0.704 for EGFR/HER2, AChE/BChE, PIK3CA/PIK3CB, and PIK3CA/mTOR, respectively. Adding docking to scaffold-grouped ECFP4 models changed cross-validated AUROC by at most 0.020; alternative PIK3CA receptors shifted discrimination in opposite directions across two related pairs; and all four primary-seed `summary_min` confidence intervals included 0.5. A pre-frozen BindingDB-native gate yielded no eligible external pair. Thus, selective hard negatives, ligand-only controls, and receptor and seed sensitivities are necessary for interpreting dual-pocket docking, while the present four-pair panel does not establish target-general performance.

**Keywords:** dual-target docking; benchmark formulation; selectivity hard negatives; chemical confounding; receptor realization; virtual screening

## 1. Introduction

Multitarget drug design aims to modulate two or more biological targets with a single small molecule, in order to address pathway redundancy, compensatory signaling, and drug resistance in complex disease. Relative to a classical single-target agent, a rationally designed multitarget ligand may act on connected nodes of a disease network and thereby achieve a more adequate pharmacological effect; this idea is now a central theme of polypharmacology.[1] Over the past decade the field has moved from largely serendipitous multi-pharmacology toward structure-guided design that combines structural biology, computational chemistry, and, increasingly, generative models.[2] Molecular docking remains one of the most widely used tools in structure-based virtual screening (SBVS): a docking engine places the ligand in a protein binding site and a scoring function ranks ligand–receptor complementarity.[3,4] A natural computational tactic in dual-target discovery is therefore to dock each candidate into both pockets and to treat favorable scores on both targets as evidence of dual-target potential. How such scores should be interpreted depends on benchmark construction. DUD and DUD-E use property-matched decoys because unmatched negatives can reduce screening to separation by coarse ligand properties.[5,6] LIT-PCBA instead uses experimental assay labels and systematically controls known decoy and chemical biases.[7] CASF-2016 evaluates scoring, ranking, docking, and screening power on protein–ligand complexes, but still poses a single-complex problem.[8] None of these resources formulates dual-target discrimination over an experimentally labeled four-state ligand space.

A strict dual-target evaluation must distinguish four ligand states: **dual-active**, **A-selective**, **B-selective**, and **neither** (Figure 1A). A-only and B-only ligands are the **selectivity hard negatives** of the task: already potent on one target, they can produce plausible scores in that pocket while lacking corresponding activity on the other. The computational endpoint therefore asks whether docking discriminates dual-actives from the matching single-target selectives in both directions. Zhou, Li, and Hou benchmarked dual-kinase docking against noninhibitors.[9] Building on that setting, we introduce experimentally defined directional hard negatives and compare apparent discrimination under different formulations on the same scores. Dual versus neither, scored on experimental inactives, is reported as a formulation contrast. A balanced four-state panel is limited by the need for paired measurements and enough selectives on both arms.

Recent dual-target generative methods still rely on docking-based success metrics relative to reference ligands.[10,11] These metrics evaluate computational dual-target design relative to reference ligands, whereas our benchmark tests discrimination against experimentally defined selective hard negatives.

Here, we ask whether benchmark formulation changes the apparent evidence for dual-target recognition. We construct DualFourClass-Bench as a four-state panel with two directional primary tasks—dual versus A-only scored in pocket B, and dual versus B-only scored in pocket A (Figure 1B)—and a conservative worst-direction discrimination summary (`summary_min`). We then test whether the observed discrimination is retained under alternative ligand panels, activity-aggregation rules, and receptor structures.

## 2. Methods

### 2.0 Analysis hierarchy

Analyses were classified before manuscript claims were finalized as **primary** (frozen K = 4 pairs, θ = 6.0 labels, Vina mode-1 pocket scores, and the two pocket-matched directional AUROCs), **prespecified sensitivity** (relabeling, missing-data bounds, receptor swaps, grouped resampling, formulation contrasts, holdout and literature-year splits, and the BindingDB-native supply audit), or **post-hoc exploratory** (operating-point filters, label-supply censuses, full-map ligand-only models, geometric snapshots, and archived exploratory docking). Post-hoc outputs are diagnostic only and are not phrased as confirmatory primary evidence.

### 2.1 Data and experimental-state definition

Ligand activities used to define operational experimental states were retrieved from the public ChEMBL Web API activity endpoint. The target-pair supply audit was frozen on 2026-07-23. pChEMBL converts standardized quantitative potency or affinity measurements (e.g. IC50, EC50, Ki, Kd, and Potency) to an approximate −log10 activity scale. Assay types, conditions, experimental systems, species, and target-confidence levels are not necessarily equivalent; pChEMBL is therefore a heterogeneous curation scale here, not an assay-equivalent ground truth.

As a post-hoc label robustness analysis, the scored-panel molecule–target records were re-fetched from the current ChEMBL API on 2026-08-26 under Homo sapiens `SINGLE PROTEIN` targets, assay confidence score ≥8, exact relations, IC50/Ki/Kd/EC50/Potency endpoints, and no data-validity or duplicate flags; maximum pChEMBL was reapplied at θ = 6.0 without changing panel membership or docking scores (Table S36).

When several pChEMBL values existed for the same ligand–target pair, the **maximum** was used as the one-to-one representative for primary curation. Activity-aggregation sensitivity replaced the maximum with the **median** under the same θ = 6.0 rule (Table S29; class-wise record counts in Table S35). Ligands missing a usable pChEMBL value on either target were excluded from analyses requiring paired labels. Consequently, the benchmark represents the subset of ChEMBL compounds measured at both targets and may enrich profiled chemotypes, known polypharmacology, and well-studied series. Absence of a measurement was not interpreted as inactivity. Salt, solvate, and multicomponent records were split by connected component, retaining the organic fragment with the most heavy atoms.

Complete-case map coverage and source-document concentration were profiled separately (Tables S37, S36); they diagnose coverage and provenance concentration and do not identify the activity distribution of unmeasured compounds.

For each pair A/B, ligands were assigned one of four experimental states: **dual** (strong on both), **A-only** (strong on A, weak on B), **B-only** (strong on B, weak on A), and **neither** (insufficient on both). A-only and B-only ligands are selectivity hard negatives.

The **strict 6.5/5.5 criterion** was: dual, both pChEMBL ≥ 6.5; A-only, A ≥ 6.5 and B ≤ 5.5; B-only symmetric; neither, both ≤ 5.5. The 5.5–6.5 gray zone was excluded from this audit. Metal-dependent systems (e.g. HDACs) were excluded in advance. A strict 6.5/5.5 criterion was used only for target-pair supply qualification. All primary benchmark labels were then defined uniformly using θ = 6.0 (dual, both ends ≥ θ; A-only, A ≥ θ and B < θ; B-only symmetric; neither, both < θ), which was fixed before docking results were inspected. Construction rules were frozen from the supply audit before sampling (Table 1). Relabeling at θ ∈ {5.5, 6.5} and under the strict 6.5/5.5 rule is reported as sensitivity (Table S4). Underpowered cells are flagged in Results.

As a count-only check of the ChEMBL supply gate, the frozen pairs were recounted in BindingDB and PubChem without docking or panel rebuild (Table S12).

A BindingDB-native external slice was rebuilt from versioned 202608 article and patent TSV archives under rules frozen in `external_slice_contract.yaml` before independent class counts could change them: curated article or patent sources; human wild-type single-chain UniProt mapping; exact IC50/Ki/Kd; both ends measured; median within ligand–target–endpoint; θ = 6.0 four-state labels; exclusion of development-panel documents, InChIKeys, ChEMBL IDs, and ECFP4-similar molecules (Tanimoto < 0.70 to development ligands). The prespecified gate required dual, A-only, and B-only each n ≥ 20, at least three sources per class, and top-document ligand fraction ≤ 50%; gate failure halted the pipeline before docking, and no external pair qualified (Supporting Information).

### 2.2 Benchmark construction

DualFourClass-Bench retains four experimental states. The primary analysis comprises two directional pairwise tasks. The neither class is retained for descriptive formulation analysis.

Candidate pairs were screened with the strict audit in Section 2.1. The frozen evaluation set comprises PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB, and EGFR/HER2. EGFR/HER2 is retained as a supply-limited case. Ligands were drawn under frozen class quotas and random seed 20260729. Where structures were available at sampling, a per-class Bemis–Murcko scaffold cap limited series over-representation: at most two molecules per scaffold in PIK3CA/mTOR (PM48) and at most five in EGFR/HER2. AChE/BChE and PIK3CA/PIK3CB used class quotas and a deterministic shuffle. Panels were not redrawn after docking scores were seen.

AChE/BChE and PIK3CA/PIK3CB were sampled under the strict supply gate (target 28 / 28 / 28 / 16; n_panel = 100). EGFR/HER2 (n_panel = 110) and PIK3CA/mTOR PM48 (n_panel = 48; constructed 18 / 14 / 12 / 4) were assembled under the primary θ = 6.0 labels (Table 1). Cross-pair AUROCs therefore mix target-pair biology with panel-construction differences. Ligand–receptor jobs that failed to yield a score were dropped; n_scored can fall below n_panel (Table 1; Table S27). Failed compounds were audited by class and RDKit properties on the retained largest fragment. For each directional arm, an arm-available sensitivity used any valid score in the required pocket, and rank-extreme bounds assigned every comparison involving a missing required-pocket score either against or in favor of the claimed direction. These bounds are deterministic missing-data stress tests, not an imputation model. An expanded PIK3CA/mTOR panel (PM110) keeps all 48 PM48 ligands as a nested size check.

**Table 1.** DualFourClass-Bench composition and docking settings. Construction labels record the supply/panel-building rule for each pair. All primary AUROCs in Tables 2–3 use unified θ = 6.0 experimental-state labels. n_panel is the frozen panel membership, including neither; n_scored is the dual / A-only / B-only count with valid both-end Vina scores entering the primary directional AUROCs.

| Pair | Construction labels | PDB (A / B) | Resolution (Å) | n_panel | n_scored (dual / A-only / B-only) | Vina exhaustiveness |
|------|---------------------|-------------|----------------|-------:|------------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | strict 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | strict 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

### 2.3 Receptor preparation and docking protocol

Receptors were PDB entries with a small-molecule cognate ligand: PIK3CA/mTOR, 4L23 / 4JT6 (X6K / PI-103); AChE/BChE, 4EY7 / 4BDS (E20 / THA); PIK3CA/PIK3CB, 4L23 / 2WXF (X6K / 039); EGFR/HER2, 3POZ / 3RCD (03P / TAK-285). The site was defined from the cognate ligand. An axis-aligned bounding box on cognate heavy atoms was expanded by 5 Å on each axis; any edge shorter than 20 Å was set to at least 20 Å (Table S2). Water and the cognate ligand were removed and Meeko wrote PDBQT. PIK3CA, mTOR, EGFR, and HER2 used hydrogen-containing protein coordinates already in the frozen directories (`mk_prepare_receptor.py --read_pdb`). AChE, BChE, and PIK3CB were extracted from deposited ATOM/TER records and prepared with `mk_prepare_receptor` (default alternate location A). Docking treated noncovalent small-molecule sites only.

Before production docking, each frozen receptor was redocked with its cognate ligand. Nine poses were requested and heavy-atom RMSD to the crystal ligand was computed in the docking frame. The prespecified pass criterion was \(\mathrm{RMSD}_{\mathrm{best9}} < 2.0\) Å, i.e. whether any retained pose lies within 2.0 Å of the crystal ligand. If default exhaustiveness failed the gate, search effort was raised to a prespecified fallback. Production docking therefore used exhaustiveness 16 for PIK3CA/mTOR and 8 for the other main panels. A post-hoc ranked-pose re-audit reported top-1, top-3, and all-deposited-pose RMSD where the underlying SDF/PDBQT artifacts permit topology-checked recomputation (Table S3). It used RDKit symmetry-aware `CalcRMS` without coordinate superposition; legacy PDBQT atoms were first mapped to the reference SDF by element-constrained crystal-coordinate assignment. The prespecified best-of-nine gate is a pose-generation/search-coverage check, not a top-ranked-pose validation. Original EGFR/HER2 nine-mode production PDBQTs were not recovered; they were re-redocked under the frozen Vina protocol (seed 20260727, exhaustiveness 8, nine modes) and labeled reconstructed QC rather than historical production artifacts (Table S3).

Ligands started from frozen ChEMBL SMILES: desalt to the largest organic fragment, add explicit hydrogens in RDKit, embed with ETKDGv3 (seed 20260727), locally optimize with MMFF (at most 200 steps), and convert with default Meeko to PDBQT. Protonation states, tautomers, and conformational ensembles were not systematically enumerated. Docking used AutoDock Vina 1.2.7 with the default `vina` scoring function, nine retained poses, `energy_range = 3` kcal mol\(^{-1}\), and random seed 20260727 (Table S1). To test scoring-function dependence, the same Vina poses were rescored with RTMScore (`rtmscore_model1`, best of nine) and GNINA 1.3.2 CNN (`--cnn_scoring rescore --minimize`, best of nine after Open Babel SDF conversion). Vina’s primary readout is the mode-1 energy; RTM and GNINA CNN are best-of-9 rescores. The primary endpoint remains Vina.

Separately, GNINA 1.3.2 was run in docking-search mode on EGFR/HER2 and PIK3CA/mTOR, generating new poses from the same frozen Meeko ligand PDBQT files, receptor coordinates, boxes, exhaustiveness (8 and 16, respectively), nine retained poses, and seed 20260727. The readout is mode-1 `minimizedAffinity`. Two ligands failed on both pockets (EGFR/HER2 neither ligand EH120_109; PIK3CA/mTOR A-only ligand PM48_19) and were omitted from analyses requiring complete scores. This protocol asks whether the formulation effect persists when pose generation is changed; it is not a multi-engine bake-off (Table S32).

Vina seed sensitivity repeated the four primary panels with four additional seeds (20260811–20260814) while holding prepared ligands, receptors, boxes, exhaustiveness, mode count, energy range, and analysis rules fixed. The production seed (20260727) remained the primary analysis. Across the five frozen seeds, we report the median, interquartile range, and range of each directional AUROC, `summary_min`, Dual-versus-neither AUROC, and their difference; no seed was selected after observing performance (Table S54).

### 2.4 Primary endpoint and statistics

Throughout this Article, “dual-target recognition” names this computational discrimination task. Two binary AUROCs were computed per pair. Dual versus A-only used the pocket B score, \( \mathrm{AUC}_{D/A} = \mathrm{AUROC}(\text{dual},\;\text{A-only};\;S_B) \). Dual versus B-only used the pocket A score, \( \mathrm{AUC}_{D/B} = \mathrm{AUROC}(\text{dual},\;\text{B-only};\;S_A) \). Dual is always the positive class. Vina reports \(E_{\mathrm{Vina}}\) (kcal mol\(^{-1}\); more negative is more favorable); \(S_{\mathrm{Vina}} = -E_{\mathrm{Vina}}\).

The worst-direction discrimination summary is \( \mathrm{summary}_{\min} = \min(\mathrm{AUC}_{D/A},\;\mathrm{AUC}_{D/B}) \). It conservatively summarizes the two directional AUROCs but is not a scoring function, calibrated probability, or biological-activity estimate. Taking the minimum of two noisy estimates introduces downward selection bias; the two component AUROCs therefore remain the inferential quantities and are always reported. Arithmetic, geometric, and harmonic means are aggregation sensitivities (Table S26). The designated descriptive endpoint is pocket-matched Vina `summary_min` under unified θ = 6.0 (Table 2; PIK3CA/mTOR uses PM48). A prespecified RDKit panel (heavy-atom count, molecular weight, cLogP, and TPSA) was evaluated with the same directional workflow; the highest AUROC among them is a best single-descriptor reference (Tables 2, S28, S19). Dual versus neither (experimental inactives; `vina_mean`) and Dual versus all non-duals are formulation contrasts on the same frozen scores (Table 3; Table S22). Because these contrasts change both negative-set composition and score aggregation, their differences are descriptive and are not paired tests of a single estimand. PIK3CA/mTOR neither n = 4 is flagged underpowered.

Primary uncertainty reports ligand-level bootstrap on the two directional arms (and on `summary_min` where reported): ligands were resampled with replacement, preserving class structure (\(B = 2000\), seed 20260729, percentile 95% CI). Document-cluster and scaffold-cluster bootstrap on the same directional AUROCs are prespecified sensitivity estimators that resample literature-connected groups and Bemis–Murcko scaffold groups rather than individual ligands (Supporting Information: `document_cluster_bootstrap_v1.csv`, `scaffold_cluster_bootstrap_v1.csv`). Paired contrasts used the same resample (Tables S17, S19). Arms or cells that cannot be stably estimated are reported as not stably estimable rather than imputed, using the same rule applied to document-blocked cross-validation. To isolate negative-class choice from score aggregation, Dual-versus-selective and Dual-versus-neither AUROCs were also compared with the same pocket score under a shared resample of dual observations and independent resamples of the two negative classes (Table S34). A simulation-based detectable-effect analysis used the observed class sizes under the same bootstrap, a binormal score model, and a grid of true AUROCs (Table S31; Figure S6).

### 2.5 Confounder, holdout, and receptor-sensitivity analyses

Scores for targets A and B were swapped as a falsification control, leaving ligands, receptors, and all other settings unchanged. Directional AUROCs were also recomputed after ligand-efficiency normalization (\(S_{\mathrm{dock}}/N_{\mathrm{heavy}}\)) and on potency- or size-constrained subsets (\(|\Delta\mathrm{pChEMBL}| \leq 0.5\); \(|\Delta N_{\mathrm{heavy}}| \leq 2\)). Logistic models compared docking alone with docking plus heavy-atom count and TPSA. Morgan/ECFP4 fingerprints (radius 2, 2048 bits) with logistic regression provided a ligand-only chemical baseline under Bemis–Murcko scaffold `GroupKFold` (Tables S5, S20, S23, S24). Nearest-neighbor Tanimoto subsets were diagnostic. A coarse contact count and whole-chain sequence identity were exploratory only (Tables S7, S11).

As a post-hoc chemical-space audit, each experimental state was summarized by median and interquartile range for molecular weight, heavy-atom count, cLogP, TPSA, formal charge, and rotatable-bond count; Bemis–Murcko scaffold count and singleton-scaffold fractions were also reported (Table S38). Maximum within-pair ECFP4 Tanimoto to the nearest dual ligand was calculated for non-dual ligands. Because the descriptors are correlated and the audit was post hoc, these summaries were not converted into a battery of unadjusted significance tests.

To test dependence on exact panel membership, an unused-pool resample was drawn after excluding all ChEMBL entries used in the main panels and PM110. Ligands still come from the same ChEMBL harvest, target pairs, and label rules; this is an internal panel-membership sensitivity, not external validation. The resample was built for PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB (20 dual / 20 A-only / 20 B-only; `HOLDOUT_SEED = 20260731`); EGFR/HER2 was not eligible. Receptor, box, ligand preparation, exhaustiveness, scoring, and statistics matched the main benchmark (Tables S8, S13).

A document-blocked analysis used the same frozen scores. Ligands that share any retained high-confidence ChEMBL `document_id` were connected into one group, so compounds from the same paper cannot appear in both training and test folds. `GroupKFold` used those groups; ECFP4, physicochemical-descriptor, and docking logistic models shared the identical folds (Tables S39, S40). Folds lacking both classes were dropped. If fewer than two valid folds remained, the arm was reported as not stably estimable; the grouping rule was not changed after AUROC was seen.

A literature-year split was frozen before AUROC was computed (`docs/TIME_SPLIT_PROTOCOL_FREEZE.md`). A ligand’s year is the earliest `document.year` among its retained high-confidence records. The primary cutoff is 2018 (train: first year < 2018; test: first year ≥ 2018); 2015 and 2020 are pre-specified sensitivities. Late ligands were not used to choose thresholds, receptors, or endpoints. Directional AUROC on already-scored test ligands is reported only if dual, A-only, and B-only each have n ≥ 10; smaller cells are counts only (Table S41).

Priority assay-context fields were extracted for 186 of 352 scored ligands (`assay_context_audit.csv`); a subsequent metadata review filled include/exclude for all 186 priority ligands (179 include / 7 uncertain / 0 exclude) from organism and assay-type fields while ChEMBL assay free-text was unavailable. No frozen DualFourClass label was changed, so Table 2 was not recomputed. This pass is not paper-level assay harmonization.

Receptor-structure sensitivity used alternate crystals that, before scores were seen, (i) matched the true target protein, (ii) contained a small-molecule cognate in the ATP or target site, (iii) had acceptable resolution, and (iv) passed the same cognate redocking QC. Structures docked were PIK3CA 4JPS and 5DXT and mTOR 4JSX. Replacement was one pocket at a time on PIK3CA/mTOR (PM48) and PIK3CA/PIK3CB (Table S10). Rigid Cα superposition and a coarse pocket-residue occupancy snapshot on PM48 PIK3CA crystals were exploratory geometric controls (Tables S10, S33).

As a post-hoc exploratory supply census, the frozen J0 candidate-pair list was recounted under θ = 6.0 four-state rules without additional docking (Table S44). Post-hoc AND-like dual pocket filters and full-map ligand-only ECFP4 logistic models on the four frozen ChEMBL maps are reported in Tables S46 and S47, respectively. Multivariate property matching on the frozen scored panels used 1:1 greedy assignment in z-scored MW/cLogP/TPSA/heavy-atom space with Euclidean calipers of 0.5 and 1.0 SD (Table S45).

MCL1/Bcl-xL was formally demoted from the primary evaluation set after prespecified LC6 topology-aware pose-gold validation was not established for the nominated receptors. Panel docking outputs are retained only as an exploratory archive in the repository; they are not reported as a fifth main pair, do not extend the disparate-fold domain claim (AChE/BChE already occupies the non-kinase role), and are not used as confirmatory evidence.

Analyses ran under Python 3 with RDKit 2026.3.1, meeko 0.7.1, AutoDock Vina 1.2.7, GNINA 1.3.2, and RTMScore. Panels, scores, scripts, and parameter tables are available in the public repository (Data and Software Availability). The evaluation contract is `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`.

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

## 4. Discussion

### 4.1 Benchmark formulation changes the evidentiary standard for dual-target docking

Negative-class definition materially changed the EGFR/HER2 result: Dual versus neither gave AUROC 0.756, whereas directional `summary_min` was 0.430 (Table 3). When aggregation was held fixed by using the same pocket A score, replacing B-only with neither negatives increased AUROC by 0.378 [0.205, 0.547] (Table S34). Independent GNINA pose generation and all five Vina seeds preserved a positive formulation gap for this pair (Tables S32, S54). The other pairs did not reproduce a gap of comparable magnitude, making the finding a pair-specific failure mode rather than a general property of dual-target docking.

Relative to Zhou et al.,[9] the comparison asks whether docking separates dual-actives from experimentally defined selectives, not only from inactives. Existing docking benchmarks have shown that decoy construction and chemical bias change virtual-screening interpretation;[5–7,12,13] the same concern applies when dual-target conclusions depend on the experimental negative class. Post-hoc AND-filter and full-map ECFP4 diagnostics in the Supporting Information reinforce that Dual versus neither is chemically easier than Dual versus selectives; they are not confirmatory primary results and do not expand docking to undocked pairs.

### 4.2 What docking adds—and does not add—beyond ligand chemistry

Physicochemical descriptors and chemotype carried substantial experimental-label information. On AChE/BChE, TPSA alone exceeded docking on the corresponding contrast, and scaffold-grouped ECFP4 exceeded docking on several arms (Results 3.3). Adding the pocket-matched docking score to ECFP4 changed scaffold-grouped CV AUROC by at most 0.020. Within these panels, ligand series, and receptor realizations, docking therefore supplied little measurable incremental discrimination beyond 2D chemistry. Ligand-only controls are consequently needed before attributing an apparent dual-target signal to receptor-specific complementarity.[7,12]

Document- and scaffold-cluster bootstrap keep the EGFR/HER2 weak arm near chance with intervals that span 0.5, so the formulation interpretation does not rest on ligand-independence assumptions alone. PIK3CA/mTOR’s four-member neither class is a single-document sample, and some document-blocked arms are not stably estimable; those cells are reported as such rather than imputed.

### 4.3 Receptor realization and evaluation conditions

Receptor realization was also a performance variable. Holding one pocket fixed and replacing the PIK3CA crystal raised discrimination for PIK3CA/PIK3CB but lowered it for PIK3CA/mTOR (Figure 5), consistent with kinase cross-docking studies.[14] A single receptor structure is therefore insufficient to support a robustness claim.

### 4.4 Implications for dual-target virtual screening

Favorable scores in both pockets do not automatically establish experimentally defined dual activity. After a dual-pocket score looks favorable, four practical checks remain: (i) directional discrimination against A-only and B-only hard negatives; (ii) whether a ligand-only ECFP or property model recovers a similar signal under a leak-resistant split; (iii) an unused ligand pool or document-blocked split; (iv) at least one alternate receptor realization (Figure 8). A failure at any step marks the claim as formulation-, chemistry-, panel-, or receptor-dependent computational evidence. This concern is consistent with JCIM work showing that docking rescoring performance varies across experimentally grounded screening sets.[15]

### 4.5 Limitations

The four target pairs are a data-constrained case panel, not a representative suite, and differences in panel construction prevent interpreting `summary_min` as a target ranking. All four primary-seed intervals include 0.5. Requiring measurements at both targets enriches jointly profiled chemistry: complete-case fractions were 14.5%–34.0%. Neither the unused-pool resample nor the failed 2018 time split provides external validation, and the BindingDB-native rebuild yielded no pair meeting the frozen external gate (Tables S41, S48–S49).[16]

The activity labels also combine heterogeneous assays, use maximum pChEMBL for primary aggregation, and lack resolved construct or mutation annotations for the audited priority set. Cognate best-of-nine QC establishes search coverage rather than top-ranked-pose accuracy. MCL1/Bcl-xL was therefore retained only as an exploratory stress test after its topology-aware pose-gold gate could not be established. These limitations restrict inference to the processable compounds, receptors, and engines examined here.

## 5. Conclusions

Evaluating dual-target docking only against neither-active compounds can conceal failure to reject single-target selectives. The EGFR/HER2 contrast persisted under independent GNINA pose generation and five Vina seeds, whereas its magnitude was not reproduced across the other pairs. Docking added at most 0.020 AUROC beyond scaffold-grouped ECFP4, and receptor substitutions shifted two PIK3CA-related results in opposite directions.

Accordingly, dual-target virtual-screening studies should report both directional selective controls, a ligand-only baseline, grouped or holdout sensitivity, and alternative receptor realizations. The present four-pair analysis identifies evaluation failure modes but, with all primary `summary_min` intervals spanning 0.5 and no eligible BindingDB external pair, does not establish target-general performance or external transfer.

## Data and Software Availability

Benchmark membership, experimental-state labels, receptor and docking-box definitions, per-ligand docking scores, analysis tables, and all scripts used to regenerate the reported statistics and figures are available in the `Dual_Target_Docking` directory of the public repository at https://github.com/1280602962-debug/gwj260531. `data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` indexes the principal results through Table S53; five-seed Vina scores and summaries for Table S54 are under `data/jcim_multiseed_v0/tables/`. The archive also includes the evaluation and external-slice contracts, reconstructed EGFR/HER2 cognate QC, independent GNINA pose-generation scores, receptor and grouped-resampling sensitivities, the failed BindingDB-native gate, and the demoted MCL1/Bcl-xL stress test. SHA-256 checksums of manuscript-facing tables are in `REVISION_CHECKSUM_MANIFEST_v1.csv`. The ChEMBL supply audit was frozen on 2026-07-23; the high-confidence activity view was fetched on 2026-08-26; the BindingDB-native archive lock is release 202608. A versioned GitHub Release and Zenodo DOI will be issued from a tagged snapshot. The BindingDB TSV archives themselves are not redistributed.

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
