# Benchmark Formulation and Chemical Confounding in Docking-Based Dual-Target Recognition

## Abstract

Whether favorable docking scores at two targets constitute evidence of dual-target recognition has not been adequately tested against experimentally defined single-target selectives. We constructed DualFourClass-Bench, a curated four-pair, four-state panel with two directional primary tasks: dual-actives versus A-only selectives scored in pocket B, and dual-actives versus B-only selectives scored in pocket A; the weaker arm is the pair-level summary (`summary_min`). On the same frozen AutoDock Vina scores, a Dual-versus-neither comparator produced a substantially stronger impression on EGFR/HER2 (AUROC 0.756 versus directional `summary_min` 0.430), while 9 of the Top-10 mixed-library ligands were experimental selectives. AChE/BChE and PIK3CA/PIK3CB showed only small, overlapping formulation increments, and the PIK3CA/mTOR Dual-versus-neither comparator was underpowered (neither n = 4). Under scaffold-grouped cross-validation, adding docking changed AUROC by at most approximately 0.02 in absolute value beyond ECFP4. Replacing maximum pChEMBL with the median among repeated measurements produced minimal pair-level changes. In contrast, alternative receptor realizations moved PIK3CA/mTOR from 0.692 to 0.486/0.505 but PIK3CA/PIK3CB from 0.500 to 0.691/0.685. Thus, apparent dual-target discrimination depended on benchmark formulation, target pair, ligand chemistry, and receptor realization, with limited incremental information from docking beyond ligand-level chemical baselines in this scaffold-aware evaluation. These results support experimentally defined selectivity hard negatives and confounder-aware controls; they do not establish a universal overestimation law or prove that docking lacks pocket-specific information.

**Keywords:** dual-target docking; benchmark formulation; selectivity hard negatives; chemical confounding; receptor realization; virtual screening

## 1. Dual-target design and the role of structure-based virtual screening

Multitarget drug design aims to modulate two or more biological targets with a single small molecule, in order to address pathway redundancy, compensatory signaling, and drug resistance in complex disease. Relative to a classical single-target agent, a rationally designed multitarget ligand may act on connected nodes of a disease network and thereby achieve a more adequate pharmacological effect; this idea is now a central theme of polypharmacology.[1] Over the past decade the field has moved from largely serendipitous multi-pharmacology toward structure-guided design that combines structural biology, computational chemistry, and, increasingly, generative models.[2]

Molecular docking remains one of the most widely used tools in structure-based virtual screening (SBVS). A docking engine places the ligand in a protein binding site and a scoring function ranks ligand–receptor complementarity, providing a rapid, structure-level filter of large compound libraries.[3,4] A natural computational tactic in dual-target discovery is therefore to dock each candidate into both pockets and to treat favorable scores on both targets as evidence of dual-target potential.

How such scores should be interpreted, however, is already known to depend on how the data set is built: the definition of negatives, hidden chemical bias, and the choice of metric all change apparent performance. The Directory of Useful Decoys (DUD) and its successor DUD-E frame single-target docking as enrichment of actives against property-matched decoys, precisely because unmatched decoys can make enrichment look like separation of crude ligand properties.[5,6] LIT-PCBA showed that artificially constructed active/decoy collections, including DUD, DUD-E, and MUV, can contain obvious and hidden chemical biases that overestimate the true accuracy of virtual screening, and replaced them with experimental dose–response labels under physicochemical-range control.[7] Structure-scoring benchmarks such as CASF-2016 ask a different, still single-complex question: scoring, ranking, docking, and screening power of scoring functions on high-quality protein–ligand complexes drawn from PDBbind-quality data.[8] None of these resources formulates dual-target recognition over an experimentally labeled four-state ligand space.

**Extending the single-target docking evaluation logic to a dual-target task is therefore not sufficient.**

## 2. Dual-target recognition is a different task from single-target virtual screening

A conventional SBVS problem can be written as Active versus Decoy. Negatives exist mainly to set a decision boundary against the actives of one protein.

For dual-target ligands the task changes. A strict dual-target benchmark must distinguish at least four ligand states with distinct biological meaning (a four-state data set, not a four-class classifier): **dual-active** ligands that act on both targets, **A-selective** ligands that act only on target A, **B-selective** ligands that act only on target B, and **neither** ligands that lack sufficient activity on either target (Figure 1A):

|  | *B*<sup>+</sup> | *B*<sup>−</sup> |
|--|:--:|:--:|
| *A*<sup>+</sup> | Dual | A-only |
| *A*<sup>−</sup> | B-only | Neither |

A-only and B-only ligands are not ordinary decoys. They are the **selectivity hard negatives** of the dual-recognition problem: already potent on one target, they can produce plausible docking scores in that pocket, yet they lack the corresponding activity on the other. For an evaluation that claims dual-target recognition, the question is therefore not whether a candidate can obtain a favorable docking score in each pocket separately, but whether it can **distinguish true dual-actives from the matching single-target selectives in both directions at once**.

That question is stricter than prior dual-target docking evaluations. Zhou, Li, and Hou assessed docking-based virtual screening on four kinase pairs, first as inhibitor versus noninhibitor on each target and then as dual-target identification, and reported structure dependence together with a high false-positive rate among predicted duals.[9] That work established that dual-target docking can be benchmarked and that docking versus inactives is not sufficient for a clean dual hit list. It did not treat experimentally labeled A-only and B-only ligands as directional hard negatives, and it did not ask whether a Dual-versus-neither (inactive) comparator on the same scores would change the interpretation of directional Dual-versus-selective discrimination. Dual-versus-neither is used here as a **nonselectivity-controlled comparator**, not as “the conventional dual-target benchmark.”

It follows that averaging, summing, or otherwise pooling the two docking scores does not describe dual-target recognition. A ligand may score very favorably on target A and poorly on target B; the mean can still look strong, but that number does not support dual activity. Likewise, beating a reference ligand’s docking score on both targets can define a computational dual-success criterion, yet it does not answer the stricter experimental question: **can the computational score separate true dual-actives from potent, target-selective hard negatives?**

The primary benchmarking problem is therefore not to find a better docking score, but to build a benchmark and a readout that match this biological state space.

## 3. Dual-target benchmarks are jointly limited by experimental data supply and chemical confounding

Instantiating that four-state data set is hard in public data. DualFourClass-Bench retains all four experimental states, but the prespecified primary evaluation is two directional pairwise discriminations (dual versus A-only and dual versus B-only), not a four-class classifier. DualFourClass-Bench is a **four-state curated benchmark with two directional primary tasks**. Single-target virtual-screening benchmarks can rely on mature active/decoy construction.[5,6] A strict dual-target benchmark instead requires the *same* compound to carry comparable experimental measurements on both targets, and further requires A-only and B-only ligands with a clear selectivity gap. A usable experimental panel therefore needs not only enough dual-actives, but also selective hard negatives that are adequate in count and potency range on **both** arms.

That requirement raises the construction bar. Assay types, activity endpoints, experimental conditions, and coverage differ across databases, and dual-target labeling adds the extra constraint that both ends must be measured. **The number of target pairs that can support a strict four-state evaluation is itself a data-supply question that must be quantified; it cannot be assumed that every pharmacologically related pair yields a balanced benchmark.** Compounds that are quantitatively measured on both ends and separated by a selectivity gap are not automatically abundant. How many pairs pass that completeness gate is a methodological bottleneck in benchmark construction, not a post-hoc apology for the size of the evaluation set.

Dual-target panels also inherit a chemical-confounding problem that is easy to miss. Molecular weight, polarity, hydrogen-bonding features, lipophilicity, and scaffold membership can jointly shape experimental activity and docking scores. If dual-active, A-only, and B-only ligands differ systematically in these ligand-level properties, a high AUROC may reflect statistical differences among ligands rather than structure complementarity in the two pockets. LIT-PCBA established that chemical bias in constructed active/decoy sets can inflate virtual-screening performance; a dual-target docking benchmark therefore needs explicit chemical and physicochemical controls of the same kind.[7]

## 4. Dual-target generative methods make stricter docking evaluation a practical need

Structure-based generative methods have made the same evaluation gap operational. Zhou, Guan, et al. formulated dual-target drug design as a generative task and introduced DualDiff: two pockets are aligned in 3D and SE(3)-equivariant messages are composed over shared ligand nodes so that a diffusion model pretrained on single-target complexes can be transferred to the dual-target setting.[10] Their docking evaluation uses AutoDock Vina redocking and reports Vina Dock on each pocket, Max Vina Dock (the worse of the two Vina Dock values), and Dual High Affinity—the fraction of generated molecules whose estimated affinities **exceed those of the reference ligands on both targets**.[10] Dual High Affinity is computational dual success relative to references, not mean pooling of the two scores; Max Vina Dock already attends to the weaker arm.

FuseDiff jointly generates a shared ligand graph and two pocket-specific binding poses and evaluates generation on an independent DualDiff (DDF) test set, again reporting Vina Dock, Max Vina Dock, and Dual High Affinity.[11] These studies show that **deciding whether a molecule satisfies a dual-target structural requirement is now a practical problem in computational dual-target design**.

That docking-based criterion is not the same as an assay-label hard-negative task. Requiring both pocket scores to beat a reference ligand measures computational dual docking success; it does not test whether the molecule can be distinguished from ligands that are potent on only one of the two targets. Generative docking-based evaluation and the **experimentally labeled dual-versus-selective discrimination** studied here are therefore complementary, not interchangeable, benchmarks. This work does not redock DualDiff or FuseDiff molecules. DualFourClass-Bench is intended as a possible downstream check of molecules that pass a Dual High Affinity filter, not as an empirical bake-off of those generators.

## 5. Aim and contribution of this study

Existing evaluations of dual-target molecular design generally assess whether a molecule simultaneously obtains favorable scores at both targets, but such criteria do not directly test discrimination against ligands that are experimentally active at only one target. This distinction is important because a favorable score at one target can coexist with poor recognition at the other, while selective ligands may remain deceptively favorable in both docking pockets.

**Here, we ask whether the formulation of the benchmark itself changes the apparent evidence for dual-target recognition.** We construct an experimentally defined four-state ligand panel, use directional pocket-matched discrimination against A-selective and B-selective hard negatives as the primary task, and compare it with a nonselectivity-controlled dual-versus-neither comparator. We then test whether the resulting signal persists after chemical, physicochemical, ligand-pool, activity-aggregation, and receptor-structure controls.

The contribution is an evaluation protocol and a curated benchmark resource, not a new docking algorithm or scoring function. DualFourClass-Bench is a **four-state curated benchmark with two directional primary tasks**: dual versus A-only scored in pocket B, and dual versus B-only scored in pocket A (Figure 1B). Neither is retained to describe the full experimental state space; it does not enter the primary AUROCs. The pair-level summary is the weaker arm (`summary_min`), so a strong score on one target cannot hide failure on the other. Dual-versus-neither on the same scores is a comparator, not “the conventional dual-target benchmark.”

The study separates three aspects of benchmark validity within one evidence chain: task formulation (Dual versus selective rather than Dual versus neither), confounder-aware evaluation (docking alongside ligand-only and wrong-pocket controls), and evaluation-condition sensitivity (activity aggregation, unused ligand pools, and receptor realization). These analyses test whether a benchmark result behaves as a fixed property of docking or as a conditional result of the evaluation setting.

A public-data supply audit first asks how many candidate pairs can support that four-state construction (Methods 2.1–2.3). The size of the evaluation set is a result of that audit, not a design target frozen in this Introduction. Pooled docking scores, a wrong-pocket control, and two-dimensional chemical and physicochemical baselines are reported as auxiliary contrasts, to separate pocket-specific signal from ligand-level confounding.

The nested experimental question remains: to what extent can existing docking scores distinguish dual-active ligands from target-selective hard negatives on both arms, and how far does that discrimination depend on the target pair, the receptor structure, or ligand chemistry? The protocol is meant as a stricter downstream check for dual-target virtual screening and for generative dual-target design—not as a bake-off of those generators, and not as a comprehensive dual-target suite.

## 2. Methods

### 2.1 Data sources and activity curation

Ligand activities used as **experimentally derived activity labels** were retrieved from the public ChEMBL Web API activity endpoint. The target-pair supply audit was frozen on 2026-07-23. pChEMBL converts molar concentration–response measurements (IC50, Ki, Kd, EC50, and related endpoints) to an approximate −log10 activity scale for large-scale integration. Assay types, conditions, and experimental systems are not equivalent; pChEMBL is used here as a curation convenience, not as an absolute affinity measured under one protocol.

When several pChEMBL values existed for the same ligand–target pair, the **maximum** was used as the one-to-one representative for primary curation. Because assay types, conditions, and experimental systems are not equivalent, activity-aggregation sensitivity was assessed in a prespecified analysis by re-fetching assay-level records and replacing the maximum with the **median** of repeated pChEMBL measurements (Table S29). This alternative aggregation was applied across all frozen benchmark panels without changing panel membership, docking parameters, or Vina scores. Class assignment was compared under the same θ = 6.0 rule. Ligands missing a usable pChEMBL value on either target were excluded from analyses requiring paired labels.

Salt, solvate, and multicomponent ChEMBL records were split by connected component; the organic fragment with the most heavy atoms was retained.

To test whether the ChEMBL supply gate is an artifact of one database, the frozen target pairs were recounted in BindingDB and PubChem (counts only; no docking and no panel rebuild). BindingDB was queried with REST `getLigandsByUniprots` (cutoff = 1 mM, so weak-end measurements are not truncated); PubChem with PUG REST `protein/accession/…/concise`. Endpoints were restricted to IC50/Ki/Kd/EC50; the representative value was the maximum converted p-activity; class rules matched the strict supply gate in Section 2.2. Identifiers were BindingDB monomerid and PubChem CID, with **no** cross-database structure merge. The primary count used equal-relation measurements (censored `>`/`<` values dropped); treating inequalities as point estimates was a sensitivity only. Counts are reported in Supporting Information Table S12.

### 2.2 Target-pair supply audit and experimental ligand-state definition

Candidate pairs were audited for whether public data can support a strict dual-target benchmark. For each pair A/B, ligands were assigned one of four **experimental states**:

- **dual**: strong activity on both targets;
- **A-only**: strong on A, weak on B;
- **B-only**: strong on B, weak on A;
- **neither**: insufficient activity on both.

A-only and B-only ligands are selectivity hard negatives, not DUD/DUD-E-style assumed decoys.

**Strict supply-audit rule (construction gate, not the sole label for every comparison).** Dual: both pChEMBL values ≥ 6.5. A-only: A ≥ 6.5 and B ≤ 5.5. B-only is symmetric. Neither: both ≤ 5.5. The 5.5–6.5 gray zone was excluded from the strict audit. The rule asks whether both arms have enough selective hard negatives for a reasonably balanced panel. How many pairs pass, and which pairs entered the frozen set, are reported in Results 3.1. Metal-dependent systems (e.g. HDACs) were excluded in advance as unsuitable primary objects for this noncovalent protocol.

**The primary manuscript comparison uses one prespecified θ = 6.0 rule.** Dual: both ends ≥ θ; A-only: A ≥ θ and B < θ; B-only is symmetric; neither: both < θ. Construction may use this single threshold when the strict rule leaves too few selectives to fill quotas. Construction rules were frozen from the supply audit before sampling and are recorded in Table 1. Thresholds were chosen to assemble analyzable quotas, not after inspecting docking scores. As a supporting sensitivity, states were relabeled at θ ∈ {5.5, 6.5} and under the strict 6.5/5.5 rule, and pocket-matched summary_min was recomputed (Table S4). That grid is not a second primary standard. Underpowered cells are flagged in Results.

### 2.3 DualFourClass-Bench panel construction

The resource is a **four-state curated benchmark with two directional primary tasks**. Dual, A-only, B-only, and neither are all retained to describe the experimental space; the prespecified primary endpoint is dual versus A-only and dual versus B-only. Neither does not enter the primary directional AUROCs. This is not a four-class classification benchmark.

Candidate pairs were screened with the strict audit in Section 2.2. The frozen evaluation set comprises PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB, and EGFR/HER2. EGFR/HER2 is retained as a **supply-limited case** and is not treated as equivalent in supply to the other pairs.

For each pair, ligands were drawn from the labeled pool under frozen class quotas and random seed 20260729. Where structures were available at sampling, a per-class Bemis–Murcko scaffold cap limited series over-representation: at most two molecules per scaffold in PIK3CA/mTOR (PM48) and at most five in EGFR/HER2. Structures were unavailable at sampling for AChE/BChE and PIK3CA/PIK3CB, so these panels used class quotas and a deterministic shuffle without an additional diversity constraint. Post-construction Murcko scaffolds are reported with the deposited tables. Final membership, state labels, ChEMBL identifiers, SMILES, and sampling scripts are deposited. Panels were not redrawn after docking scores were seen.

Construction rules were not identical across pairs. AChE/BChE and PIK3CA/PIK3CB were sampled under the strict 6.5/5.5 gate; EGFR/HER2 and PIK3CA/mTOR used θ = 6.0 because the strict gate left too few B-only ligands. Cross-pair AUROCs therefore mix target-pair biology with panel-construction differences (sample size, threshold, series composition, receptor) and are not interpreted as purely intrinsic docking performance.

Quotas and construction labels were as follows. AChE/BChE and PIK3CA/PIK3CB: strict 6.5/5.5, target dual / A_only / B_only / neither = 28 / 28 / 28 / 16 (panel n = 100). EGFR/HER2: θ = 6.0 construction rule (n = 110). PIK3CA/mTOR: θ = 6.0, main comparison panel PM48 (n = 48; constructed 18 / 14 / 12 / 4), on which receptors and the docking protocol were frozen.

Ligand–receptor jobs that failed to yield a score were dropped for that receptor; ligands missing a usable score on either end were omitted from pocket-matched AUROCs that require both scores, so analysis counts can fall below construction quotas (Table 1). AUROC tables are therefore **conditional on compounds the docking engine can process**. Attempted / successful / failed counts, including chemical-coverage failures such as unsupported AutoDock atom type `B`, are reported in Table S27.

An expanded PIK3CA/mTOR panel (PM110) keeps all 48 PM48 ligands and adds molecules under the strict rule, targeting 30 / 30 / 30 / 25. PM110 is a superset of PM48, used to check whether point estimates stay in the same direction after increasing panel size. It is not an independent primary benchmark or replicate. Cross-pair comparison in the main text uses PM48.

Ligand-side generalization was assessed with one unused-pool holdout (Section 2.11), because the remaining hard-negative supply could not support a distribution of non-overlapping balanced panels. Ligand-level bootstrap (Section 2.8) instead quantifies uncertainty within a fixed panel and is not unused-pool resampling.

**Table 1.** DualFourClass-Bench composition and docking settings. Construction labels record the supply/panel-building rule for each pair. All primary AUROCs in Tables 2–3 use unified θ = 6.0 experimental-state labels; strict 6.5/5.5 is a supply and construction gate, with relabeling reported as sensitivity in Table S4.

| Pair | Construction labels | PDB (A / B) | Resolution (Å) | Panel n | Analysis n (dual / A_only / B_only) | Vina exhaustiveness |
|------|---------------------|-------------|----------------|-------:|------------------------------------:|--------------------:|
| PIK3CA/mTOR | θ = 6.0 | 4L23 / 4JT6 | 2.50 / 3.60 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | strict 6.5/5.5 | 4EY7 / 4BDS | 2.35 / 2.10 | 100 | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | strict 6.5/5.5 | 4L23 / 2WXF | 2.50 / 1.90 | 100 | 28 / 27 / 28 | 8 |
| EGFR/HER2 | θ = 6.0 | 3POZ / 3RCD | 1.50 / 3.21 | 110 | 28 / 38 / 32 | 8 |

### 2.4 Protein structures and binding-site definition

Receptors were PDB entries with experimental structures and a small-molecule cognate ligand: PIK3CA/mTOR, 4L23 / 4JT6 (X6K / PI-103); AChE/BChE, 4EY7 / 4BDS (E20 / THA); PIK3CA/PIK3CB, 4L23 / 2WXF (X6K / 039); EGFR/HER2, 3POZ / 3RCD (03P / TAK-285). Resolutions are in Table 1.

The site was defined from the cognate ligand. An axis-aligned bounding box on cognate heavy atoms was expanded by 5 Å on each axis; any edge shorter than 20 Å was set to at least 20 Å. Box centers and sizes are frozen in JSON and listed in Supporting Information Table S2.

Water and the cognate ligand were removed and Meeko wrote PDBQT. PIK3CA, mTOR, EGFR, and HER2 used hydrogen-containing protein coordinates already in the frozen directories (`mk_prepare_receptor.py --read_pdb`). AChE, BChE, and PIK3CB were extracted from deposited ATOM/TER records (waters and hetero atoms removed) and prepared with `mk_prepare_receptor` (default alternate location A). PDBFixer was not used to rebuild missing atoms, nor Reduce for independent pH-dependent protonation or histidine tautomer enumeration; protonation is part of the frozen preparation. Docking treated noncovalent small-molecule sites only; metals and other cofactors were not extra dockable components in the box.

### 2.5 Cognate redocking quality control

Before production docking, each frozen receptor was redocked with its cognate ligand to test whether the box, receptor preparation, and search settings can generate a near-native pose **within the retained pose ensemble**.

Nine poses were generated and heavy-atom RMSD to the crystal ligand was computed in the docking frame (no protein superposition). PIK3CA/mTOR and EGFR/HER2 used meeko `REMARK SMILES IDX` mapping and the minimum CalcRMS over graph automorphisms; AChE/BChE and PIK3CB used Hungarian matching on heavy atoms. Define

\[
\mathrm{RMSD}_{\mathrm{best9}} = \min_{i=1,\ldots,9} \mathrm{RMSD}_i.
\]

The prespecified pass criterion was \(\mathrm{RMSD}_{\mathrm{best9}} < 2.0\) Å.

This QC tests **pose-generation capability**: whether a near-native pose appears among the retained poses. It does **not** require the top-ranked Vina pose (mode 1) to be near-native. Best-of-9 QC and mode-1 scoring are different evaluations.

If the default exhaustiveness failed the gate, search effort was raised to a prespecified fallback without changing the box, receptor, or random seed, and QC was repeated. Production docking therefore used receptor-specific frozen exhaustiveness: 16 for PIK3CA/mTOR and 8 for the other main panels. QC values and mode-1 versus best-of-9 comparisons are in Supporting Information Table S3.

### 2.6 Ligand preparation and molecular docking

Ligands started from frozen ChEMBL SMILES: desalt to the largest organic fragment, add explicit hydrogens in RDKit, embed with ETKDGv3 (seed 20260727), locally optimize with MMFF (at most 200 steps), and convert with default Meeko to PDBQT. Protonation states, tautomers, and conformational ensembles were not systematically enumerated. Schrödinger LigPrep was not used.

Docking used AutoDock Vina 1.2.7 with the default `vina` scoring function. Nine poses were retained per ligand–receptor pair, `energy_range = 3` kcal mol\(^{-1}\), random seed 20260727 (same as ETKDG). Exhaustiveness followed Table 1. Ligand preparation, box rules, and the scoring function were shared across main panels; only receptor coordinates, box numbers, and the prespecified exhaustiveness varied. Full parameters are in Supporting Information Table S1.

### 2.7 Alternative scoring channels

To test dependence on one scoring function, the **same Vina-generated poses** were rescored with RTMScore and GNINA CNN.

RTMScore used public weights `rtmscore_model1` and scored all nine Vina poses; the highest RTMScore in that pocket was kept.

GNINA 1.3.2 ran CNN rescoring on CPU (`--cnn_scoring rescore --minimize`). The final protocol converted each of the nine Vina poses to SDF (Open Babel) and took the highest CNNscore per pocket, matching RTM pose coverage. Mode-1-only GNINA scores are retained as a historical sensitivity control, not the channel readout.

Vina’s primary readout is the mode-1 energy; RTM and GNINA are best-of-9 rescores. The three channels do not aggregate the nine poses in the same way and are **not** a head-to-head docking-engine competition. They are a scoring-channel sensitivity analysis. The primary endpoint is always Vina.

### 2.8 Primary endpoint and statistical analysis

#### 2.8.1 Pocket-matched directional AUROC

Two binary AUROCs were computed per pair. Dual versus A-only used the **pocket B** score:

\[
\mathrm{AUC}_{D/A} = \mathrm{AUROC}(\text{dual},\;\text{A-only};\;S_B),
\]

asking whether docking can use structural information at the nonselective target B to separate dual-actives from A-only ligands that are already potent at A. Dual versus B-only used the pocket A score:

\[
\mathrm{AUC}_{D/B} = \mathrm{AUROC}(\text{dual},\;\text{B-only};\;S_A).
\]

Dual is always the positive class. Neither does not enter these contrasts.

Vina reports \(E_{\mathrm{Vina}}\) (kcal mol\(^{-1}\); more negative is more favorable). Define

\[
S_{\mathrm{Vina}} = -E_{\mathrm{Vina}}
\]

so that larger primary scores mean stronger predicted binding. RTMScore and GNINA CNN scores are already higher-better.

#### 2.8.2 summary_min

The pair summary is the weaker arm:

\[
\mathrm{summary}_{\min} = \min(\mathrm{AUC}_{D/A},\;\mathrm{AUC}_{D/B}).
\]

This is a task-aligned **worst-arm aggregation**, not a new scoring function. The minimum prevents a strong arm from hiding failure on the other; it is not the unique statistically natural aggregator. Arithmetic, geometric, and harmonic means of the two directional AUROCs are reported as sensitivity aggregators (Table S26). Pair ranking and the direction of the EGFR Dual-versus-neither contrast are unchanged under all four aggregators. The single primary endpoint is pocket-matched Vina `summary_min` under unified θ = 6.0 (Table 2; PIK3CA/mTOR uses PM48). Prespecified secondary endpoints are the two directional arms, pocket-matched RTMScore, pocket-matched GNINA CNN best-of-9, and the descriptor panel in Section 2.8.3. Sensitivity / falsification endpoints are the θ grid, PM110, E = 8, unused-pool holdout, receptor replacement, and wrong-pocket (including paired Δ). Exploratory endpoints are ECFP4, contact_count (not a PLIF), and Top-10 hard-negative counts on pooled `vina_mean`. The hierarchy is in Supporting Information Table S16. Pooled `vina_mean` directional AUROC is **not** Table 2.

#### 2.8.3 Physicochemical descriptor controls

A **prespecified** RDKit panel was computed: heavy-atom count, molecular weight, cLogP, and TPSA. Each descriptor was evaluated with the same directional AUROC workflow; **all four are reported** (Table 2; Table S28). The highest AUROC among them is a **best single-descriptor reference** — a descriptive post-hoc maximum, not a confirmatory competitor. Paired Δ between docking and that reference is not a confirmatory test of “beats the selected best descriptor” (Table S19).

#### 2.8.4 Score-aggregation controls

Pooled means of the two pocket scores, wrong-pocket assignment (Section 2.9.1), and worst-pocket aggregation were computed as auxiliaries, not as the primary endpoint (Table S6).

#### 2.8.5 Bootstrap uncertainty

AUROC and summary_min uncertainty used ligand-level bootstrap: ligands were resampled with replacement, preserving class structure, and both directional AUROCs and summary_min were recomputed. \(B = 2000\), seed 20260729, percentile 95% CI \([P_{2.5}, P_{97.5}]\). Paired contrasts used the **same** resample to form \(\Delta = \mathrm{Metric}_1 - \mathrm{Metric}_2\) (Tables S17, S19). Murcko-scaffold resampling is reported as a control; the text uses ligand-level intervals. Intervals are descriptive. Outside the prespecified primary endpoint, this work does not treat “whether the CI crosses 0.5” as a formal significance test across many pairs and controls.

#### 2.8.6 Benchmark-formulation comparison

As an auxiliary contrast on the **same** frozen Vina scores, a **Dual-versus-neither comparator** (experimental inactives; `vina_mean` and `vina_worst`) and Dual versus all non-duals were computed beside the directional primary endpoint. Dual-versus-neither is a **nonselectivity-controlled comparator** on this panel, not a claim that established dual-target benchmarks use Dual versus neither as their official task. Neither ligands are used here; they still do not enter Table 2. PIK3CA/mTOR neither n = 4 is flagged underpowered. The comparison asks whether omitting selective hard negatives can change the apparent evidence for dual-target recognition; it is not a second primary endpoint and is not a paired significance test (different negative sets; Table 3; Table S22). Single-target-style analogues—(dual + A-only) versus (B-only + neither) in pocket A, and the symmetric B contrast—are reported in Table S22.

### 2.9 Confounder and falsification analyses

#### 2.9.1 Wrong-pocket falsification control

Scores for targets A and B were swapped; ligands, receptors, and all other settings were unchanged. Directional AUROCs and summary_min were recomputed. This is a **falsification control**, not a positive control designed to prove pocket specificity. Matched > wrong on a fixed panel is **not** taken as evidence of pocket-specific signal. Wrong-pocket performance near or above matched-pocket performance counts against a pocket-specific reading. Holdout reversal further means wrong-pocket is **not a reliable universal negative control under panel shift**.

#### 2.9.2 Ligand-efficiency normalization

Each pocket score was divided by heavy-atom count, \(S_{\mathrm{LE}} = S_{\mathrm{dock}} / N_{\mathrm{heavy}}\), and directional AUROCs and summary_min were recomputed.

#### 2.9.3 Potency- and size-matched subsets

Subsets with \(|\Delta\mathrm{pChEMBL}| \leq 0.5\) or \(|\Delta N_{\mathrm{heavy}}| \leq 2\) were formed and directional AUROCs recomputed. Matching reduces n; the analysis asks whether the direction changes, not whether a small subset is independent strong evidence (Table S5).

#### 2.9.4 Covariate-adjusted analysis

Logistic regression compared

\[
\mathrm{Model}_1:\ Y \sim S_{\mathrm{dock}}, \qquad
\mathrm{Model}_2:\ Y \sim S_{\mathrm{dock}} + N_{\mathrm{heavy}} + \mathrm{TPSA},
\]

where \(Y\) is the dual versus selective-hard-negative label. scikit-learn `LogisticRegression` was used (\(C = 1.0\), `max_iter = 2000`). Model AUROC, the docking coefficient, and its odds ratio are reported. The question is residual discrimination after size and polarity, not a new primary predictor.

#### 2.9.5 Two-dimensional chemical baseline

Morgan/ECFP4 fingerprints (radius 2, 2048 bits) with the same logistic settings provided a ligand-only chemical baseline. Evaluation used Bemis–Murcko scaffold `GroupKFold` with \(K = \min(5, N_{+}, N_{-}, N_{\mathrm{scaffold}})\) and at least two folds, so the same scaffold does not span train and test. High CV AUROC therefore means discrimination remains when molecules from the same Murcko scaffold are not shared between folds; it is **not** target-external generalization. On PIK3CA/mTOR, \(n_{\mathrm{scaffolds}} \approx n\), so the split is nearly leave-one-scaffold. Random `StratifiedKFold` is a leakage check only (Table S20). Incremental models (physchem, ECFP4, docking, and combinations) use the same split; logistic docking AUROC is not the rank AUROC in Table 2 (Table S24). Nearest-neighbor ECFP4 Tanimoto matching of A-only/B-only ligands to duals is reported at T ≥ 0.3 / 0.4 / 0.5 because T ≥ 0.7 matching is empty on these panels (Table S23). T ≥ 0.3 is a **similarity-constrained subset**, not a chemically matched analogue set.

#### 2.9.6 Scoring-independent contact count

On frozen Vina **mode-1** poses, a scoring-free geometric descriptor was the number of ligand heavy atoms within 4.0 Å of any receptor heavy atom:

\[
N_{\mathrm{contact}} = \#\{i:\ \min_j d_{ij} \le 4.0\,\text{Å}\}.
\]

This descriptor does not use the docking energy function. \(N_{\mathrm{contact}}\) was used for dual versus A-only in pocket A and dual versus B-only in pocket B, isomorphic to the wrong-pocket pocket-wise comparison, as a geometric confounder control (size/burial). The 4.0 Å cutoff is coarse and **not** a PLIF. Magnitude agreement with Vina wrong-pocket is not assumed (Table S11).

#### 2.9.7 Cross-pair sequence identity (exploratory)

The longest protein chain was read from each frozen `*_protein.pdb` (standard amino-acid ATOM records). Pairwise global identity used Biopython `PairwiseAligner` (BLOSUM62, gap open = −11, extend = −1), normalized by alignment length or by the shorter chain (Table S7). The quantity is a coarse whole-chain similarity proxy, not pocket residue correspondence, pocket RMSD, or a PLIF.

### 2.10 Single-target enrichment reference

On PIK3CA 4L23 and mTOR 4JT6, single-target active versus weak-active sets were built. Actives: pChEMBL ≥ 6.5. Weak actives: measured on the same target with pChEMBL ≤ 5.5 and property-matched to actives within ±50 Da (MW), ±1.5 (cLogP), and ±25 Å² (TPSA). MW and logP windows follow common property-matched decoy practice (Mysinger et al., *J. Med. Chem.* **2012**, *55*, 6582–6594); TPSA is an added polarity match. Target size was about 50 actives and 150 weak actives. Preparation, receptor, box, and Vina settings matched the PIK3CA/mTOR main panel (exhaustiveness = 16). AUROC, EF1%, and EF5% are reported as a single-target enrichment backdrop, not as a substitute for dual-target summary_min.

### 2.11 Unused-pool holdout

To test dependence on the exact frozen panel members, all ChEMBL entries used in the main panels and in PM110 were excluded from the strict-label pool, and an **unused-pool, panel-external holdout** was drawn from what remained. This is not cross-database or cross-assay external validation: ligands still come from the same ChEMBL harvest, the same target pairs, and the same label rules.

Holdout was built only for pairs whose unused pool could supply 20 dual, 20 A-only, and 20 B-only ligands: PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB. EGFR/HER2 was not eligible for an equivalent draw and was not patched with a non-equivalent sample. PIK3CA/mTOR exclusion used the PM110 superset, which covers PM48. The draw used `HOLDOUT_SEED = 20260731` (distinct from the construction seed) and a Murcko cap of three members per state class. The list was frozen before docking.

The holdout did not enter panel construction, protocol tuning, or primary-endpoint choice. Receptor, box, ligand preparation, exhaustiveness, scoring, and statistics matched the main benchmark, including `summary_min` and ligand-level bootstrap. Jobs that produced no Vina score were dropped as in Section 2.3. Descriptor, wrong-pocket, matched-subset, and contact-count analyses were repeated on the holdout (Tables S8, S13). Matched-subset diagnostics do not rewrite Table S8.

### 2.12 Receptor-structure sensitivity analysis

To test sensitivity of the benchmark conclusion to receptor choice, alternate crystals were required, **before scores were seen**, to (i) have a polymer entity that is the true target protein (no chimeras or off-target scaffolds), (ii) contain a small-molecule cognate in the ATP or target site, (iii) have acceptable resolution, and (iv) pass the same cognate redocking QC as Section 2.5. Structures actually docked were PIK3CA 4JPS and 5DXT and mTOR 4JSX. This is a **receptor-structure sensitivity analysis** (a receptor-realization effect), not a robustness check and not a certification of a “more correct” crystal. The purpose was to quantify sensitivity of the dual-target discrimination endpoint to receptor realization rather than to identify a superior receptor structure. PIK3CA/mTOR is not treated as a prespecified structure-invariant positive case.

Replacement was **one pocket at a time**. On PIK3CA/mTOR (PM48), 4JPS/5DXT replaced pocket A while pocket B kept frozen 4JT6 scores; 4JSX replaced pocket B while pocket A kept frozen 4L23 scores. Exhaustiveness was 16, matching the PM48 main panel. On PIK3CA/PIK3CB, the same prepared 4JPS and 5DXT receptors replaced pocket A while pocket B kept frozen 2WXF scores; exhaustiveness was 8, matching that main panel. New boxes used that crystal’s own cognate ligand and the Section 2.4 AABB rule. Ligand preparation, seed (20260727), scoring function, and the primary endpoint matched the corresponding main analysis. Jobs that produced no Vina score were dropped as in Section 2.3; attempted / successful / failed counts are reported with the swap tables.

As an exploratory, docking-free geometric control, rigid superposition used Biopython `Superimposer` on residue-matched Cα atoms of the longest chain (Kabsch). Pocket residues were those within 5 Å of the reference cognate heavy atoms; local pocket RMSD used the **same** transform, with no second local fit. Alternate cognate centroids were projected in that frame. Matched Cα counts can differ, so global RMSDs are not equal-coverage comparisons. The set of alternates is small; Cα RMSD is not assumed to explain AUROC changes quantitatively (Table S10).

### 2.13 Software and data availability

Analyses ran under Python 3 with RDKit 2026.3.1, meeko 0.7.1, AutoDock Vina 1.2.7, GNINA 1.3.2, and RTMScore (`rtmscore_model1`); Open Babel converted Vina poses to SDF. Superposition and chain alignment used Biopython. AUROC, logistic regression, and cross-validation used NumPy, SciPy, scikit-learn, and pandas. Panels, scores, scripts, and parameter tables will be released with the public pack (Data and Software Availability).

## 3. Results

### 3.1 Experimental Data Supply Limits Strict Dual-Target Benchmark Construction

To determine whether public bioactivity data can support a strict dual-target recognition evaluation, we first audited ligand supply for 49 ChEMBL-cached candidate target pairs (Figure 2). Dual-target docking evaluation requires four experimentally labeled ligand states: dual, A-only, B-only, and neither (Figure 1A). The benchmark is four-state with two directional primary tasks; neither is curated but does not enter the primary AUROCs. Ligands that meet the activity threshold on one target and are explicitly inactive on the other are directional selective hard negatives, used to test whether a docking score can suppress both single-target arms.

Under the strict labeling rule (dual: both ends pChEMBL ≥ 6.5; selective: active end ≥ 6.5 and opposite end ≤ 5.5), pairs that simultaneously supply enough A-only and B-only hard negatives were scarce despite the large candidate list. Only four pairs met the thick-panel gate of ≥50 strict hard negatives on **both** ends. After excluding metal-dependent HDAC1/HDAC6, which is unsuitable as a routine small-molecule docking benchmark, PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB formed three relatively well-supplied pairs. EGFR/HER2 retained only 7 strict B-only ligands and was therefore kept as a supply-limited case, not as a thick panel equivalent to the first three (Table 1).

This supply constraint is not a ChEMBL-only counting artifact. A zero-docking BindingDB / PubChem count check on the same four pairs (Supporting Information Table S12) left the thick-panel gate intact under an equal-relation rule that more closely matches pChEMBL (`equal_only`): min hard-negative counts for the three frozen thick pairs were 76 / 92 / 58 in BindingDB and 86 / 97 / 61 in PubChem (ChEMBL cache: 80 / 78 / 56), all still ≥ 50. EGFR/HER2 rose to about 30 B-end hard negatives in the other databases — enough for a thin (≥ 20) pool, not a thick (≥ 50) panel. Treating censored inequality records as point estimates (`as_is`) would inflate EGFR/HER2 supply (BindingDB min HN = 85), but 49 of 92 as-is B_only ligands have **only** `>` records on EGFR. That changes the definition of “both ends quantitatively measured with equal-relation assays” and was not used to freeze the benchmark. PubChem tracks BindingDB closely (deposition overlap) and is not a second independent census.

The size of the final benchmark was therefore constrained by the availability of experimentally defined directional hard negatives, not by post-hoc selection of docking-favorable pairs. DualFourClass-Bench is a constrained but experimentally labeled four-pair evaluation panel plus protocol, not a complete sample of dual-target tasks and not a comprehensive benchmark suite (Methods 2.1–2.3). The strict 6.5/5.5 rule quantifies supply and records panel construction, whereas θ = 6.0 defines the experimental-state labels for all primary AUROCs (Methods 2.2). The following analyses proceed from benchmark formulation (Section 3.2), through ligand-level chemical baselines (Section 3.3), to evaluation-condition sensitivity (Section 3.4) and falsification controls (Section 3.5).

### 3.2 Benchmark Formulation Changes Apparent Dual-Target Recognition

On the frozen four pairs, AutoDock Vina scores were evaluated under one unified θ = 6.0 label rule using pocket-matched directional AUROC (Figure 1B; Methods 2.8). Scores are \(S=-E_{\mathrm{Vina}}\) (higher better); dual is the positive class. The prespecified pair summary is `summary_min`, the smaller of the two directional arms, so a strong arm cannot hide a weak one. Arithmetic, geometric, and harmonic means of the two arms are sensitivity aggregators only; pair ranking is unchanged under all four aggregators (Table S26). For AChE/BChE and PIK3CA/PIK3CB, construction used the stricter 6.5/5.5 rule, but θ = 6.0 gives identical ligand classification and AUROC on this data (Table S4). EGFR/HER2 and PIK3CA/mTOR are more threshold-sensitive; under the strict rule both become underpowered on B_only, so the strict rule is a supporting sensitivity analysis, not a second primary standard. Ranking trends held across the full threshold grid (Figure S1A).

These four `summary_min` values are therefore **not** four interchangeable estimates of intrinsic docking performance. AChE/BChE and PIK3CA/PIK3CB were built under the strict supply gate; EGFR/HER2 and PIK3CA/mTOR used θ = 6.0; the panels also differ in n, series composition, and receptor. Cross-pair differences mix those construction factors with target-pair biology.

EGFR/HER2, AChE/BChE, PIK3CA/PIK3CB, and PIK3CA/mTOR gave directional `summary_min` values of 0.430, 0.606, 0.500, and 0.692, respectively (Table 2; Figure 4A; Figure S4). Different pairs are limited by different weak arms: dual-versus-B-only AUROC is 0.430 on EGFR/HER2 and 0.500 on PIK3CA/PIK3CB, whereas PIK3CA/mTOR reaches 0.714 and 0.692 on the two directions (Figure 4A). Relative to pooling, pocket matching raised point estimates without changing rank order (Table S6).

The same frozen scores were then scored under a Dual-versus-neither comparator and under Dual versus all non-duals (Table 3; Figure 3). Dual-versus-neither is a **nonselectivity-controlled comparator** on this panel (experimental inactives; `vina_mean`), not a claim that Dual versus neither is the official task of prior dual-target benchmarks. The two AUROCs use different negative sets and are a **descriptive formulation contrast**, not a paired significance test.

EGFR/HER2 provides the clearest formulation example. Dual versus neither yielded AUROC 0.756 [0.562, 0.920] (n_neg = 12), whereas directional `summary_min` remained 0.430 [0.284, 0.576]. Dual versus all non-duals collapsed to 0.551 [0.443, 0.666], showing that the extra difficulty is the selectives. In a mixed-library ranking of all 110 EGFR/HER2 ligands by `vina_mean`, the Top-10 contained 1 dual, 5 A-only, 4 B-only, and 0 neither (EF10 = 0.393; hard-negative fraction = 0.90); EF5 was also below random (Table S25). A Dual-versus-neither readout would therefore have supported docking-based dual recognition on this pair, while the directional task and the screening-facing Top-10 both show preferential enrichment of selectives.

That formulation gap is **pair-dependent**, not a four-pair overestimation law. AChE/BChE and PIK3CA/PIK3CB showed only small Dual-versus-neither increments (0.649 and 0.559) whose intervals overlap the directional arms. PIK3CA/mTOR Dual versus neither is underpowered (neither n = 4) and is not interpreted as a reverse effect; Dual versus all non-duals on that pair was 0.674, close to `summary_min` 0.692. PIK3CA/mTOR is therefore treated as a **conditional directional signal**, not as the paper’s central success case (Results 3.4).

**Table 2.** Pocket-matched directional AUROC on the frozen K = 4 set (Vina; unified θ = 6.0), with all four prespecified descriptor `summary_min` values. The highest descriptor is a best single-descriptor reference, not a confirmatory competitor. Wrong-pocket and ligand-efficiency controls are in Table S6; full descriptor arms are in Table S28.

| Pair | n (dual / A_only / B_only) | dual vs A_only (pocket B) | dual vs B_only (pocket A) | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.284, 0.576] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.440, 0.740] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.347, 0.648] | 0.622 | 0.620 | 0.595 | 0.418 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.464, 0.802] | 0.463 | 0.448 | 0.310 | 0.260 |

**Table 3.** Same Vina scores under Dual-versus-neither versus directional formulations (unified θ = 6.0). Dual-versus-neither uses experimental inactives (`vina_mean`); Dual versus all non-duals counts A-only, B-only, and neither as negatives. Directional CIs are from Table 2. The contrast is descriptive; negative sets differ. PIK3CA/mTOR Dual versus neither is underpowered (n_neg = 4).

| Pair | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|---|---:|---:|---:|---:|
| EGFR/HER2 | 0.430 [0.284, 0.576] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.440, 0.740] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/PIK3CB | 0.500 [0.347, 0.648] | 0.559 [0.373, 0.746] | 16 | 0.556 [0.437, 0.672] |
| PIK3CA/mTOR | 0.692 [0.464, 0.802] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |

Docking discrimination was therefore not a consistent cross-pair capability, and no pair clearly exceeded chance after uncertainty accounting. Under the primary frozen-receptor protocol, PIK3CA/mTOR had the highest `summary_min` point estimate (0.692), but its 95% bootstrap CI (0.464–0.802) included 0.5 and its paired difference from the best single-descriptor reference included 0 (Table S19). This point-estimate advantage was not invariant to receptor realization (Results 3.4). AChE/BChE (0.606) remained below TPSA (0.733), while EGFR/HER2 (0.430) and PIK3CA/PIK3CB (0.500) showed no clear advantage over their descriptor references.

Docking coverage is not complete. On the main panels, both-end scores were obtained for 110/110 EGFR/HER2 ligands, 95/100 AChE/BChE ligands, 99/100 PIK3CA/PIK3CB ligands, and 48/48 PIK3CA/mTOR ligands (Table S27). The single PIK3CA/PIK3CB failure is `PAB_034` (A-only; CHEMBL5089694), a docking timeout on 4L23 (`timeout_900s`, 23 rotatable bonds), not a label filter; PIK3CB 2WXF succeeded. AUROCs are therefore conditional on compounds AutoDock Vina can process. Substituting RTMScore or GNINA as an alternative scoring channel did not change the overall ranking. After unified best-of-9 pose coverage, GNINA pocket-matched `summary_min` remained at or below same-panel Vina on EGFR/HER2, AChE/BChE, and PIK3CA/mTOR. On PIK3CA/PIK3CB, GNINA best-of-9 was 0.533 versus Vina 0.500; both are near chance with overlapping intervals (Tables S14–S15; Figure S1B). GNINA is a single CNN channel, not a three-engine competition. The protocol passed cognate pose-generation QC; that QC is not screening-performance validation.

### 3.3 Ligand Properties and Chemotype Explain Much of the Apparent Signal

To test whether docking discrimination exceeds simple ligand-level signal, pocket-matched docking was first compared with four prespecified physicochemical descriptors (Figure 4B; Table 2). Relative to the **best single-descriptor reference** on each pair, the paired difference in summary_min was −0.052, −0.128, −0.122, and +0.229 for EGFR/HER2, AChE/BChE, PIK3CA/PIK3CB, and PIK3CA/mTOR; all four 95% confidence intervals include 0 (Table S19; Figure S3C). Even the largest positive point difference, on PIK3CA/mTOR, is therefore not distinguishable from the ligand-property reference with the present sample. This comparison uses pocket-matched summary_min, not the pooled `vina_mean` gate (EGFR/HER2 `vina_mean` 0.2824 ≠ Table 2 0.4297).

AChE/BChE is a direct confounding case. Mean TPSA was ≈ 75 for dual ligands versus ≈ 51 for selective hard negatives (Figure 4C). TPSA alone gave AUROC ≈ 0.769, above Vina under the same contrast (≈ 0.56). Adding heavy-atom count and TPSA raised dual-versus-B-only AUROC from 0.606 to 0.807, while the docking-score odds ratio was only ≈ 1.18 (Figure 7C). On this arm, docking discrimination therefore depends largely on ligand physicochemical properties and cannot be read as independent pocket-specific information.

PIK3CA/mTOR differs in degree. Adding heavy-atom count and TPSA shifted AUROC by about +0.07 to +0.11, with docking odds ratios ≈ 2.19 and 3.08, suggesting some residual pocket-related signal. The paired difference versus the descriptor baseline still includes 0, so that residual cannot be claimed as a confirmed independent advantage. After ligand-efficiency normalization, only PIK3CA/mTOR remained above the heavy-atom baseline (0.657 versus 0.463).

A two-dimensional chemical baseline makes the same point (Figure 7A). ECFP4 logistic regression under Bemis–Murcko scaffold GroupKFold yielded fold AUROCs of about 0.78–0.91 on several arms, well above the corresponding docking contrasts — for example 0.85 versus 0.43 for EGFR/HER2 dual-versus-B-only. That result means discrimination remains when molecules from the same Murcko scaffold are not shared between folds; it is not target-external generalization. On PIK3CA/mTOR, \(n_{\mathrm{scaffolds}} \approx n\), so the split is nearly leave-one-scaffold. A random `StratifiedKFold` check sits on average +0.011 above the scaffold split across eight directional contrasts (Table S20; Figure S3D); leakage is small. Dual/selective labels are systematically associated with chemotype, so an AUROC on docking scores alone does not establish pocket-specific physical recognition.

Under the present scaffold-grouped benchmark, adding the pocket-matched docking score to ECFP4 changed CV AUROC by at most approximately 0.02 in absolute value (largest: −0.0198 on PIK3CA/mTOR dual versus A-only), and the change was negative on several arms (Table S24). This is not a claim that docking encodes no structural information in general: the logistic architecture is simple, K = 4, and there is no nested model comparison. Logistic docking AUROC is not the rank AUROC in Table 2 and is often lower. Chemotype-constrained A-only/B-only subsets at ECFP4 Tanimoto ≥ 0.7 were empty. At T ≥ 0.3, the strongest unmatched arm (PIK3CA/PIK3CB dual versus A-only, 0.691) fell to 0.503 (n_neg = 11), whereas distant hard negatives (T < 0.3) rose to 0.819 (Table S23). T ≥ 0.3 is a similarity-constrained subset, not a chemically matched analogue set. T ≥ 0.4/0.5 cells are often n_neg ≤ 7 and are not interpreted as a second primary result. The scarcity of close analogues is itself a second data-supply bottleneck beyond the four-state label requirement.

On potency- or size-matched subsets, dual-versus-B-only remained weak or near chance on EGFR/HER2 and PIK3CA/PIK3CB (about 0.45–0.52). The PIK3CA/mTOR ranking trend was unchanged, but per-arm n was often < 15 with wide intervals (Table S5; Figure 7D). All four descriptors are shown in Figure 7B; none is treated as a confirmatory competitor.

### 3.4 Evaluation-Condition Sensitivity: Activity Aggregation, Ligand Panels, and Receptor Realization

Primary labels use the maximum available pChEMBL value. Replacing that aggregation with the median among repeated measurements, after re-fetching assay-level records for every scored ligand, changed four-state class assignment at θ = 6.0 for 7/110 EGFR/HER2 ligands (label agreement 103/110 = 93.6%), 1/95 AChE/BChE ligands (94/95 = 98.9%), 1/99 PIK3CA/PIK3CB ligands (98/99 = 99.0%), and 0/48 PIK3CA/mTOR ligands (48/48 = 100%) (Table S29). Numeric max ≠ median was more common than class flips (40/110, 13/95, 25/99, 27/48). On API-refetched labels, `summary_min` moved from 0.417 to 0.424 (EGFR/HER2), 0.606 to 0.629 (AChE/BChE), 0.500 to 0.500 (PIK3CA/PIK3CB), and 0.692 to 0.692 (PIK3CA/mTOR). Frozen Table 2 EGFR/HER2 is 0.430 rather than 0.417 because one cache/API mismatch (`EH120_060` / CHEMBL24828) reclassifies that ligand as dual under API max; relative to the frozen table, median aggregation still leaves EGFR/HER2 at 0.424. Pair ranking and the directional conclusions are therefore insensitive to this aggregation choice. Assay-level heterogeneity remains, because pChEMBL values are not assay-equivalent.

Ligand-panel and protocol-level sensitivity analyses next tested whether the higher PIK3CA/mTOR summary_min was an artifact of one panel composition or one search setting (Figure S5). Lowering exhaustiveness from 16 to 8 moved summary_min from 0.692 to 0.660 (Δ ≈ 0.03), much smaller than between-pair differences (Figure S1D).

On PM110, which retains every PM48 ligand and expands to actual n = 115 (30 / 30 / 30 dual / A_only / B_only for analysis), Vina summary_min was 0.648 [0.51, 0.76], about 0.04 below PM48 (0.692), with the same ranking trend (Figure S1C). The directional signal is therefore not driven solely by particular PM48 members. PM110 is not an independent validation set; it is a stability check. Same-panel RTMScore was 0.576; GNINA best-of-9 was 0.613 [0.46, 0.74] on PM110 and 0.655 [0.43, 0.81] on PM48, still not above same-panel Vina.

On the unused-pool holdout — ligands that did not enter main-panel construction or protocol tuning (20 / 20 / 20 per pair; seed 20260731; EGFR/HER2 not eligible) — PIK3CA/mTOR summary_min was 0.765 [0.603, 0.891], above the main-panel 0.692. AChE/BChE was 0.618 [0.422, 0.759], close to the main panel but with a confidence interval that spans 0.5. PIK3CA/PIK3CB fell to 0.425 [0.241, 0.618] (Tables S8 / S16). Of 60 PIK3CA/PIK3CB holdout ligands attempted, 59 yielded both-end scores; HOAP_028 failed on both pockets because AutoDock atom type `B` is unsupported (boron; Table S27). AChE and PIK3CA/mTOR holdouts were 60/60 successful. The boron failure is a chemical-coverage limit of the engine, not silent missingness; AUROC is conditional on processable compounds. The holdout shares the same ChEMBL extraction batch and is not independent cross-database validation. It supports persistence of the observed signal in an unused ligand pool.

The PIK3CA/mTOR directional signal is therefore still visible among same-ecosystem ligands unseen at panel construction, whereas the PIK3CA/PIK3CB signal is not. Docking performance remains a target-pair property rather than a transferable attribute across pairs.

We next tested whether directional discrimination depends on a particular receptor realization, holding one pocket frozen and replacing the other (Figure 5; Table S9; Table S30). Three alternate crystals passed cognate redocking QC, with best-of-9 RMSD 0.607 Å (4JPS), 0.624 Å (5DXT), and 0.515 Å (4JSX); chimeric 3T8M remains excluded.

On PIK3CA/mTOR, replacing PIK3CA 4L23 with 4JPS or 5DXT, while holding mTOR at 4JT6, dropped PM48 summary_min from 0.692 to 0.486 [0.259, 0.692] and 0.505 [0.292, 0.696] (Figure 5A). The change concentrated in the D/B direction that depends on the alternate PIK3CA structure; the D/A direction, which still uses original 4JT6, stayed at 0.714. Replacing mTOR 4JT6 with 4JSX gave summary_min 0.639 [0.418, 0.776]. After the mTOR swap the point estimate remains above 0.5, but the 95% CI includes 0.5.

The same PIK3CA crystals were then used on the PIK3CA/PIK3CB panel, holding 2WXF scores frozen (exhaustiveness 8, matching the main panel; Figure 5B). Replacement **raised** summary_min from 0.500 to 0.691 [0.516, 0.779] (4JPS) and 0.685 [0.506, 0.768] (5DXT). Dual versus A-only, which still uses frozen 2WXF, stayed at 0.691. Dual versus B-only, which uses the alternate PIK3CA score, rose from 0.500 to 0.707 (4JPS) and 0.685 (5DXT). The weak arm therefore switches: originally D/B on 4L23 (0.500); after 4JPS the bottleneck is the frozen 2WXF arm (0.691); after 5DXT the two arms are nearly balanced. Both alternate jobs attempted 100 ligands and succeeded on 99; the missing ligand is again `PAB_034` (timeout at 600 s on both 4JPS and 5DXT). The same ligand already timed out on original 4L23, so the 99-ligand set is identical to Table 2. A 100-ligand AUROC is not available under any PIK3CA crystal in this protocol; the failure is a docking timeout, not an experimental label filter.

The two pairs therefore move in opposite directions under the same PIK3CA perturbation. Receptor choice can both create and attenuate apparent dual-target discrimination. This is a receptor-realization effect, not a robustness demonstration and not a unidirectional collapse. The design holds one receptor fixed, so the Δ is attributable to the replaced pocket rather than a simultaneous two-structure swap. The two pairs share PIK3CA; the pattern is not a universal law over K = 4.

Cα comparison shows a local pocket Cα RMSD of only 0.343 Å between 5DXT and 4L23, yet PIK3CA/mTOR summary_min still fell to 0.505, so backbone similarity is not sufficient to preserve discrimination (Table S10). Global Cα RMSD among these deposited PIK3CA structures (1.44–1.49 Å) exceeds that among these mTOR structures (0.45 Å), consistent in direction with greater PIK3CA-end movement on PIK3CA/mTOR but not a quantitative causal explanation of the opposite PIK3CA/PIK3CB shift: 5DXT matched 862 Cα atoms versus 982 for 4JPS, and each target has only one or two alternates. Cognate-ligand centroid distances of 2.1–2.6 Å indicate the same general ATP-competitive site. Passing pose-generation QC is not the same as transferable screening discrimination. The protocol passed cognate pose-generation QC; it was not “validated” as a virtual-screening method.

### 3.5 Wrong-Pocket Controls Show That a Benchmark Control Can Fail under Panel Shift

On the main panels, pocket-matched summary_min exceeded the wrong-pocket control on all four pairs; matched-minus-wrong differences were 0.170, 0.161, 0.151, and 0.090. The EGFR/HER2 and AChE/BChE intervals exclude 0; the PIK3CA/PIK3CB and PIK3CA/mTOR intervals include 0 (Tables S6, S17; Figure 6A; Figure S3A). Wrong-pocket summary_min values were 0.260, 0.444, 0.349, and 0.602. That main-panel pattern is **not** taken as evidence of pocket-specific signal.

That relationship reversed on the unused-pool holdout (Figure 6B). Wrong-pocket summary_min was 0.788, 0.643, and 0.520 for PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB, versus matched-pocket 0.765, 0.618, and 0.425. All three matched-minus-wrong point differences were negative (−0.023 / −0.025 / −0.095), and every 95% confidence interval includes 0 (Table S17; Figure S3B). The reversal is a point-estimate pattern, not an interval that excludes zero. Wrong-pocket is therefore **not a reliable universal negative control under panel shift**.

Potency- and size-matched comparisons tested whether ligand potency or size differences in the holdout could explain the reversal (Figure 6C; Table S13). Wrong-pocket remained ≥ matched-pocket (after potency matching: AChE/BChE 0.642 versus 0.593, n_min = 18; PIK3CA/PIK3CB 0.562 versus 0.363, n_min = 11; PIK3CA/mTOR 0.734 versus 0.715, n_min = 12). Sampling shift is real — most clearly on PIK3CA/mTOR, where holdout dual / A_only mean pA is about 1.1–1.3 lower than the main panel and B_only mean pB is about 1.8 lower — but matching does not remove the paradox.

Scoring-independent contact_count reached AUROC 0.698–0.714 on the B direction, indicating a real size/burial contribution on that arm, but the magnitude does not explain the full Vina wrong-pocket signal (Figure 6D; Table S11). On PIK3CA/mTOR, Vina wrong-pocket summary_min was 0.788, whereas the weaker contact_count arm was only 0.552. On the A arm, dual versus A_only size gaps are small and contact_count AUROC is closer to chance (0.552–0.622).

The holdout wrong-pocket reversal is therefore an unresolved failure mode exposed by the benchmark, not a phenomenon that a single size or potency factor currently explains.

### 3.6 Structural Context Provides Only Exploratory Clues

As an exploratory analysis, we compared within-pair whole-chain sequence identity with summary_min (Table S7). Among the four pairs, PIK3CA/mTOR has the lowest whole-chain identity (18.1% over alignment length) and the highest summary_min, whereas EGFR/HER2 is the reverse (71.4%; ErbB-family kinase domains are highly homologous). That pattern is inconsistent with a simple “more similar targets are harder to distinguish” rule. With only four pairs, and with whole-chain identity not a direct measure of binding-pocket similarity, the observation is a structural background clue, not correlation evidence. PIK3CA and mTOR belong to the PIKK-related superfamily and share known local homology at the ATP-competitive site; low whole-chain identity must not be read as dissimilar pockets.

Existing pose-level diagnostics on PIK3CA/mTOR show two representative failure typologies (not a panel-wide PLIF). T2: a selective hard negative forms geometrically clean, hinge-positive ATP-like poses in both pockets (for example amino-triazine / morpholine–ATP chemotypes that remain high-occupancy and hinge-positive on the weak mTOR end), inflating both scores. T5: some classic duals (e.g., Torin1, omipalisib) rank well under Vina, but alternative rescoring prefers poses off the PIK3CA hinge / cognate site. Cognate ligands PI-103 / X6K recover near-native poses under protocol checks (Table S3). These are observed pose patterns, not a residue-level mechanism, and pose-generation QC is not screening validation.

## 4. Discussion

### 4.1 Benchmark formulation changes the evidentiary standard for dual-target docking

The primary finding is not that one docking scoring function attained the highest dual-target AUROC. It is that benchmark formulation can change what “dual-target docking success” appears to mean. Standard docking benchmarks typically separate actives from decoys. Here a model must distinguish dual-active ligands from single-target selective hard negatives in both directions. Those hard negatives are experimentally active on one target and therefore cannot be treated as ordinary decoys. Dual-versus-neither on the same scores is a **nonselectivity-controlled comparator**, not “the conventional dual-target benchmark.” The supply audit showed that public bioactivity data rarely provide enough such ligands on both arms at once; only a few of 49 candidate pairs met the strict thick-panel gate (Results 3.1; Figure 2).

Zhou, Li, and Hou already showed that docking-based dual-kinase screening can look useful against noninhibitors, is structure-dependent, and still admits a high false-positive rate among predicted duals.[9] DualFourClass-Bench asks a narrower follow-up on the same scores: whether a Dual-versus-neither (inactive) readout and a directional Dual-versus-selective readout agree. They do not agree on EGFR/HER2 (Results 3.2, Table 3, Figure 3): Dual versus neither was 0.756, directional `summary_min` was 0.430, and mixed-library Top-10 ranking enriched selectives (9/10). **EGFR/HER2 provides a clear example, not a four-pair law.** AChE/BChE and PIK3CA/PIK3CB showed only small, overlapping increments; PIK3CA/mTOR Dual versus neither is underpowered (neither n = 4). The contrast is descriptive rather than a paired significance test. The increment relative to 2013 is that formulation gap—not another four-pair docking survey.

That data constraint is itself methodological. DUD, DUD-E, and LIT-PCBA already showed that decoy construction, chemical bias, and assay-derived labels change virtual-screening conclusions.[5–7] Simple methods and some unbiasing procedures can also overestimate structure-based virtual screening by learning ligand distributions.[12] Recent bioassay-derived evaluation further emphasizes that real assay data can expose limits that constructed ligand/decoy sets conceal.[13] DualFourClass-Bench does not use those single-target collections and does not evaluate DiffDock-Pocket. It extends the same concern to a dual-target setting: the evaluation depends on how hard negatives are defined experimentally, not on how long the candidate-pair list is.

The main value of DualFourClass-Bench is therefore not dataset size. It is the conversion of dual-target recognition into an experimentally labeled hard-negative discrimination task that must hold in both directions. The resource is a curated four-pair panel plus evaluation protocol, not a comprehensive dual-target suite.

A remaining out-of-panel failure mode belongs with this evidentiary standard. On the main panels, pocket-matched scores exceeded the wrong-pocket control; on the unused-pool holdout the inequality reversed as a point-estimate pattern whose paired intervals still include zero (Results 3.5; Figure 6). Wrong-pocket is therefore not a reliable universal negative control under panel shift, and main-panel matched > wrong is not taken as pocket-specific proof.

### 4.2 Chemical information can substitute for apparent docking signal

The difficulty of dual-target docking is not the sum of two single-target docking problems. A favorable score for an A-active ligand in pocket A does not imply activity at target B. Even under the task-aligned directional metric, no pair yielded discrimination clearly above chance after uncertainty accounting. Under the primary frozen-receptor protocol, PIK3CA/mTOR had the highest directional point estimate (0.692), but its 95% CI included 0.5 and its paired difference from the best single-descriptor reference included 0 (Results 3.2; Table S19). It is therefore a **conditional directional signal**, not a generalizable success case.

On AChE/BChE, TPSA alone discriminated better than docking, and scaffold-grouped ECFP4 exceeded docking on several arms (Results 3.3). The strong ligand-only baseline demonstrates that the experimental labels are associated with chemical-space differences that can be exploited without receptor information; docking performance must therefore be interpreted relative to ligand-only baselines rather than in isolation. Under the present scaffold-grouped benchmark, adding the docking score to ECFP4 changed CV AUROC by at most approximately 0.02 in absolute value, and the change was negative on several arms. Dual-target discrimination was strongly target-pair dependent, and docking provided limited incremental information beyond ligand-level chemical baselines under scaffold-aware evaluation. That is a statement about this panel, not a proof that docking encodes no pocket-specific information.

Chemotype-constrained hard negatives make the same point at a coarser similarity cutoff. T ≥ 0.7 matched A-only/B-only subsets were empty. At T ≥ 0.3, the strongest unmatched arm (PIK3CA/PIK3CB dual versus A-only, 0.691) fell to 0.503 (n_neg = 11), whereas distant hard negatives (T < 0.3) rose to 0.819. T ≥ 0.3 is a similarity-constrained subset, not a chemically matched analogue set.

That observation is consistent with recent attention to chemical bias in virtual-screening benchmarks: simple models or poorly constructed decoys can look strong by learning ligand composition rather than target-specific recognition.[7,12] Without A-only/B-only hard negatives and ligand-property / chemical baselines, an apparently strong dual-target docking result may only be recovering molecular properties associated with the dual label.

### 4.3 Receptor realization is an independent source of uncertainty

PIK3CA/mTOR remains the case that most warrants further study and the case that most requires caution. Main-panel summary_min was 0.692, PM110 was 0.648, and the unused-pool holdout was 0.765, so the directional signal is not driven solely by a few PM48 members. It is not receptor-invariant: replacing the PIK3CA structure dropped summary_min to 0.486 and 0.505, whereas replacing the mTOR structure left 0.639 (Results 3.4; Figure 5A). The accurate claim is not that PIK3CA/mTOR docking reliably identifies dual-target ligands. It is that this pair shows a limited directional signal under a particular receptor realization — persistent under ligand-panel replacement, not assumed invariant under receptor replacement.

The same PIK3CA replacement on PIK3CA/PIK3CB moved the estimate in the opposite direction, from 0.500 to 0.691 and 0.685, while the B-end receptor was held frozen (Figure 5B). Receptor choice can therefore alter not only the magnitude but the direction of apparent dual-target discrimination. The contrasting effects argue against interpreting a receptor swap simply as a loss of docking accuracy: receptor realization is part of the evaluation condition and an independent source of variance that can create or attenuate apparent discrimination. The two-pair, one-pocket-at-a-time design is stronger evidence than a single collapse anecdote, but it is not a universal law (K = 4; both pairs share PIK3CA). No molecular mechanism is claimed.

A coupled-task reading is useful and remains a hypothesis. `summary_min` tracks the weaker arm, so replacing \(S_A\) can change which arm is the bottleneck. On PIK3CA/PIK3CB the original weak arm was D/B on 4L23 (0.500); after 4JPS that arm rose to 0.707 and the frozen 2WXF arm (0.691) became limiting. Other non-exclusive hypotheses include local side-chain or pocket-geometry changes that reorder dual versus selective scores, and different ligand chemical distributions on the two panels facing the same PIK3CA crystals. Without a panel-wide residue-level PLIF analysis, these alternatives remain unresolved.

Local pocket Cα RMSD between 5DXT and 4L23 is only 0.343 Å, yet PIK3CA/mTOR summary_min still fell to 0.505. Structural similarity and transferable screening discrimination are therefore not the same question. Passing pose-generation QC is also not screening-performance invariance. That distinction is consistent with recent kinase cross-docking benchmarks that treat receptor representation as an independent performance variable; those studies used different docking engines and are not an extrapolation of the present protocol.[14]

### 4.4 Implications for dual-target virtual screening and generative design

The results bear directly on dual-target virtual screening and generative design. Favorable docking scores in both pockets do not automatically mean that a generated molecule is an experimentally plausible dual-active ligand. If dual-target generators use docking as a downstream filter, they should not report two-pocket scores alone; they should also be evaluated against selective hard negatives and ligand-only chemical controls. Even after single-target ultralarge docking, postprocessing and rescoring have been shown not to separate known binders from inactives robustly across assays.[15] A dual-target setting additionally requires suppressing both experimental hard-negative arms, so two favorable scores — or their simple average — are not sufficient evidence.

This study does not show that existing dual-target generators fail, and it does not evaluate DualDiff, FuseDiff, or other generative models.[10,11] DualDiff’s Dual High Affinity is dual success versus reference-ligand dock scores, not mean pooling; FuseDiff’s independent test set is the DualDiff benchmark (DDF). Those papers ask whether generated structures can obtain favorable docking scores. DualFourClass-Bench can serve as a downstream evaluation layer for such methods: whether generated molecules outrank experimentally defined single-target hard negatives, rather than merely optimizing a docking score.

### 4.5 Limitations

Five limitations define the scope of interpretation.

First, the benchmark contains only four target pairs because experimentally defined dual-target hard negatives are scarce. K = 4 is a data-constrained case panel, not a comprehensive suite. The four `summary_min` values also mix panel-construction differences (strict 6.5/5.5 versus θ = 6.0; unequal n) with target-pair biology and should not be read as a population-level ranking of intrinsic docking performance.

Second, ground truth is ChEMBL-derived. The unused-pool holdout remains within the same extraction batch and is not independent cross-database validation. BindingDB/PubChem checks were count-level only.

Third, assay heterogeneity remains after activity-aggregation control. Primary curation uses maximum pChEMBL; replacing that rule with the median among repeated measurements produced minimal pair-level changes (Results 3.4; Table S29), so max aggregation is a controlled limitation rather than an unresolved fatal threat. pChEMBL measurements are still not fully assay-equivalent. Confidence≥8 and Homo sapiens filters were not rebuilt.

Fourth, receptor realization can raise or lower pair-level discrimination, but the experiments do not identify a molecular origin. Pocket-local Cα RMSD alone could not explain the observed performance change, and residue-level PLIF/side-chain analyses were not systematically performed. The one PIK3CA/PIK3CB docking timeout (`PAB_034`) is reported as 100 attempted / 99 successful / 1 failed on the original and both alternate PIK3CA crystals; it was not excluded because of its label.

Fifth, this study evaluates computational discrimination rather than experimentally testing newly predicted dual-target compounds. The benchmark addresses the reliability of docking-based ranking, not the prospective biological efficacy of selected molecules. The study does not aim to prove that docking is universally good or bad for dual-target discovery; it asks whether the evaluation formulation itself changes the apparent evidence for dual-target recognition.

## 5. Conclusions

We established DualFourClass-Bench as an experimentally grounded evaluation setting for docking-based dual-target recognition, explicitly testing whether docking can distinguish dual-active ligands from A-selective and B-selective hard negatives. Across four frozen target pairs, dual-target discrimination was strongly target-pair dependent (`summary_min` AUROCs ranging from 0.430 to 0.692), and docking provided limited incremental information beyond ligand-level chemical baselines under scaffold-aware evaluation. Primary pair-level estimates were largely insensitive to replacing maximum pChEMBL with the median among repeated measurements. PIK3CA/mTOR exhibited the strongest point estimate and retained a positive directional signal in an unused ligand-pool holdout; however, the uncertainty of the primary estimate and its sensitivity to receptor structure preclude interpreting this result as a generalizable dual-target decision rule.

The broader analyses indicate that the apparent performance of dual-target docking is jointly determined by task formulation, ligand chemical composition, and receptor realization. On EGFR/HER2, a Dual-versus-neither comparator (AUROC 0.756) would have supported docking-based dual recognition, whereas the directional worst arm remained 0.430; that formulation gap is pair-dependent, not a four-pair reversal, and is a descriptive contrast rather than a paired significance test. In several target pairs, ligand-property or chemotype-based references matched or exceeded docking discrimination, while the unused-pool holdout revealed an unresolved wrong-pocket reversal that was not eliminated by potency or size matching, although the corresponding paired confidence intervals included zero. Receptor realization can alter the magnitude and even the direction of apparent discrimination; it is a realization effect, not a robustness certificate. Collectively, these findings argue against interpreting a favorable score in two pockets as sufficient evidence of dual-target activity. Instead, dual-target virtual screening should incorporate experimentally defined single-target hard negatives, ligand-level confounder controls, out-of-panel ligand evaluation, and receptor-sensitivity analysis. The principal contribution of DualFourClass-Bench is therefore not a universal docking winner, but a systematic protocol for defining the evidentiary and reliability boundaries of docking-based dual-target recognition.

## Data and Software Availability

Benchmark membership, experimental-state labels, receptor and docking-box definitions, per-ligand docking scores, analysis tables, and all scripts used to regenerate the reported statistics and figures are available in the `Dual_Target_Docking` directory of the public repository at https://github.com/1280602962-debug/gwj260531. `data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` indexes the principal numerical results and their source tables. The analysis environment and zero-docking reproduction commands are documented in the repository README. A versioned Zenodo archive, including the submission data package and its DOI, will be deposited before journal submission; the DOI is therefore not claimed in this draft.

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

(11) Wu, J.; Qiao, A.; Wang, Z.; Wei, Z.; Chen, S. FuseDiff: Symmetry-Preserving Joint Diffusion for Dual-Target Structure-Based Drug Design. arXiv:2603.05567, 2026. (preprint)

(12) Tran-Nguyen, V.-K.; Ballester, P. J. Beware of Simple Methods for Structure-Based Virtual Screening: The Critical Importance of Broader Comparisons. *J. Chem. Inf. Model.* **2023**, *63*, 1401–1405. DOI: 10.1021/acs.jcim.3c00218.

(13) Ahmed, F.; Soellner, M. B.; Brooks, C. L., III. Real-World Assessment of Machine-Learned Docking Using Bioassay-Derived Benchmarks. *J. Chem. Inf. Model.* **2026**, *66*, 8752–8759. DOI: 10.1021/acs.jcim.5c03020.

(14) Schaller, D. A.; Christ, C. D.; Chodera, J. D.; Volkamer, A. Benchmarking Cross-Docking Strategies in Kinase Drug Discovery. *J. Chem. Inf. Model.* **2024**, *64*, 8848–8858. DOI: 10.1021/acs.jcim.4c00905.

(15) Sindt, F.; Bret, G.; Rognan, D. On the Difficulty to Rescore Hits from Ultralarge Docking Screens. *J. Chem. Inf. Model.* **2025**, *65*, 5553–5566. DOI: 10.1021/acs.jcim.5c00730.
