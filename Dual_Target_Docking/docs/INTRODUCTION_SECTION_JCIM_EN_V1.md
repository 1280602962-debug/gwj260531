# Introduction (JCIM Articles draft, English)

> Companion to [`INTRODUCTION_DRAFT_ZH_JCIM_V1.md`](INTRODUCTION_DRAFT_ZH_JCIM_V1.md) (Chinese working draft).  
> Numbered citations and verification bounds: [`INTRODUCTION_REFS_JCIM_V1.md`](INTRODUCTION_REFS_JCIM_V1.md).  
> Positioning: [`POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`](POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md).  
> **Do not freeze *K* = 4 in this section**; the size of the audited evaluation set belongs in Methods / Results. Numbered subsections map onto Methods and Figure 1 and may be dropped in the typeset manuscript.

---

## 1. Dual-target design and the role of structure-based virtual screening

Multitarget drug design aims to modulate two or more biological targets with a single small molecule, in order to address pathway redundancy, compensatory signaling, and drug resistance in complex disease. Relative to a classical single-target agent, a rationally designed multitarget ligand may act on connected nodes of a disease network and thereby achieve a more adequate pharmacological effect; this idea is now a central theme of polypharmacology.[1] Over the past decade the field has moved from largely serendipitous multi-pharmacology toward structure-guided design that combines structural biology, computational chemistry, and, increasingly, generative models.[2]

Molecular docking remains one of the most widely used tools in structure-based virtual screening (SBVS). A docking engine places the ligand in a protein binding site and a scoring function ranks ligand–receptor complementarity, providing a rapid, structure-level filter of large compound libraries.[3,4] A natural computational tactic in dual-target discovery is therefore to dock each candidate into both pockets and to treat favorable scores on both targets as evidence of dual-target potential.

How such scores should be interpreted, however, is already known to depend on how the data set is built: the definition of negatives, hidden chemical bias, and the choice of metric all change apparent performance. The Directory of Useful Decoys (DUD) and its successor DUD-E frame single-target docking as enrichment of actives against property-matched decoys, precisely because unmatched decoys can make enrichment look like separation of crude ligand properties.[5,6] LIT-PCBA showed that artificially constructed active/decoy collections, including DUD, DUD-E, and MUV, can contain obvious and hidden chemical biases that overestimate the true accuracy of virtual screening, and replaced them with experimental dose–response labels under physicochemical-range control.[7] Structure-scoring benchmarks such as CASF-2016 ask a different, still single-complex question: scoring, ranking, docking, and screening power of scoring functions on high-quality protein–ligand complexes drawn from PDBbind-quality data.[8] None of these resources formulates dual-target recognition over an experimentally labeled four-class state space.

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

A public-data supply audit first asks how many candidate pairs can support that four-state construction (Methods 2.1–2.3). The size of the evaluation set is a result of that audit, not a design target frozen in this Introduction. Pooled docking scores, a wrong-pocket control, and two-dimensional chemical and physicochemical baselines are reported as auxiliary contrasts, to separate pocket-specific signal from ligand-level confounding.

The nested experimental question remains: to what extent can existing docking scores distinguish dual-active ligands from target-selective hard negatives on both arms, and how far does that discrimination depend on the target pair, the receptor structure, or ligand chemistry? The protocol is meant as a stricter downstream check for dual-target virtual screening and for generative dual-target design—not as a bake-off of those generators, and not as a comprehensive dual-target suite.

---

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
