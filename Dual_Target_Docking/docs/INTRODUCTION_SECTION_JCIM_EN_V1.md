# Introduction (JCIM Articles draft, English)

## 1. Introduction

### 1.1 Dual-target design and why conventional docking benchmarks are insufficient

Multitarget drug design aims to modulate two or more biological targets with a single small molecule, in order to address pathway redundancy, compensatory signaling, and drug resistance in complex disease. Relative to a classical single-target agent, a rationally designed multitarget ligand may act on connected nodes of a disease network and thereby achieve a more adequate pharmacological effect; this idea is now a central theme of polypharmacology.[1] Over the past decade the field has moved from largely serendipitous multi-pharmacology toward structure-guided design that combines structural biology, computational chemistry, and, increasingly, generative models.[2]

Molecular docking remains one of the most widely used tools in structure-based virtual screening (SBVS). A docking engine places the ligand in a protein binding site and a scoring function ranks ligand–receptor complementarity, providing a rapid, structure-level filter of large compound libraries.[3,4] A natural computational tactic in dual-target discovery is therefore to dock each candidate into both pockets and to treat favorable scores on both targets as evidence of dual-target potential.

How such scores should be interpreted depends on benchmark construction: negative definitions, chemical bias, and metric choice all change apparent performance. DUD and DUD-E use property-matched decoys because unmatched negatives can reduce screening to separation by coarse ligand properties.[5,6] LIT-PCBA instead uses experimental assay labels and systematically controls known decoy and chemical biases to make virtual-screening evaluation more realistic.[7] CASF-2016 evaluates scoring, ranking, docking, and screening power on protein–ligand complexes, but still poses a single-complex problem.[8] None of these resources formulates dual-target discrimination over an experimentally labeled four-state ligand space.

**Extending the single-target docking evaluation logic to a dual-target task is therefore not sufficient.**

For dual-target ligands the task changes. A strict dual-target benchmark must distinguish at least four ligand states with distinct biological meaning (a four-state data set, not a four-class classifier): **dual-active** ligands that act on both targets, **A-selective** ligands that act only on target A, **B-selective** ligands that act only on target B, and **neither** ligands that lack sufficient activity on either target (Figure 1A):

|  | *B*<sup>+</sup> | *B*<sup>−</sup> |
|--|:--:|:--:|
| *A*<sup>+</sup> | Dual | A-only |
| *A*<sup>−</sup> | B-only | Neither |

A-only and B-only ligands are not ordinary decoys. They are the **selectivity hard negatives** of the task: already potent on one target, they can produce plausible scores in that pocket while lacking corresponding activity on the other. The computational endpoint therefore asks whether docking discriminates dual-actives from the matching single-target selectives in both directions; it does not by itself establish independent pocket-specific recognition or biological dual activity.

Zhou, Li, and Hou benchmarked dual-kinase docking against noninhibitors and reported structure dependence and false positives among predicted duals.[9] Building on that evaluation setting, we introduce experimentally defined A-only and B-only directional hard negatives and directly compare apparent discrimination under different formulations. Dual-versus-neither is used as a **nonselectivity-controlled comparator**, not as “the conventional dual-target benchmark.”

Pooling two pocket scores can hide a weak direction, and beating reference-ligand scores on both targets defines computational success rather than experimental dual activity. The benchmark and readout must therefore match the four-state biological space.

Constructing such a panel from public data requires the same compound to have comparable measurements on both targets and enough selective hard negatives on both arms. Assay types, conditions, and coverage differ, so the number of target pairs that support a balanced four-state evaluation is itself a data-supply question rather than an assumed resource.

Dual-target panels also inherit chemical confounding. If dual-active and selective ligands differ in molecular properties or scaffolds, AUROC can reflect label-associated ligand distributions rather than pocket complementarity. A dual-target benchmark therefore requires explicit physicochemical and ligand-only controls.[7]

### 1.2 Docking-based dual-target design makes stricter evaluation practical

DualDiff and FuseDiff illustrate why this distinction is operational: both use docking-based criteria, including the worse pocket score and the fraction of generated molecules outperforming reference ligands at both targets, to evaluate dual-target designs.[10,11] Those metrics assess computational success relative to references, not discrimination against experimentally selective ligands. The reported generative metrics and the present hard-negative endpoint are therefore complementary rather than competing benchmarks. We do not redock their generated molecules or reinterpret their reported metrics as competitors to the present primary endpoint.

### 1.3 Objective and contribution

**Here, we ask whether benchmark formulation changes the apparent evidence for dual-target recognition.** We construct an experimentally defined four-state panel, use pocket-matched directional discrimination against A-selective and B-selective hard negatives as the primary task, and compare it with a nonselectivity-controlled Dual-versus-neither comparator. We then stress-test the observed discrimination with ligand-only, physicochemical, wrong-pocket, ligand-pool, activity-aggregation, and receptor-realization controls.

DualFourClass-Bench is a curated benchmark with two directional primary tasks: dual versus A-only scored in pocket B, and dual versus B-only scored in pocket A (Figure 1B). Neither describes the full experimental state space but does not enter the primary AUROCs. `summary_min` is the conservative worst-direction discrimination summary, not a new docking score. A public-data audit determines how many candidate pairs can support this construction; evaluation-set size is an audit result rather than a design target.

The contribution is an evaluation protocol and resource, not a new docking algorithm. It connects task formulation, confounder-aware evaluation, and evaluation-condition sensitivity to ask whether apparent discrimination behaves as a fixed property of docking or a conditional result of the benchmark setting. The protocol is intended as a downstream check for dual-target virtual screening and generative design, not as a competition with existing generators or a comprehensive dual-target suite.

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

(11) Wu, J.; Qiao, A.; Wang, Z.; Wei, Z.; Chen, S. FuseDiff: Symmetry-Preserving Joint Diffusion for Dual-Target Structure-Based Drug Design. In *Proceedings of the 32nd ACM SIGKDD Conference on Knowledge Discovery and Data Mining, Vol. 2*; ACM: New York, 2026; pp 12432–12443. DOI: 10.1145/3770855.3819050.
