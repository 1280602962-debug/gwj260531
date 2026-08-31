# Introduction (JCIM Articles draft, English)

## 1. Introduction

Multitarget drug design aims to modulate two or more biological targets with a single small molecule, in order to address pathway redundancy, compensatory signaling, and drug resistance in complex disease. Relative to a classical single-target agent, a rationally designed multitarget ligand may act on connected nodes of a disease network and thereby achieve a more adequate pharmacological effect; this idea is now a central theme of polypharmacology.[1] Over the past decade the field has moved from largely serendipitous multi-pharmacology toward structure-guided design that combines structural biology, computational chemistry, and, increasingly, generative models.[2] Molecular docking remains one of the most widely used tools in structure-based virtual screening (SBVS): a docking engine places the ligand in a protein binding site and a scoring function ranks ligand–receptor complementarity.[3,4] A natural computational tactic in dual-target discovery is therefore to dock each candidate into both pockets and to treat favorable scores on both targets as evidence of dual-target potential. How such scores should be interpreted depends on benchmark construction. DUD and DUD-E use property-matched decoys because unmatched negatives can reduce screening to separation by coarse ligand properties.[5,6] LIT-PCBA instead uses experimental assay labels and systematically controls known decoy and chemical biases.[7] CASF-2016 evaluates scoring, ranking, docking, and screening power on protein–ligand complexes, but still poses a single-complex problem.[8] None of these resources formulates dual-target discrimination over an experimentally labeled four-state ligand space.

A strict dual-target evaluation must distinguish four ligand states: **dual-active**, **A-selective**, **B-selective**, and **neither** (Figure 1A). A-only and B-only ligands are the **selectivity hard negatives** of the task: already potent on one target, they can produce plausible scores in that pocket while lacking corresponding activity on the other. The computational endpoint therefore asks whether docking discriminates dual-actives from the matching single-target selectives in both directions. Zhou, Li, and Hou benchmarked dual-kinase docking against noninhibitors.[9] Building on that setting, we introduce experimentally defined directional hard negatives and compare apparent discrimination under different formulations on the same scores. Dual versus neither, scored on experimental inactives, is reported as a formulation contrast. A balanced four-state panel is limited by the need for paired measurements and enough selectives on both arms.

Recent work has made the dual-target setting experimentally concrete without answering the same evaluation question. Wu et al. showed that large-library docking can prospectively yield joint binders for selected target pairs, while also reporting that subsequent optimization remained challenging.[19] That study asks whether docking can find dual actives; the present study asks whether retrospective evidence for dual-target recognition changes when experimentally selective ligands, rather than nonbinders, define the negative class. POLYGON generated dual-target chemical space from binding data for more than 100,000 compounds and synthesized 32 MEK1/mTOR compounds with experimental functional validation;[20] generative polypharmacology should therefore not be summarized as docking-derived success metrics alone, even though some dual-target generators still report docking-relative design metrics.[10,11] A 147-target AI docking benchmark further showed that method ranking depends on whether negatives are experimentally determined low-activity TrueDecoys or random commercial-library decoys.[21] Kinase-Bench assembled 6875 selective ligands, 75 kinases, and 422,799 decoys to test selective enrichment against kinase-specific decoys and interaction filters,[22] whereas DualFourClass constructs experimentally measured dual, A-only, B-only, and neither states on each pair.

Here, we ask whether benchmark formulation changes the apparent evidence for dual-target recognition. We construct DualFourClass-Bench as a four-state panel with two directional primary tasks—dual versus A-only scored in pocket B, and dual versus B-only scored in pocket A (Figure 1B)—and summarize the weaker of those two arms as worst-arm AUROC (`summary_min`). We then test whether the observed discrimination is retained under alternative ligand panels, activity-aggregation rules, and receptor structures.

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

(19) Wu, Y.; et al. Large Library Docking for Polypharmacology. *J. Med. Chem.* **2026**, *69*, 6210–6229.

(20) Munson, B. P.; et al. POLYGON. *Nat. Commun.* **2024**, *15*, 3636.

(21) Wei, et al. Kinase-Bench. *J. Chem. Inf. Model.* **2024**, *64*, 9528–9550.

(22) Gu, Shen, Zhang, et al. AI docking benchmark. *Nat. Mach. Intell.* **2025**, *7*, 509–520.
