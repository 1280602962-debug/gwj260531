# Introduction (JCIM Articles draft, English)

## 1. Introduction

Multitarget drug design uses one molecule to modulate two or more disease-related targets. Relative to combination therapy, a single multitarget ligand can reduce pharmacokinetic mismatch and formulation complexity.[1,2] Generative design and ultra-large library docking have already produced experimentally confirmed dual-target candidates. POLYGON tested 32 MEK1/mTOR compounds.[20] Ultra-large library docking has also found multitarget ligands among hundreds of millions of molecules, and experimental structures support some predicted binding modes.[19] Such work shows that computational methods can find dual-active molecules in a candidate library. Structure-based dual-target generative models can likewise design against two pockets at once.[10,11] How docking scores change meaning once ligands already carry both-end experimental labels still needs a separate evaluation.

Molecular docking searches ligand poses in a binding site and scores them. It remains a common ranking method in virtual screening.[3,4] Many benchmarks treat known actives versus artificial decoys as a binary task, as in DUD and DUD-E. LIT-PCBA instead uses high-throughput assay actives and inactives. DOCKSTRING provides a reproducible docking workflow for comparing ligand-design methods.[5–7,18] The source and definition of control compounds can change virtual-screening results and even reorder methods.[12,13,21] Rescoring hits from ultra-large docking screens is also not always stable.[15] Control-set composition is therefore part of virtual-screening evaluation itself.

When both-end experimental activities are available, dual-target screening can define four experimental states: dual, A-only, B-only, and neither. Dual differs from A-only or B-only at only one target, and from neither at both targets. Dual-versus-selective and dual-versus-neither discrimination therefore do not mean the same thing. Structure-based dual-target screening usually docks each candidate into both pockets and prefers molecules that score favorably on both. Zhou and coworkers evaluated docking-based dual-kinase virtual screening and found that single-target inhibitors were an important source of predicted dual-inhibitor false positives.[9] Their analysis covered several false-positive types, but it did not treat dual-versus-A-only and dual-versus-B-only as two directional tasks tied to the corresponding targets. Kinase-Bench further showed that experimental selectives can test structural discrimination between related kinases, with the aim of identifying inhibitors selective for one kinase. Distinguishing dual from A-only and B-only ligands among compounds already measured at both ends is a different problem.[22]

Prior work has shown that single-target selectives affect dual-target screening. Less systematically examined is how evaluation-task definition changes the interpretation of docking performance when both-end experimental activities are already available. Does changing the experimental-state comparison, at a fixed score channel, change the performance judgment? How far can ligand chemistry without receptor structure explain the apparent discrimination? Is remaining discrimination consistent with pocket-level structural information, and is that correspondence stable? Answers to these questions clarify what dual-target docking performance evaluation means and where the observed discrimination comes from.

---

## References

(1) Anighoro, Bajorath, Rastelli, *J. Med. Chem.* **2014**, *57*, 7874–7887.
(2) Proschak, Stark, Merk, *J. Med. Chem.* **2019**, *62*, 420–444.
(3) Kitchen et al., *Nat. Rev. Drug Discov.* **2004**, *3*, 935–949.
(4) Eberhardt et al., *J. Chem. Inf. Model.* **2021**, *61*, 3891–3898.
(5) Huang, Shoichet, Irwin, *J. Med. Chem.* **2006**, *49*, 6789–6801.
(6) Mysinger et al., *J. Med. Chem.* **2012**, *55*, 6582–6594.
(7) Tran-Nguyen, Jacquemard, Rognan, *J. Chem. Inf. Model.* **2020**, *60*, 4263–4273.
(9) Zhou, Li, Hou, *J. Chem. Inf. Model.* **2013**, *53*, 982–996.
(12) Tran-Nguyen, Ballester, *J. Chem. Inf. Model.* **2023**, *63*, 1401–1405.
(13) Ahmed, Soellner, Brooks, *J. Chem. Inf. Model.* **2026**, *66*, 8752–8759.
(19) Wu, Vigneron, Braz et al., *J. Med. Chem.* **2026**, *69*, 6210–6229.
(20) Munson, Chen, Bogosian et al., *Nat. Commun.* **2024**, *15*, 3636.
(21) Gu, Shen, Zhang et al., *Nat. Mach. Intell.* **2025**, *7*, 509–520.
(22) Wei, Zhou, Jing et al., *J. Chem. Inf. Model.* **2024**, *64*, 9528–9550.
