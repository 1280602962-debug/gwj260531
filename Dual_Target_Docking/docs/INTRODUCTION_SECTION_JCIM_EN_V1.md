# Introduction (JCIM Articles draft, English)

## 1. Introduction

Multitarget drug design uses one molecule to modulate two or more disease-related targets. Relative to combination therapy, a single multitarget ligand can reduce pharmacokinetic mismatch and formulation complexity.[1,2] Generative design and ultra-large library docking have already produced experimentally confirmed dual-target candidates. POLYGON tested 32 MEK1/mTOR compounds.[19] Ultra-large library docking has also found multitarget ligands among hundreds of millions of molecules, and experimental structures support some predicted binding modes.[18] Such work shows that computational methods can find dual-active molecules in a candidate library. Structure-based dual-target generative models can likewise design against two pockets at once.[10,11] How docking scores change meaning once ligands already carry both-end experimental labels still needs a separate evaluation.

Molecular docking searches ligand poses in a binding site and scores them. It remains a common ranking method in virtual screening.[3,4] Many benchmarks treat known actives versus artificial decoys as a binary task, as in DUD and DUD-E. LIT-PCBA instead uses high-throughput assay actives and inactives. DOCKSTRING provides a reproducible docking workflow for comparing ligand-design methods.[5–7,17] The source and definition of control compounds can change virtual-screening results and even reorder methods.[12,13,20] Rescoring hits from ultra-large docking screens is also not always stable.[15] Control-set composition is therefore part of virtual-screening evaluation itself.

When the same ligand has usable measurements at both targets, an activity threshold assigns it to one of four states: dual-active at both ends; selective for target A or target B; or below threshold at both ends. These states are labeled dual, A-only, B-only, and neither. Dual differs from A-only or B-only at only one target, and from neither at both, so dual-versus-A-only, dual-versus-B-only, and dual-versus-neither answer different questions. Structure-based dual-target screening usually docks each candidate into both pockets and prefers molecules that score favorably on both. Zhou and coworkers evaluated docking-based dual-kinase virtual screening and found that single-target inhibitors were an important source of predicted dual-inhibitor false positives.[9] Their analysis covered several false-positive types, but it did not treat the two selective states as separate directional tasks tied to the corresponding targets. Kinase-Bench uses experimental selectives to test structural discrimination between related kinases, with the aim of identifying inhibitors selective for one kinase. Distinguishing dual from A-only and B-only ligands among compounds already measured at both ends is a different problem.[21]

Prior work has shown that single-target selectives affect dual-target screening. What remains less clear is how the evaluation task itself changes the interpretation of docking performance when both-end experimental activities are already available. We asked three questions. First, after the score channel is held fixed, does changing the experimental-state comparison change the performance judgment? Second, what discrimination can ligand chemistry obtain without receptor structure? Third, is the observed docking discrimination consistent with the corresponding pocket, and is that correspondence stable under the available controls?

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
(18) Wu, Vigneron, Braz et al., *J. Med. Chem.* **2026**, *69*, 6210–6229.
(19) Munson, Chen, Bogosian et al., *Nat. Commun.* **2024**, *15*, 3636.
(20) Gu, Shen, Zhang et al., *Nat. Mach. Intell.* **2025**, *7*, 509–520.
(21) Wei, Zhou, Jing et al., *J. Chem. Inf. Model.* **2024**, *64*, 9528–9550.
