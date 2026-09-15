# Title and abstract (English working draft)

## Title

**Evaluating Dual-Target Molecular Docking: Experimental-State Comparisons across Multiple Target Pairs and Sources of Discrimination**

## Abstract

Molecular docking is often used to rank candidates in dual-target virtual screening. In retrospective evaluation, the control set may be compounds below the activity threshold at both targets, or single-target selectives that remain active at one target. We defined a four-state evaluation on eight human target pairs (dual, A-only, B-only, and neither) and compared those states after holding the same target score channel fixed. Changing the control class and its compounds altered apparent docking discrimination, and the change varied by pair and direction. On EGFR/HER2, the EGFR-pocket AUROC for dual-versus-B-only was 0.324 and rose to 0.786 against neither (difference 0.462 [0.262, 0.651]); JAK1/TYK2 showed a similar difference (0.444 [0.263, 0.620]). Effects on the other pairs were smaller or less certain. Receptor-free ECFP4 carried substantial class information on several directions. Under the same scaffold-grouped logistic-regression setting, adding the matched-pocket docking score changed AUROC by at most 0.023. Only one of eight main-panel matched-minus-mismatched 95% intervals excluded zero (AChE/BChE; EGFR/HER2 included zero), and none of the seven available holdout intervals did; these intervals were not multiplicity-adjusted. These results support a narrower interpretation: among compounds with experimental measurements at both targets, strong dual-versus-neither discrimination does not establish strong discrimination against single-target selectives. Retrospective dual-target evaluation should report both directional selectivity comparisons, together with ligand-chemistry and pocket-correspondence controls.

## Keywords

dual-target docking; experimental state; selectivity; AutoDock Vina; virtual screening evaluation
