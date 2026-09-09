# Title and abstract (English working draft)

## Title

**Docking-Based Dual-Target Recognition: Multi-Pair Evaluation Design and Sources of Discrimination**

## Abstract

Molecular docking is often used to rank and prioritize candidate ligands in dual-target virtual screening. In retrospective evaluation, the control set may be compounds below the activity threshold on both targets, or experimentally selective single-target ligands. We built a four-state evaluation on eight human target pairs, comprising dual, A-only, B-only, and neither ligands. Holding the same target score channel fixed, changing the control class altered apparent docking discrimination, and the change was not uniform across pairs or evaluation directions. On EGFR/HER2, the EGFR-pocket AUROC for dual-versus-B-only was 0.430 and rose to 0.808 against neither (difference 0.378 [0.205, 0.547]). JAK1/TYK2 showed a similar difference (0.444 [0.263, 0.620]), and independent GNINA pose generation reproduced both task differences. A scaffold-grouped ligand-only ECFP4 baseline matched docking on several directions. Adding the matched-pocket docking score changed AUROC by at most 0.023. Only the EGFR/HER2 and AChE/BChE main panels had higher directional AUROC summaries for matched pockets. On the seven internal holdouts, matched-minus-mismatched intervals all included 0. A matched-pocket advantage was not stably recovered. EGFR/HER2 has no unused-pool holdout. Among compounds with both-end experimental activities, strong dual-versus-neither discrimination does not necessarily imply strong discrimination against single-target selectives. Retrospective dual-target evaluation should report both directional selectivity comparisons, together with ligand-chemistry and pocket-correspondence controls.

## Keywords

dual-target docking; experimental state; selectivity; AutoDock Vina; virtual screening evaluation
