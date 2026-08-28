# Title and Abstract (JCIM Articles draft, English)

## Title

**A Four-Pair Formulation Audit of Docking-Based Dual-Target Recognition**

## Abstract

Docking is often used to interpret dual-target recognition, but the negative class is usually a nonbinder or unmatched decoy rather than an experimentally selective ligand. We freeze DualFourClass as a four-state, four-pair formulation audit of that choice. Under AutoDock Vina, experimentally defined dual ligands were ranked against A-only and B-only selectives. Worst-arm AUROC values were 0.430 (EGFR/HER2), 0.606 (AChE/BChE), 0.500 (PIK3CA/PIK3CB), and 0.692 (PIK3CA/mTOR); all 95% confidence intervals included 0.5. On EGFR/HER2 the same scores looked much stronger when dual ligands were compared with neither-class negatives (AUROC 0.756) than in the directional comparison (0.430). Independent GNINA pose generation reproduced that formulation gap (0.783 versus 0.220). Scaffold-grouped ligand-only models already captured most of the apparent ranking, and substituting alternative PIK3CA crystal structures reversed the sign of the PI3K-family contrast. Two-pocket score filters therefore need selectivity-aware negatives; favorable occupancy of two pockets is not, by itself, evidence of dual-target recognition.

## Keywords

dual-target docking; selectivity; hard negatives; AutoDock Vina; GNINA; virtual screening
