# Title and Abstract (JCIM Articles draft, English)

## Title

**A Multi-Pair Formulation Audit of Docking-Based Dual-Target Recognition**

## Abstract

Docking is often used to interpret dual-target recognition, but the negative class is usually a nonbinder or unmatched decoy rather than an experimentally selective ligand. We construct DualFourClass as a four-state formulation audit of that choice. The original primary set was the 2026-07-23 four-pair freeze. After a post-hoc ChEMBL 37 universe census, five ordinary noncovalent pairs that passed pre-declared gates were added to the same analysis stack—collection completion, not a pre-registered eight-pair freeze. A docked fourth original candidate, PIK3CA/PIK3CB, was withdrawn after a receptor-identity audit found that its "PIK3CB" receptor (PDB 2WXF) is murine PIK3CD, not human PIK3CB, and is reported only as a documented receptor-identity failure. The primary table therefore has eight rows: PIK3CA/mTOR, AChE/BChE, supply-limited EGFR/HER2, F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, and PPARA/PPARD. Under AutoDock Vina, experimentally defined dual ligands were ranked against A-only and B-only selectives. Worst-arm AUROC values were 0.430 (EGFR/HER2), 0.606 (AChE/BChE), 0.692 (PIK3CA/mTOR), 0.345 (F2/F10), 0.365 (JAK1/TYK2), 0.588 (JAK1/JAK2), 0.649 (PPARG/PPARA), and 0.446 (PPARA/PPARD). Only PPARG/PPARA has a 95% confidence interval entirely above 0.5, and that interval is membership- and formulation-sensitive. On EGFR/HER2 the same scores looked much stronger when dual ligands were compared with neither-class negatives (AUROC 0.756) than in the directional comparison (0.430); JAK1/TYK2 showed the same formulation gap (0.770 versus 0.365). Independent GNINA pose generation reproduced both gaps. Scaffold-grouped ligand-only models already captured most of the apparent ranking, and substituting alternative PIK3CA crystal structures dropped the PIK3CA/mTOR contrast to near chance. Two-pocket score filters therefore need selectivity-aware negatives; favorable occupancy of two pockets is not, by itself, evidence of dual-target recognition.

## Keywords

dual-target docking; selectivity; hard negatives; AutoDock Vina; GNINA; virtual screening
