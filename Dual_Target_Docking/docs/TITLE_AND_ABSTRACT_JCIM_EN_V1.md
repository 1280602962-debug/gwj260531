# Title and Abstract (JCIM Articles draft, English)

## Title

**A Four-Pair Formulation Audit of Docking-Based Dual-Target Recognition**

## Abstract

Whether favorable docking scores at two targets constitute evidence of dual-target activity depends on the negative class used for evaluation. We performed a four-pair formulation audit using ChEMBL-derived operational activity states: dual, A-only, B-only, and neither. Two pocket-matched directional AUROCs test dual ligands against the corresponding single-target selectives; their minimum (`summary_min`) is reported as a conservative descriptive summary. On EGFR/HER2, Vina Dual versus neither yielded AUROC 0.756, whereas directional `summary_min` was 0.430; 9 of the mixed-library Top-10 compounds were experimental selectives. An independent GNINA pose-generation protocol reproduced this descriptive gap (0.783 versus 0.220). The other pairs did not show the same separation, and the PIK3CA/mTOR neither comparator contained only four compounds. In scaffold-grouped models, adding docking to ECFP4 changed cross-validated AUROC by at most 0.020 in absolute value. Alternative receptors lowered PIK3CA/mTOR `summary_min` from 0.692 to 0.486/0.505 but raised PIK3CA/PIK3CB, demonstrating receptor dependence rather than structural robustness. All four primary `summary_min` confidence intervals included 0.5. Blocking ligands that share a ChEMBL document left the EGFR/HER2 directional weak arm at 0.430; a pre-frozen 2018 literature-year split did not yield two evaluable pairs and is not claimed as external validation. Thus, this data-constrained case panel supports selective hard negatives and confounder-aware controls as evaluation requirements, but does not establish target-general docking performance or biological recognition.

## Keywords

dual-target docking; benchmark formulation; selectivity hard negatives; chemical confounding; receptor realization; virtual screening
