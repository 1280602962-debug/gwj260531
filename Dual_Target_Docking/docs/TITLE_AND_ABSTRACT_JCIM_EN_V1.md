# Title and Abstract (JCIM Articles draft, English)

## Title

**A Four-Pair Formulation Audit of Docking-Based Dual-Target Recognition**

## Abstract

Whether favorable docking scores at two targets constitute evidence of dual-target activity depends on how the negative class is defined for evaluation. Using ChEMBL-derived four-state labels (dual, A-only, B-only, neither) on four frozen target pairs, we evaluated AutoDock Vina with two pocket-matched directional AUROCs that compare dual ligands with the corresponding single-target selectives; their minimum (`summary_min`) is reported only as a conservative descriptive summary, with uncertainty from ligand, document-cluster, and scaffold-cluster bootstrap resampling of the same primary contrasts. On EGFR/HER2, Dual versus neither gave AUROC 0.756 while directional `summary_min` was 0.430, and an independent GNINA pose-generation protocol reproduced that descriptive gap (0.783 versus 0.220); the other pairs did not show the same separation, scaffold-grouped models gained at most 0.020 AUROC when docking was added to ECFP4, alternative receptors changed `summary_min` in opposite directions across pairs, all four primary `summary_min` intervals included 0.5, and a pre-frozen BindingDB-native external gate yielded zero dockable pairs. These data-constrained results, taken as a whole, support selective hard negatives and confounder-aware controls as evaluation requirements but do not establish target-general docking performance or external transfer.

## Keywords

dual-target docking; benchmark formulation; selectivity hard negatives; chemical confounding; receptor realization; virtual screening
