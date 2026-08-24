# Title and Abstract (JCIM Articles draft, English)

## Title

**Benchmark Formulation and Chemical Confounding in Docking-Based Dual-Target Recognition**

## Abstract

Whether favorable docking scores at two targets constitute evidence of dual-target recognition has not been adequately tested against experimentally defined single-target selectives. We constructed DualFourClass-Bench, a curated four-pair, four-state panel with two directional primary tasks: dual-actives versus A-only selectives scored in pocket B, and dual-actives versus B-only selectives scored in pocket A; the weaker arm is the pair-level summary (`summary_min`). On the same frozen AutoDock Vina scores, a Dual-versus-neither comparator produced a substantially stronger impression on EGFR/HER2 (AUROC 0.756 versus directional `summary_min` 0.430), while 9 of the Top-10 mixed-library ligands were experimental selectives. AChE/BChE and PIK3CA/PIK3CB showed only small, overlapping formulation increments, and the PIK3CA/mTOR Dual-versus-neither comparator was underpowered (neither n = 4). Under scaffold-grouped cross-validation, adding docking changed AUROC by at most approximately 0.02 in absolute value beyond ECFP4. Replacing maximum pChEMBL with the median among repeated measurements produced minimal pair-level changes. In contrast, alternative receptor realizations moved PIK3CA/mTOR from 0.692 to 0.486/0.505 but PIK3CA/PIK3CB from 0.500 to 0.691/0.685. Thus, apparent dual-target discrimination depended on benchmark formulation, target pair, ligand chemistry, and receptor realization, with limited incremental information from docking beyond ligand-level chemical baselines in this scaffold-aware evaluation. These results support experimentally defined selectivity hard negatives and confounder-aware controls; they do not establish a universal overestimation law or prove that docking lacks pocket-specific information.

## Keywords

dual-target docking; benchmark formulation; selectivity hard negatives; chemical confounding; receptor realization; virtual screening
