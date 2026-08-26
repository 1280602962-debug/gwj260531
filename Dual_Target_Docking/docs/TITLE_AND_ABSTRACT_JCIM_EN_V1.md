# Title and Abstract (JCIM Articles draft, English)

## Title

**Benchmark Formulation and Chemical Confounding in Docking-Based Dual-Target Recognition**

## Abstract

Whether favorable docking scores at two targets constitute evidence of dual-target recognition has not been adequately tested against experimentally defined single-target selectives. We therefore evaluated whether benchmark formulation changes the apparent evidence for dual-target recognition. To this end, we constructed DualFourClass-Bench, a curated four-pair, four-state panel with two directional primary tasks and a conservative worst-direction discrimination summary (`summary_min`). Within this four-pair, Vina-based benchmark, EGFR/HER2 Dual versus neither yielded AUROC 0.756, whereas directional `summary_min` was 0.430; in the EGFR/HER2 mixed-library Top-10, 9 of 10 compounds were experimental selectives. An independent GNINA pose-generation protocol left that gap intact (0.783 versus 0.220). The other pairs did not show the same gap, and the PIK3CA/mTOR Dual-versus-neither comparator was underpowered. In scaffold-grouped models, the largest absolute AUROC change after adding docking to ECFP4 was 0.020, while replacing maximum pChEMBL with the median of repeated measurements produced minimal pair-level changes. Alternative receptors moved PIK3CA/mTOR from 0.692 to 0.486/0.505 and raised PIK3CA/PIK3CB. These findings define a reliability boundary within the experimentally supported target pairs and protocols evaluated here, and support selective hard negatives and confounder-aware controls as complementary requirements for dual-target docking evaluation.

## Keywords

dual-target docking; benchmark formulation; selectivity hard negatives; chemical confounding; receptor realization; virtual screening
