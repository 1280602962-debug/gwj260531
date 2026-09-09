# Title and abstract (English working draft)

## Title

A Multi-Pair Formulation Audit of Docking-Based Dual-Target Recognition

## Abstract

Docking is often used to interpret dual-target recognition, but the negative class is usually a nonbinder or unmatched decoy rather than an experimentally selective ligand. We construct DualFourClass as a four-state evaluation of that choice on eight pairs: PIK3CA/mTOR, AChE/BChE, supply-limited EGFR/HER2, F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, and PPARA/PPARD. All pairs share one ChEMBL extract and four-state logic; candidate-pool rules, class quotas, and scaffold caps differ by pair because strict bidirectional supply and structure feasibility differ (Table 1). Holding the score channel fixed, Dual versus neither versus Dual versus B-only differed by 0.378 [0.205, 0.547] on EGFR/HER2 and 0.444 [0.263, 0.620] on JAK1/TYK2. Scaffold-cluster intervals excluded 0 for both; the EGFR/HER2 document-cluster interval also excluded 0, whereas the JAK1/TYK2 document-cluster interval included 0. Independent GNINA pose generation reproduced both task gaps. Worst-arm AUROC values ranged from 0.345 to 0.692; only PPARG/PPARA has a 95% interval entirely above 0.5, and that interval is membership- and formulation-sensitive. Adding a pocket-matched docking score to scaffold-grouped ECFP4 changed AUROC by at most 0.023, which does not imply that all structural information is uninformative. Matched-minus-mismatched `summary_min` excluded 0 on two main panels but included 0 on all scored holdouts, so a stable matched-pocket advantage was not obtained. No pair was packaged as an external evaluation set after the BindingDB-native remainder. On ligands with both-end measurements, neither-class controls can hide a failure to recognize single-target actives; directional selectivity comparisons and accompanying controls bound how dual-target screening results should be interpreted.

## Keywords

dual-target docking; selectivity; experimental-state comparison; AutoDock Vina; GNINA; virtual screening
