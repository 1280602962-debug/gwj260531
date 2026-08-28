# Title and Abstract (JCIM Articles draft, English)

## Title

**A Four-Pair Formulation Audit of Docking-Based Dual-Target Recognition**

## Abstract

Benchmark construction determines what favorable scores at two targets can establish about dual activity. We assigned ChEMBL ligands for four frozen target pairs to dual, A-only, B-only, or neither states and evaluated AutoDock Vina in two pocket-matched directions: dual versus A-only in pocket B and dual versus B-only in pocket A. The smaller directional AUROC (`summary_min`) served as a conservative descriptive summary. On EGFR/HER2, Dual versus neither produced an AUROC of 0.756, whereas `summary_min` was 0.430; independent GNINA pose generation preserved the contrast (0.783 versus 0.220). The other pairs did not show the same gap. Across five frozen Vina seeds, the sign of the Dual-versus-neither minus `summary_min` gap was unchanged for every pair, with median `summary_min` values of 0.373, 0.599, 0.478, and 0.704 for EGFR/HER2, AChE/BChE, PIK3CA/PIK3CB, and PIK3CA/mTOR, respectively. Adding docking to scaffold-grouped ECFP4 models changed cross-validated AUROC by at most 0.020; alternative PIK3CA receptors shifted discrimination in opposite directions across two related pairs; and all four primary-seed `summary_min` confidence intervals included 0.5. A pre-frozen BindingDB-native gate yielded no eligible external pair. Thus, selective hard negatives, ligand-only controls, and receptor and seed sensitivities are necessary for interpreting dual-pocket docking, while the present four-pair panel does not establish target-general performance.

## Keywords

dual-target docking; benchmark formulation; selectivity hard negatives; chemical confounding; receptor realization; virtual screening
