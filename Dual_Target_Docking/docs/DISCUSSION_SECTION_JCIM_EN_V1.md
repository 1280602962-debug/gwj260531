# Discussion (JCIM Articles draft, English)

## 4. Discussion

### 4.1 Benchmark formulation changes the evidentiary standard for dual-target docking

Negative-class definition materially changed the EGFR/HER2 result: Dual versus neither gave AUROC 0.756, whereas directional `summary_min` was 0.430 (Table 3). When aggregation was held fixed by using the same pocket A score, replacing B-only with neither negatives increased AUROC by 0.378 [0.205, 0.547] (Table S34). Independent GNINA pose generation and all five Vina seeds preserved a positive formulation gap for this pair (Tables S32, S54). The other pairs did not reproduce a gap of comparable magnitude, making the finding a pair-specific failure mode rather than a general property of dual-target docking.

Relative to Zhou et al.,[9] the comparison asks whether docking separates dual-actives from experimentally defined selectives, not only from inactives. Existing docking benchmarks have shown that decoy construction and chemical bias change virtual-screening interpretation;[5–7,12,13] the same concern applies when dual-target conclusions depend on the experimental negative class. Post-hoc AND-filter and full-map ECFP4 diagnostics in the Supporting Information reinforce that Dual versus neither is chemically easier than Dual versus selectives; they are not confirmatory primary results and do not expand docking to undocked pairs.

### 4.2 What docking adds—and does not add—beyond ligand chemistry

Physicochemical descriptors and chemotype carried substantial experimental-label information. On AChE/BChE, TPSA alone exceeded docking on the corresponding contrast, and scaffold-grouped ECFP4 exceeded docking on several arms (Results 3.3). Adding the pocket-matched docking score to ECFP4 changed scaffold-grouped CV AUROC by at most 0.020. Within these panels, ligand series, and receptor realizations, docking therefore supplied little measurable incremental discrimination beyond 2D chemistry. Ligand-only controls are consequently needed before attributing an apparent dual-target signal to receptor-specific complementarity.[7,12]

Document- and scaffold-cluster bootstrap keep the EGFR/HER2 weak arm near chance with intervals that span 0.5, so the formulation interpretation does not rest on ligand-independence assumptions alone. PIK3CA/mTOR’s four-member neither class is a single-document sample, and some document-blocked arms are not stably estimable; those cells are reported as such rather than imputed.

### 4.3 Receptor realization and evaluation conditions

Receptor realization was also a performance variable. Holding one pocket fixed and replacing the PIK3CA crystal raised discrimination for PIK3CA/PIK3CB but lowered it for PIK3CA/mTOR (Figure 5), consistent with kinase cross-docking studies.[14] A single receptor structure is therefore insufficient to support a robustness claim.

### 4.4 Implications for dual-target virtual screening

Favorable scores in both pockets do not automatically establish experimentally defined dual activity. After a dual-pocket score looks favorable, four practical checks remain: (i) directional discrimination against A-only and B-only hard negatives; (ii) whether a ligand-only ECFP or property model recovers a similar signal under a leak-resistant split; (iii) an unused ligand pool or document-blocked split; (iv) at least one alternate receptor realization (Figure 8). A failure at any step marks the claim as formulation-, chemistry-, panel-, or receptor-dependent computational evidence. This concern is consistent with JCIM work showing that docking rescoring performance varies across experimentally grounded screening sets.[15]

### 4.5 Limitations

The four target pairs are a data-constrained case panel, not a representative suite, and differences in panel construction prevent interpreting `summary_min` as a target ranking. All four primary-seed intervals include 0.5. Requiring measurements at both targets enriches jointly profiled chemistry: complete-case fractions were 14.5%–34.0%. Neither the unused-pool resample nor the failed 2018 time split provides external validation, and the BindingDB-native rebuild yielded no pair meeting the frozen external gate (Tables S41, S48–S49).[16]

The activity labels also combine heterogeneous assays, use maximum pChEMBL for primary aggregation, and lack resolved construct or mutation annotations for the audited priority set. Cognate best-of-nine QC establishes search coverage rather than top-ranked-pose accuracy. MCL1/Bcl-xL was therefore retained only as an exploratory stress test after its topology-aware pose-gold gate could not be established. These limitations restrict inference to the processable compounds, receptors, and engines examined here.
