# Discussion (JCIM Articles draft, English)

## 4. Discussion

### 4.1 Benchmark formulation changes the evidentiary standard for dual-target docking

The primary question is whether negative-class definition changes apparent dual-target docking evidence. On the same frozen EGFR/HER2 scores, Dual versus neither gave AUROC 0.756, whereas directional `summary_min` was 0.430 (Results 3.2; Table 3). Holding the pocket A score fixed retained a descriptive difference of 0.378 [0.205, 0.547] between neither and B-only negatives (Table S34). A dual-versus-inactive formulation can therefore present a more favorable picture than the directional hard-negative task in this case; the remaining pairs do not support a universal effect. Under independent GNINA pose generation, Dual versus neither remained 0.783 and directional `summary_min` 0.220, so the EGFR/HER2 finding is not Vina-specific here (Results 3.2; Table S32).

Relative to Zhou et al.,[9] the comparison asks whether docking separates dual-actives from experimentally defined selectives, not only from inactives. Existing docking benchmarks have shown that decoy construction and chemical bias change virtual-screening interpretation;[5–7,12,13] the same concern applies when dual-target conclusions depend on the experimental negative class. Post-hoc AND-filter and full-map ECFP4 diagnostics in the Supporting Information reinforce that Dual versus neither is chemically easier than Dual versus selectives; they are not confirmatory primary results and do not expand docking to undocked pairs.

### 4.2 What docking adds—and does not add—beyond ligand chemistry

Physicochemical descriptors and chemotype already carry much of the experimental-label information. On AChE/BChE, TPSA alone discriminated better than docking on the corresponding contrast, and scaffold-grouped ECFP4 exceeded docking on several arms (Results 3.3). Adding the pocket-matched docking score to ECFP4 changed scaffold-grouped CV AUROC by at most 0.020. That result should not be read as “docking has no value” or “ECFP4 explains everything.” It shows that, on these frozen panels and scaffold-grouped tasks, docking did not provide a stable large incremental separation beyond 2D chemistry; the finding is limited to the present labels, ligand series, and receptor realizations. Without ligand-level controls, an apparently strong dual-target docking result may only recover molecular properties associated with the dual label.[7,12]

Document- and scaffold-cluster bootstrap keep the EGFR/HER2 weak arm near chance with intervals that span 0.5, so the formulation interpretation does not rest on ligand-independence assumptions alone. PIK3CA/mTOR’s four-member neither class is a single-document sample, and some document-blocked arms are not stably estimable; those cells are reported as such rather than imputed.

### 4.3 Receptor realization and evaluation conditions

Apparent discrimination also depends on how the evaluation condition is specified. Holding one pocket frozen and replacing the other raised discrimination on one PIK3CA-related pair and lowered it on the other (Figure 5), consistent with kinase cross-docking work that treats receptor representation as a performance variable.[14] Receptor sensitivity is therefore a boundary of the claim, not evidence of structural robustness.

### 4.4 Implications for dual-target virtual screening

Favorable scores in both pockets do not automatically establish experimentally defined dual activity. After a dual-pocket score looks favorable, four practical checks remain: (i) directional discrimination against A-only and B-only hard negatives; (ii) whether a ligand-only ECFP or property model recovers a similar signal under a leak-resistant split; (iii) an unused ligand pool or document-blocked split; (iv) at least one alternate receptor realization (Figure 8). A failure at any step marks the claim as formulation-, chemistry-, panel-, or receptor-dependent computational evidence. This concern is consistent with JCIM work showing that docking rescoring performance varies across experimentally grounded screening sets.[15]

### 4.5 Limitations

The docking evaluation contains only four target pairs because experimentally defined dual-target hard negatives are scarce under the construction gate. K = 4 is a data-constrained case panel, not a representative dual-target suite; primary `summary_min` values mix panel-construction differences with biology and should not be read as a target ranking. All four primary intervals include 0.5 (Table S31).

Operational labels require usable measurements at both targets, enriching jointly profiled chemistry (complete-case fractions 14.5%–34.0%). The unused-pool resample is internal, not external. BindingDB/PubChem Table S12 is count-level only. A BindingDB-native 202608 rebuild with literature, structure, and ECFP4 filters yielded zero pairs meeting the pre-frozen primary external gate and was not docked (Tables S48–S49).[16] A pre-frozen 2018 literature-year split was not evaluable as temporal external validation. Assay heterogeneity, max-pChEMBL aggregation, and incomplete construct/mutation metadata remain. Cognate best-of-nine QC establishes search coverage, not necessarily top-ranked pose correctness. MCL1/Bcl-xL was formally demoted and is not a fifth main pair. The claimed domain is the protocol-processable chemical subset under the frozen receptors and engines used here.
