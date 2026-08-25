# Discussion (JCIM Articles draft, English)

## 4. Discussion

### 4.1 Benchmark formulation changes the evidentiary standard for dual-target docking

The benchmark operationalizes dual-target recognition as discrimination against experimentally defined A-selective and B-selective ligands in both directions. On the same frozen EGFR/HER2 scores, Dual versus neither gave AUROC 0.756, whereas directional `summary_min` was 0.430 (Results 3.2; Table 3). A dual-versus-inactive formulation can therefore present a more favorable picture than the directional hard-negative task.

Relative to Zhou et al.,[9] the present comparison asks whether docking separates dual-actives from experimentally defined selectives, not only from inactives. Existing docking benchmarks have shown that decoy construction, chemical bias, and assay-derived labels change virtual-screening interpretation;[5–7,12,13] the same concern applies here, because dual-target conclusions depend on how the negative class is defined experimentally.

### 4.2 Apparent docking signal cannot be attributed to docking scores alone

Physicochemical descriptors and chemotype already carry much of the experimental-label information. On AChE/BChE, TPSA alone discriminated better than docking, and adding heavy-atom count and TPSA raised dual-versus-B-only AUROC while the docking-score odds ratio remained near one (Results 3.3). Scaffold-grouped ECFP4 exceeded docking on several arms—for example 0.85 versus 0.43 on EGFR/HER2 dual-versus-B-only—so the labels contain ligand-structure-associated information that can be used without receptor information. Adding the pocket-matched docking score to ECFP4 changed CV AUROC by at most 0.020 in absolute value. At T ≥ 0.3, the strongest unmatched arm (PIK3CA/PIK3CB dual versus A-only) fell from 0.691 to 0.503, whereas distant hard negatives (T < 0.3) rose to 0.819. Without these ligand-level controls, an apparently strong dual-target docking result may only be recovering molecular properties associated with the dual label.[7,12]

### 4.3 Receptor realization is another important dimension of the evaluation condition and can raise or lower apparent discrimination

Apparent discrimination also depends on how the evaluation condition is specified. Replacing maximum pChEMBL with the median of repeated measurements left pair-level conclusions essentially unchanged (Results 3.4). An unused-pool holdout from the same ChEMBL harvest preserved some ranking trends and reversed others. Holding one pocket frozen and replacing the other raised apparent discrimination on one PIK3CA-related pair and lowered it on the other (Figure 5). These shifts are consistent with kinase cross-docking work that treats receptor representation as a performance variable.[14]

### 4.4 Implications for dual-target virtual screening and generative design

These results have direct methodological implications for dual-target virtual screening and may also apply to generative design workflows that use docking as a downstream filter. Favorable scores in both pockets do not automatically establish experimentally defined dual activity. This concern is consistent with recent JCIM studies showing that docking rescoring performance can vary substantially across experimentally grounded screening sets.[15] In the present dual-target task, docking scores therefore need to be read together with experimentally defined selective hard negatives and ligand-only chemical controls.

The same four checks can be applied as a practical diagnostic (Figure 8). After a dual-pocket score looks favorable: (i) require directional discrimination against A-only and B-only hard negatives; (ii) ask whether a ligand-only ECFP or property model recovers a similar signal; (iii) test an unused ligand pool; (iv) replace at least one receptor realization. A failure at any step marks the claim as formulation-, chemistry-, panel-, or receptor-dependent computational evidence.

### 4.5 Limitations

First, the benchmark contains only four target pairs because experimentally defined dual-target hard negatives are scarce. K = 4 is a data-constrained case panel rather than a comprehensive dual-target suite. The four `summary_min` values also mix panel-construction differences (strict 6.5/5.5 versus θ = 6.0; unequal n) with target-pair biology. Current class sizes resolve large directional effects more readily than moderate ones (Table S31).

Second, experimental labels are ChEMBL-derived. The unused-pool holdout remains within the same extraction batch and is therefore not independent external validation. BindingDB/PubChem checks were count-level only.

Third, assay heterogeneity remains. Primary curation uses maximum pChEMBL. Confidence≥8 and Homo sapiens filters were not rebuilt.

Fourth, receptor replacement can raise or lower pair-level discrimination, but the experiments do not identify a molecular origin. The two receptor-sensitivity examples share PIK3CA.

Fifth, the primary protocol is AutoDock Vina; GNINA and RTMScore rescored the same poses. An independent pose-generation test is specified for local execution on EGFR/HER2 and PIK3CA/mTOR. This study does not include prospective testing of newly predicted dual-target compounds.
