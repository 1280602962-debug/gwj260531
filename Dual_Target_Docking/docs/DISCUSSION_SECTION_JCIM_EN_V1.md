# Discussion (JCIM Articles draft, English)

## 4. Discussion

### 4.1 Benchmark formulation changes the evidentiary standard

The primary question is whether negative-class definition changes apparent dual-target docking evidence. On frozen EGFR/HER2 scores, Dual versus neither looked substantially stronger than the directional hard-negative task (Results 3.2; Table 3). Independent GNINA pose generation reproduced that gap, so the observation is not Vina-specific here. The EGFR/HER2 formulation gap remained positive across all five seeds (Table S54). The remaining pairs did not support a universal effect. Relative to Zhou et al.,[9] the comparison therefore asks whether docking separates dual-actives from experimentally defined selectives, not only from inactives. Existing docking benchmarks have shown that decoy construction and chemical bias change virtual-screening interpretation;[5–7,12,13,21] the same concern applies when dual-target conclusions depend on the experimental negative class. Wu et al. demonstrated that large-library docking can prospectively yield joint binders;[19] DualFourClass addresses a different question—whether retrospective evidence for dual-target recognition changes when experimentally selective ligands rather than nonbinders define the negative class. Kinase-Bench tests selective enrichment against kinase-specific decoys,[22] whereas DualFourClass constructs four experimentally measured states on each pair.

### 4.2 Ligand chemistry and receptor realization

Physicochemical descriptors and chemotype already carry much of the experimental-label information. Docking contributed little incremental discrimination beyond scaffold-grouped ligand-only models (Results 3.3; Table S24). Without ligand-level controls, an apparently strong dual-target docking result may only recover molecular properties associated with the dual label.[7,12]

Apparent discrimination also depends on receptor realization. Holding one pocket frozen and replacing the other raised discrimination on one PIK3CA-related pair and lowered it on the other (Results 3.4; Figure 5), consistent with kinase cross-docking work that treats receptor representation as a performance variable.[14]

### 4.3 Implications for dual-target virtual screening

Favorable scores in both pockets do not automatically establish experimentally defined dual activity. After a dual-pocket score looks favorable, four practical checks remain: (i) directional discrimination against A-only and B-only hard negatives; (ii) whether a ligand-only ECFP or property model recovers a similar signal under a leak-resistant split; (iii) an unused ligand pool or document-blocked split; (iv) at least one alternate receptor realization (Figure 8). POLYGON and related generative methods show that dual-target design can be experimentally validated,[20] but docking-based dual claims still need selectivity-aware evaluation before two-pocket occupancy is treated as recognition evidence. This concern is consistent with JCIM work showing that docking rescoring performance varies across experimentally grounded screening sets.[15]

### 4.4 Limitations

The analysis covers four data-rich target pairs rather than a target-general benchmark. Labels are derived from heterogeneous public bioactivity records and require complete measurements at both targets (complete-case fractions 14.5%–34.0%). Receptor dependence was evaluated only for selected structures, and no prospective experimental validation of newly predicted dual ligands was performed. No target pair met the prespecified BindingDB external-slice gate, so the slice is not claimed as external validation. Cognate best-of-nine QC establishes search coverage, not a top-ranked-pose validation. MCL1/Bcl-xL was formally demoted to Supporting Information.
