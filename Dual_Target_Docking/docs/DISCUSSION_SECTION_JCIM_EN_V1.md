# Discussion (JCIM Articles draft, English)

## 4. Discussion

### 4.1 Benchmark formulation changes the evidentiary standard for dual-target docking

The primary finding is that benchmark formulation can change what dual-target docking success appears to mean. Standard docking benchmarks typically separate actives from decoys. Here a model must distinguish dual-active ligands from single-target selective hard negatives in both directions. Dual versus neither on the same scores is used as a nonselectivity-controlled comparator. Public bioactivity data rarely provide enough such ligands on both arms at once (Results 3.1; Figure 2).

Relative to the dual-target evaluation setting of Zhou et al.,[9] this study adds experimentally defined A-only/B-only directional hard negatives and compares formulations on the same scores. EGFR/HER2 provides the clearest example: a nonselectivity-controlled comparator suggests favorable discrimination, whereas the directional hard-negative task does not. The corresponding mixed-library ranking likewise favored selective ligands. AChE/BChE and PIK3CA/PIK3CB showed only small, overlapping increments, and the PIK3CA/mTOR comparator was underpowered (neither n = 4).

Existing docking benchmarks have shown that decoy construction, chemical bias, and assay-derived labels change virtual-screening interpretation;[5–7,12,13] the present results extend that concern to a dual-target task whose conclusions depend on how hard negatives are defined experimentally.

### 4.2 Ligand-only chemical baselines reveal confounding in apparent docking discrimination

Apparent dual-target docking performance must be interpreted together with ligand-only baselines. All four `summary_min` bootstrap intervals included 0.5, so no pair yielded clear evidence excluding chance-level discrimination at the present sample sizes. Under the primary frozen-receptor protocol, PIK3CA/mTOR had the highest point estimate, but its paired difference from the best single-descriptor reference included 0 (Results 3.2; Table S19).

On AChE/BChE, TPSA alone discriminated better than docking, and scaffold-grouped ECFP4 exceeded docking on several arms (Results 3.3). The experimental labels therefore contain ligand-structure-associated information that can be exploited without receptor information. Adding docking to ECFP4 produced little incremental improvement in CV AUROC (largest absolute change 0.020). At T ≥ 0.3, the strongest unmatched arm (PIK3CA/PIK3CB dual versus A-only, 0.691) fell to 0.503 (n_neg = 11), whereas distant hard negatives (T < 0.3) rose to 0.819. Without A-only/B-only hard negatives and ligand-property or chemical baselines, an apparently strong dual-target docking result may only be recovering molecular properties associated with the dual label.[7,12]

### 4.3 Receptor realization is another dimension of the evaluation condition

On PIK3CA/mTOR, replacing PIK3CA lowered summary_min from 0.692 to 0.486 and 0.505, whereas replacing mTOR left 0.639 (Results 3.4; Figure 5A). The same PIK3CA replacement on PIK3CA/PIK3CB raised the estimate from 0.500 to 0.691 and 0.685 (Figure 5B). Receptor realization can therefore raise or lower apparent discrimination depending on the target pair. Both examples share PIK3CA, so additional pairs would be needed before generalizing. The contrasting shifts may reflect changes in the limiting arm, local pocket geometry, or ligand-distribution effects; the present data cannot distinguish among these possibilities. Simple backbone similarity also failed to predict the PIK3CA/mTOR change (Results 3.6), consistent with kinase cross-docking work that treats receptor representation as a performance variable.[14]

### 4.4 Implications for dual-target virtual screening and generative design

These findings have implications for dual-target virtual screening and may be relevant to generative design workflows that use docking as a downstream filter. Favorable scores in both pockets do not automatically establish experimentally defined dual activity. If dual-target generators use docking as a filter, they should also be evaluated against selective hard negatives and ligand-only chemical controls. DualFourClass-Bench can serve as a downstream evaluation layer for such methods, testing whether generated molecules outrank experimentally defined single-target hard negatives.[10,11]

### 4.5 Limitations

First, the benchmark contains only four target pairs because experimentally defined dual-target hard negatives are scarce. K = 4 is a data-constrained case panel rather than a comprehensive dual-target suite. The four `summary_min` values also mix panel-construction differences (strict 6.5/5.5 versus θ = 6.0; unequal n) with target-pair biology.

Second, experimental labels are ChEMBL-derived. The unused-pool holdout remains within the same extraction batch and is therefore not independent external validation. BindingDB/PubChem checks were count-level only.

Third, assay heterogeneity remains after activity-aggregation control. Primary curation uses maximum pChEMBL; replacing that rule with the median among repeated measurements produced minimal pair-level changes (Results 3.4; Table S29). Confidence≥8 and Homo sapiens filters were not rebuilt.

Fourth, receptor realization can raise or lower pair-level discrimination, but the experiments do not identify a molecular origin. The two receptor-sensitivity examples share PIK3CA.

Fifth, this study evaluates computational discrimination and does not include prospective testing of newly predicted dual-target compounds.
