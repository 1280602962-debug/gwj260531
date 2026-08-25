# Discussion (JCIM Articles draft, English)

## 4. Discussion

### 4.1 Benchmark formulation changes the evidentiary standard for dual-target docking

The primary finding is that benchmark formulation can change what dual-target docking success appears to mean. Standard docking benchmarks typically separate actives from decoys. Here a model must distinguish dual-active ligands from single-target selective hard negatives in both directions. Those hard negatives are experimentally active on one target. Dual versus neither on the same scores is used as a nonselectivity-controlled comparator. Public bioactivity data rarely provide enough such ligands on both arms at once (Results 3.1; Figure 2).

Relative to the dual-target evaluation setting of Zhou et al.,[9] this study adds experimentally defined A-only/B-only directional hard negatives and compares formulations on the same scores. On EGFR/HER2, Dual versus neither yielded 0.756, directional `summary_min` yielded 0.430, and mixed-library Top-10 ranking enriched selectives (9/10). AChE/BChE and PIK3CA/PIK3CB showed only small, overlapping increments, and the PIK3CA/mTOR comparator was underpowered (neither n = 4).

That data constraint is itself methodological. DUD, DUD-E, and LIT-PCBA already showed that decoy construction, chemical bias, and assay-derived labels change virtual-screening conclusions,[5–7] and simple methods can overestimate structure-based virtual screening by learning ligand distributions.[12] Recent bioassay-derived evaluation further emphasizes that real assay data can expose limits that constructed ligand/decoy sets conceal.[13] DualFourClass-Bench extends the same concern to a dual-target setting: evaluation depends on how hard negatives are defined experimentally.

### 4.2 Ligand-only chemical baselines reveal confounding in apparent docking discrimination

The difficulty of dual-target docking is not the sum of two single-target docking problems. All four `summary_min` bootstrap intervals included 0.5, so no pair yielded clear evidence excluding chance-level discrimination at the present sample sizes. Under the primary frozen-receptor protocol, PIK3CA/mTOR had the highest point estimate (0.692), but its paired difference from the best single-descriptor reference included 0 (Results 3.2; Table S19).

On AChE/BChE, TPSA alone discriminated better than docking, and scaffold-grouped ECFP4 exceeded docking on several arms (Results 3.3). The ligand-only baseline shows that the experimental labels contain ligand-structure-associated information that can be exploited without receptor information; docking performance must therefore be interpreted relative to ligand-only baselines. Under the present scaffold-grouped benchmark, adding docking to ECFP4 produced little incremental improvement in CV AUROC (largest absolute change 0.020). This panel-specific result does not establish that docking provides no additional information.

Chemotype-constrained hard negatives make the same point at a coarser similarity cutoff. T ≥ 0.7 matched subsets were empty. At T ≥ 0.3, the strongest unmatched arm (PIK3CA/PIK3CB dual versus A-only, 0.691) fell to 0.503 (n_neg = 11), whereas distant hard negatives (T < 0.3) rose to 0.819. Without A-only/B-only hard negatives and ligand-property or chemical baselines, an apparently strong dual-target docking result may only be recovering molecular properties associated with the dual label.[7,12]

### 4.3 Receptor realization is another dimension of the evaluation condition

PIK3CA/mTOR is the pair with the highest primary point estimate and the largest receptor-swap movement. Main-panel summary_min was 0.692, PM110 was 0.648, and the unused-pool holdout was 0.765. Their common direction suggests that the result is not determined solely by the exact PM48 membership, but both checks remain within the same ChEMBL data ecosystem. Replacing PIK3CA dropped summary_min to 0.486 and 0.505, whereas replacing mTOR left 0.639 (Results 3.4; Figure 5A).

The same PIK3CA replacement on PIK3CA/PIK3CB moved the estimate from 0.500 to 0.691 and 0.685, while the B-end receptor was held frozen (Figure 5B). Thus, replacing the same PIK3CA receptor raised or lowered apparent discrimination depending on the target pair. Receptor realization is another dimension of the evaluation condition. Both examples share PIK3CA, so additional pairs would be needed before generalizing. Simple backbone similarity also failed to predict the PIK3CA/mTOR change (Results 3.6), consistent with recent kinase cross-docking work that treats receptor representation as an independent performance variable.[14]

### 4.4 Implications for dual-target virtual screening and generative design

These results have potential implications for dual-target virtual screening and for generative methods that use docking as a downstream evaluation. Favorable docking scores in both pockets do not automatically mean that a molecule is an experimentally plausible dual-active ligand. If dual-target generators use docking as a filter, they should also be evaluated against selective hard negatives and ligand-only chemical controls. DualFourClass-Bench can serve as a downstream evaluation layer for such methods, testing whether generated molecules outrank experimentally defined single-target hard negatives.[10,11]

### 4.5 Limitations

Five limitations define the scope of interpretation.

First, the benchmark contains only four target pairs because experimentally defined dual-target hard negatives are scarce. K = 4 is a data-constrained case panel. The four `summary_min` values also mix panel-construction differences (strict 6.5/5.5 versus θ = 6.0; unequal n) with target-pair biology and should not be read as a ranking of intrinsic docking performance.

Second, experimental labels are ChEMBL-derived. The unused-pool holdout remains within the same extraction batch. BindingDB/PubChem checks were count-level only.

Third, assay heterogeneity remains after activity-aggregation control. Primary curation uses maximum pChEMBL; replacing that rule with the median among repeated measurements produced minimal pair-level changes (Results 3.4; Table S29). pChEMBL measurements are still not fully assay-equivalent. Confidence≥8 and Homo sapiens filters were not rebuilt.

Fourth, receptor realization can raise or lower pair-level discrimination, but the experiments do not identify a molecular origin. The two receptor-sensitivity examples share PIK3CA. Residue-level PLIF/side-chain analyses were not systematically performed.

Fifth, this study evaluates computational discrimination and does not include prospective testing of newly predicted dual-target compounds. The benchmark addresses the reliability of docking-based ranking.
