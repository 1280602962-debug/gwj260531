# Discussion (JCIM Articles draft, English)

> Companion to [`DISCUSSION_DRAFT_ZH_JCIM_V1.md`](DISCUSSION_DRAFT_ZH_JCIM_V1.md) (Chinese authoritative for this rewrite cycle).  
> Citation audit: [`DISCUSSION_REFS_JCIM_V1.md`](DISCUSSION_REFS_JCIM_V1.md); numbering continues from [`INTRODUCTION_REFS_JCIM_V1.md`](INTRODUCTION_REFS_JCIM_V1.md).  
> Finding → interpretation → alternative explanation → evidence → implication. Open results stay open.

## 4. Discussion

### 4.1 A Strict Dual-Target Benchmark Exposes a Task That Conventional Docking Evaluation Does Not Directly Test

The primary finding is not that one docking scoring function attained the highest dual-target AUROC. It is that a strict dual-target evaluation is a different task from conventional single-target virtual screening. Standard docking benchmarks typically separate actives from decoys. Here a model must distinguish dual-active ligands from single-target selective hard negatives in both directions. Those hard negatives are experimentally active on one target and therefore cannot be treated as ordinary decoys. The supply audit showed that public bioactivity data rarely provide enough such ligands on both arms at once; only a few of 49 candidate pairs met the strict thick-panel gate (Results 3.1; Figure 2).

Zhou, Li, and Hou already showed that docking-based dual-kinase screening can look useful against noninhibitors, is structure-dependent, and still admits a high false-positive rate among predicted duals.^(9) DualFourClass-Bench asks a narrower follow-up on the same scores: whether a Dual-versus-neither (inactive) readout and a directional Dual-versus-selective readout agree. They do not agree on EGFR/HER2 (Results 3.2, Table 3). That formulation gap—not another four-pair docking survey—is the increment relative to 2013.

That data constraint is itself methodological. DUD, DUD-E, and LIT-PCBA already showed that decoy construction, chemical bias, and assay-derived labels change virtual-screening conclusions.^(5–7) Simple methods and some unbiasing procedures can also overestimate structure-based virtual screening by learning ligand distributions.^(12) Recent bioassay-derived evaluation further emphasizes that real assay data can expose limits that constructed ligand/decoy sets conceal.^(13) DualFourClass-Bench does not use those single-target collections and does not evaluate DiffDock-Pocket. It extends the same concern to a dual-target setting: the evaluation depends on how hard negatives are defined experimentally, not on how long the candidate-pair list is.

The main value of DualFourClass-Bench is therefore not dataset size. It is the conversion of dual-target recognition into an experimentally labeled hard-negative discrimination task that must hold in both directions.

### 4.2 Why Dual-Target Docking Is Harder Than Two Single-Target Docking Tasks

The difficulty of dual-target docking is not the sum of two single-target docking problems. A favorable score for an A-active ligand in pocket A does not imply activity at target B. Strict dual recognition still requires suppressing A-only ligands in pocket B and B-only ligands in pocket A. A pooled score can hide the weaker arm; pocket-matched directional AUROC tests the extra recognition on the non-selective target more directly.

Even under that task-aligned metric, the four pairs were heterogeneous. Three pairs sat near chance or below the strongest trivial physicochemical baseline. Only PIK3CA/mTOR gave a higher point estimate, and its confidence interval remained compatible with chance (Results 3.2). On AChE/BChE, TPSA alone discriminated better than docking, and scaffold-grouped ECFP4 exceeded docking on several arms (Results 3.3). Dual/selective labels therefore carry ligand-level information that can produce a strong apparent signal.

That observation is consistent with recent attention to chemical bias in virtual-screening benchmarks: simple models or poorly constructed decoys can look strong by learning ligand composition rather than target-specific recognition.^(7,12) The present work extends the same requirement to dual-target evaluation. Without A-only/B-only hard negatives and ligand-property / chemical baselines, an apparently strong dual-target docking result may only be recovering molecular properties associated with the dual label.

### 4.3 The PIK3CA/mTOR Case: Limited Directional Signal Rather Than a Generalizable Rule

PIK3CA/mTOR is the case that most warrants further study and the case that most requires caution. Main-panel summary_min was 0.692, PM110 was 0.648, and the unused-pool holdout was 0.765, so the directional signal is not driven solely by a few PM48 members. It is not receptor-invariant: replacing the PIK3CA structure dropped summary_min to 0.486 and 0.505, whereas replacing the mTOR structure left 0.639 (Results 3.4). The accurate claim is not that PIK3CA/mTOR docking reliably identifies dual-target ligands. It is that this pair shows a limited directional signal under a particular receptor realization — persistent under ligand-panel replacement, not assumed invariant under receptor replacement.

Both PIK3CA and mTOR present ATP-site chemotypes with accessible binding modes, so some dual ligands can form reasonable hinge-oriented poses in both pockets. The same ATP-site compatibility can also give a selective hard negative a geometrically plausible pose in the second pocket, producing false dual recognition. Pose-level observations are compatible with that possibility (T2 / T5 in Results 3.6), but residue-level PLIF analysis was not completed, so the account is not a confirmed structural mechanism.

Local pocket Cα RMSD between 5DXT and 4L23 is only 0.343 Å, yet summary_min still fell to 0.505. Structural similarity and transferable screening discrimination are therefore not the same question. Passing pose-generation QC is also not screening-performance robustness. That distinction is consistent with recent kinase cross-docking benchmarks that treat receptor representation as an independent performance variable; those studies used different docking engines and are not an extrapolation of the present protocol.^(14)

### 4.4 Implications for Dual-Target Virtual Screening and Generative Design

The results bear directly on dual-target virtual screening and generative design. Favorable docking scores in both pockets do not automatically mean that a generated molecule is an experimentally plausible dual-active ligand. If the scoring function is rewarded by size, polarity, or chemotype, a generator can optimize those attributes into a high dual docking score without independent binding advantage at both targets. Even after single-target ultralarge docking, postprocessing and rescoring have been shown not to separate known binders from inactives robustly across assays.^(15) A dual-target setting additionally requires suppressing both experimental hard-negative arms, so two favorable scores — or their simple average — are not sufficient evidence.

Downstream evaluation of dual-target generative design should therefore include at least three layers: experimentally labeled discrimination of dual-actives from A-only/B-only selective hard negatives; ligand-property and ligand-only chemical baselines; and receptor-structure sensitivity. Reporting the two pocket scores, or their mean, does not answer those questions.

This study does not show that existing dual-target generators fail, and it does not evaluate DualDiff, FuseDiff, or other generative models.^(10,11) DualDiff’s Dual High Affinity is dual success versus reference-ligand dock scores, not mean pooling; FuseDiff’s independent test set is the DualDiff benchmark (DDF). Those papers ask whether generated structures can obtain favorable docking scores. DualFourClass-Bench can serve as a downstream evaluation layer for such methods: whether generated molecules outrank experimentally defined single-target hard negatives, rather than merely optimizing a docking score.

### 4.5 Wrong-Pocket Reversal Is an Unresolved Benchmark Failure Mode

The benchmark also exposes an unresolved failure mode: the wrong-pocket control reversed between the main panels and the unused-pool holdout. On the main panels, pocket-matched scores exceeded the wrong-pocket control. On the holdout, wrong-pocket scores were at least as high as matched-pocket scores on all three eligible pairs. Potency and size matching did not remove the reversal, and scoring-independent contact_count explained only part of the B-arm signal (Results 3.5).

The pattern should not be attributed to one docking engine’s scoring artifact: it appears under the same Vina protocol, and a coarse geometric metric also shows size/burial-related signal. Contact_count still cannot reproduce Vina wrong-pocket discrimination in magnitude, so the unique source remains unidentified. Ligand-distribution shift, pose selection, receptor-specific interaction patterns, and nonlinear size dependence of the scoring function are all possible; they were not isolated here.

For a benchmark paper, the unresolved result is itself useful. A pocket-specificity control that looks reasonable on a fixed panel need not hold on an unseen ligand pool. Future dual-target docking benchmarks should therefore report wrong-pocket, chemical-property controls, and a panel-external holdout together with matched-pocket performance.

### 4.6 Limitations

Only the five highest-priority limits are stated here; the full inventory is in the Limitations draft. Closing claims are in Conclusions and are not repeated here.

First, the benchmark contains only four target pairs because experimentally defined dual-target hard negatives are scarce. The observed target-pair heterogeneity should not be interpreted as a population-level estimate across all possible target pairs. The ligand-level bootstrap (B = 2000) describes uncertainty inside a fixed panel; leftover strict hard negatives after the main panels and holdout cannot support anything close to 1000 non-overlapping balanced panels.

Second, activity labels were aggregated using the maximum available pChEMBL value. This may inflate apparent activity when measurements vary across assays. The current frozen data package does not retain sufficient assay-level metadata to perform a complete max-versus-median or assay-confidence sensitivity analysis.

Third, the unused-pool holdout remains within the ChEMBL-derived data ecosystem and therefore should not be considered an independent cross-database validation.

Fourth, receptor-swap experiments demonstrate structure dependence but do not identify its molecular origin. Pocket-local Cα RMSD alone could not explain the observed performance change, and residue-level PLIF/side-chain conformational analyses were not systematically performed.

Finally, this study evaluates computational discrimination rather than experimentally validating newly predicted dual-target compounds. The benchmark addresses the reliability of docking-based ranking, not the prospective biological efficacy of selected molecules.
