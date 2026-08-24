# Title and Abstract (JCIM Articles draft, English)

> Working draft. Not yet typeset. Claim ceiling: [`CLAIM_CEILING.md`](../data/jcim_bench_v0/CLAIM_CEILING.md).  
> Dual-versus-neither is a **nonselectivity-controlled comparator**, not “the conventional dual-target benchmark.”  
> Dual-vs-neither vs directional AUROC is a **descriptive formulation contrast**, not a paired significance test.

---

## Title (preferred)

**Benchmarking Docking-Based Dual-Target Recognition with Directional Selectivity Hard Negatives**

Alternate (more pointed; slightly perspective-like; not preferred for a first JCIM submission):

*Benchmark Formulation Matters: Evaluating Docking-Based Dual-Target Recognition with Directional Selectivity Hard Negatives*

Do **not** use: “Docking can/cannot identify dual-target ligands”; DualFourClass as a comprehensive suite; Dual-versus-neither as “the conventional benchmark.”

---

## Abstract (draft)

Whether favorable docking scores at two targets constitute evidence of dual-target recognition has not been adequately tested against directional single-target hard negatives. We constructed DualFourClass-Bench, a curated four-pair, four-state evaluation panel with two directional primary tasks: dual-actives versus A-only selectives scored in pocket B, and dual-actives versus B-only selectives scored in pocket A. The pair-level summary is the weaker arm (`summary_min`). On the same frozen AutoDock Vina scores, a Dual-versus-neither comparator that omits selectives can give an overly favorable impression of dual-target recognition in some pair contexts. EGFR/HER2 is the proof-of-principle case: Dual versus neither yielded AUROC 0.756, whereas directional `summary_min` was 0.430; mixed-library ranking further placed 9 of the Top-10 ligands among experimental selectives. AChE/BChE and PIK3CA/PIK3CB showed only small, overlapping increments; PIK3CA/mTOR Dual versus neither is underpowered (neither n = 4) and is not interpreted as a reverse effect. Under scaffold-grouped cross-validation, adding the docking score produced little incremental AUROC beyond ECFP4. These results do not establish a four-pair overestimation law, and they do not prove that docking encodes no pocket-specific information. They show that benchmark formulation, chemotype composition, and receptor realization jointly determine what “dual-target docking success” appears to mean. DualFourClass-Bench is a data-constrained evaluation protocol, not a comprehensive dual-target suite.
