# Results (JCIM English manuscript draft)

> DualFourClass-Bench · Evaluation / benchmark Article  
> Venue target: *J. Chem. Inf. Model.*  
> Skills applied from GitHub (read-only; not vendored into this repo):  
> - [Boom5426/Nature-Paper-Skills](https://github.com/Boom5426/Nature-Paper-Skills): `results-section-revision`, `scientific-prose-style`, `scientific-writing`  
> - [Yuan1z0825/nature-skills](https://github.com/Yuan1z0825/nature-skills): `nature-polishing` (paper_type=methods/evaluation, section=results, journal=generic)  
> - [stephenturner/skill-deslop](https://github.com/stephenturner/skill-deslop): phrase/structure deslop checklist  
> Numbers: `jcim_j0j1_v0`, `jcim_bench_v0/tables/pocket_matched_*`, `jcim_strengthen_t0t1_v0`.  
> Claim ceiling: evaluation article; primary metric = pocket-matched directional AUROC; no universal scorer claim.

---

## 3. Results

### 3.1 Strict hard-negative supply on public target pairs

Dual-target docking evaluation needs ligands in four activity classes: dual actives, A-only hard negatives, B-only hard negatives, and inactive/neither. We audited 49 ChEMBL-backed target pairs under a strict label rule (dual: both pChEMBL ≥ 6.5; A_only or B_only: active end ≥ 6.5 and opposite end ≤ 5.5). Only four pairs supplied at least 50 strict hard negatives on both arms. After excluding the metal-dependent HDAC1/HDAC6 pair, three pairs remained usable for conventional thick docking panels: PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB (Table 1). As a contrast, the literature-common EGFR/HER2 pair yielded only seven strict B_only ligands under the same rule and did not meet the thick four-class panel threshold. The frozen K=4 evaluation set was fixed from this audit (construction rules in Methods).

### 3.2 Pooled scores versus pocket-matched directional AUROC

A dual-target score must suppress both single-target hard-negative arms. When both contrasts reuse one pooled score such as `vina_mean`, a strong arm can mask a weak arm. On EGFR/HER2, the pooled summary AUROC was near 0.50, whereas the weaker directional arm (dual versus B_only) fell to about 0.28 (Figure 1). We therefore report pocket-matched directional AUROC as the primary metric: dual versus A_only uses the pocket-B score, dual versus B_only uses the pocket-A score, and `summary_min` is the minimum of the two arm AUROCs. Pooled scores are retained only as controls. Pocket matching raised point estimates on all four pairs relative to pooling, but the ranking did not change: only PIK3CA/mTOR stayed clearly above chance, and the other three pairs remained at `summary_min` ≤ 0.61 (Table 2). As an aggregation control, worst-pocket `summary_min` values were 0.271 (EGFR/HER2), 0.579 (AChE/BChE), 0.439 (PIK3CA/PIK3CB), and 0.627 (PIK3CA/mTOR); the full aggregation contrast is in Supporting Information Table S6.

### 3.3 Directional discrimination on the frozen K=4 set

We scored each panel with AutoDock Vina, RTMScore best-of-K, and GNINA CNN (mode_01 rescore) under a common RDKit/meeko ligand protocol. Vina pocket-matched `summary_min` is the primary report; RTM and GNINA are channel controls. Bootstrap 95% confidence intervals used B = 2000 resamples with seed 20260729 (Table 2; Figure 2).

**Table 2. Pocket-matched directional AUROC (Vina) on the K=4 set.**

| Target pair | n (D / A / B) | D vs A (pocket B) | D vs B (pocket A) | summary_min [95% CI] | Wrong-pocket min | LE-PM min | Best trivial baseline |
|---|---:|---:|---:|---|---:|---:|---|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.281, 0.576] | 0.260 | 0.311 | cLogP 0.482 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.442, 0.737] | 0.444 | 0.413 | TPSA 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.340, 0.648] | 0.349 | 0.332 | heavy atoms 0.622 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.457, 0.813] | 0.602 | 0.657 | heavy atoms 0.463 |

PIK3CA/mTOR was the only pair whose Vina `summary_min` point estimate exceeded both 0.5 and the heavy-atom baseline (0.463; Δ ≈ +0.23). The lower CI bound remained near 0.5. EGFR/HER2 and PIK3CA/PIK3CB gave `summary_min` of 0.430 and 0.500, both below their best trivial baselines. AChE/BChE docking (0.606) was below TPSA (0.733). RTM and GNINA did not reverse this pattern: both were weaker than Vina on PIK3CA/mTOR and also failed the baseline gate on the negative pairs.

Rank readout agreed with the weak-arm AUROC. Among the ten ligands with the highest Vina pooled scores on EGFR/HER2, nine were hard negatives (bootstrap mean ≈ 8.9; CI ≈ 7–10).

### 3.4 Baseline gate against trivial descriptors

For each pair we subtracted the best trivial descriptor `summary_min` (heavy atoms, molecular weight, cLogP, or TPSA) from the docking `summary_min`. Under the pooled-score bootstrap gate, the Δ confidence intervals for EGFR/HER2 and PIK3CA/PIK3CB lay entirely below zero. AChE/BChE failed on the point estimate, with a Δ interval that just crossed zero. PIK3CA/mTOR beat heavy atoms on the point estimate, but the Δ 95% CI still included zero (Figure 3). After pocket matching, the same qualitative gaps remained: docking trailed cLogP on EGFR/HER2, TPSA on AChE/BChE, and heavy atoms on PIK3CA/PIK3CB. These gates limit claims to evaluation diagnostics rather than a general dual-target decision score.

### 3.5 Wrong-pocket, ligand-efficiency, polarity, and matched-subset controls

If docking scores mainly tracked ligand properties, wrong-pocket controls should depart from 0.5. Wrong-pocket `summary_min` values were 0.260 (EGFR/HER2), 0.444 (AChE/BChE), 0.349 (PIK3CA/PIK3CB), and 0.602 (PIK3CA/mTOR). The gap between pocket-matched and wrong-pocket minima exceeded 0.09 for every pair (Table 2). After ligand-efficiency normalization, only PIK3CA/mTOR remained above the heavy-atom baseline (0.657 versus 0.463). The other three pairs lost directional support under that normalization.

On AChE/BChE, mean TPSA was ≈ 75 for dual actives and ≈ 51 for hard negatives. TPSA alone separated dual from hard-negative ligands with AUROC ≈ 0.769, above Vina under the same contrast (≈ 0.56). Adding heavy-atom count and TPSA as covariates raised the pocket-matched dual-versus-B AUROC from 0.606 to 0.807 (Δ ≈ +0.20), with an odds ratio of 1.18 for the docking score. On PIK3CA/mTOR the corresponding AUROC shifts were smaller (≈ +0.07 to +0.11), with odds ratios of 2.19 and 3.08 for dual-versus-A and dual-versus-B, indicating residual directional information after size and polarity adjustment.

Potency-matched (|ΔpChEMBL| ≤ 0.5) and size-matched (|Δheavy atoms| ≤ 2) subsets still gave weak or near-chance dual-versus-B AUROCs on EGFR/HER2 and PIK3CA/PIK3CB (≈ 0.45–0.52). Directional signal on PIK3CA/mTOR persisted in matched subsets, but arm sizes were often below 15 and intervals were wide; those strata are reported in the Supporting Information.

### 3.6 Scaffold-grouped ligand fingerprint baseline

As a ligand-only baseline we trained ECFP4 logistic regression with Murcko-scaffold GroupKFold so that no scaffold spanned training and test folds. Scaffold-fold AUROCs typically fell between 0.78 and 0.91 and exceeded the matching pocket-matched docking arms. On EGFR/HER2 dual versus B_only, the fingerprint model reached 0.85 while docking reached 0.43. Random ligand folds averaged only about 0.01 higher than scaffold folds (maximum ≈ 0.10), consistent with many near-singleton scaffolds and limited leakage control. We therefore report scaffold-fold numbers in the main text and treat random-fold scores as a leakage diagnostic only. The fingerprint baseline shows that labels track chemotype; it does not by itself prove pocket-specific physics.

### 3.7 Exhaustiveness, single-target enrichment, and PM110 expansion

Repeating PIK3CA/mTOR docking at exhaustiveness 8 instead of 16 changed Vina `summary_min` from 0.692 to 0.660 (Δ = +0.032 for E16 versus E8). That gap is too small to explain why this pair outranks the others, so both settings are reported.

Single-target enrichment used property-matched ChEMBL weak binders (pChEMBL ≤ 5.5) as decoys rather than random unrelated molecules. Enrichment AUROC was 0.603 on 4L23 (PIK3CA) and 0.629 on 4JT6 (mTOR), with EF1% of 2.04 and 2.00 and EF5% of 1.22 and 3.20. Both ends showed weak but non-trivial enrichment, consistent with docking retaining some single-target signal without acting as a strong virtual-screening engine.

PM110 expands PM48 by keeping all 48 ligands and adding quota-sampled molecules; it is not an independent replicate. Vina pocket-matched `summary_min` on PM110 was 0.648 [0.51, 0.76], compared with 0.692 on PM48 (Δ ≈ −0.04), with a narrower interval and the same direction. RTM and GNINA on PM110 gave 0.576 and 0.522.

### 3.8 Robustness under a unified label rule

To address the concern that the four panels were built under different label rules, we relabeled each existing panel under θ ∈ {5.5, 6.0, 6.5} and under the strict 6.5/5.5 rule, then recomputed pocket-matched `summary_min` on the same Vina scores (Supporting Information Table S4). Under the strict rule, AChE/BChE and PIK3CA/PIK3CB remained at 0.606 and 0.500; PIK3CA/mTOR was 0.639 (versus 0.692 at θ = 6.0 on the main table); EGFR/HER2 was 0.324 (versus 0.430). EGFR/HER2 and PIK3CA/mTOR were underpowered on the strict rule (7 and 4 B_only ligands). The ranking was unchanged: PIK3CA/mTOR remained highest, and the other three pairs stayed ≤ 0.61.

---

## Revision notes

- Applied `results-section-revision`: topic titles, first-paragraph bridges, one job per paragraph, consequence sentences that set up the next subsection.
- Applied `scientific-prose-style` / `nature-polishing` stance: no em dashes, past-tense observation language, hedge once, openers varied, sentences mostly ≤ 30 words.
- Applied `deslop`: removed meta-advice about what not to claim, throat-clearing, and editorial defense; left flat negative reporting.
- Adapted methods-paper Results job to an evaluation article: fair baselines and failure anatomy, not method-win framing.
- Claim ceiling unchanged: PIK3CA/mTOR is an exploratory positive control; Δ versus trivial baseline still does not exclude zero.
