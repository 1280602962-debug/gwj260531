# Results (JCIM Articles draft, English)

## 3. Results

### 3.1 Experimental data supply limits strict dual-target benchmark construction

To determine whether public bioactivity data can support a strict dual-target recognition evaluation, we first audited ligand supply for 49 ChEMBL-cached candidate target pairs (Figure 2). Ligands that meet the activity threshold on one target and are explicitly inactive on the other are directional selective hard negatives.

Under the strict labeling rule (dual: both ends pChEMBL ≥ 6.5; selective: active end ≥ 6.5 and opposite end ≤ 5.5), pairs that simultaneously supply enough A-only and B-only hard negatives were scarce. Only four pairs met the thick-panel gate of ≥50 strict hard negatives on **both** ends. After excluding metal-dependent HDAC1/HDAC6, PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB formed three relatively well-supplied pairs. EGFR/HER2 retained only 7 strict B-only ligands and was therefore kept as a supply-limited case (Table 1).

This supply constraint is not a ChEMBL-only counting artifact. A zero-docking BindingDB / PubChem count check on the same four pairs (Table S12) left the thick-panel gate intact under an equal-relation rule that more closely matches pChEMBL (`equal_only`): min hard-negative counts for the three frozen thick pairs were 76 / 92 / 58 in BindingDB and 86 / 97 / 61 in PubChem (ChEMBL cache: 80 / 78 / 56), all still ≥ 50. EGFR/HER2 rose to about 30 B-end hard negatives in the other databases, enough for a thin (≥ 20) pool but not a thick (≥ 50) panel. Treating censored inequality records as point estimates (`as_is`) would inflate EGFR/HER2 supply (BindingDB min HN = 85), but 49 of 92 as-is B_only ligands have **only** `>` records on EGFR and were not used to freeze the benchmark.

The size of the final benchmark was therefore constrained by the availability of experimentally defined directional hard negatives. The strict 6.5/5.5 rule quantifies supply and records panel construction, whereas θ = 6.0 defines the experimental-state labels for all primary AUROCs (Methods 2.1). Subsequent analyses examine benchmark formulation (Section 3.2), ligand-level chemical baselines (Section 3.3), evaluation-condition sensitivity (Section 3.4), and falsification controls (Section 3.5).

### 3.2 Benchmark formulation changes apparent dual-target discrimination

On the frozen four pairs, AutoDock Vina scores were evaluated under one unified θ = 6.0 label rule using pocket-matched directional AUROC (Figure 1B; Methods 2.4). Scores are \(S=-E_{\mathrm{Vina}}\) (higher better); dual is the positive class. The prespecified worst-direction discrimination summary is `summary_min`, the smaller directional AUROC. Arithmetic, geometric, and harmonic means are aggregation sensitivities only; pair ranking is unchanged under all four summaries (Table S26). For AChE/BChE and PIK3CA/PIK3CB, construction used the stricter 6.5/5.5 rule, but θ = 6.0 gives identical ligand classification and AUROC on this data (Table S4). EGFR/HER2 and PIK3CA/mTOR become underpowered on B_only under the strict rule. Ranking trends held across the full threshold grid (Figure S1A).

EGFR/HER2, AChE/BChE, PIK3CA/PIK3CB, and PIK3CA/mTOR gave directional `summary_min` values of 0.430, 0.606, 0.500, and 0.692, respectively (Table 2; Figure 4A). Dual-versus-B-only AUROC is 0.430 on EGFR/HER2 and 0.500 on PIK3CA/PIK3CB, whereas PIK3CA/mTOR reaches 0.714 and 0.692 on the two directions. Relative to pooling, pocket matching raised point estimates without changing rank order (Table S6).

The same frozen scores were then scored under Dual versus neither and Dual versus all non-duals (Table 3; Figure 3). Dual versus neither uses experimental inactives (`vina_mean`). EGFR/HER2 provides the clearest formulation example. Dual versus neither yielded AUROC 0.756 [0.562, 0.920] (n_neg = 12), whereas directional `summary_min` remained 0.430 [0.284, 0.576]. Dual versus all non-duals collapsed to 0.551 [0.443, 0.666]. In a mixed-library ranking of all 110 EGFR/HER2 ligands by `vina_mean`, the Top-10 contained 1 dual, 5 A-only, 4 B-only, and 0 neither (EF10 = 0.393; hard-negative fraction = 0.90); EF5 was also below random (Table S25). AChE/BChE and PIK3CA/PIK3CB showed only small Dual-versus-neither increments (0.649 and 0.559) whose intervals overlap the directional arms. PIK3CA/mTOR Dual versus neither is underpowered (neither n = 4); Dual versus all non-duals on that pair was 0.674, close to `summary_min` 0.692.

**Table 2.** Pocket-matched directional AUROC on the frozen K = 4 set (Vina; unified θ = 6.0), with all four prespecified descriptor `summary_min` values. `n` gives the scored dual / A_only / B_only class sizes entering the primary AUROCs after requiring both-end scores; neither ligands are excluded from these AUROCs. Constructed panel size and both-end docking coverage are reported in Table 1 and Table S27. The highest descriptor is a best single-descriptor reference. Wrong-pocket and ligand-efficiency controls are in Table S6; full descriptor arms are in Table S28.

| Pair | n (dual / A_only / B_only) | dual vs A_only (pocket B) | dual vs B_only (pocket A) | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.284, 0.576] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.440, 0.740] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.347, 0.648] | 0.622 | 0.620 | 0.595 | 0.418 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.464, 0.802] | 0.463 | 0.448 | 0.310 | 0.260 |

**Table 3.** Same Vina scores under Dual-versus-neither versus directional formulations (unified θ = 6.0). Dual-versus-neither uses experimental inactives (`vina_mean`); Dual versus all non-duals counts A-only, B-only, and neither as negatives. Directional CIs are from Table 2. Negative sets differ. PIK3CA/mTOR Dual versus neither is underpowered (n_neg = 4).

| Pair | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|---|---:|---:|---:|---:|
| EGFR/HER2 | 0.430 [0.284, 0.576] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.440, 0.740] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/PIK3CB | 0.500 [0.347, 0.648] | 0.559 [0.373, 0.746] | 16 | 0.556 [0.437, 0.672] |
| PIK3CA/mTOR | 0.692 [0.464, 0.802] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |

All four 95% bootstrap intervals for `summary_min` included 0.5; thus, at the present sample sizes, no target pair yielded clear evidence excluding chance-level discrimination. PIK3CA/mTOR had the highest point estimate (0.692; 95% CI 0.464–0.802), but its paired difference from the best single-descriptor reference also included 0 (Table S19). AChE/BChE (0.606) remained below TPSA (0.733), while EGFR/HER2 (0.430) and PIK3CA/PIK3CB (0.500) showed no clear advantage over their descriptor references.

Both-end scores were obtained for 110/110 EGFR/HER2 ligands, 95/100 AChE/BChE ligands, 99/100 PIK3CA/PIK3CB ligands, and 48/48 PIK3CA/mTOR ligands (Table S27). One A-only ligand failed PIK3CA docking because of computational timeout and was omitted from analyses requiring that score (Tables S27, S30). AUROCs are therefore conditional on compounds AutoDock Vina can process. Alternative scorers on the same poses did not change the overall ranking (Tables S14–S15; Figure S1B).

### 3.3 Ligand properties and chemotype explain much of the apparent signal

Pocket-matched docking was first compared with four prespecified physicochemical descriptors (Figure 4B; Table 2). Relative to the best single-descriptor reference on each pair, the paired difference in summary_min was −0.052, −0.128, −0.122, and +0.229 for EGFR/HER2, AChE/BChE, PIK3CA/PIK3CB, and PIK3CA/mTOR; all four 95% confidence intervals include 0 (Table S19; Figure S3C). Even the largest positive point difference, on PIK3CA/mTOR, is therefore not distinguishable from the ligand-property reference with the present sample.

AChE/BChE is a direct confounding case. Mean TPSA was ≈ 75 for dual ligands versus ≈ 51 for selective hard negatives (Figure 4C). TPSA alone gave AUROC ≈ 0.769, above Vina under the same contrast (≈ 0.56). Adding heavy-atom count and TPSA raised dual-versus-B-only AUROC from 0.606 to 0.807, while the docking-score odds ratio was only ≈ 1.18 (Figure 7C). Thus, part of the apparent docking discrimination on this arm can be accounted for by ligand physicochemical information.

PIK3CA/mTOR differs in degree. Adding heavy-atom count and TPSA shifted AUROC by about +0.07 to +0.11, with docking odds ratios ≈ 2.19 and 3.08. The paired difference versus the descriptor baseline still includes 0. After ligand-efficiency normalization, only PIK3CA/mTOR remained above the heavy-atom baseline (0.657 versus 0.463).

A two-dimensional chemical baseline makes the same point (Figure 7A). ECFP4 logistic regression under Bemis–Murcko scaffold GroupKFold yielded fold AUROCs of about 0.78–0.91 on several arms, well above the corresponding docking contrasts — for example 0.85 versus 0.43 for EGFR/HER2 dual-versus-B-only. The experimental labels therefore contain ligand-structure-associated information that can be exploited without receptor information. On PIK3CA/mTOR, \(n_{\mathrm{scaffolds}} \approx n\), so the split is nearly leave-one-scaffold. A random `StratifiedKFold` check sits on average +0.011 above the scaffold split across eight directional contrasts (Table S20; Figure S3D).

Under the present scaffold-grouped benchmark, the largest absolute change after adding the pocket-matched docking score to ECFP4 was 0.020 (unrounded −0.0198 on PIK3CA/mTOR dual versus A-only), and several changes were negative (Table S24). At T ≥ 0.3, PIK3CA/PIK3CB dual versus A-only fell from 0.691 to 0.503 (n_neg = 11), whereas distant hard negatives (T < 0.3) yielded 0.819; T ≥ 0.4/0.5 cells were often n_neg ≤ 7, and T ≥ 0.7 was empty (Table S23). On potency- or size-matched subsets, dual-versus-B-only remained weak or near chance on EGFR/HER2 and PIK3CA/PIK3CB (about 0.45–0.52), with per-arm n often < 15 (Table S5; Figure 7D).

### 3.4 Activity aggregation, ligand panels, and receptor realization affect evaluation outcomes

Primary labels use the maximum available pChEMBL value. Assay-level records were re-fetched and re-aggregated with both maximum and median rules, without changing panel membership, docking parameters, or Vina scores. Median aggregation changed 7/110, 1/95, 1/99, and 0/48 ligand-state assignments (label agreements 93.6%, 98.9%, 99.0%, and 100%). On the API-refetched labels, pair-level `summary_min` moved only modestly (0.417→0.424, 0.606→0.629, 0.500→0.500, and 0.692→0.692). These API-refetched estimates are a max-versus-median sensitivity reported alongside Table 2 (Table S29). Assay-level heterogeneity remains because pChEMBL values are not assay-equivalent.

Lowering exhaustiveness from 16 to 8 moved PIK3CA/mTOR summary_min from 0.692 to 0.660 (Figure S1D). On PM110 (analysis n = 115; 30 / 30 / 30 dual / A_only / B_only), Vina summary_min was 0.648 [0.51, 0.76], about 0.04 below PM48, with the same ranking trend (Figure S1C). On the unused-pool holdout (20 / 20 / 20 per pair; seed 20260731; EGFR/HER2 not eligible), PIK3CA/mTOR summary_min was 0.765 [0.603, 0.891], AChE/BChE was 0.618 [0.422, 0.759], and PIK3CA/PIK3CB fell to 0.425 [0.241, 0.618] (Tables S8, S16). The holdout shares the same ChEMBL extraction batch.

We next tested receptor realization, holding one pocket frozen and replacing the other (Figure 5; Tables S9, S30). Three alternate crystals passed cognate redocking QC (best-of-9 RMSD 0.607 Å for 4JPS, 0.624 Å for 5DXT, and 0.515 Å for 4JSX). On PIK3CA/mTOR, replacing PIK3CA 4L23 with 4JPS or 5DXT, while holding mTOR at 4JT6, dropped PM48 summary_min from 0.692 to 0.486 [0.259, 0.692] and 0.505 [0.292, 0.696] (Figure 5A). The change concentrated in the D/B direction that depends on the alternate PIK3CA structure; D/A stayed at 0.714. Replacing mTOR 4JT6 with 4JSX gave summary_min 0.639 [0.418, 0.776].

The same PIK3CA crystals on the PIK3CA/PIK3CB panel, with 2WXF frozen, **raised** summary_min from 0.500 to 0.691 [0.516, 0.779] (4JPS) and 0.685 [0.506, 0.768] (5DXT) (Figure 5B). Dual versus A-only stayed at 0.691; dual versus B-only rose from 0.500 to 0.707 and 0.685. All three PIK3CA conditions used the same 99 scored ligands (Table S30). Thus, replacing the same PIK3CA receptor altered the apparent discrimination in opposite directions across the two target pairs, demonstrating receptor-realization sensitivity rather than a uniform loss of docking performance.

### 3.5 Wrong-pocket controls reveal an unresolved out-of-panel failure

On the main panels, pocket-matched summary_min exceeded the wrong-pocket control on all four pairs; matched-minus-wrong differences were 0.170, 0.161, 0.151, and 0.090. The EGFR/HER2 and AChE/BChE intervals exclude 0; the PIK3CA/PIK3CB and PIK3CA/mTOR intervals include 0 (Tables S6, S17; Figure 6A; Figure S3A). Wrong-pocket summary_min values were 0.260, 0.444, 0.349, and 0.602.

The point-estimate relationship reversed on the unused-pool holdout (Figure 6B). Wrong-pocket summary_min was 0.788, 0.643, and 0.520 for PIK3CA/mTOR, AChE/BChE, and PIK3CA/PIK3CB, versus matched-pocket 0.765, 0.618, and 0.425. All three matched-minus-wrong point differences were negative (−0.023 / −0.025 / −0.095), and every 95% confidence interval included 0 (Table S17; Figure S3B). After potency or size matching, wrong-pocket remained ≥ matched-pocket (Table S13). Scoring-independent contact_count reached AUROC 0.698–0.714 on the B direction, but this coarse surrogate cannot account for the magnitude of the observed Vina wrong-pocket discrimination (Figure 6D; Table S11). The holdout wrong-pocket reversal is therefore an unresolved failure mode exposed by the benchmark.

### 3.6 Exploratory structural context

Local pocket Cα RMSD between 5DXT and 4L23 is only 0.343 Å, yet PIK3CA/mTOR summary_min still fell to 0.505, so simple backbone similarity does not predict benchmark transferability (Table S10). Whole-chain identity and representative pose diagnostics are retained in the Supporting Information (Table S7; Note S1).
