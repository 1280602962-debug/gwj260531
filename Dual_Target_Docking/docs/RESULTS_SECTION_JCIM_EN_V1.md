# Results (JCIM Articles draft, English)

> Companion to [`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) (Chinese authoritative for this rewrite cycle).  
> Numbers from `PRIMARY_METRIC_V2.md`, `PRIMARY_METRIC_CLAIM_GATE.md`, SI Tables S4–S6. No fabricated experiments.

## 3. Results

### 3.1 Construction of a dual-target recognition benchmark: public-data limits on hard-negative supply

Dual-target docking evaluation requires four ligand classes: dual, A-selective, B-selective, and neither. Experimentally defined selective ligands on each arm serve as hard-negative selective ligands for testing whether a score can suppress both single-target arms.

Across 49 audited ChEMBL target pairs under the primary strict rule (dual: both ends ≥ 6.5; selective: active end ≥ 6.5 and opposite ≤ 5.5), only four pairs retained ≥50 hard-negative selective ligands on **both** ends. Despite the large number of ChEMBL target pairs, **balanced dual-target benchmarking was severely constrained by the scarcity of experimentally characterized hard-negative selective ligands**. After excluding metal-dependent HDAC1/HDAC6, three pairs supported reasonably balanced strict panels; EGFR/HER2 entered as a supply-limited case (few strict B-selective ligands), not as a thick panel. The frozen K = 4 set follows this audit, not post-hoc selection of docking-favorable pairs (Methods 2.1–2.3).

### 3.2 Cross-pair robustness under a unified label rule

Construction used the strict rule for AChE/BChE and PIK3CA/PIK3CB, and θ = 6.0 for EGFR/HER2 and PIK3CA/mTOR (PM48) when strict selectives were scarce (Table 1). To address threshold cherry-picking, all four panels were relabeled under a single rule with θ ∈ {5.5, 6.0, 6.5} and the strict 6.5/5.5 rule, and pocket-matched summary_min was recomputed on existing Vina scores (Table S4). This is the primary cross-pair robustness analysis; construction-time estimates are a construction readout only.

Under the strict rule, AChE/BChE and PIK3CA/PIK3CB remained at 0.606 and 0.500; PIK3CA/mTOR was 0.639 (slightly below the construction θ = 6.0 value of 0.692); EGFR/HER2 was 0.324. EGFR/HER2 and PIK3CA/mTOR had only 7 and 4 B-selective ligands under the strict rule and are marked underpowered. Ranking trends remained consistent: PIK3CA/mTOR highest; the other three pairs ≤ 0.61.

### 3.3 Docking shows limited, pair-dependent ability to discriminate true dual-target ligands

A dual-target score must suppress both selective arms. Pooling the two pocket scores (e.g., mean) can let the stronger arm mask failure on the weaker arm. The primary metric is therefore pocket-matched directional AUROC: dual versus A_only scored in pocket B, dual versus B_only scored in pocket A, with summary_min as the smaller arm (Methods 2.6). Scores are \(S=-E_{\mathrm{Vina}}\) (higher better); dual is the positive class.

Two points must be kept distinct. First, **direction-specific discrimination failure**: on EGFR/HER2, pocket-matched dual-versus-B_only AUROC is 0.430, and the weaker arm under a pooled protocol can read near ~0.28; that low value reflects failure of that direction itself, not a pooling arithmetic that “creates” 0.28 from 0.50. Second, **pooling can mask the weak arm**: on the same pair, a pooled summary can approach ~0.50 and appear merely mediocre. Relative to pooling, pocket matching raised point estimates across pairs without changing rank order (Table S6).

**Table 2.** Pocket-matched directional AUROC on the frozen K = 4 set (Vina; construction-time labels). Wrong-pocket, ligand-efficiency, and descriptor baselines are in Supporting Information Table S6.

| Pair | n (dual / A_only / B_only) | dual vs A_only (pocket B) | dual vs B_only (pocket A) | summary_min [95% CI] |
|---|---:|---:|---:|---|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.281, 0.576] |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.442, 0.737] |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.340, 0.648] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.457, 0.813] |

Docking signal is generally weak and highly target-pair dependent. Only PIK3CA/mTOR has summary_min point estimate above both 0.5 and the best trivial descriptor baseline (heavy-atom count 0.463), yet its 95% CI lower bound still approaches 0.5. EGFR/HER2 and PIK3CA/PIK3CB do not beat their best descriptor baselines; AChE/BChE (0.606) remains below the TPSA baseline (0.733). RTMScore and GNINA did not change this pattern.

### 3.4 Physicochemical and structural confounding dominates several apparent dual-target signals

The central finding is not that “docking is bad,” but that **many apparent dual-target signals are explained by ligand properties**. Relative to the best trivial descriptor, docking–descriptor Δ intervals lie below 0 for EGFR/HER2 and PIK3CA/PIK3CB; AChE/BChE does not clear the gate; PIK3CA/mTOR exceeds heavy-atom count at the point estimate, but the Δ 95% CI still includes 0 (Table S6). Wrong-pocket summary_min values were 0.260, 0.444, 0.349, and 0.602; pocket-matched gains over wrong-pocket exceeded 0.09 on all pairs, underscoring ligand-level confounding.

On AChE/BChE, mean TPSA was ≈ 75 (dual) versus ≈ 51 (hard-negative selective ligands). TPSA alone separated dual from hard negatives with AUROC ≈ 0.769, above Vina under the same contrast (≈ 0.56). Adding heavy-atom count and TPSA raised dual-versus-B AUROC from 0.606 to 0.807. That increase indicates that **the apparent docking contribution was largely dependent on physicochemical covariates**; the docking odds ratio near 1 (OR ≈ 1.18) should not be read as retained independent directional information. On PIK3CA/mTOR, covariate-adjusted AUROC shifts were smaller (≈ +0.07 to +0.11) with OR ≈ 2.19 and 3.08, suggesting only residual pocket-specific signal, to be read with Δ intervals that still include 0.

Scaffold-grouped ECFP4 baselines further elevate the confounding narrative: fold AUROCs often 0.78–0.91 and frequently above the corresponding docking arms (e.g., EGFR/HER2 dual-versus-B_only: fingerprint 0.85 vs docking 0.43). Labels correlate with chemotype; this alone does not prove pocket-physical specificity of docking scores.

### 3.5 Robustness checks and case-dependent success

Changing PIK3CA/mTOR exhaustiveness from 16 to 8 moved Vina summary_min from 0.692 to 0.660 (Δ ≈ +0.03)—far smaller than between-pair differences. Single-target enrichment against property-matched weak binders (pChEMBL ≤ 5.5) gave AUROC 0.603/0.629 and EF1% 2.04/2.00 on 4L23/4JT6: docking retained limited enrichment capability, not a strong VS engine.

Expanding PIK3CA/mTOR to PM110 is a **stability check**, not independent validation and not an attempt to “rescue” the estimate with a larger panel. PM48 itself is small (18/14/12). On PM110, Vina summary_min was 0.648 [0.51, 0.76] versus 0.692 on PM48 (Δ ≈ −0.04); ranking trend remained consistent.

Across §3.2–3.5, **only PIK3CA/mTOR showed reproducible but modest pocket-related discrimination**; apparent signals on the other three pairs are largely explained by ligand properties or 2D chemotype.

### 3.6 Structural clues for the only reproducible dual-target signal (case-level)

PIK3CA (4L23) and mTOR (4JT6) are ATP-competitive kinase-related pockets; cognate ligands recover near-native poses under protocol checks (Table S3). Existing pose-level failure typology shows hinge/ATP-like duals that can rank well on both ends (type T2) versus cases where rescoring prefers off-hinge poses (type T5). Even on this best pair, ATP-site cross-chemotypes can be misread as dual when both pockets yield geometrically clean, hinge-positive poses. We do not claim a completed shared-residue / PLIF campaign; a dedicated structural-determinants analysis remains future work. The PIK3CA/mTOR advantage should be read as limited directional signal for some chemotypes under a shared ATP recognition framework, not as a validated general dual-target decision rule.
