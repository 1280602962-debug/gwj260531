# Results (JCIM Articles draft, English)

> Companion to [`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) (Chinese authoritative for this rewrite cycle).  
> Numbers from `PRIMARY_METRIC_V2.md`, `PRIMARY_METRIC_CLAIM_GATE.md`, SI Tables S4–S6. No fabricated experiments.

## 3. Results

### 3.1 Construction of a dual-target recognition benchmark: public-data limits on hard-negative supply

Starting from 49 audited ChEMBL target pairs, only four pairs retained ≥50 experimentally characterized hard-negative selective ligands on **both** ends under the primary strict rule (A_only: A ≥ 6.5 and B ≤ 5.5; B_only symmetrically). This is the first major finding of the study: despite the large number of ChEMBL target pairs, **balanced dual-target benchmarking was severely constrained by the scarcity of experimentally characterized hard-negative selective ligands**. The frozen evaluation set therefore comprises four pairs (Table 1): PIK3CA/mTOR (PM48), AChE/BChE, PIK3CB/mTOR, and EGFR/HER2 as a supply-limited case (B_only = 8). Cognate redocking and pocket-matched scoring definitions follow Methods 2.4–2.5.

### 3.2 Cross-pair robustness under a unified label rule

To address threshold cherry-picking, all four panels were relabeled under a single rule with θ ∈ {5.5, 6.0, 6.5} and the strict 6.5/5.5 rule, and pocket-matched summary_min was recomputed on existing Vina scores (Table S4). This is the primary cross-pair robustness analysis. Under the strict rule, PIK3CA/mTOR remained the highest among thick panels (summary_min = 0.639); EGFR/HER2 remained lowest (0.324) with underpowered B_only. Construction-time labels in Table 2 are a construction readout and must not be mixed with Table S4 as competing primary standards.

### 3.3 Docking shows limited, pair-dependent ability to discriminate true dual-target ligands

Pooled dual-end scores (min/mean/max of the two pocket scores) collapsed directional information and often approached chance (e.g., EGFR/HER2 pooled AUROC near 0.50). The primary metric is therefore pocket-matched directional AUROC (Methods 2.5). On EGFR/HER2, the weaker pocket-matched arm fell to about 0.28; this is a **direction-specific discrimination failure**, not an artifact produced by pooling. Table 2 reports only the two directional AUROCs, summary_min, and 95% CIs. Wrong-pocket, ligand-efficiency, and descriptor baselines are deferred to Supporting Information (Tables S5–S6 and related notes).

**Table 2.** Pocket-matched directional AUROC (Vina; construction-time labels).

| Pair | AUROC A (dual vs A_only) | AUROC B (dual vs B_only) | summary_min [95% CI] |
|---|---:|---:|---:|
| PIK3CA/mTOR | 0.692 | 0.690 | 0.690 [0.52, 0.84] |
| AChE/BChE | 0.553 | 0.606 | 0.553 [0.45, 0.66] |
| PIK3CB/mTOR | 0.439 | 0.674 | 0.439 [0.32, 0.56] |
| EGFR/HER2 | 0.281 | 0.545 | 0.281 [0.15, 0.43] |

Only PIK3CA/mTOR had both directional point estimates above 0.5; its summary_min CI still spans 0.5. Worst-pocket aggregation (Table S6) preserved the same ranking: PIK3CA/mTOR highest (0.627), EGFR/HER2 lowest (0.271).

### 3.4 Physicochemical and structural confounding dominates several apparent dual-target signals

The central scientific finding is not that “docking is bad,” but that **many apparent dual-target signals are explained by ligand properties**. Heavy-atom count and TPSA often matched or exceeded Vina under the same four-class contrasts. Scaffold-grouped ECFP4 baselines frequently outperformed docking (e.g., AChE/BChE summary_min ≈ 0.78; EGFR/HER2 ≈ 0.72; PIK3CA/mTOR ≈ 0.65), indicating that dual-target recognition is largely chemotype-determined in these panels. Wrong-pocket controls and matched-subset analyses are consistent with this interpretation (Supporting Information).

On AChE/BChE, mean TPSA was ≈ 75 (dual) versus ≈ 51 (hard-negative selective ligands). TPSA alone separated dual from hard negatives with AUROC ≈ 0.769, above Vina under the same contrast (≈ 0.56). Adding heavy-atom count and TPSA raised dual-versus-B AUROC from 0.606 to 0.807. That increase indicates that **the apparent docking contribution was largely dependent on physicochemical covariates**; the docking odds ratio near 1 (OR ≈ 1.18) should not be read as retained independent directional information. On PIK3CA/mTOR, AUROC shifts after covariates were smaller (≈ +0.07 to +0.11) with OR ≈ 2.19 and 3.08, suggesting only residual pocket-specific signal, to be read together with Δ intervals that still include 0.

### 3.5 Robustness checks and case-dependent success

Raising Vina exhaustiveness from 8 to 16 left ranking trends consistent (ranking trend remained consistent). Single-target enrichment against property-matched weak binders (pChEMBL ≤ 5.5) gave AUROC 0.603/0.629 and EF1% 2.04/2.00 on 4L23/4JT6—docking retained limited enrichment capability, not a strong VS engine. Expanding PIK3CA/mTOR to PM110 is a **stability check**, not independent validation: summary_min moved from 0.690 to 0.648 [0.51, 0.76] with the same ranking direction.

Across §3.2–3.5, **only PIK3CA/mTOR showed reproducible but modest pocket-related discrimination**; apparent signals on the other three pairs are largely explained by ligand properties or 2D chemotype.

### 3.6 Structural clues for the only reproducible dual-target signal (case-level)

PIK3CA and mTOR are ATP-competitive kinase-related pockets. Existing failure typology on this pair shows hinge/ATP-like duals that rank well on both ends (type T2) versus large fused aromatics that fail on one end (type T5). We do not claim a completed shared-residue / PLIF campaign here; a dedicated structural-determinants analysis (pocket RMSD, conserved hinge contacts, interaction fingerprints) is the natural next step to explain why this pair alone retains residual signal.
