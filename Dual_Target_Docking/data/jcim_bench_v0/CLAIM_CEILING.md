# CLAIM_CEILING — DualFourClass-Bench (JCIM)

## Allowed claims
1. This is an **evaluation / benchmark** article for dual-target docking decision metrics.
2. Primary metric = **pocket-matched directional AUROC** (D vs A_only uses pocket-B score; D vs B_only uses pocket-A score) with a **prespecified four-descriptor panel** (heavy atoms, MW, cLogP, TPSA) reported in full. The strongest descriptor is a **descriptive** baseline, not a confirmatory “best-of-four” hypothesis test. Pooled/wrong-pocket controls are reported in parallel. DualFourClass-Bench is a **four-state** curated resource; the primary endpoint is **two directional pairwise** tasks, not a four-class classifier.
3. K=4 pairs are a **frozen evaluation set**, not a claim that the metric generalizes to all target pairs.
4. EGFR/HER2 is a **supply-limited case study** (existing unified RDKit EH110); no claim from new EGFR docking. A later BindingDB/PubChem **count-level** check (Table S12) does not rebuild this panel: under equal-relation measurements EGFR still fails the ≥50 thick-panel gate (min HN ≈ 30).
5. Prep protocol is frozen: **RDKit ETKDG + meeko**. Do **not** mention Schrodinger LigPrep in the manuscript (no formal license; early borrow was internal-only).
6. Contribution language may say we **established a systematic benchmarking framework / evaluation protocol** and released **DualFourClass-Bench**. That means an evaluation system, **not** a new docking algorithm.

## Forbidden / over-claim
1. Do **not** claim a universal “decision arm” or that `rtm_min_z` (or any pooled score) is validated as a general master score.
2. Do **not** claim score-function invariance beyond the three channels actually run (Vina / RTM / GNINA mode_01 CNN).
3. Do **not** present Track B method competition results as the paper’s core.
4. Do **not** fold shortfall/clash flags or architecture tags into the primary score.
5. Do **not** mix LigPrep and RDKit poses in primary tables.
6. Do **not** invent a grand method acronym (e.g. **D-DRAF**, “Dual-target Docking Reliability Assessment Framework”) or write “we developed a **novel** dual-target docking framework named …”.
7. Do **not** sell ordinary analysis steps (data audit, hard-negative construction, AUROC, confounder baselines, robustness) as numbered **Framework Step 1–5** in Intro/Abstract/TOC as if they were a new algorithm. Methods section numbering is fine; packaging them as a named invention is not.
8. Do **not** use absolute titles of the form “Docking can/cannot identify dual-target ligands.” Prefer *evaluating the reliability and limitations of docking-based dual-target recognition*.
9. Do **not** write that public-data hard-neg supply is “ChEMBL-invariant” without Table S12. BindingDB/PubChem `as_is` counts can pass ≥50 on EGFR/HER2 because of `>` censored values; the matched `equal_only` rule does **not** recast EGFR as a thick panel.
10. Do **not** claim a BindingDB- or PubChem-derived docking panel that was not actually docked.
11. Do **not** write that docking the leftover unused ChEMBL pool would yield **1000 independent non-overlapping** balanced panels. Strict hard-neg leftover after main+holdout is 37/39 (PM), 141/30 (AChE), 8/19 (PIK3CB), 22/0 (EGFR). Holdout is one unused-pool draw, not a 1000-panel distribution. Do not relabel the existing ligand bootstrap as unused-pool resampling.
12. Do **not** promote the 27-ligand max-vs-median API diagnostic (`max_vs_median_diagnostic_sample_v1.csv`) to a completed SI sensitivity table.
13. Do **not** call the AChE/BChE or PIK3CA/PIK3CB ChEMBL-id prefix filter a chemical-diversity constraint.
14. Do **not** call unused-pool holdout “external validation.”
15. Do **not** present wrong-pocket as a positive control that proves pocket specificity; it is a **falsification** control.
16. Do **not** label receptor replacement as “structure robustness that confirms stability”; it is **receptor-structure sensitivity**.
17. Do **not** call `summary_min` a novel scoring function; it is a worst-arm aggregation of two AUROCs.
18. Do **not** treat max(heavy, MW, cLogP, TPSA) as a prespecified confirmatory baseline.
19. Do **not** write that conventional Dual-vs-neither evaluation **systematically overestimates** dual-target docking on all K = 4 pairs. EGFR/HER2 is the pair where Dual-vs-neither (0.756) and directional summary_min (0.430) diverge; AChE and PIK3CA/PIK3CB show small overlapping increments; PIK3CA/mTOR Dual-vs-neither is underpowered (neither n = 4) and must not be used as a reverse-overestimation story.
20. Do **not** promote Tanimoto ≥ 0.7 chemotype-matched AUROCs (the matched sets are empty). Modest T ≥ 0.3 drops are allowed with n_neg reported.
21. Do **not** treat logistic docking AUROC as Table 2 rank AUROC. Incremental ECFP+docking ≈ ECFP is allowed; “docking is information-free on every dual-target pair” is not.
22. Do **not** promote the 27-ligand max-vs-median diagnostic, or a failed live ChEMBL activity fetch, to a completed SI table. Full-panel median relabel is still outstanding (`scripts/assay_aggregation_max_vs_median_v1.py`).

## Conclusions claim ceiling

Allowed closing claims:
- DualFourClass-Bench is an **experimentally grounded evaluation setting** (once in Conclusions; not a new algorithm).
- Discrimination was **limited and strongly target-pair-dependent** (summary_min 0.430–0.692).
- On EGFR/HER2, a conventional Dual-versus-neither readout looked substantially stronger than the directional worst arm; this is **pair-dependent**, not a four-pair overestimation law.
- PIK3CA/mTOR had the strongest point estimate and a positive directional signal in an unused ligand pool; uncertainty + receptor sensitivity **preclude a generalizable dual-target decision rule**.
- Apparent signals can be substantially influenced by ligand properties, chemotype, and receptor realization.
- Unused-pool wrong-pocket reversal is **unresolved**; paired CIs included zero.
- Contribution = **reliability / evidentiary boundary**, not a docking winner.

Forbidden in Conclusions:
- “Docking is ineffective for dual-target recognition.”
- “Docking can identify dual-target ligands.”
- *validated* / *robust performance* / *generalizable dual-target docking strategy*.
- PDB IDs, RTM/GNINA numbers, 4JPS/5DXT/4JSX AUROCs.
- Repeating *experimentally grounded* more than once.

## Mechanism-analysis claim ceiling (Results 3.4–3.6)

Allowed:
- PIK3CA inter-crystal global Cα RMSD (1.44–1.49 Å) is larger than mTOR (0.45 Å) **in this structure set**, consistent in **direction** with greater PIK3CA-end sensitivity.
- Local pocket Cα conservation on 5DXT (0.343 Å) with collapsed summary_min (0.505) means Cα pocket conservation is **not sufficient**.
- Scoring-free `contact_count` shows a real size/burial confound, **especially on the B arm** (0.698–0.714).
- Holdout wrong-pocket ≥ pocket-matched **survives** potency matching (|Δp|≤0.5) and size matching (|Δheavy|≤2) on all three holdout pairs (Table S13). Sampling shift exists (PM holdout weaker than the main panel) but is not a sufficient explanation.

Forbidden:
- Do **not** write that Cα RMSD “quantitatively explains” or causes the AUROC collapse (n = 2 / n = 1; 5DXT matched 862 vs 982 Cα).
- Do **not** write that contact-count matches Vina wrong-pocket in **magnitude** (PM: Vina 0.788 vs contact min 0.552). Dual vs A_only size gap on AChE/BChE is 35.1 vs 34.0.
- Do **not** write that potency matching “solved” the holdout wrong-pocket paradox.
- Do **not** write PLIF, rotamer, or “mechanism solved.”
- Do **not** write that DualDiff/FuseDiff were re-scored on DualFourClass-Bench. The Intro/Discussion sentence is a **use-case** for the benchmark, not a generative bake-off.
- Do **not** write that DualDiff’s Dual High Affinity is mean-pooling; it is dual success vs a **reference ligand’s dock scores**. The gap is missing experimental hard-negatives, not the algebraic form of the mean.
- GNINA best-of-9 fair rescore (2026-08-24, real, not fabricated): pocket-matched summary_min moves by only −0.04 to +0.08 across K=4. On three of four pairs (EGFR/HER2, AChE/BChE, PIK3CA/mTOR) it stays at or below the same-panel Vina pocket-matched value. On **PIK3CA/PIK3CB it is marginally above** the Vina reference under both mode01 (0.554) and best9 (0.533) vs Vina 0.500 — this predates the best9 push (mode01 was already above) and the margin is well inside the bootstrap CI overlap, so it is **not** a "GNINA beats Vina" finding; both are near chance. Do **not** write "GNINA never exceeds Vina on any pair" (false for PIK3CA/PIK3CB) or that GNINA best-of-9 "resolves" the weak-signal finding on EGFR/HER2 or AChE/BChE (both remain ≤ 0.41, below chance). Do **not** call the `min(score_A, score_B)`-for-both-contrasts number in `GNINA_BEST9_STATUS.md` "pocket-matched" — it is **worst-pocket**; use `GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md` for the true directional number comparable to Table 2.
- Do **not** write that PIK3CA/mTOR is a structure-invariant reproducible success; holdout is ligand-set same-direction, crystal swap is receptor-dependent, and mTOR-swap CI includes 0.5.

## Engine stack actually run
| Channel | Status |
|---------|--------|
| AutoDock Vina | primary sampler (done) |
| RTMScore best-of-K | primary rescorer (done) |
| GNINA CNN | **DONE** — best-of-9 `--cnn_scoring rescore --minimize` (v1.3.2 CPU), pose-symmetric with RTM as of 2026-08-24; mode_01-only results retained as historical backup, not primary |
