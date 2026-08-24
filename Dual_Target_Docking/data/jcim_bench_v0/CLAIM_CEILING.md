# CLAIM_CEILING — DualFourClass-Bench (JCIM)

## Allowed claims
1. This is an **evaluation / benchmark** article for dual-target docking decision metrics.
2. Primary metric = **pocket-matched directional AUROC** (D vs A_only uses pocket-B score; D vs B_only uses pocket-A score) with **trivial baselines** and pooled/wrong-pocket controls reported in parallel.
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

See also: `docs/POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`.

## Mechanism-analysis claim ceiling (Results 3.9 / 3.11)

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
