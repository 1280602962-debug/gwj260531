# CLAIM_CEILING — DualFourClass-Bench (JCIM)

## Allowed claims
1. This is an **evaluation / benchmark** article for dual-target docking decision metrics.
2. Primary metric = **pocket-matched directional AUROC** (D vs A_only uses pocket-B score; D vs B_only uses pocket-A score) with **trivial baselines** and pooled/wrong-pocket controls reported in parallel.
3. K=4 pairs are a **frozen evaluation set**, not a claim that the metric generalizes to all target pairs.
4. EGFR/HER2 is a **supply-limited case study** (existing unified RDKit EH110); no claim from new EGFR docking.
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

See also: `docs/POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`.

## Mechanism-analysis claim ceiling (Results 3.9 / 3.11)

Allowed:
- PIK3CA inter-crystal global Cα RMSD (1.44–1.49 Å) is larger than mTOR (0.45 Å) **in this structure set**, consistent in **direction** with greater PIK3CA-end sensitivity.
- Local pocket Cα conservation on 5DXT (0.343 Å) with collapsed summary_min (0.505) means Cα pocket conservation is **not sufficient**.
- Scoring-free `contact_count` shows a real size/burial confound, **especially on the B arm** (0.698–0.714).

Forbidden:
- Do **not** write that Cα RMSD “quantitatively explains” or causes the AUROC collapse (n = 2 / n = 1; 5DXT matched 862 vs 982 Cα).
- Do **not** write that contact-count matches Vina wrong-pocket in **magnitude** (PM: Vina 0.788 vs contact min 0.552). Dual vs A_only size gap on AChE/BChE is 35.1 vs 34.0.
- Do **not** write PLIF, rotamer, or “mechanism solved.”
- Do **not** write that PIK3CA/mTOR is a structure-invariant reproducible success; holdout is ligand-set same-direction, crystal swap is receptor-dependent, and mTOR-swap CI includes 0.5.

## Engine stack actually run
| Channel | Status |
|---------|--------|
| AutoDock Vina | primary sampler (done) |
| RTMScore best-of-K | primary rescorer (done) |
| GNINA CNN | **DONE** — mode_01 `--cnn_scoring rescore --minimize` (v1.3.2 CPU) |
