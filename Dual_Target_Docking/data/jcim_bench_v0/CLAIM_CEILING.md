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

## Engine stack actually run
| Channel | Status |
|---------|--------|
| AutoDock Vina | primary sampler (done) |
| RTMScore best-of-K | primary rescorer (done) |
| GNINA CNN | **DONE** — mode_01 `--cnn_scoring rescore --minimize` (v1.3.2 CPU) |
