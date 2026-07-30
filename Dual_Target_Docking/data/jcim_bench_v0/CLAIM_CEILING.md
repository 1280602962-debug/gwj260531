# CLAIM_CEILING — DualFourClass-Bench (JCIM)

## Allowed claims
1. This is an **evaluation / benchmark** article for dual-target docking decision metrics.
2. Primary metric = **pocket-matched directional AUROC** (D vs A_only uses pocket-B score; D vs B_only uses pocket-A score) with **trivial baselines** and pooled/wrong-pocket controls reported in parallel.
3. K=4 pairs are a **frozen evaluation set**, not a claim that the metric generalizes to all target pairs.
4. EGFR/HER2 is a **supply-limited case study** (existing unified RDKit EH110); no claim from new EGFR docking.
5. Prep protocol is frozen: **RDKit ETKDG + meeko**. Do **not** mention Schrodinger LigPrep in the manuscript (no formal license; early borrow was internal-only).

## Forbidden / over-claim
1. Do **not** claim a universal “decision arm” or that `rtm_min_z` (or any pooled score) is validated as a general master score.
2. Do **not** claim score-function invariance beyond the three channels actually run (Vina / RTM / GNINA mode_01 CNN).
3. Do **not** present Track B method competition results as the paper’s core.
4. Do **not** fold shortfall/clash flags or architecture tags into the primary score.
5. Do **not** mix LigPrep and RDKit poses in primary tables.

## Engine stack actually run
| Channel | Status |
|---------|--------|
| AutoDock Vina | primary sampler (done) |
| RTMScore best-of-K | primary rescorer (done) |
| GNINA CNN | **DONE** — mode_01 `--cnn_scoring rescore --minimize` (v1.3.2 CPU) |
