# CLAIM_CEILING — DualFourClass-Bench (JCIM)

## Allowed claims
1. This is an **evaluation / benchmark** article for dual-target docking decision metrics.
2. Primary metric = **directional AUROC** (dual vs A_only; dual vs B_only) with **trivial baselines** reported.
3. K=4 pairs are a **frozen evaluation set**, not a claim that the metric generalizes to all target pairs.
4. EGFR/HER2 is a **supply-limited case study** (existing unified RDKit EH110); no claim from new EGFR docking.
5. Prep protocol is frozen: **RDKit ETKDG + meeko** (LigPrep used only as sensitivity delta on PM48).

## Forbidden / over-claim
1. Do **not** claim a universal “decision arm” or that `rtm_min_z` (or any pooled score) is validated as a general master score.
2. Do **not** claim score-function invariance: **GNINA CNN rescore was SKIPPED** (no binary on host).
3. Do **not** present Track B method competition results as the paper’s core.
4. Do **not** fold shortfall/clash flags or architecture tags into the primary score.
5. Do **not** mix LigPrep and RDKit poses in primary tables.

## Engine stack actually run
| Channel | Status |
|---------|--------|
| AutoDock Vina | primary sampler |
| RTMScore best-of-K | primary rescorer |
| GNINA CNN | SKIPPED — see `analysis/GNINA_STATUS.md` |
