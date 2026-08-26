# Post-review experiment plan (six suggestions)

Judgment for a JCIM evaluation article. No new docking was run in the cloud
environment (no `vina` / `gnina` / `smina` on PATH; git has no full pose archive).

| Suggestion | Do it? | Where | Status |
|---|---|---|---|
| Unify all work to 6.5/5.5 **or** θ=6.0 | **No.** Keep dual-purpose rules. | Methods 2.1 wording | Done in this revision |
| Independent pose-generation engine | **Yes, P0.** One engine, EGFR/HER2 + PM48 | Local, now deposited | Done: `data/jcim_independent_dock_v0/`; EGFR gap remains (0.783 vs 0.220) |
| Receptor/wrong-pocket structural mechanism | **Partial P1.** PLIF as hypothesis, not cause | Local poses, now deposited | Done as geometric occupancy: `analysis/plif_v1/`; Cα RMSD already in 3.6 |
| Chemical confounding → diagnostic workflow | **Yes, P1.** Existing experiments, no docking | Figure 8, Discussion 4.4 | Done |
| Lower abstract/conclusion scope | **Yes, P0** | Abstract, Conclusions | Done: four-pair Vina-based benchmark; reliability boundary |
| Power analysis | **Yes, as detectable-effect simulation, not observed power** | Table S31, Figure S6, Results 3.2 | Done: N_MC=1000, B=2000, seed 20260729 |

## What the simulation says (summary_min)

P(95% CI excludes 0.5) at observed n:

| Pair | n | true 0.60 | true 0.70 | true 0.75 |
|---|---|---:|---:|---:|
| EGFR/HER2 | 28/38/32 | 0.065 | 0.621 | 0.907 |
| AChE/BChE | 27/25/28 | 0.049 | 0.504 | 0.828 |
| PIK3CA/PIK3CB | 28/27/28 | 0.041 | 0.564 | 0.849 |
| PIK3CA/mTOR | 18/14/12 | 0.025 | 0.219 | 0.452 |

Allowed sentence: the benchmark resolves large directional effects more readily
than moderate ones; CI including 0.5 does not establish chance equivalence.

## Local docking that would still help most

Independent GNINA docking search and the PIK3CA occupancy snapshot are now
deposited (Tables S32–S33). Do **not** add five scoring functions or unify
the two thresholds.
