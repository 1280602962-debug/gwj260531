# Claim-hardening round (2026-08-24)

Internal note. Cloud-doable fixes only; no new docking. Scripts: `data/jcim_novelty_v0/scripts/claim_hardening_v1.py`.

## What was closed in this round

| Reviewer issue | Action |
|---|---|
| Dual-vs-neither called “conventional benchmark” | Renamed **dual-versus-neither / nonselectivity-controlled comparator**. Descriptive contrast only. |
| `summary_min` as unique aggregator | Table S26: min / mean / harmonic. Pair ranking unchanged (PM > AChE > PIK3CB > EGFR). EGFR gap remains under all three. |
| Panel-construction difference | Methods 2.3 + Results 3.2: four AUROCs mix construction rules, n, series, receptor. |
| ECFP incremental overclaim | “Under the present scaffold-grouped benchmark, adding docking produced little incremental AUROC beyond ECFP4.” |
| Scaffold CV as generalization | Explicitly not target-external; PM nearly leave-one-scaffold. |
| Best-of-4 descriptor competitor | Table 2 + Table S28 list all four; **best single-descriptor reference**. |
| Wrong-pocket as specificity proof | Unresolved; **not a reliable universal negative control under panel shift.** |
| HOAP_028 silent missingness | Table S27 census; chemical-coverage failure. Main-panel AChE 95/100, PIK3CB 99/100. |
| Identifier prefix as diversity | Deleted. Class quotas + deterministic shuffle; no extra diversity constraint. |
| 0.756 vs 0.430 as significance | Forbidden. Descriptive formulation contrast. |
| T ≥ 0.3 as matched analogues | Similarity-constrained subset only. |
| Four-class classification | Four-state benchmark with two directional primary tasks. |
| Story A vs B | Formulation (EGFR) is the spine; chemotype / receptor are mechanism. PM is a conditional exception. |

## Not closed here (still local / API)

- Regenerating numbered main-text figures so Fig 3 is the formulation bar chart (plan updated; plot CSVs exist as Fig S4). A4 and B5 are now in the manuscript and Fig 5B.

## Closed after A4 / B5 (separate round)

- Full-panel max vs median: Table S29; label agreement + metric shift.
- Second-pair receptor swap: Table S30; opposite-direction realization effect; PAB_034 100/99/1.

## Claim that remains defendable

EGFR/HER2 provides a direct example in which omitting selective hard negatives gives an overly favorable impression (0.756 vs 0.430) that also appears as mixed-library enrichment of selectives. Other pairs show weaker or inconclusive formulation effects. K = 4 remains a data-constrained case panel.
