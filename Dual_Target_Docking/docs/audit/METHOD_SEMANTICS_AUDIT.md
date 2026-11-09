# Method-semantics audit

Date: 2026-09-20

Checks: common score scale; train/test isolation; pair/arm/fold keys; no cross-arm lookup; point and CI from the same statistic; label vs source field; no borrowed constituent CI; descriptive vs predictive wording.

| Analysis | Classification | Notes |
|----------|----------------|-------|
| Full-panel single-descriptor screen | VALID_WITH_LIMITATION | Same raw-descriptor scale within each named column. Best-descriptor Δ CI is conditional on full-panel selection. Labeled descriptive in Methods / SI S5. |
| Descriptor nested CV (pre-correction) | INVALID_RECOMPUTED | Heterogeneous raw pooling; 21 fold-level D-vs-B test cells used the other arm’s descriptor; no inner CV; no `summary_min` CI. Replaced by probability-scale nested CV. |
| Descriptor nested CV (corrected) | VALID_WITH_LIMITATION | Common P(dual) scale. Inner CV on outer-train only. Keys are pair/arm/outer_fold. `summary_min` CI is ligand-level bootstrap of fixed OOF probabilities, not selection-adjusted model uncertainty. |
| ECFP4 / ECFP4+docking | VALID | Per-arm GroupKFold; unscaled logistic P(dual); same splits for incremental Δ. Does not share descriptor selection or the `directional_smin(val_da)` bug. |
| ECFP4 scaler sensitivity | VALID | StandardScaler fit on each training fold only. Not a replacement for the unscaled 0.023 primary increment. |
| Matched / mismatched pocket | VALID | Same Vina scores; pocket assignment swapped. Δ CI from paired ligand bootstrap of Δ = matched − mismatched `summary_min`. |
| Fixed-score negative-class Δ | VALID | Shared dual resample; selective and neither resampled independently. Point and CI are the Δ statistic. |
| Five-seed Vina | VALID | Same boxes and eligible labels; seed is the only intended change. Not a CI for a min-of-seeds statistic. |
| Receptor substitution | VALID | One pocket replaced; `summary_min` CI from `summary_min_stratified`. Existing point/CI unchanged this round. |
| Independent GNINA | VALID_WITH_LIMITATION | `summary_min` CI was already min-inside-replicate. Weaker-arm display previously borrowed that CI (INVALID display). Corrected: directional CI stored and shown. Complete-case n now taken from the same rows. |
| RTM / CNN same-pose rescore | VALID | Same Vina poses, different score channel. `summary_min` CI from `summary_min_stratified`. |
| Max vs median | VALID | Aggregation sensitivity on stored poses; not a new docking run. |
| Cluster bootstrap | VALID_WITH_LIMITATION | Ligand-level primary Δ remains official. Cluster intervals are sensitivity; JAK1/TYK2 document cluster is not recomputed. |
| Detectable-effect simulation | VALID_WITH_LIMITATION | Binormal Monte Carlo P(CI excludes 0.5). Not observed power. Uses the same shared-dual bootstrap protocol as Table 2. |
| Table 2 primary directional / `summary_min` | VALID | Point from complete sample; CI from shared-dual min-inside-replicate. PIK3CA Vina shows the distinction ([0.480, 0.802] vs weaker-arm [0.491, 0.868]). |
| Table 3 two-pocket mean / Top-10 / AND | VALID | Ranking on two-pocket mean; Top-10 and AND are descriptive operating points. |

## Explicit answers

1. Corrected descriptor OOF scores share a probability scale. Old raw pooling did not.
2. Different folds may select different descriptors; the pooled score is now P(dual).
3. Outer-test IDs are disjoint from selection, scaler, and logistic-fit IDs (leakage audit = 0).
4. Lookup keys are pair + arm + outer_fold. Duplicate keys fail the generator.
5. D-vs-A and D-vs-B no longer share a two-arm `directional_smin` call on one array.
6–8. CI semantic audit: `same_statistic = 1` on current tables (`results/qa/metric_ci_semantic_audit.csv`).
9. Descriptive screen and nested predictive baseline are named separately.
10. Full-panel best descriptor is not treated as prespecified predictive performance.

## Table 2 / Table 3 / ECFP4

`TABLE2_AFFECTED = NO`  
`TABLE3_AFFECTED = NO`  
`ECFP4 NOT AFFECTED` (max |Δ AUROC| = 0.0234)
