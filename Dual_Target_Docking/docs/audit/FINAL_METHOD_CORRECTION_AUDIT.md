# FINAL METHOD CORRECTION AUDIT

Date: 2026-09-20  
Previous `READY_FOR_FINAL_MANUSCRIPT_DRAFTING` marks are superseded by this audit.

## Verdict

**READY_FOR_FINAL_MANUSCRIPT_DRAFTING**

P0 = 0  
Actionable P1 = 0  
P2 = 0  

## CLOSED

| Item | Status |
|------|--------|
| raw heterogeneous descriptor pooling | **CLOSED** — pooled score is P(dual) ∈ [0,1] |
| 20/21 D-vs-B cross-arm mapping bug | **CLOSED** — 80/80 pair×arm×fold units CORRECT; generator fails if `selection_source_arm != arm` |
| descriptor test leakage | **CLOSED** — leakage = 0 (IDs and scaffolds) |
| corrected OOF score scale | **CLOSED** — probability |
| Table S7a `summary_min` CI mismatch | **CLOSED** — nested `summary_min` CI is joint shared-dual min bootstrap; GNINA weaker-arm CI is now directional ([0.419, 0.825] vs `summary_min` [0.410, 0.769] on PIK3CA) |
| metric-to-CI semantic audit | **PASS** — `results/qa/metric_ci_semantic_audit.csv`, all `same_statistic = 1` |
| Table 2 unaffected | **confirmed** — file identical to pre-correction snapshot; 128/128 lock MATCH |
| Table 3 unaffected | **confirmed** — `two_pocket_mean_ranking.csv` and Top-10 identical |
| ECFP4 unaffected / revalidated | **confirmed** — `ecfp4_incremental_information.csv` identical; max \|Δ AUROC\| = 0.0234 |
| Figure 2 overlap | **CLOSED** — layout-only fix; n=4 outside CI |
| Figure S4 overlap | **CLOSED** — layout-only fix; dagger under colorbar |
| 11-figure numerical trace | **PASS** — 347/347 MATCH |

## Rebuild

Pinned analysis environment (Python 3.12.3, `requirements-analysis.txt`):

- `rebuild_freeze.py` → `/tmp/dual_target_freeze_rebuild` complete
- `compare_rebuild_to_canonical.py` PASS (91224 exact cells; 0 tolerated)
- `verify_freeze_rebuild.py` PASS (Table 2, Table 3, five-seed, ECFP, descriptor nested assertions)
- `check_current_chain.py` PASS
- descriptor mapping / leakage / independent replay / canonical diff PASS
- `UNEXPECTED_CHANGE = 0`

`TABLE2_AFFECTED = NO`  
`TABLE3_AFFECTED = NO`  
`ECFP4 NOT AFFECTED`

## Corrected nested descriptor `summary_min`

| Pair | Old invalid raw-pool | Corrected P(dual) OOF [95% CI] |
|------|---------------------:|--------------------------------|
| EGFR/HER2 | 0.4171 | 0.6325 [0.4666, 0.6926] |
| JAK1/JAK2 | 0.5996 | 0.4380 [0.2939, 0.5913] |
| JAK1/TYK2 | 0.5630 | 0.4587 [0.3019, 0.5393] |
| PIK3CA/mTOR | 0.3542 | 0.3889 [0.1943, 0.5927] |
| AChE/BChE | 0.7422 | 0.6496 [0.4972, 0.7964] |
| F2/F10 | 0.4768 | 0.3538 [0.2278, 0.4925] |
| PPARG/PPARA | 0.4595 | 0.3579 [0.2266, 0.4990] |
| PPARA/PPARD | 0.4756 | 0.4438 [0.3008, 0.5410] |

The chemistry-control claim is weakened. Full-panel TPSA on AChE/BChE remains a descriptive observation only.

## Reports

- `docs/audit/DESCRIPTOR_CV_ROOT_CAUSE_AUDIT.md`
- `docs/audit/DESCRIPTOR_RESULT_IMPACT.md`
- `docs/audit/DESCRIPTOR_METHOD_CORRECTION_REPORT.md`
- `docs/audit/METHOD_SEMANTICS_AUDIT.md`
- `docs/audit/FIGURE_VISUAL_QA.md`
- `results/qa/descriptor_cv_fold_mapping_audit.csv`
- `results/qa/descriptor_nested_cv_leakage_audit.csv`
- `results/qa/metric_ci_semantic_audit.csv`
- `results/qa/method_correction_canonical_diff.csv`
- `results/qa/figure_value_trace.csv`

## Stop condition

All listed stop conditions are met. Core docking conclusions are unchanged. The auxiliary single-descriptor nested-CV statistic is now defined as a probability-scale OOF AUROC with a matching `summary_min` bootstrap CI.
