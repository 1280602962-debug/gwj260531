# CI rebuild comparison policy

This is repository reproducibility / CI policy. It is not manuscript Methods.

GitHub-hosted `revision-validate` runners may produce machine-level floating-point differences in logistic-regression OOF probabilities at approximately `1e-6`. Observed case: `ecfp4_oof_predictions.csv` row JAK1/JAK2 | J1J2_095 | D_vs_A | ECFP4+docking, column `oof_prob`, canonical `0.340862` vs CI rebuild `0.340861`.

`scripts/qa/compare_rebuild_to_canonical.py` therefore uses:

- exact comparison for identifiers, memberships, labels, counts, fold assignments, and all canonical summary outputs (Table 2 / Table 3 / AUROC / CI / class counts / ECFP ΔAUROC);
- tolerance `abs(rebuilt − canonical) <= 1e-5` **only** for `ecfp4_oof_predictions.csv:oof_prob` and `descriptor_nested_oof_predictions.csv:oof_probability`;
- a printed `TOLERATED:` line for every in-tolerance cell; differences larger than `1e-5` still FAIL.

Independent AUROC / ΔAUROC replay in `scripts/qa/verify_freeze_rebuild.py` remains mandatory and is not relaxed.

Do not widen this tolerance to other files or columns. A new unexpected difference is a stop-and-report event, not an automatic tolerance increase.

Python on CI is 3.12.x. Scientific packages stay exact-pinned (`requirements-analysis.txt`). The 1e-6 OOF difference is not treated as a Python-patch failure.
