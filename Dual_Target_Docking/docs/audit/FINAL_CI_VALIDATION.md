# FINAL CI VALIDATION

Date: 2026-09-20  
Branch: `cursor/scientific-freeze-c7cc`  
PR: https://github.com/1280602962-debug/gwj260531/pull/39

This is repository CI policy evidence. It is not manuscript Methods.

## Verdict

**FINAL_CI_VALIDATION = PASS**

Failed run `#348` was not rewritten. A new `revision-validate` pair was triggered by commit `a088db7e`.

| Run | Event | Result | URL |
|-----|-------|--------|-----|
| #348 | pull_request (pre-fix) | failure at `compare_rebuild_to_canonical.py` | historical |
| #349 | push `a088db7e` | **success** | https://github.com/1280602962-debug/gwj260531/actions/runs/35482799463 |
| #350 | pull_request `a088db7e` | **success** | https://github.com/1280602962-debug/gwj260531/actions/runs/35482799840 |

## Commit

`a088db7e` Tolerate machine-level 1e-6 OOF probability drift only in the CI CSV compare.

Comparator change: exact text for identifiers, memberships, labels, counts, and all canonical summary tables; `abs(rebuilt − canonical) <= 1e-5` **only** for `ecfp4_oof_predictions.csv:oof_prob`. Policy: `docs/CI_REBUILD_COMPARE.md`.

`results/canonical/*.csv` were not modified.

## Environment (#349)

- Python 3.12.14
- numpy 2.5.2
- pandas 3.0.5
- scipy 1.18.1
- scikit-learn 1.9.0
- RDKit 2026.03.5
- matplotlib 3.11.1
- Meeko 0.7.1
- gemmi 0.7.5

`scripts/check_analysis_env.py`: PIN_MATCH.

## Workflow steps (#349 / #350)

| Step | Result |
|------|--------|
| Analysis environment | PASS |
| Fresh zero-dock rebuild | PASS |
| `verify_freeze_rebuild.py` | PASS (`result: PASS`; ECFP max \|Δ\| = 0.0234) |
| `compare_rebuild_to_canonical.py` | PASS |
| Current scientific chain | PASS |
| Adjudication tables unchanged | PASS |

## Compare details (#349)

- exact-match cells: 73222
- tolerated numeric cells: **0**
- maximum tolerated absolute difference: **0**

The `#348` mismatch (`JAK1/JAK2 | J1J2_095 | D_vs_A | ECFP4+docking`, `oof_prob` 0.340862 vs 0.340861, abs 1e-6) did not recur on `#349`/`#350`. A local synthetic check of that exact cell is `TOLERATED` at `1e-6` and FAIL at `>1e-5`. Identifier fields remain exact.

## Local confirmation (same freeze env, before push)

- verify = PASS
- Table 2 MATCH (64/64 EN+ZH)
- Table 3 MATCH (64/64 EN+ZH)
- stale errors = 0
- current chain = PASS
- canonical CSV SHA-256 unchanged 31/31

## What was not done

No docking, no θ change, no PRIMARY_PAIRS change, no rewrite of Table 2 / Table 3 / AUROC / CI / ECFP ΔAUROC / fold assignments, and no widening of tolerance beyond `ecfp4_oof_predictions.csv:oof_prob`.
