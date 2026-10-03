# THRESHOLD_PROVENANCE_GATE

Two gates. They are independent. This file does not relabel any ligand.

Primary θ remains 6.0.

## GATE_THRESHOLD_PRESPEC = CONFIRMED

Question: were 5.5, 6.5, and `strict_6.5_5.5` specified as sensitivities before V4.2 PRIMARY directional results existed?

### Evidence on the current orphan branch

`FORMAL_METRICS_ANALYSIS_FREEZE_FINAL.yaml` first appears in commit `af328e3e52a6bdd4e1c49575c1f6f515e227ce10` (2026-09-30), the same commit that also adds `results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv`.

On this branch alone, git cannot separate the YAML from the PRIMARY tables.

The YAML text (also present at `af328e3e`) is:

```yaml
labels:
  theta_primary: 6.0
  sensitivity_only: [5.5, 6.5, strict_6.5_5.5]
```

Header: written by `pre_metrics_readiness_audit_v42.py`; `metrics_execution_unlocked: false`.

### Evidence from other git objects (not `/tmp/pr39_fiveseed` files)

`git log --all -S 'strict_6.5_5.5'` hits `cursor/scientific-freeze-c7cc` commit `a0ef2726036e35f5bf2b78412f9a3ba43506e8e9` (2026-09-20).

That commit does **not** contain `PRIMARY_DIRECTIONAL_METRICS.csv`.

`git show a0ef2726:Dual_Target_Docking/scripts/analysis/compute_canonical_results.py` defines:

```python
rules = [("theta_5.5", 5.5, False), ("theta_6.0", 6.0, False), ("theta_6.5", 6.5, False), ("strict_6.5_5.5", None, True)]
```

inside `compute_label_sensitivity`.

Those `scripts/analysis/*` files are **forbidden as future compute imports**. They are used here only as historical provenance that the four-rule list existed before V4.2 PRIMARY results.

The current formal YAML restates that same list as `sensitivity_only`.

### Gate decision

`GATE_THRESHOLD_PRESPEC = CONFIRMED`

The sensitivity **definitions** predate V4.2 PRIMARY results. This is not a new cutoff invented after seeing M0–M3 AUROC.

## GATE_ACTIVITY_SOURCE_READ = NOT_GRANTED

Question: may this project now read continuous pA/pB and relabel under the pre-specified thresholds?

No.

Reasons:

1. `FORMAL_AUTHORITY_PATHS.yaml` lists `data/processed/activity_adjudication/ligand_activity_aggregate_v1.csv` under `forbidden_activity_regeneration`.
2. That file is **absent** from the current release tree.
3. It exists as a git object on `cursor/scientific-freeze-c7cc` (`a0ef2726`) with columns `historical_pA`, `historical_pB`, `max_A`, `median_A`, `max_B`, `median_B`, `class_max`, `class_median`. Seeing those column names is not a read-permission for formal relabel.
4. This Stage 1 file does not copy that CSV into the release tree.

`data_readiness` for Module G: `NEED_FROZEN_SOURCE_PERMISSION`

A later human decision may grant `SENSITIVITY_ONLY_READ_PERMISSION` limited to:

- read-only frozen continuous values;
- relabel only at 5.5 / 6.0 / 6.5 / strict 6.5/5.5;
- no live ChEMBL;
- no re-adjudication;
- no change to primary θ=6.0 or the published primary classes.

Until that permission exists, threshold AUROC and transition tables must not run.

CONFIRMED on gate 1 does not pass gate 2.
