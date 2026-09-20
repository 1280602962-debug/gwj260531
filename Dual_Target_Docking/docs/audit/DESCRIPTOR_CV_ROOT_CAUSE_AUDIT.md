# Descriptor CV root-cause audit

Date: 2026-09-20  
Code inspected: `scripts/analysis/compute_descriptor_baselines.py` (pre-correction), `scripts/analysis/fit_ecfp4_models.py`, `scripts/analysis/bootstrap_metrics.py`, `scripts/analysis/patch_publication_text.py`, `results/canonical/descriptor_nested_scaffold_cv.csv`, `results/canonical/descriptor_baselines.csv`, SI Table S5 / S7a.

No production code was changed before this report was written.

## Locators

| Item | Location |
|------|----------|
| Single-descriptor generator | `scripts/analysis/compute_descriptor_baselines.py` (`nested_oof`, `main`) |
| Fold table | `results/canonical/descriptor_nested_scaffold_cv.csv` (40 pair×fold rows) |
| Pair-level nested numbers | `results/canonical/descriptor_baselines.csv` columns `nested_scaffold_cv_oof_*` |
| Descriptor selection table | same fold table, columns `selected_D_vs_A`, `selected_D_vs_B` |
| OOF prediction table | **absent** (no ligand-level descriptor OOF file) |
| Fold assignment table | **absent** for descriptors; ECFP4 uses `model_fold_assignments.csv` (separate pipeline) |
| Table S5 source | `descriptor_baselines.csv` + `ecfp4_incremental_information.csv` via `patch_publication_text.py` `si_s5_desc` / `si_s5_ecfp` |
| Table S7a source | `patch_publication_text.py` lines 948–974; GNINA cells from `computational_robustness.csv` |
| Manuscript / SI citations | Methods (descriptor vs nested OOF), Results (AChE TPSA 0.742 / 0.801; PIK3CA heavy 0.463); SI S5 prose |
| Figures using descriptor results | Figure S1 (descriptive screen, not nested OOF). Figure 3 is ECFP4, not single-descriptor nested CV |
| `summary_min` CI functions | `bootstrap_metrics.summary_min_stratified` (correct min-inside-replicate). Descriptor nested OOF has **no CI**. S7a GNINA weaker-arm CI copies `summary_min_ci_*` |

---

## A–L answers

### A. How are outer folds defined?

`nested_oof` builds one `sklearn.model_selection.GroupKFold` on **all ligands of the pair** (dual, A-only, B-only, and neither). Dummy `X` is zeros; grouping variable is Bemis–Murcko scaffold (`scaffold_of`). `n_splits = min(GROUPKFOLD_MAX_SPLITS=5, n_scaffolds)`. Both directional arms reuse this same split.

### B. Is scaffold GroupKFold used?

Yes, for the outer split only.

### C. Is each pair × arm partitioned independently?

**No.** One pair-level GroupKFold is shared by D-vs-A and D-vs-B. ECFP4 (`fit_ecfp4_models.py`) *does* split each arm independently on contrast ligands only. The two pipelines are not the same.

### D. Where is the descriptor selected?

Inside the outer-fold loop of `nested_oof`, lines 131–139. For each of `{tpsa, clogp, heavy, mw}`, training-set raw values are scored with `directional_smin`. The name with the highest training D-vs-A AUROC is `best_da_name`; the name with the highest training D-vs-B AUROC is `best_db_name`.

### E. Does selection use only outer-training data?

**Yes** for the *choice* of name. Selection AUROC is `directional_smin(val[tr], tr_cls)`.

### F. Is there an inner CV?

**No.** There is no inner GroupKFold. The analysis is named “nested scaffold CV” but selection is a single training-set AUROC screen. That name is not justified.

### G. Does the test fold participate in descriptor selection?

**No.** Test indices are not used in the `da > best_da` / `db > best_db` loop.

### H. How is descriptor orientation set?

There is **no model**. Raw descriptor values are used as ranking scores. `auroc(pos, neg)` treats **higher value as more dual**. There is no training-learned coefficient, no sign flip, and no StandardScaler. If the true association is “lower TPSA → dual”, the raw AUROC is < 0.5 and another descriptor usually wins. This is not a learned orientation; it is a fixed “higher-is-dual” ranking.

### I. What is written for the test fold?

**Raw descriptor values**, not rank, z-score, or probability:

```python
oof["selected_DvsA"][te] = val_da[te]
oof["selected_DvsB"][te] = val_db[te]
```

### J. Which field is used for pooled AUROC?

Pooled D-vs-A AUROC uses `oof["selected_DvsA"]` restricted to dual vs A-only. Pooled D-vs-B uses `oof["selected_DvsB"]` restricted to dual vs B-only. Those arrays are concatenations of **heterogeneous raw descriptor columns**.

### K. Why were different raw descriptors allowed to be concatenated?

The docstring (lines 107–109) states this was accepted as “the nested-selection readout”. There is no scale unification. TPSA (~0–200 Å²), cLogP (~−2–8), MW (~200–600), and heavy-atom count (~15–40) are pooled as if they were one score. That pooled AUROC has no common scale and is statistically invalid.

### L. Root cause of the D-vs-B mapping errors

**Not** a missing `arm` dictionary key, **not** a pair/fold-only merge, **not** a dataframe join, **not** column-name transposition, and **not** loop-variable overwrite of `best_db_name`.

The stored *selections* `selected_D_vs_A` / `selected_D_vs_B` are computed separately and are internally consistent.

The bug is the **test-fold evaluation call**:

```146:146:scripts/analysis/compute_descriptor_baselines.py
        te_da, te_db, te_sm = directional_smin(val_da[te], te_cls) if (te_cls == "dual").any() else (np.nan, np.nan, np.nan)
```

`directional_smin(vals, cls)` always returns **both** D-vs-A and D-vs-B AUROCs from the **same** `vals` array. Here `vals = val_da[te]`, i.e. the D-vs-A selected descriptor. Therefore:

* `test_D_vs_A` uses the D-vs-A selection (correct when the score is that raw column).
* `test_D_vs_B` uses the **D-vs-A** selection (wrong whenever `selected_D_vs_A != selected_D_vs_B`).
* `test_summary_min` is `min` of those two, so it is also wrong on those folds.

`val_db` is computed and written into `oof["selected_DvsB"]`, but it is never passed to `directional_smin` for the fold-level test metrics.

Pre-correction count: **21 of 40** pair×fold rows have `selected_D_vs_A != selected_D_vs_B`. Those 21 D-vs-B test cells are cross-arm. The user-facing “20” count is this set (21 in the frozen CSV). AChE/BChE selected TPSA on both arms in every fold, so it is not in the error set.

Pooled pair-level `nested_scaffold_cv_oof_D_vs_B` uses `oof["selected_DvsB"]` and is **not** cross-arm. It is still invalid because of heterogeneous raw pooling (issue A).

Mechanism class:

| Hypothesis | Verdict |
|------------|---------|
| Dictionary key missing `arm` | No |
| Merge key only pair/fold | No (no merge) |
| D-vs-A selection reused for D-vs-B **test AUROC** | **Yes** — `directional_smin(val_da[te], …)` |
| Loop-variable overwrite of selection names | No |
| Dataframe join missing contrast | No |
| Column naming swap in the CSV writer | No |

---

## Table S7a / `summary_min` CI

Two separate CI problems:

1. **Descriptor nested `summary_min` has no CI at all.** `descriptor_baselines.csv` stores point OOF AUROCs only. Nothing in SI Table S5 displays a nested `summary_min` CI. The chemistry-control number is therefore a point estimate without a matching interval.

2. **SI Table S7a (independent GNINA) weaker-arm CI is copied from `summary_min` CI.** In `patch_publication_text.py` lines 969–972 the GNINA weaker-arm cell uses `gp['summary_min_ci_lo/hi']`, not a directional AUROC CI. `computational_robustness.csv` does not store directional CIs. `summary_min_stratified` *does* compute them (`auroc_D_vs_A_ci_*`, `auroc_D_vs_B_ci_*`) but `gnina_rows` discards them.

   The GNINA `summary_min` interval itself is produced by `summary_min_stratified` (min inside each shared-dual replicate). That point/CI pair is semantically matched. The **weaker-arm** column is not: it prints a directional point with the min-of-two interval.

   Vina-primary S7a rows are correct: `summary_min` CI from `primary_summary_min.csv`, weaker-arm CI from `primary_directional_auroc.csv`. PIK3CA/mTOR Vina shows the distinction: `0.692 [0.480, 0.802]` vs weaker-arm `0.692 [0.491, 0.868]`.

S7a GNINA class counts in the generator are hardcoded (`28 / 37 / 32 / 11` for EGFR/HER2) and do not match `computational_robustness` complete-case counts (EGFR/HER2 independent GNINA is 20 / 33 / 26 / 10 in the current SI). That is a display-source bug, not a docking-score bug.

---

## What is not contaminated

* Table 2 / Table 3 generators (`compute_canonical_results.py`) do not call `nested_oof` or `compute_descriptor_baselines.py`.
* ECFP4 (`fit_ecfp4_models.py`) has its own GroupKFold, its own per-arm loop, and writes probabilities via `LogisticRegression.predict_proba`. It does not import descriptor selection or the `directional_smin(val_da)` call.
* Figure 3 uses ECFP4 OOF AUROC, not single-descriptor nested OOF.
* Figure S1 uses the full-panel descriptive screen (`best_descriptor_summary_min`), not the invalid pooled raw OOF.

---

## Required remediation (not done in this report)

1. Replace raw-value pooling with train-only inner scaffold CV → outer-train StandardScaler + single-variable logistic regression → held-out P(dual).
2. Evaluate each arm with that arm’s selected descriptor / model. Lookup key: `pair`, `arm`, `outer_fold`.
3. Write ligand-level OOF provenance.
4. Compute `summary_min` CI from shared-dual bootstrap of the **min of the two directional OOF AUROCs**, labeled as ligand-level bootstrap of fixed OOF predictions.
5. Store and display directional CIs separately from `summary_min` CIs (descriptor table and Table S7a GNINA).
6. Keep the full-panel univariate screen as a descriptive analysis; do not call it selection-adjusted prediction.
