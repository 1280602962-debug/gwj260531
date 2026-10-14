# NLRP3 classifier: read-only audit and corrected pipeline (2026-10-14)

Branch `cursor/nlrp3-model-calibration-audit-20261014`, started from
`cursor/urat1-nlrp3-dualtarget-aidd-e43d` at commit `bc6d555`.

Scope: the NLRP3 functional-activity classifier only. No frozen artefact was
modified: the 1588-molecule shrink pool, the 12 nominated and 21 reserve
candidates, all docking output, the URAT1 analyses, the PF-03882845 MD plan and
`MANUSCRIPT_FACTS.yaml` are untouched. The 8319-compound clinical library was not
re-scored.

All numbers below were recomputed from `data/processed/nlrp3_records.csv` on this
branch; none is quoted from an earlier report.

---

## 1. Data and labels (verified, not changed)

| Item | Verified value |
| --- | --- |
| Records | 609 |
| Unique compounds | 513 |
| Unique assays | 39 |
| Unique Bemis-Murcko scaffolds | 255 |
| Records per molecule | mean 1.19, max 15; 55 molecules have >1 record |
| Record-level active rate | 0.591 |
| Molecule-level active rate | 0.604 (310 active / 203 inactive) |

Curation (`utils_ml.curate_nlrp3_records`): ChEMBL export filtered to assay
descriptions matching `IL-1`, `Assay Type == "B"`, `Standard Relation == "="`,
pChEMBL between 4 and 10, assays with at least 5 records; SMILES canonicalised
with RDKit (`MolFromSmiles` → `SanitizeMol` → `MolToSmiles`); `active =
pActivity >= 6.0`; scaffolds from `MurckoScaffoldSmiles(includeChirality=False)`.

Molecule-level aggregation is `max` over records, i.e. a molecule counts as
active if any assay called it active. 23 molecules carry discordant records
(both an active and an inactive measurement); they are listed in
`results/nlrp3_model_v2/data_label_audit.json`. **This label definition was kept
unchanged**, as instructed. It is a real methodological weakness worth a
sensitivity analysis later — `max` aggregation is optimistic for compounds
measured in several assays of unequal stringency — but changing it at the same
time as fixing the calibration would have made the two effects inseparable.

Per-assay active counts are in the same JSON. Several assays are strongly
unbalanced (single-class assays exist), which is the original motivation for the
assay one-hot feature.

No defect was found in featurisation itself: `featurize_smiles` is a pure
per-molecule function with no fitted state, so it cannot leak across splits. The
leaks in the historical pipeline are in the *split logic*, below.

---

## 2. Defects found in the historical pipeline

### 2.1 The production calibrator was fitted in-sample (critical)

`scripts/02_train_asymmetric_models.py`, `fit_nlrp3_final()`:

```python
scaler = StandardScaler()
x_s = scaler.fit_transform(x)          # line 354 — fitted on ALL 609 records

base = _xgb_classifier(seed=seed)
base.fit(x_s, y, sample_weight=sample_w)   # line 357 — trained on ALL 609 records

cal_mask = df["canonical_smiles"].isin(cal_smiles).values   # line 360
raw_prob = base.predict_proba(x_s[cal_mask])[:, 1]          # line 363
calibrator = IsotonicRegression(out_of_bounds="clip")
calibrator.fit(raw_prob, y[cal_mask])                       # line 365
```

`scaffold_holdout_split` does carve out a scaffold-disjoint calibration set, but
that set is never withheld from `base.fit`. Reproducing the split with seed 42
gives **66 calibration molecules / 71 records (41 active, 30 inactive), 100 % of
which were also in the base model's training set**. The StandardScaler was
likewise fitted on all 609 records, including the calibration and the nominal
holdout molecules.

Consequence: the isotonic regression was fitted to probabilities the model had
already memorised, so it saw a far sharper separation than it will ever see on
new compounds. The deposited production calibrator
(`data/models/nlrp3_model.joblib`) is a 10-knot step function emitting only five
distinct values — {0, 0.25, 0.333, 0.917, 1.0} — with the jump between 0.333 and
0.917 straddling the decision threshold. Scoring the grid shows that
"calibrated probability ≥ 0.5" is exactly equivalent to "raw probability ≥
0.4505". The published cutoff is therefore not a probability statement at all; it
is an arbitrary rank cutoff produced by an over-fitted monotone transform, and it
is sensitive to the single calibration split.

What this does **not** mean: the historical scaffold-CV AUROC is not invalidated
by this. `train_nlrp3_cv` never touches the calibrator — it reports raw
probabilities from fold models fitted on training records only. Base-model
*discrimination* and final *probability reliability* are separate claims, and only
the second one is broken here. Section 4 confirms the first one reproduces.

Fix: nested design. The base model is fitted on inner-training scaffolds, the
calibrator on inner-calibration scaffolds the base model never saw, and the outer
test fold sees neither. Scalers are fitted on base-fit data only.

### 2.2 Cross-validation scored molecules with information unavailable at deployment

`train_nlrp3_cv` (lines 286–288) builds the test matrix from
`te_records["Assay ChEMBL ID"]`, the *true* assay each test molecule was measured
in, and takes the max over that molecule's real records. `screen_repurposing_library.py`
(`predict_nlrp3_ensemble`, lines 39–53) instead expands every clinical molecule
into the five most frequent training assays and takes `np.max` over those five
synthetic contexts. The two routines are not the same function, so the reported
AUROC did not evaluate the deployed rule.

A second, smaller issue: `top_assay_ids` is computed from the full dataset
(line 268) before the CV loop, so the assay vocabulary is chosen with test records
in view.

Fix: the deployment rule is evaluated directly — each outer-test molecule is
scored under five fixed assay contexts chosen from outer-training records only,
aggregated with `max`, and compared against its true molecule-level label. The
historical "true assay ID" rule is retained as a labelled reference variant so the
two can be compared on identical folds.

### 2.3 Related observations (reported, not acted on)

- `results/training/` does not exist on this branch; the production bundle lives
  at `data/models/nlrp3_model.joblib`, which is what `screen_repurposing_library.py`
  falls back to. The historical report and manifest are consistent on this point.
- `data/si/nlrp3_threshold_sensitivity/summary.json` already noted that per-molecule
  probabilities differ by <1e-6 across the five assay contexts for clinical
  molecules, so max/mean/median aggregation gave an identical pool. Section 4
  explains why: the assay one-hot carries almost no signal for molecules outside
  the training assays.
- The 1588 vs 1377 pool-size gap between the production run and the SI retrain is
  consistent with the two runs differing in calibration (production applies the
  in-sample isotonic, i.e. an effective raw cutoff of 0.4505; the SI retrain used
  raw probabilities at 0.5) plus retraining variance. This explanation is
  consistent with the artefacts but was not verified by re-scoring the library,
  which is out of scope here.

---

## 3. Corrected design

Implemented in `scripts/nlrp3_model_v2_validate.py`, with the model and
calibrator classes in `scripts/nlrp3_v2_lib.py` and the prediction entry point in
`scripts/nlrp3_v2_predict.py`.

- **Outer**: 5-fold `GroupKFold` on Murcko scaffold at the molecule level (the
  same `scaffold_cv_indices` helper and the same molecule ordering as the
  historical CV, so folds are comparable).
- **Inner**: within each outer-training set, a scaffold-grouped
  `GroupShuffleSplit` reserves 25 % of molecules for calibration; the base model
  is fitted on the remaining 75 %.
- Every record of a molecule stays with that molecule; scaffolds never cross a
  boundary that is required to be independent. Verified per fold and written to
  `validation_report.json → leakage_checks`: all SMILES and scaffold overlaps
  between fit, calibration and test are 0 in all 15 fold × repeat combinations.
- StandardScaler, assay vocabulary and the five deployment contexts are all
  derived from base-fit data only.
- The inner split is repeated with 3 seeds per outer fold to measure how much the
  calibration depends on which scaffolds land in the calibration set.
- Calibration methods compared: none, Platt (logistic regression on the logit),
  isotonic. Calibrators are only fitted when the calibration set has ≥30
  molecules and both classes present (satisfied everywhere: 87–121 molecules).
- Three scoring variants on identical outer folds: `A_structure_only`,
  `B_assay_deployment`, `B_assay_legacy_true_assay`.
- Model comparison by paired bootstrap over molecules (2000 resamples) on the
  pooled outer-test predictions.
- Thresholds 0.3/0.4/0.5/0.6/0.7 were fixed before looking at results and are
  reported for all variants.

Nothing in the historical scripts was edited; the corrected pipeline is additive.

---

## 4. Results

Full tables: `results/nlrp3_model_v2/performance_summary.csv` and
`validation_report.json`. Pooled outer-test (513 molecules, 60.4 % active):

| Variant | Calib. | AUROC | AUPRC | Acc | Prec | Rec | Spec | F1 | EF@10% | Brier |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Structure-only | none | **0.890** | **0.913** | 0.850 | 0.877 | 0.874 | 0.813 | 0.876 | 1.56 | **0.126** |
| Structure-only | Platt | 0.877 | 0.895 | 0.836 | 0.872 | 0.855 | 0.808 | 0.863 | 1.56 | 0.131 |
| Structure-only | isotonic | 0.869 | 0.873 | 0.864 | 0.882 | 0.894 | 0.818 | 0.888 | 1.46 | 0.132 |
| Assay-cond. (deployment) | none | 0.881 | 0.907 | 0.832 | 0.866 | 0.855 | 0.798 | 0.860 | 1.59 | 0.133 |
| Assay-cond. (deployment) | Platt | 0.867 | 0.894 | 0.827 | 0.855 | 0.858 | 0.778 | 0.857 | 1.62 | 0.136 |
| Assay-cond. (deployment) | isotonic | 0.856 | 0.873 | 0.832 | 0.850 | 0.877 | 0.764 | 0.864 | 1.53 | 0.142 |
| Assay-cond. (legacy true assay) | none | 0.882 | 0.909 | 0.832 | 0.866 | 0.855 | 0.798 | 0.860 | 1.59 | 0.132 |

AUPRC must be read against the 60.4 % active rate: a random classifier scores
0.604, so 0.91 is a real but not spectacular margin. EF@10 % is likewise capped
at 1/0.604 = 1.66, so the observed 1.46–1.62 is close to the ceiling and carries
little information.

**Deployment vs historical evaluation.** Legacy minus deployment AUROC is
+0.0015 (95 % CI −0.0001 to +0.0031). The historical evaluation mode is
therefore *not* shown to have materially overstated deployed performance. The
assay-information mismatch is a genuine methodological defect, but on this
dataset its numerical impact is negligible, and I am not claiming an
overestimate that the data do not support.

**Structure-only vs assay-conditioned.** On identical outer folds, assay-conditioned
minus structure-only AUROC is −0.0089 (95 % CI −0.020 to +0.0013, p = 0.088).
The assay one-hot does not improve outer-test prediction; the point estimate
slightly favours dropping it. This is consistent with the earlier SI observation
that the five assay contexts give near-identical probabilities for clinical
molecules — the encoding mostly describes the training corpus, not the compound.

**Fold variability.** Structure-only uncalibrated per-fold AUROC: 0.906, 0.913,
0.915, 0.824, 0.908 (mean 0.893, SD 0.039). Fold 3 is consistently the hardest
for every variant, so the spread reflects scaffold difficulty rather than one
model's instability. Fold sizes are 102–103 molecules with active rates 0.52–0.68.

**Calibration.** Neither post-hoc method improves the Brier score: 0.126
uncalibrated, 0.131 Platt, 0.132 isotonic for the structure-only model, and the
same ordering for the assay-conditioned model. Across the 3 inner repeats the
Brier SD is ±0.025, i.e. larger than the differences between methods, so no
calibration method can be declared superior here. The uncalibrated reliability
curve is mildly over-confident at both tails (bin 0.0–0.1: predicted 0.025,
observed 0.123; bin 0.9–1.0: predicted 0.972, observed 0.902) and well behaved in
between. Isotonic on ~120 calibration molecules collapses the output to ~55
distinct values and produces a non-monotone reliability curve in the 0.3–0.5
region — exactly the instability that made the historical in-sample calibrator
look so convincing.

Conclusion on calibration: with 513 molecules there is no evidence that post-hoc
calibration helps, so the uncalibrated model is retained and the limitation is
reported rather than papered over.

**Threshold 0.5.** For the structure-only uncalibrated model the prespecified
sweep gives accuracy 0.842 / 0.850 / 0.850 / 0.834 / 0.823 and F1 0.873 / 0.877 /
0.876 / 0.860 / 0.847 at 0.3 / 0.4 / 0.5 / 0.6 / 0.7. The 0.4–0.5 region is a flat
optimum, so keeping 0.5 is defensible — but as a conventional operating point
with a known 81 % specificity and 87 % recall, not as a calibrated probability.

**Comparison with the archived numbers.** Historical mean-of-folds AUROC 0.8934
and AUPRC 0.9143 came from the legacy evaluation mode; their honest counterparts
here are the legacy variant's 0.8831 and 0.9087 (mean of folds), and the
structure-only mean of folds is 0.8932. Base-model discrimination therefore
reproduces almost exactly — the historical AUROC claim survives. The historical
Brier score of 0.1284 has no fair counterpart: it was computed through the
in-sample isotonic calibrator, so it is a resubstitution figure. The honest
value for the recommended model is 0.126 at pooled outer-test level, which
happens to be similar, but the two are not the same quantity.

---

## 5. Recommendation for the next stage

Use `results/nlrp3_model_v2/nlrp3_v2_structure_only.joblib`, uncalibrated,
threshold 0.5, via `scripts/nlrp3_v2_predict.py`. Basis: equal-or-better AUROC,
AUPRC, accuracy and Brier than the assay-conditioned model on identical folds;
no evidence the assay encoding helps; one prediction per molecule instead of a
five-context maximum, which removes an arbitrary aggregation step and the
question of which assay contexts to pick; and the deployment rule is then
identical to the validated rule.

The assay-conditioned bundle is kept alongside it so that the historical pipeline
remains reproducible for comparison, not because it is recommended.

Open items deliberately left for a later stage: sensitivity of the `max`
molecule-level label to a stricter aggregation; applicability-domain analysis for
clinical compounds, none of which appear in any training assay; and re-scoring of
the 8319-compound library, which must be reported as a separate pool and must not
overwrite the frozen 1588.
