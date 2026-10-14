# NLRP3 classifier v2 — corrected training and validation

Everything here was produced by

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project/scripts
python3 nlrp3_model_v2_validate.py          # ~2 min on 8 cores
```

The audit that motivated this rebuild, and the interpretation of the numbers, are
in `docs/NLRP3_MODEL_V2_AUDIT_20261014.md`. This directory is additive: no frozen
production artefact, screening pool, candidate list, docking result or manuscript
fact table was modified.

## Files

| File | Contents |
| --- | --- |
| `data_label_audit.json` | Recomputed record/compound/assay/scaffold counts, label rule, per-assay class balance, discordant molecules |
| `validation_report.json` | Full design description, per-fold leakage checks, all metrics, reliability curves, bootstrap comparisons, software versions |
| `performance_summary.csv` | Pooled outer-test metrics, one row per (variant, calibration) |
| `outer_test_predictions.csv` | Per-molecule outer-test predictions (fold, label, raw and calibrated probability) for every variant |
| `fold_assignments.csv` | Outer fold membership of each molecule |
| `inner_repeat_sensitivity.csv` | Metrics for each outer fold × inner calibration split × variant × calibration method |
| `inner_repeat_summary.csv` | Mean and SD of AUROC / Brier / F1 across inner repeats |
| `final_model_split_membership.csv` | Which molecules trained the final base models and which were reserved for calibration |
| `model_card.json` | Recommended model, hyper-parameters, prediction rule, threshold, held-out performance, software versions |
| `nlrp3_v2_structure_only.joblib` | **Recommended model.** Molecule-level XGBoost on Morgan2048 + 12 RDKit descriptors |
| `nlrp3_v2_assay_conditioned.joblib` | Assay-conditioned model (historical feature set), kept for comparison |
| `calibration_curves.png` | Reliability curves, three variants × three calibration methods |
| `probability_distributions.png` | Outer-test probability histograms split by true label |
| `fold_auroc_comparison.png` | Per-fold outer-test AUROC against the historical 0.8934 reference line |

Both `.joblib` bundles are dicts with keys `model`, `calibrators`
(`none` / `platt` / `isotonic`) and `prediction_rule`; the classes live in
`scripts/nlrp3_v2_lib.py`, which must be importable when loading them.

## Stage 2

`stage2/` and `nlrp3_v2_structure_only_n513_production.joblib` are the
follow-up. The 392-molecule bundle above is unchanged. The production model is
the 513-molecule file: same hyperparameters, no calibration, historical max
label. It scores 1525 of 8319 clinical compounds at probability 0.5 (Jaccard
0.670 with the frozen 1588). Interpretation, the 80/32 recount, candidate
retention and the 39-molecule docking gap are in
`docs/NLRP3_MODEL_V2_STAGE2_20261014.md`.

## Headline result

On 5-fold scaffold-grouped outer folds with calibration fitted on scaffolds the
base model never saw, the structure-only model reaches pooled outer-test AUROC
0.890, AUPRC 0.913 (active rate 0.604), Brier 0.126. The assay-conditioned model
is not better (ΔAUROC −0.009, 95 % CI −0.020 to +0.001), and neither Platt nor
isotonic calibration improves the Brier score at this sample size.

## Next stage: re-scoring the clinical library

Not run here, deliberately. When it is run:

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project/scripts
python3 nlrp3_v2_predict.py \
  --input ../data/repurposing/repurposing_manifest.csv \
  --model ../results/nlrp3_model_v2/nlrp3_v2_structure_only.joblib \
  --calibration none \
  --threshold 0.5 \
  --output ../results/nlrp3_model_v2/predictions/clinical_library_v2.csv
```

Rules for that run:

1. Write to a new path. `data/repurposing/screening/docking_pool_p05.csv`
   (n = 1588), the 12 nominated and 21 reserve candidates, all docking output and
   the PF-03882845 MD files stay exactly as they are.
2. Report the v2 pool as a separate, additional pool and report its overlap with
   the frozen 1588 as an observation. Do not reconcile one to the other.
3. Keep threshold 0.5 and `--calibration none`, the combination validated here.
   If a different operating point is wanted, justify it from the prespecified
   sweep in `validation_report.json`, not from which candidates survive.
4. Before interpreting the new pool, check applicability domain: no clinical
   compound appears in any training assay, and max-Tanimoto to the 513 training
   molecules should be reported alongside the probabilities.
