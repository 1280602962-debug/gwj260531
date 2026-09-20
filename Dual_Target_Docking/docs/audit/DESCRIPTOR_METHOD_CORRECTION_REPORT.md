# Descriptor method correction report

Date: 2026-09-20

## 1. Old method

Outer scaffold GroupKFold on all ligands of a pair. Each fold selected the training-set best raw descriptor for D-vs-A and for D-vs-B. Test scores were the raw selected columns. Pooled AUROC concatenated those raw values across folds.

## 2. Why invalid

TPSA, cLogP, MW, and heavy-atom count do not share a scale. Concatenating them and computing one AUROC is not a probability or calibrated decision score. The pooled statistic had no common meaning.

## 3. Implementation bug

`nested_oof` called `directional_smin(val_da[te], te_cls)`. That function returns both directional AUROCs from one array. Fold-level `test_D_vs_B` therefore used the D-vs-A selected descriptor whenever the two arms chose different names (21 of 40 pair×fold rows). Pooled `oof["selected_DvsB"]` was assigned correctly; the fold table was not. There was no inner CV.

Table S7a GNINA weaker-arm CI was copied from `summary_min_ci_*`. PIK3CA/mTOR independent GNINA therefore printed 0.633 [0.410, 0.769] for the weaker arm instead of the directional interval [0.419, 0.825].

## 4. Corrected method

For each pair × arm:

1. Independent outer scaffold GroupKFold on contrast ligands only.
2. Inner scaffold GroupKFold on outer-training data selects one prespecified descriptor `{tpsa, clogp, heavy, mw}`.
3. StandardScaler + univariate logistic regression are fit on the outer-training fold only.
4. Held-out output is P(dual) ∈ [0, 1].
5. Pooled AUROC uses those probabilities.
6. `summary_min` CI is a ligand-level shared-dual bootstrap of min(AUROC_DvsA*, AUROC_DvsB*) on the fixed OOF probabilities (B = 2000, seed = 20260729). It is not a re-run of nested selection.

## 5–6. Old vs corrected values

See `docs/audit/DESCRIPTOR_RESULT_IMPACT.md`. JAK1/JAK2 nested `summary_min` fell from 0.600 to 0.438. EGFR rose from 0.417 to 0.632. AChE/BChE fell from 0.742 to 0.650 [0.497, 0.796]. PPARG/PPARA nested `summary_min` is 0.358 [0.227, 0.499].

Descriptive full-panel screen is unchanged under the pinned RDKit 2026.3.5 environment.

## 7. Tables changed

- `descriptor_baselines.csv` — nested OOF columns and CIs; two cLogP descriptive cells at 4 decimals
- `descriptor_nested_scaffold_cv.csv` — now 80 pair×arm×fold units with explicit `pair/arm/outer_fold` keys
- `descriptor_nested_oof_predictions.csv` — new
- SI Table S5 — nested probability-scale table added
- SI Table S7a — GNINA class counts from complete cases; weaker-arm CI is directional
- `computational_robustness.csv` — added directional CIs and class counts; existing `summary_min` point/CI unchanged

## 8. Figures changed

Figure 2 and Figure S4: layout only (spacing, n=4 placement, S4 dagger). Scientific values unchanged. Figure 3 is ECFP4 and was not numerically changed.

## 9. Manuscript claims changed

Methods 2.6.1 and Results 3.3 now separate the descriptive screen from the nested P(dual) baseline and weaken the “simple properties explain docking” reading. See impact report.

## 10. Core paper conclusion

Unchanged. The paper remains a negative-class / ranking evaluation of dual-target docking. Receptor-free ECFP4 still carries substantial class information; docking still adds at most 0.023 AUROC.

## 11. Table 2 / Table 3

Unchanged. `TABLE2_AFFECTED = NO`. `TABLE3_AFFECTED = NO`. Those tables are produced by `compute_canonical_results.py` from docking scores and labels only.

## 12. ECFP4 conclusion

`ECFP4 NOT AFFECTED`. Separate script, separate per-arm GroupKFold, probability OOF. Max |Δ AUROC| remains 0.0234.
