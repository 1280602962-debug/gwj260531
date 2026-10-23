# Current analysis baseline (2026-10-23)

This is the current entry for statistics, figures, tables and manuscript writing.
Branch: `cursor/v4-2-final-audited-results-20260930`.
Source commit: `dc028b3350e804cbf8f40acea08261ed4f89c4fa`.
Paths below are relative to `Dual_Target_Docking/`.

## Current population and authorities

PIK3CA/mTOR uses accepted Q110: 110 observations (dual/A_only/B_only/neither =
32/32/32/14). Other seven pairs retain frozen V4.2 rows.

- `CURRENT_RESULT_AUTHORITIES.json`: current machine-readable authority.
- `results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv`: 869 rows,
  including 867 eligible observations (248/256/252/111).
- `results/current/pair_ligand_mapping.csv`: 870 rows including AB_040 alias.
- `results/canonical/model_fold_assignments.csv`: other seven frozen folds plus Q110 folds.
- `results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv`: 120 current rows.
- `results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv`: 32 current paired deltas.
- `results/jcim_stage2/`: current B/C/D/D1/E/H/I tables; F has three pairs/six conditions.
- `analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv`: 1004/1000/994 directional records.
- `results/current/RESULT_SYNC_REPORT.json`: merge checks and population totals.
- `manuscript/`: synchronized text, tables, review PDFs and source-data package.

PIK3CA/mTOR Q110 primary and supporting source files are in
`extension_jcim_validation/pik3ca_mtor_q110_v42/`.
The strict 6.5/5.5 subset has 56 members and is SI only.
Q110 alternative-receptor calculations were not performed. The old PM48 F rows
are excluded from the current F table, figures and manuscript.

## Definitions

M0: Vina seed42 mode1, score = minus affinity. M1: GNINA CNNscore on the same pose.
M1b: maximum CNNscore on saved Vina poses (reranking). M2: mode1 RTMScore model1,
not maximum saved RTMScore. M3: independent GNINA seed42 mode1 CNNscore.
`summary_min = min(AUC_B, AUC_A)`, not AUROC of per-ligand minimum scores.
D_vs_A uses score_B; D_vs_B uses score_A. Activity-ineligible rows are not neither.
AB_040 is an alias of AB_046, not an observation.
Bootstrap definitions and docking parameters are unchanged.

## Reproduction and historical material

`python Dual_Target_Docking/scripts/rerun_v4_2/sync_current_q110_results.py`
rebuilds current CSVs from the immutable source commit and accepted Q110 files.
It checks that other-seven-pair rows are unchanged and does not run docking,
AUROC, bootstrap or model fitting. Run in a Git checkout containing the source commit.

Raw V4.2 receptor, box, ligand, score and protocol files under `reruns/` remain
provenance and score sources for the other seven pairs and the nine reused Q110
ligands. Their historical counts and old PM48 rows do not describe the current
population. Old freeze/acceptance/audit documents remain dated records, not
current-result entry points. Re-run the historical PM48 regression only at the
source commit. Its accepted report is retained under `regression_pm48/`.

Superseded PM48 derived rows, redundant regression result copies and obsolete
non-formal score summaries are removed from the current tree. Full historical
files remain available through Git at the source commit. Do not restore them as
current results. Do not import old AUROC files as current inputs.

## Evidence scope

The old ChEMBL source audit covers its original 805 observations, not all 867
current observations. Q110 membership provides two-end pChEMBL and ChEMBL IDs;
a full per-assay source export for Q110 is not present. Seven pairs use their
original frozen folds; Q110 uses its documented pre-score Murcko-group fold rule.
These are not external independent validation results. G/J/K remain unexecuted.
