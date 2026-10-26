# Current analysis baseline (2026-10-26)

This is the current entry for statistics, figures, tables and manuscript writing.
Branch: `cursor/v4-2-final-audited-results-20260930`.
Integration input: `7ef3b25fddb32b3d7089db90974fdfccb405847f`.
Primary score/statistics source: `dc028b3350e804cbf8f40acea08261ed4f89c4fa`.
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

`python Dual_Target_Docking/scripts/rerun_v4_2/jcim_manuscript_integrate.py --supporting-update`
updates supplementary tables and review copies from already accepted files.
The historical `sync_current_q110_results.py` is the original Q110 promotion
helper, not the current evidence-integration entry; run it only at its source
commit when reproducing that historical stage.

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

## Current provenance and CAP5 evidence

- `results/current/provenance_audit_20261023/`: 867 observations have both-end raw
  activity records; 855 have complete specified document identifiers, 12 do not.
  Frozen maxima and theta6 classes have zero mismatches. Source traceability
  does not establish a common assay platform or external independent validation.
- `results/current/chemical_diversity_cap5_20261023/`: supporting sensitivity,
  retaining 861 and excluding six observations. It does not replace the primary
  population, scores, folds, PRIMARY tables, Stage 2 tables, or original sampling rules.
- Primary population remains 867 eligible observations. Q110 remains 110.
- Seven historical selection procedures remain DOCUMENTED_NOT_REPLAYED.
- Figures and main tables retain the primary population. CAP5 is reported in
  Results R8, SI sections 10–11 and Tables S23–S31. Source evidence is in S19–S22.
- G/J/K remain unexecuted. Q110 alternative-receptor results remain unavailable.

The old 805-observation activity export remains a provenance source for reused
records, joined through the current 867-observation crosswalk. The current source
audit supersedes earlier statements that Q110 has no per-assay source export.
