# ANALYSIS_FREEZE_JCIM

Document status: STAGE1_RULE_FREEZE
Written on branch: `cursor/v4-2-final-audited-results-20260930`
Worktree: `/tmp/v4_2_final_audited_release_20260930`
Baseline HEAD at Stage 0: `2503f61570305e98faf135d167976b079eb6023c`
This file freezes subsequent JCIM analyses. It does not compute them.

Two statuses are separate:

- `ANALYSIS_FREEZE_JCIM` judges whether rules, populations, metrics, and interpretation boundaries are locked.
- `DATA_READINESS` is per-module and is filled only after Stage 4 technical audit.

This document does not treat `analysis_freeze.yaml` as method authority. That file is `STALE_PROVENANCE_ONLY`.

## 1. Central scientific question

Molecular docking's ability to separate experimentally dual-active ligands from single-target-active ligands; whether that directional discrimination appears as different single-target contamination in dual-target joint ranking; and whether that judgment is stable to scoring method, ligand chemistry, stochastic search, receptor conformation, and experimental activity definition.

The manuscript must not be rewritten as:

- which scoring method is best;
- GNINA is better than Vina;
- RTMScore is better than Vina;
- which receptor is best;
- which seed is best;
- an ordinary new-library virtual-screening application paper.

## 2. Evidence layers

1. Eight-pair directional benchmark.
2. Why some systems succeed and others fail.
3. Whether directional differences translate into joint-ranking differences.
4. A later, separately frozen new-library application.

This freeze covers layers 1–3 only. Layer 4 is not executed, does not select a pair, does not list a shortlist, and does not rank candidates.

## 3. Formal pairs

Exactly these eight, including known failures:

- EGFR/HER2
- JAK1/JAK2
- JAK1/TYK2
- PIK3CA/mTOR
- AChE/BChE
- F2/F10
- PPARG/PPARA
- PPARA/PPARD

EGFR/HER2, AChE/BChE, and PPARG/PPARA must remain. Low AUROC, AUROC below 0.5, or wide CI is not a deletion criterion.

## 4. Formal authority

In this order:

1. `00_protocol/FORMAL_AUTHORITY_PATHS.yaml`
2. `00_protocol/FORMAL_METRICS_ANALYSIS_FREEZE_FINAL.yaml`
3. `00_protocol/SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml`
4. `00_protocol/ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml`

`analysis_freeze.yaml` is historical provenance only. Its old receptors, score masters, activity paths, and analysis scripts must not override the four files above.

Formal statistics resolve `PROJECT_ROOT` as the `Dual_Target_Docking` directory. They must not use `/tmp/pr39_fiveseed` as an input root.

## 5. Path-dependency classes

Absolute paths are classified as:

| Class | Meaning | Stage-1 consequence |
|---|---|---|
| `DENYLIST_MENTION` | String exists only to block a bad read | not a blocker |
| `RUNTIME_BINARY_LOCAL_ONLY` | Historical GNINA/RTM/Vina binary path | not a blocker for zero-docking analyses |
| `PROVENANCE_FILE_READ_NO_SCIENTIFIC_USE` | File opened, contents unused for scores/population/statistics | record only |
| `ANALYSIS_INPUT_DEPENDS_ON_HISTORICAL_ABS_PATH` | Formal statistics actually consume scientific data via a historical absolute path | `BLOCKED_BASELINE_DEPENDENCY` |

Recorded at Stage 1 (read-only):

- `scripts/rerun_v4_2/run_primary_metrics_v42.py`: `/tmp/pr39_fiveseed` is in `FORBIDDEN_READ_SUBSTR` → `DENYLIST_MENTION`.
- Same runner loads `13_qa/FORMAL_ANALYSIS_INPUT_MANIFEST.csv` and discards it (`_ = load_csv(...)`) → `PROVENANCE_FILE_READ_NO_SCIENTIFIC_USE`.
- `scripts/rerun_v4_2/ablation_config.py` `GNINA_*` / `RTM_*` (`/mnt/d/...`, `/home/gwj/...`) → `RUNTIME_BINARY_LOCAL_ONLY`.
- `scripts/rerun_v4_2/formal_metrics_lib.py` and the primary runner bind `PROJECT_ROOT` from `--project-root` or `Path(__file__).resolve().parents[2]`. No Stage-1 evidence of `ANALYSIS_INPUT_DEPENDS_ON_HISTORICAL_ABS_PATH`.

This freeze does not modify runners.

## 6. Methods M0–M3

| Method | Definition | Role |
|---|---|---|
| M0 | Vina seed=42 mode=1; score = `-vina affinity` / `-minimizedAffinity`; higher-better | full-text primary method |
| M1 | same Vina seed42 mode1 pose; GNINA CNNscore | fixed-pose scoring change |
| M1b | maximum GNINA CNNscore on the same saved Vina pose set | saved-pose reranking only; not independent GNINA docking |
| M2 | same Vina seed42 mode1 pose; RTMScore model1 column `M2_RTMScore` | fixed-pose scoring change |
| M3 | independent GNINA docking, seed42, formal mode1 CNNscore | whole docking-workflow difference vs M0 |

Forbidden: `max_saved_RTMScore` as formal M2; promoting M1/M1b/M2/M3 to a new primary after seeing results.

Ablation wording:

- M0 → M1: scoring change at a fixed pose.
- M1 → M1b: reranking within the saved pose set.
- M0 → M3: overall docking-workflow difference.
- M1b → M3 must not be called a pure search-algorithm contribution.

## 7. Directional definitions

Inherited TYPE A:

- dual vs A_only uses score_B.
- dual vs B_only uses score_A.
- `summary_min = min(AUC_B, AUC_A)`.
- Forbidden: `AUROC(min(score_A, score_B))`.
- AUROC is Mann–Whitney; ties contribute 0.5.
- `activity_eligible=0` is not `neither`.
- AB_040 is an alias of AB_046 and is not an independent observation.

M0 shared-dual universe (inherited):

- Dual: `class=dual`, `activity_eligible=1`, M0 score_A finite, M0 score_B finite.
- A_only (primary directional): `class=A_only`, `activity_eligible=1`, M0 score_B finite.
- B_only (primary directional): `class=B_only`, `activity_eligible=1`, M0 score_A finite.

## 8. Bootstrap

Inherited TYPE A:

- B = 10000
- seed = 271828
- `numpy.random.default_rng`
- reset the RNG at the start of every analysis unit
- CI = `numpy.percentile(..., [2.5, 97.5], method="linear")`
- class-stratified paired identical indices
- point estimate from the original sample, not the bootstrap mean
- no silent NaN/Inf drop
- independent CI-overlap tests are forbidden

## 9. Scope contraction

Only two analyses use all of M0–M3:

1. the already completed formal directional benchmark;
2. the dual-target joint-ranking challenge.

All other new sensitivities are M0-only: negative-class, wrong-pocket, chemistry, five-seed, alternative receptor, threshold.

Forbidden: a 5-method × negatives × seeds × receptors × thresholds result matrix.

## 10. Module A — existing M0–M3 benchmark

Status: `ALREADY_COMPLETED`.

Authority:

- `results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv`
- `results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv`

Do not recompute, redock, or edit these tables. Report all eight pairs.

Manuscript role: main-text layer 1.

## 11. Module B — negative-class sensitivity (M0 only)

Question: if the negative class changes from single-target-active ligands to `neither`, how much does the docking judgment change?

Direction B:

- primary: dual vs A_only using score_B
- support: dual vs neither using score_B
- `delta_negative_B = AUC(dual vs neither | B) − AUC(dual vs A_only | B)`

Direction A:

- primary: dual vs B_only using score_A
- support: dual vs neither using score_A
- `delta_negative_A = AUC(dual vs neither | A) − AUC(dual vs B_only | A)`

Do not invent `summary_min_neither` unless a later freeze adds it.

Dual population: the formal M0 shared-dual universe. Both directions share that dual set.

Negatives:

- dual vs A_only: A_only with M0 score_B finite
- dual vs neither on B: neither with `activity_eligible=1` and M0 score_B finite
- dual vs B_only: B_only with M0 score_A finite
- dual vs neither on A: neither with `activity_eligible=1` and M0 score_A finite

Bootstrap (TYPE B detail, aligned to the existing paired rule):

- one dual resample index shared by the two AUROCs being subtracted;
- A_only / B_only resampled independently;
- neither resampled independently;
- do not use `n < 10 → NA` or `n < 10 → drop`.
- if a class is empty, write NA and do not start bootstrap (inherited empty-class rule).
- report `n_dual`, `n_negative`, point estimate, and CI for every non-empty comparison, including PIK3CA/mTOR neither if n is very small. Small n lowers interpretive strength only.

Manuscript role: main-text layer 2.

## 12. Module C — correct vs wrong pocket (M0 only)

Question: is directional discrimination aligned with the pocket that carries the experimental difference?

| Task | correct pocket | wrong pocket |
|---|---|---|
| dual vs A_only | B | A |
| dual vs B_only | A | B |

Population is stricter than the primary directional A_only/B_only rule because both pocket scores are required:

- dual: formal M0 shared-dual universe (both scores finite);
- relevant negatives: `activity_eligible=1` and **both** M0 score_A and M0 score_B finite;
- correct and wrong AUCs use this identical ligand set;
- if this drops ligands relative to the primary directional table, record the n change and do not revise the primary table.

`delta_pocket = AUC_correct − AUC_wrong`.

Bootstrap: one dual index and one negative-class index per replicate; both pocket AUCs use those same draws.

Interpretation: supports or weakens pocket-specific discrimination. Forbidden: claiming a true binding mechanism, or claiming chemistry bias is excluded.

Manuscript role: main-text layer 2, with full tables in SI.

## 13. Module D — ligand chemistry (M0 directional tasks only)

See `CHEMISTRY_BASELINE_FREEZE.md`.

Blocking granularity is `pair × arm`, not the whole module:

- `CHEMISTRY_SMILES_PARSE_FAILURE`, `TRAIN_SINGLE_CLASS`, or `OOF_SINGLE_CLASS` block only that `pair × arm`;
- other executable units remain eligible;
- the module is not globally BLOCKED if at least one `pair × arm` is executable.

`TEST_SINGLE_CLASS` is a warning only. Formal chemistry AUROC, when later permitted, is one OOF-merged AUROC per `pair × arm`, never a training AUROC and never a per-fold AUROC.

Frozen fold coverage is incomplete for 10 formal `pair × arm` records (7 unique ligands). Do not resplit. `fold-unassigned` is only the technical absence of a frozen CV fold; it is not poor data quality, unreliable activity, or a reason to drop the ligand from the primary docking evaluation set. Those ligands are excluded only from chemistry OOF models. D1 description still includes them. All chemistry models, and any chemistry-vs-docking supporting comparison, use the same `CHEMISTRY_OOF_POPULATION` (fold-assigned, parsable SMILES, not alias, required M0 score finite). Primary docking populations stay unchanged. A chemistry-matched M0 AUROC, if computed, is `SUPPORTING_CHEMISTRY_MATCHED_M0` only.

Manuscript role: main-text layer 2.

## 14. Module E — five-seed sensitivity (M0 only)

Seeds: 17, 29, 42, 71, 101.

Count names (TYPE A):

- `n_pair_expanded_seed42 = 1614`
- `n_unique_physical_jobs = 1592`
- `n_pair_expanded_job_seed = 8070`
- `n_physical_job_seed = 7960`

`7960` is not “unique physical jobs”.

`13_qa/phase7_fiveseed_master.csv` is a candidate master, not authority, until Stage 4 closure. `primary_score` must not be treated as a docking score until its historical meaning is documented. The candidate raw affinity field is `mode1_affinity`. Future M0 score, if the master closes, is `-mode1_affinity` after SUCCESS/TIMEOUT provenance is confirmed. This Stage 1 file does not compute that AUROC.

Future pair×seed outputs: `AUC_B`, `AUC_A`, `summary_min`. Per pair: seed42, median, min, max, IQR. Forbidden: best seed, mean-of-five-scores as a new primary, t-tests, treating five seeds as five independent experiments.

Manuscript role: SI tables; short main-text robustness summary.

## 15. Module F — alternative receptor (M0 only)

Authority: `ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml`.

Pairs: PIK3CA/mTOR, AChE/BChE, F2/F10, PPARA/PPARD.

Substitutions only: 4L2Y, 4JT5, 4EY6, 1P0M, 3SHC, 2Y5F, 6KAX, 5U46.

One receptor at a time (`A2/B1` or `A1/B2`). Future outputs must include `delta_AUC_affected_direction`, `delta_AUC_unreplaced_direction` (internal negative control), and may include `delta_summary_min`. Do not report `summary_min` alone.

Forbidden: new PDBs, result-driven receptor swap, RMSD-driven swap, picking the better alternative.

Manuscript role: SI tables; short main-text robustness summary.

## 16. Module G — activity-threshold sensitivity

Primary θ remains 6.0.

Candidate sensitivities already listed in `FORMAL_METRICS_ANALYSIS_FREEZE_FINAL.yaml`: 5.5, 6.5, `strict_6.5_5.5`.

Two independent gates:

1. `GATE_THRESHOLD_PRESPEC` — were these sensitivities written before PRIMARY results?
2. `GATE_ACTIVITY_SOURCE_READ` — may a later phase grant `SENSITIVITY_ONLY_READ_PERMISSION` for frozen continuous pA/pB?

`ligand_activity_aggregate_v1.csv` remains in `forbidden_activity_regeneration`. Gate 1 CONFIRMED does not grant gate 2.

If both later pass: M0 only; AUROC plus a class-transition table. No live ChEMBL. No primary-label rewrite.

Manuscript role: SI, and only if both gates pass.

## 17. Module H — missingness

Output columns: pair, method, class, target, expected, valid, missing, missing_reason.

Prefer formal codes: `TIMEOUT`, `MISSING_NO_VINA_POSE`, `ATOM_MAPPING_FAILURE`, `TECHNICAL_FAIL`.

Imputation is forbidden.

Method comparisons keep `pairwise_common_population = required`. `all_method_intersection` remains `secondary_only`.

Manuscript role: SI / QA.

## 18. Module I — dual-target joint-ranking challenge

See `JOINT_RANKING_FREEZE.md`.

This is the only new primary analysis that uses all of M0–M3. It translates directional AUROC into candidate-list composition. It is not external validation, not independent replication, and not a natural-library prevalence simulation.

Connection to benchmark: within each pair, for M0→M1, M0→M1b, M0→M2, M0→M3, describe whether `delta_summary_min` and `delta_top10_single_target_fraction` move in the expected direction (higher `summary_min` with lower top-10 single-target fraction). Outputs are a paired/linked plot and concordant/discordant counts. Spearman or Pearson correlation is not computed.

Manuscript role: main-text layer 3.

## 19. Holdout

`holdout_membership_freeze.csv` freezes membership only. `activity_eligible=1` counts are not a complete labeled replication set.

This phase audits provenance only. Future role, if later approved: within-source reserved-sample check, M0 seed42 mode1 only. Not external validation.

## 20. Future new-library application

See `FUTURE_APPLICATION_SELECTION_RULES.md`. No pair is selected in this phase. Role = benchmark-informed application, not blind prospective validation.

## 21. Manuscript vs SI

Main text:

1. 8 pair × M0–M3 directional benchmark;
2. M0 negative-class effect;
3. M0 chemistry + wrong-pocket;
4. short M0 seed/receptor robustness summary;
5. M0–M3 joint-ranking challenge;
6. within-pair link from `summary_min` change to top-10 contamination change;
7. M0 holdout replication only if holdout later becomes ready.

SI: full five-seed and receptor tables, missingness, common-complete, descriptor notes, full wrong-pocket, threshold (if gates pass), holdout provenance, other QA.

## 22. Rule types

- TYPE A: preexisting frozen rule.
- TYPE B: new rule specified in this directory before any new scientific result is computed.
- TYPE C: conditional on provenance/readiness.

TYPE B must not be described as historical preregistration.

## 23. ANALYSIS_FREEZE_JCIM status

`ANALYSIS_FREEZE_JCIM = PASS`

This PASS means the scientific question, populations, parameters, statistics, prohibitions, and interpretation boundaries are locked. It does not mean every module is ready to compute.

Adversarial answers (Stage 5):

1. NO — no new rule drops EGFR/HER2, AChE/BChE, or PPARG/PPARA.
2. NO — 5.5 / 6.5 / strict 6.5/5.5 were already listed as `sensitivity_only` before this freeze; no new cutoff was added after seeing a new analysis.
3. NO — C=1.0, nBits=2048, liblinear, midrank, and minimax were fixed before those analyses exist.
4. NO — every module serves the frozen central question.
5. NO — joint ranking is defined as a translation, not independent validation.
6. NO — the ranking pool is named a selectivity-enriched challenge population.
7. NO — M0 remains primary; M1–M3 stay ablation plus the ranking challenge.
8. NO — no future-application pair is selected.
9. NO — the design contracted (M0-only sensitivities; Top 10% only; two chemistry models plus one increment).
10. NO — remaining execution details required by later compute are written in this directory (pair×arm chemistry blocks, wrong-pocket both-score population, ranking midrank/ties/random expectation, five-seed `mode1_affinity` vs `primary_score`). `CHEMISTRY_FOLD_UNASSIGNED_POLICY = FROZEN_BEFORE_COMPUTE` and is no longer an execution-time free parameter.

Data gaps are recorded in `DATA_READINESS_AUDIT.csv` and do not change this PASS.

## 24. Forbidden additions

No extra engines, receptors, metrics (EF1%, BEDROC, MCC, F1, accuracy, PR-AUC), consensus scores, live ChEMBL, MD-as-activity, or new compute modules justified as “more complete”. See `FORBIDDEN_POSTHOC_ACTIONS.md`.
