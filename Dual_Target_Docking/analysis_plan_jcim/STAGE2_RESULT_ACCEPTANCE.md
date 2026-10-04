# STAGE2_RESULT_ACCEPTANCE

Formal Stage 2 outputs were produced once from the audited implementation and then checked by a separate read-only script. This file does not replace `STAGE2_IMPL_ACCEPTANCE.md` or `STAGE2_PREP_ACCEPTANCE.md`.

AUDITED_IMPLEMENTATION_BASE_SHA=3f3e0d26052178c38ccb43dec637683a60f97a9d
EXECUTION_SHA=3f3e0d26052178c38ccb43dec637683a60f97a9d
SECOND_STAGE_IMPLEMENTATION_READY=YES
SECOND_STAGE_COMPUTE_EXECUTED=YES
SECOND_STAGE_COMPUTE_STATUS=COMPLETED
SECOND_STAGE_RESULTS_ACCEPTED=YES
G_J_K_EXECUTED=NO

## Run record

| item | value |
|---|---|
| entry | `python3 Dual_Target_Docking/scripts/rerun_v4_2/jcim_stage2_compute.py --project-root Dual_Target_Docking --execute-compute` |
| log | `results/jcim_stage2_run_logs/stage2_execute.log` |
| project root | `/tmp/v4_2_final_audited_release_20260930/Dual_Target_Docking` |
| start | 2026-10-04 21:37:17 CST |
| end | 2026-10-04 21:40:18 CST |
| duration | 181 s (log timestamps) |
| exit code | 0 |
| python | 3.11.15 |
| numpy | 2.4.6 |
| rdkit | 2025.09.5 |
| scikit-learn | 1.9.0 |
| bootstrap | B=10000, seed=271828, numpy percentile linear, 2.5 and 97.5 |
| status at entry | `COMPUTE_STATUS=RUNNING` |
| script end line | `SECOND_STAGE_COMPUTE_EXECUTED=YES` |

`RUNNING` is the entry marker. `COMPUTE_EXECUTED=YES` is justified because Module B's paired bootstrap and the later modules actually ran. Passing the input load is not that marker: the loader can still exit in the integrity check before any statistic. A failure after the bootstrap had started would still be executed compute, even if the B table had not been finished.

A new file only means writing has started. B is treated as written because `module_B_negative_class.csv` was re-read with 16 data rows, the fixed schema, finite CI bounds, and `n_replicates=10000` on every row. Overall completion uses the exit code, the full output set, `RUN_MANIFEST.json`, and the independent audit together.

## Process, completeness, and faults

| layer | judgment |
|---|---|
| process | exit 0; manifest lists B C D E F H I as executable and G J K as blocked; D/E/H/I `bootstrap_B` is null; B/C/F `bootstrap_B` is 10000 |
| unit completeness | B 16, C 16, D 16, OOF 928, D1 938, E 40+24, F 8, H 320, I descriptive 104, I compare 32, I ligand ranks 9058 |
| planned NA / blocks | G, J, K not run. No chemistry unit has `blocking_reason` or `incomplete_reason`. E summaries are 24/24 with five finite seeds. M1b has eight supporting rows and zero STOP rows |
| new technical faults | none found by the independent audit |

`jcim_stage2_result_audit.py` does not import the Stage 2 compute functions. Result: **39/39 PASS**.

For the inputs that were prechecked as executable, the complete-success OOF target is 928, and the written OOF table has 928 rows. That number is the check for this input. It does not override the freeze: a later, explicitly recorded train-single-class block would be accepted only when the successful OOF count equals the member sum of the remaining executable units. A new missing prediction, fit failure, or unexpected block would be PARTIAL or NO. The sample is not reduced in order to call the run complete. An unexpected M1b STOP would be reported from the member evidence and would not be filled with a delta. None of those events occurred.

Test-single-class is an observation, not a new gate. This run has `units_with_test_single_warning=0`. The frozen rule remains: a train-single-class fold blocks that pair×arm; a test-single-class fold is still predicted, warned, and omitted from that fold's AUROC, while the unit result stays the full OOF. This run did not resplit, skip predictions, or block the chemistry module because of a test-single fold.

## What the audit checked

Statistical parameters were read from the written tables and the manifest. Computable B, C, and F units have `n_replicates=10000` and bootstrap started. The manifest seed is 271828, with exactly three `bootstrap_B: 10000` entries and null bootstrap fields for the modules that do not bootstrap. CI fields are finite and the lower bound is not above the upper bound. The audit does not require a CI to contain its point or to contain zero.

Signs were checked against the frozen definitions: B is neither−single and the single-class AUROC matches PRIMARY M0; C is correct−wrong, with `extra_drop` taken from member sets (AChE D_vs_A stays 25→25); D increment is ECFP4+M0−ECFP4, recomputed from the OOF columns. The paired bootstrap implementation was not rerun at B=10000. It is accepted from the already audited code and the artificial-path evidence.

E has 40 seed rows and 24 summary rows. Seed scores used in the independent check are −mode1 affinity. Seed 42 matches formal M0. AB_046 seed 17, receptor 4BDS, is listed as TIMEOUT and is not imputed. Each summary uses the linear percentile only when five points are finite; all 24 summaries met that condition. Fewer than five would have been NA.

I was not accepted from row counts alone. The audit rebuilt common members, average midranks, worst rank, the two tie-breaks, k=ceil(0.10×n), the actual Top-k members, and the three class counts, then compared them with the ligand and descriptive tables. PRIMARY deltas and the four concordance states were checked separately. M1b is supporting only: eight rows, `bootstrap_ci_added=NO`, no new CI, and no unexpected STOP.

D's three OOF prediction columns are complete, finite, and on the same members. Supporting M0 uses that same layer-3 population. D1 remains the labeled layer (938), including AB_001 and AB_053 on AChE D_vs_A. An independent recount of the three layers is 938 / 934 / 928.

## Warnings

| class | count | action |
|---|---|---|
| RDKit MorganGenerator deprecation | 1866 log lines | recorded; run continued. Fingerprint call was not changed |
| scikit-learn `FutureWarning` on `penalty` | 240 log lines | recorded separately from the RDKit deprecation; run continued. Solver, `max_iter`, and the estimator interface were not changed |
| `ConvergenceWarning` | 0 | none to list |
| fit exception, non-finite prediction, incomplete OOF | 0 observed | would have been a technical fault, not a warning to silence |

## Not verified, and not treated as faults

- The 10000 bootstrap replicates were not repeated, so bit-level reproducibility of those draws is not re-proven here.
- The formal CSV does not contain `scaler_probe`. Train-fold scaling of the added M0 column remains the artificial `module_d` evidence from the implementation check.
- Low AUROC, negative deltas, wide intervals, and intervals that contain zero are results, not technical failures.
- G, J, and K were not entered.
