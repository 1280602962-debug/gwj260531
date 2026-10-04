# FIVE_SEED_READINESS

Stage 1 technical audit only. No AUROC.

## Count names (do not mix)

| Name | Value | Meaning |
|---|---:|---|
| `n_pair_expanded_seed42` | 1614 | pair × pocket rows at seed 42 |
| `n_unique_physical_jobs` | 1592 | unique physical Vina jobs (one receptor+ligand, seed-independent identity) |
| `n_pair_expanded_job_seed` | 8070 | 1614 × 5 |
| `n_physical_job_seed` | 7960 | 1592 × 5 |

`7960` is not unique physical jobs.

## Candidate master

Path: `reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/phase7_fiveseed_master.csv`

This file is a **candidate** until the closure below. It is not a new primary score authority for seed 42; seed-42 formal M0 remains `13_qa/official_primary_seed42_score_master_8pair.csv`.

## Rechecked counts (this audit)

- rows = **8070**
- seeds = {17, 29, 42, 71, 101}
- each seed = **1614** pair-expanded rows
- STATUS: **SUCCESS = 8039**, **TIMEOUT = 31**
- other statuses = 0

These match the expected check values. They were re-counted from the current file, not copied from an old report.

## `primary_score` field

`primary_score` is **not** a docking score.

Observed values: `{0, 1}`.

- `primary_score=1` on all 1614 seed-42 rows and on no other seed.
- `primary_score=0` on all rows with seed 17, 29, 71, or 101.

Interpretation: boolean flag “this row is the primary seed (42)”. Using it as affinity would be a scientific error.

## Candidate raw affinity field

`mode1_affinity`

- SUCCESS: 8039 / 8039 finite
- TIMEOUT: 0 / 31 finite (missing)

This matches TIMEOUT-as-missing. No imputation in the master.

Future M0 five-seed score, if a later phase is approved: `-mode1_affinity` (higher-better), same orientation as formal seed-42 M0.

## Closure against `06_vina_fiveseed/jobs`

- job directories = 8070
- master `job_id` minus directories = 0
- directories minus master = 0
- `status.json` missing = 0
- `status.json` text SUCCESS = 8039, TIMEOUT = 31

## Evidence status

`CANDIDATE_FIVE_SEED_MASTER_VERIFIED`

`data_readiness = READY_EXISTING_RAW`

No raw PDBQT reparse is required for a later M0 five-seed AUROC, provided that later phase uses `mode1_affinity` and the formal class rules.

## Seed-specific n change (document, do not impute)

Per-seed directional n is allowed to differ from seed42 when a job is TIMEOUT. This is missingness, not a reason to average seeds or drop the pair.

Documented example that later compute must report:

- pair `AChE/BChE`, ligand `AB_046`, seed17, pdb `4BDS`: `TIMEOUT` / `JOB_TIMEOUT`
- same ligand × pocket at seed42: `SUCCESS` (`mode1_affinity = -8.556`)
- the other AChE/BChE pocket `4EY7` is SUCCESS at both seed17 and seed42

Do not impute the missing seed17 4BDS score. Do not reuse the seed42 score for seed17.

## Five-seed summary quantiles (TYPE B, before compute)

Per pair, `median`, `min`, `max`, and `IQR` are computed on the five seed **point estimates** of a named metric (`AUC_B`, `AUC_A`, or `summary_min`).

Fixed definition:

- require five finite point estimates;
- `P25, P75 = numpy.percentile(values, [25, 75], method="linear")`;
- `IQR = P75 − P25`.

If fewer than five finite point estimates exist, do not report median/min/max/IQR for that metric. Write the reason. Do not use `nanpercentile` or any other silent drop of missing seeds. Do not try alternative quantile definitions after seeing results.

This five-seed summary is not a B=10000 bootstrap. Bootstrap, if any, is recorded per module.

## Still required before later compute (not done now)

- Join master rows to the formal pair-observation population / class labels.
- Keep pair-expanded vs physical-job identities distinct.
- Report per-seed n, including the AB_046 seed17/4BDS TIMEOUT effect.
- Do not average five seed scores.
- Do not replace official seed-42 M0 with a five-seed consensus.
