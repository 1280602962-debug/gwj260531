# V4.2 formal analysis baseline

This is the **only** human-readable entry for later statistics, sensitivity, figures, tables, and writing.

## 1. Branch and freeze

- Branch: `cursor/v4-2-final-audited-results-20260930`
- Cleanup-before commit: `37355a50`
- Worktree used for this cleanup: `/tmp/v4_2_final_audited_release_20260930`
- Do not use `/tmp/pr39_fiveseed` as a future execution root.

## 2. Baseline status

`V4_2_ANALYSIS_BASELINE` is recorded in `reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/V4_2_ANALYSIS_BASELINE.md` after QA.

PRIMARY seed42 formal metrics exist.
Five-seed sensitivity has **not** been run.

## 3. Formal target pairs (exactly 8)

- EGFR/HER2
- JAK1/JAK2
- JAK1/TYK2
- PIK3CA/mTOR
- AChE/BChE
- F2/F10
- PPARG/PPARA
- PPARA/PPARD

## 4. Formal receptors (exactly 14)

3POZ, 3RCD, 6N7A, 8BXH, 3LXP, 4L23, 4JT6, 4EY7, 4BDS, 4UDW, 2JKH, 9V8H, 6LXA, 5U3Q

Authority: V4.2 `01_receptors` + `03_boxes` + `00_protocol/SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml`.

`00_protocol/receptor_registry.csv` still has `run_id=UNIFORM_RERUN_V4_1_20260921`. That string is historical metadata. It does **not** mean receptors should be re-selected from V4.1.

## 5. Ligand universe

- SDF = 785
- PDBQT = 785
- ID sets identical
- `02_ligands` is a real Git tree, not a symlink

## 6–11. Unique formal score authorities

| Method | Definition | Authority file |
|---|---|---|
| M0 | Vina seed=42 mode=1; higher-better = `-vina affinity` | `reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/official_primary_seed42_score_master_8pair.csv` |
| M1 | Same Vina seed42 mode1 pose; GNINA CNNscore | `.../09_secondary_scoring/GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv` |
| M1b | max CNNscore on the **same** saved Vina poses; rerank only; not independent docking | same verified M1 master (`M1b_CNNscore`) + verified pose-long |
| M2 | Same Vina seed42 mode1 pose; RTMScore model1 column `M2_RTMScore` | `.../09_secondary_scoring/RTMSCORE_VINA_POSE_RESCORE_MASTER.csv` |
| M3 | Independent GNINA dock seed42; topology-verified official rows | `.../10_gnina_independent_docking/GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv` |

M2 must **not** use `max_saved_RTMScore`.
M0 must **not** use seeds 17/29/71/101 or best-of-9.
M1b is not independent GNINA docking.

## 12–14. Population and published metrics

- Population: `reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv` (807 rows)
- Primary metrics: `results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv` (120 rows)
- Method deltas: `results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv` (32 rows)
- `summary_min = min(AUC_B, AUC_A)` inside a method. Forbidden: `AUROC(min(score_A, score_B))`
- dual vs A_only uses score_B; dual vs B_only uses score_A
- `activity_eligible=0` is not `neither`
- AB_040 is an alias of AB_046, not an independent observation

## 15. Historical / non-formal

See `reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/HISTORICAL_NON_FORMAL_FILES.md` and `HISTORICAL_PATH_MAPPING.csv`.

## 16. Stale metadata

- `00_protocol/analysis_freeze.yaml`: `STALE_PROVENANCE_ONLY` / `NOT_FORMAL_METHOD_AUTHORITY`
- `receptor_registry.csv` V4.1 `run_id`: historical only
- `13_qa/FORMAL_AUTHORITY_PATHS_RESOLVED.json` and `FORMAL_ANALYSIS_INPUT_MANIFEST.csv`: `LOCAL_ENVIRONMENT_PROVENANCE_ONLY` (original `/tmp/pr39_fiveseed` resolution). Not future execution entry.

## 17. Formal missingness (from final tables)

Count names are fixed:

- `n_pair_expanded_seed42 = 1614`
- `n_unique_physical_jobs = 1592`
- `n_pair_expanded_job_seed = 8070`
- `n_physical_job_seed = 7960` (= 1592 × 5)

### M0

- pair-expanded seed42: 1614 accounted; 1608 SUCCESS; 6 TIMEOUT
- unique physical jobs: 1592 expected; 1586 SUCCESS; 6 TIMEOUT

### M1

- unique physical jobs expected 1592
- SUCCESS 1569
- ATOM_MAPPING_FAILURE 17
- MISSING_NO_VINA_POSE 6

### M2

- unique physical jobs expected 1592
- SUCCESS 1586
- MISSING_NO_VINA_POSE 6

### M3

- unique physical jobs expected 1592
- SUCCESS 1586
- TIMEOUT 6

Final M3 is **not** 1590/2. That number is pre-topology-repair only.

## 18–19. Five-seed data vs sensitivity

- FIVE-SEED RAW DOCKING = PRESENT (seeds 17, 29, 42, 71, 101; `n_pair_expanded_job_seed=8070`)
- FIVE-SEED SENSITIVITY ANALYSIS = **NOT YET RUN**

Raw jobs existing does not mean sensitivity is complete.

## 20. Allowed for later formal analysis

- This README and `ANALYSIS_BASELINE_MANIFEST.csv`
- Formal freeze / authority YAMLs listed in the manifest
- Formal score tables listed above
- Formal population and PRIMARY metric CSVs
- `scripts/rerun_v4_2/formal_metrics_lib.py` and the primary/post QA runners
- V4.2 `01_receptors`, `03_boxes`, `02_ligands`
- Five-seed raw jobs **only when a later approved sensitivity protocol says so**

## 21. Forbidden for later formal analysis

- Any `HISTORICAL_NON_FORMAL_*` file
- `09_secondary_scoring/gnina_rescore/`
- `scripts/analysis/*` (absent on this branch; do not restore)
- `results/canonical/current_score_master.csv` (absent; do not restore)
- `analysis_freeze.yaml` as method authority
- V4.1 `run_id` in `receptor_registry.csv` as a receptor selector
- `/tmp/pr39_fiveseed/...`
- Historical AUROC tables as inputs
- `max_saved_RTMScore` as formal M2
- Pre-repair M3 1590/2 audit as final status

## 22. Local software paths

GNINA and RTMScore paths in `SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml` and `scripts/rerun_v4_2/ablation_config.py` (`/mnt/d/...`, `/home/gwj/...`) are `LOCAL_ENVIRONMENT_ONLY`. They record the original run. Formal statistics do not call those binaries.

## Path resolution

`00_protocol/FORMAL_AUTHORITY_PATHS.yaml` stores PROJECT_ROOT-relative paths.
`PROJECT_ROOT` is the `Dual_Target_Docking` directory, from `--project-root` or from the script file location. Do not rely on process cwd.
