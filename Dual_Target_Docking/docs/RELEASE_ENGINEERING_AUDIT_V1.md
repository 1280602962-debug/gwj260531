# Release engineering / defensive-programming audit V1

Inspected scientific/content HEAD (read first; do not trust older audit commit fields):

`0642f839727d409581e20e47446990e823a5da32`

source_snapshot_commit (commit used to generate the submission pack):

`0642f839727d409581e20e47446990e823a5da32`

The packaging git commit is **not** recorded in this file or in `submission_pack/README.md`. Lock it with a Git tag / release metadata after `revision-validate` is green. Do not confuse source_snapshot_commit with current HEAD after this release-engineering change lands.

Scope: PR #38 (`cursor/methods-sentence-audit-c7cc`). No scientific-analysis definitions were changed. Vina / GNINA / RTM were not re-run.

Current CI entry (`.github/workflows/revision-validate.yml`) now runs, any step exit 1 fails the job:

```
python scripts/check_analysis_env.py
python scripts/audit/smoke_postfix_v4.py
python scripts/audit/validate_postfix_v4.py
python scripts/audit/check_canonical_checksums_v2.py
python scripts/audit/check_submission_pack_v4.py
```

Official freeze entry: `scripts/audit/freeze_submission_postfix_v4.py`
Official pack entry: `scripts/build_submission_pack_postfix_v4.py`

---

## BLOCKING

1. **GitHub Actions `revision-validate` has not yet been observed green on the packaging HEAD.**
   Local equivalents of the five CI steps are green on this worktree (see PASS). The previous failing workflow called `data/jcim_novelty_v0/scripts/validate_revision_v1.py`, which hard-codes pre-fix EGFR `summary_min` 0.4297, fixed-score delta 0.3783, 3POZ top-1 RMSD 9.505, and old AChE counts. That script is now LEGACY and exits 1. Until the updated workflow runs on GitHub, this item remains blocking.

`RELEASE PIPELINE READY` is **not** written.

---

## MAJOR

1. **Host-specific docking drivers still hard-code `/home/gwj` and `/mnt/d/...`.**
   These are as-run recipes under `data/*/scripts/` and `remediation_outputs/scripts/`. They are not CI, freeze, or pack dependencies. Official lookup is `scripts/check_docking_env.py` (`PATH` or `VINA_BIN` / `GNINA_BIN` / `RTMSCORE_PYTHON`). Heavy outputs stay locked by `data/manuscript_lock/HEAVY_OUTPUT_LOCK_V4.csv`.

2. **Figure regeneration is not an atomic directory replace.**
   `figures/jcim_article/scripts/update_figures_pr32.py` writes stems sequentially and requires Arial. CI does not regenerate artwork. Failure mid-run could leave mixed old/new PNGs; freeze therefore supports `--skip-figures` and `--dry-run`.

3. **Historical analysis scripts still use `assert` as a correctness gate.**
   Examples: `data/jcim_novelty_v0/scripts/eight_pair_ranking_operating_point_v1.py`, `review_statistics_sensitivity_v1.py`. They are not the current release pipeline. Canonical ranking/metrics are frozen CSVs checked by `smoke_postfix_v4.py` / `validate_postfix_v4.py` with explicit `fail()` + nonzero exit. The pre-fix `validate_revision_v1.py` asserts remain in the file but `main()` exits 1 before they run.

4. **Pinned analysis environment is COMPATIBLE, not PIN_MATCH, in this cloud image.**
   numpy 2.4.4 vs pin 2.5.2; matplotlib 3.11.2 vs 3.11.1; scikit-learn 1.9.1 vs 1.9.0; rdkit 2026.03.6 vs 2026.3.5; meeko/gemmi missing. Required packages import. CI installs `requirements-analysis.txt` with an unpinned numpy/pandas fallback if pins fail. `check_analysis_env.py` reports PIN_MATCH vs COMPATIBLE explicitly and exits 0 only when required packages work.

5. **Clean-room figure rebuild is skipped when Arial is absent.**
   `run_cleanroom_smoke_v4.py --skip-figures` is the operator path on ubuntu-latest / this image. Requesting figure rebuild without Arial is a hard FAIL (`update_figures_pr32.py`), not a silent PASS.

---

## MINOR

1. **Non-canonical leftover Figure 6 artwork** (`figures/jcim_article/Fig6_evidence_boundary.{png,tif,pdf}`).
   Not packed. Checksum scanner WARNs. FAIL if it appears under `submission_pack/`.

2. **Superseded lock documents remain in `docs/`** (`FIGURE_PANEL_LOCK_V2.md`, `STATISTICAL_LOCK_V1.md`, `CONTRIBUTION_LOCK_V1.md`, `MANUSCRIPT_LOCK_INVENTORY_V1.md`).
   WARN only. Current unique lock is `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`. They are not canonical checksum dependencies.

3. **`plotted_values.json` (pre-postfix schema) remains beside `plotted_values_postfix.json`.**
   Allowed as superseded companion; not a canonical dependency.

4. **`check_local_env.py` is a LEGACY pointer** that delegates to `check_analysis_env.py`. Prefer the split env scripts.

5. **Clean-room pip install of pinned `requirements-analysis.txt` may fail on this image.**
   The env gate then decides COMPATIBLE vs FAIL. That is an explicit compatibility result, not a silent PASS.

---

## PASS

1. **Pre-fix validator is abandoned, not rewritten as a historical test.**
   `data/jcim_novelty_v0/scripts/validate_revision_v1.py` is marked LEGACY / PRE-REMEDIATION. `main()` prints the refuse banner and exits 1. Pre-fix literals remain unused below the raise.

2. **New post-fix validator reads machine-readable files, not manuscript text.**
   `scripts/audit/validate_postfix_v4.py` + `scripts/audit/smoke_postfix_v4.py` ground truth: V4 lock, `post_fix_master_metrics.csv`, `unified_threshold_sensitivity_v2.csv`, ranking CSV, box JSON, RMSD table, GNINA table, `plotted_values_postfix.json`, pack hashes.

3. **Smoke (`smoke_postfix_v4.py`): 23 PASS in 0.01 s (< 1 minute).**
   - exactly 8 primary pairs
   - theta=6 primary labels
   - D/A → pocket B; D/B → pocket A
   - dual-positive at theta 6
   - mode1 primary score
   - `summary_min = min(two arms)`
   - fixed-score channel identical; EGFR fixed delta ≈ 0.4621
   - EGFR `summary_min` ≈ 0.3237
   - top-10% `k = ceil(0.1 n)`; top D+A+B+N = k; EF formula
   - ECFP scaffold leakage = 0
   - main/holdout exact ligand overlap = 0
   - AChE n_scored 27/26/28
   - EGFR matched CI includes 0; AChE matched CI excludes 0
   - corrected box JSON matches official 3POZ/3RCD
   - no active analysis input resolves under `data/_legacy_archive`
   - withdrawn PIK3CA/PIK3CB is not in the current primary pair set

4. **Canonical validation (`validate_postfix_v4.py`): 12 PASS.**
   Includes 98 PASS / 0 WARNING / 0 FAIL scientific audit pointer, 3POZ 1.019 Å / 3RCD 1.947 Å, GNINA EGFR 0.2645 / 0.7370, pack EN manuscript raw SHA256 parity.

5. **Checksum coverage (`check_canonical_checksums_v2.py`): 119 paths PASS.**
   Single allowlist: `scripts/audit/canonical_paths_v4.py`. Expected path set strictly equals current required path set. Missing required file → FAIL. Unexpected `submission_pack/` file → FAIL. Raw byte size + raw SHA256 compared. `normalized_sha256` stored separately and not mixed into `sha256_raw`. Coverage includes V4 lock, corrected EGFR box JSON, `post_fix_master_metrics.csv`, `plotted_values_postfix.json`, current manuscript/SI, current figures, submission pack, and `HEAVY_OUTPUT_LOCK_V4.csv`. Legacy/pre-fix files are not canonical dependencies.

6. **Submission pack hash parity (`check_submission_pack_v4.py`) PASS.**
   `submission_pack/RELEASE_MANIFEST.csv` fields: path, role, bytes, sha256_raw, source_path, source_sha256_raw, source_snapshot_commit. Pack copy SHA256 == source SHA256 for manuscript / SI / figures / source tables. README uses `source_snapshot_commit:` (not `commit SHA:`). Snapshot is unambiguous: `0642f839727d409581e20e47446990e823a5da32`.

7. **Legacy freeze/pack/checksum entry points refuse to run.**
   - `scripts/audit/freeze_submission_v1.py` → exit 1 (still mentions `audit_submission_five_rounds_v1.py` / `pack_submission_v1.py` only as dead code)
   - `scripts/audit/pack_submission_v1.py` → exit 1
   - `data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py` → exit 1 (WATCH still lists superseded V3 assets; unused)

8. **Current freeze does not call the old pipeline.**
   `freeze_submission_postfix_v4.py` calls: analysis env, post-fix smoke, post-fix validation, current figure generation (optional), V4 pack build, checksum `--write` then `--check`, pack hash-parity. `--source-snapshot` defaults to the existing pack README, never silent packaging HEAD. `--dry-run` / `--skip-figures` / `--skip-promote` supported. Local `--dry-run --skip-figures` completed without replacing canonical outputs.

9. **Environment split.**
   - `scripts/check_analysis_env.py`: Python >= 3.10 and pinned-or-compatible packages.
   - `scripts/check_docking_env.py`: PATH or `VINA_BIN` / `GNINA_BIN` / `RTMSCORE_PYTHON`. No `/home/gwj` or `/mnt/d/...` official lookup. Missing docking tools: exit 0 unless `--require`.

10. **Defensive programming on current publication/release scripts.**
    - Current validators use explicit `fail()` + `SystemExit(1)`, not assert-as-gate.
    - Required columns validated before access; duplicate IDs rejected; NaN/inf rejected for required numerics.
    - Allowed classes limited to dual / A_only / B_only / neither (D/A/B/N).
    - Pair set exact-match checked against the locked eight pairs.
    - Canonical required paths refuse `_legacy_archive` / `submission_pack_pre_v4_archive`.
    - Ranking tie-breaking is recorded as descending `score_mean` then ascending ligand ID.
    - Pack / checksum destructive writes support `--dry-run` or staging; pack uses staging → validate → `atomic_replace_dir`; failed freeze does not promote `submission_pack`.
    - `fail()` in pack build is `SystemExit`; staging directory is removed on `BaseException` and is gitignored.

11. **Heavy scientific outputs are SHA256-locked.**
    `HEAVY_OUTPUT_LOCK_V4.csv` records path, bytes, sha256_raw, and `do_not_redock_unless_input_or_protocol_hash_changes`. CI/freeze/clean-room do not invoke Vina/GNINA/RTM.

12. **Clean-room runner exists:** `scripts/audit/run_cleanroom_smoke_v4.py`.
    Fresh `git clone --local`, install `requirements-analysis` (optional), smoke, validate, optional figure rebuild from frozen CSVs, pack rebuild with locked `source_snapshot_commit`, checksum + pack parity, hash compare vs original pack, git-diff allowlist. Pre-commit dirty trees copy the working `Dual_Target_Docking` files and compare pre- vs post-rebuild status.

13. **No scientific-conclusion edits** in this audit. EGFR `summary_min` 0.3237, fixed delta 0.4621, AChE 27/26/28, 3POZ 1.019 Å remain the post-fix canonical values.

---

## Verdict

**NOT RELEASE PIPELINE READY**

Remaining gate: observe GitHub Actions `revision-validate` = PASS on PR #38 packaging HEAD.

Local evidence already green:

- smoke test green
- canonical checksum green (119 paths, raw SHA256)
- post-fix scientific audit green (98/0/0 pointer + validator)
- submission pack hash parity green
- no legacy file is a canonical dependency
- current source snapshot unambiguous (`0642f839727d409581e20e47446990e823a5da32`)
