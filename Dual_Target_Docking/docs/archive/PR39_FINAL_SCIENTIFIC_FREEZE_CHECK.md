<!--
============================================================
SUPERSEDED — HISTORICAL SNAPSHOT
============================================================
This file is not current scientific truth.

Current EGFR/HER2 (uniform RDKit/ETKDGv3/MMFF/Meeko, seed 20260727):
  primary n = 28/37/31
  summary_min = 0.3341
  95% CI = [0.1970, 0.4713]

The historical EGFR value 0.3237 in this snapshot must not be cited as a current result.

Current authority:
  results/canonical/*.csv
  docs/audit/FINAL_FULL_PROJECT_AUDIT.md
============================================================
-->
# PR #39 final scientific freeze check

Scope: scientific/data/code provenance only. No manuscript-language edit. No new pairs, docking, GNINA, RTMScore, MD, or receptor calculations. θ = 6.0, bootstrap B = 2000, seed = 20260729, and primary Table 2 AUROCs were not changed. ECFP incremental max |Δ AUROC| remains **0.0234** (PPARA/PPARD D vs B).

Working branch: `cursor/scientific-freeze-c7cc`  
Clean rebuild directory: `/tmp/dual_target_final_rebuild`  
Pinned venv: `/tmp/dual_target_pr39_env`

---

## Current environment

`scripts/check_analysis_env.py` → **PIN_MATCH** (exit 0)

| Package | Pin / freeze |
|---------|----------------|
| Python | 3.12.3 |
| numpy | 2.5.2 |
| pandas | 3.0.5 |
| scipy | 1.18.1 |
| scikit-learn | 1.9.0 |
| rdkit | 2026.03.5 (`importlib.metadata` 2026.3.5) |

`requirements-analysis.txt` was not edited.

---

## Clean rebuild PASS

```
python3 scripts/analysis/rebuild_freeze.py --outdir /tmp/dual_target_final_rebuild
python3 scripts/qa/verify_freeze_rebuild.py --outdir /tmp/dual_target_final_rebuild
```

`verify_freeze_rebuild.py` → **result: PASS**

Independent replay from the rebuilt master matched Table 2 directional AUROC, scheme-B CIs, ranking, ECFP OOF Δ, and deposited RMSD source copy. Primary Table 2 cells were identical to the pre-rebuild canonical file (no calculation error). ECFP max |Δ| replayed as 0.0234 at PPARA/PPARD D vs B.

---

## Canonical / rebuild equality PASS

```
python3 scripts/qa/compare_rebuild_to_canonical.py \
  --rebuilt /tmp/dual_target_final_rebuild \
  --canonical results/canonical
```

**PASS: 29/29 CSV files identical** (filename set, header, row count, cell text). Hash-free. `score_master_migration_diff.md` is not part of the CSV gate.

`python3 scripts/qa/check_current_chain.py` → **PASS**

`git diff --exit-code -- data/processed/activity_adjudication` → **unchanged**

---

## Table pack 29/29 PASS

Every `results/canonical/*.csv` matches `submission_pack/tables/canonical/` by `cmp` after pack rebuild.

Absent, as required (historical / unavailable; must be absent):

- `results/freeze_rebuild/`
- `figures/jcim_article/plotted_values.json`
- leftover publication ECFP 0.0112 CSV as a current source

---

## Figure source validation PASS

Official writer: `figures/jcim_article/scripts/update_figures_pr32.py`  
Helper only: `figures/jcim_article/scripts/plot_jcim_article_figures_v3.py` (legacy fig1–fig8 builders archived under `figures/jcim_article/scripts/archive/`)  
Plotted-value file: `figures/jcim_article/plotted_values_postfix.json` only

Figure 5C reads `results/canonical/five_seed_summary_min.csv`. EGFR/HER2 uses corrected-box five-seed Vina scores and current activity-eligible labels (`comparable_to_current_primary=1`). Production seed 20260727 matches Table 2 `summary_min` 0.3237 (n = 28/37/32/12). The five-seed range is 0.3237–0.3504. Figure 5C overlays the current primary seed for all eight pairs.

---

## Activity provenance summary

Labels were not changed. Added `activity_source_kind` / `activity_source_file` on `current_score_master`.

| pair | activity_source_kind | n_main | n_complete_case | n_primary_used |
|------|----------------------|-------:|----------------:|---------------:|
| EGFR/HER2 | chembl_assay_adjudicated | 110 | 110 | 109 |
| PIK3CA/mTOR | chembl_assay_adjudicated | 48 | 48 | 48 |
| AChE/BChE | chembl_assay_adjudicated | 95 | 95 | 94 |
| AChE/BChE | panel_pchembl_no_audit_rows | 5 | 1 | 1 |
| JAK1/JAK2 | chembl37_dump_panel | 110 | 110 | 110 |
| JAK1/TYK2 | chembl37_dump_panel | 110 | 109 | 109 |
| F2/F10 | chembl37_dump_panel | 110 | 107 | 107 |
| PPARG/PPARA | chembl37_dump_panel | 110 | 109 | 109 |
| PPARA/PPARD | chembl37_dump_panel | 110 | 110 | 110 |

EGFR/HER2 `EH120_059` and AChE/BChE `AB_087`: `unresolved_missing_arm`, activity_eligible=0, **not in primary**.

AChE fallback ligands `AB_001`, `AB_053`, `AB_054`, `AB_056`, `AB_097` are `panel_pchembl_no_audit_rows`. Only **AB_056** is complete-case and in primary.

---

## AChE AB_056 sensitivity (QA only; Table 2 not replaced)

| setting | D_vs_A | D_vs_B | summary_min |
|---------|--------|--------|-------------|
| current primary | 0.6524 | 0.6058 | 0.6058 |
| exclude AB_056 | 0.6504 | 0.6058 | 0.6058 |

The panel-pChEMBL fallback does not change `summary_min`.

---

## Active science pipeline has no SHA/hash gate

```
grep -RInE 'hashlib|sha256|checksum|git_head|rev-parse' \
  scripts figures/jcim_article/scripts
```

Hits are explanatory comments only (`Does not use checksums`, `Does not use SHA/hash/... as a PASS/FAIL gate`). `rebuild_freeze.py` no longer writes `git rev-parse HEAD`. Archive/data checksum manifests remain as archive integrity files and are not scientific PASS/FAIL.

---

## Remaining scientific limitations

- EGFR/HER2 five-seed Vina uses the corrected cognate-heavy-atom boxes and current activity-eligible labels (n = 28/37/32/12). Production seed `summary_min` 0.3237 matches Table 2. Per-ligand scores: `data/jcim_multiseed_v0/tables/scores_vina_mode1_EGFR_corrected_box_fiveseed.csv` (restored from commit `40edc431`; no new docking in this round).
- Track-B five pairs inherit ChEMBL37 dump panel pChEMBL without assay-row adjudication.
- AChE `AB_056` remains in Table 2 via panel pChEMBL with no audit rows; excluding it does not move `summary_min`.
- JAK1/TYK2 document-cluster bootstrap remains `unresolved_mapping_unavailable`.
- Detectable-effect rows are Monte Carlo simulation (n_mc = 1000), not observed power.
- Independent GNINA n can be smaller than Vina n; RTM/CNN on PPARG/PPARA are same-pose rescores, not independent docking.
- ECFP GroupKFold `fold_id` is pinned to this sklearn/rdkit environment (930 rows; no scaffold split; three models share `fold_id`). It is not required to match an older sklearn allocation. A deterministic scaffold allocator was not introduced in this round.
- No new docking, GNINA, RTM, MD, pairs, or receptor calculations were added.

Current scientific audit (design freeze): `docs/PR39_SCIENTIFIC_DATA_AUDIT.md`  
Superseded old-environment report: `docs/archive/SCIENTIFIC_FREEZE_REPORT.md`
