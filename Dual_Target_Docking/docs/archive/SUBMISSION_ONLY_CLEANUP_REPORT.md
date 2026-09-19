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
# Submission-only repository cleanup

Branch: `cursor/methods-sentence-audit-c7cc`  
Scope: `Dual_Target_Docking/`  
Pre-cleanup HEAD: `d824eb6b`

No new analysis, no redocking, no change to canonical statistics or manuscript numbers.

## 1. Deleted categories

| Category | What was removed |
|----------|------------------|
| Historical artwork | `Fig8_*`, `FigS4_pocket_matched_forest`, `FigS7_*`, `FigS8_*`, `FigS_detectable_effect`, standalone `Fig1_C_chEMBL_supply`, `ARCHIVED_FIGURES.md`, `FIGURE_AUDIT_PR32.md`, `FIGURE_ASSET_NOTES.md`, `FIGURE_SHA256.json`, `README_PR32.md`, `figures/editable/*.pptx` |
| Dual-source manuscript drafts | `assemble_manuscript_en.py` / `assemble_manuscript_zh.py` and all section/draft/title-abstract files; old Introduction/Discussion refs |
| Historical freeze/audit docs | `WRITING_FREEZE_REPORT.md`, `ISSUE_REGISTER_WRITING_FREEZE.md`, `TIME_SPLIT_PROTOCOL_FREEZE.md`, `PIK3CA_MTOR_STRUCTURE_FREEZE.md`, SOPs, `PUBLIC_TARGET_PAIR_SELECTION_REPORT.md` |
| Superseded scripts | `scripts/freeze/*`, `export_publication_tables.py`, panel-builder scripts that would re-download ChEMBL |
| Abandoned target-pair trees | `egfr_her2_panel40_v0`, `egfr_her2_panel40_reprep_rdkit_v0`, `pik3ca_mtor_panel48_v0`, `stage_m_v0`, `census/`, `protocols/` |
| Historical analysis drivers | `data/*/scripts/`, `data/*/analysis/*.md` (historical / unavailable), `jcim_novelty_v0/analysis/`, `MASTER_RESULTS_TABLE.csv` |
| Non-current bench/J0 files | `jcim_bench_v0` except `assembled_AChE_BChE.csv`; `jcim_j0j1_v0` except `j0_strict_label_supply.csv` |
| Old lock index | `manuscript_lock/` except `SI_TABLE_MERGE_MAP_v1.csv` |
| Old freeze outputs | tracked `remediation_outputs/freeze_audit/`; untracked docking dumps not committed |

Tracked deletions: **1392 files**. Production docking-score trees (universe / independent GNINA / structure-robust / holdout / panel120 / PM48 rdkit / PM110) were not removed.

Packaging-only script edits: stop writing `Fig1_C_chEMBL_supply`; drop it from `rebuild_submission_pack.py` STEMS. Figure 1 panel C is unchanged inside `Fig1_four_state_and_supply`.

`docs/LOCAL_RUN.md` and `figures/jcim_article/scripts/README.md` now describe the current zero-dock chain only.

## 2. Historical-looking files kept because the current chain reads them

- `figures/jcim_article/scripts/plot_jcim_article_figures_v3.py` — imported as `v` by `update_figures_pr32.py`
- `data/jcim_j0j1_v0/tables/j0_strict_label_supply.csv` — `v.load()`
- `data/jcim_bench_v0/tables/assembled_AChE_BChE.csv` — `v.load()` (Figure S1 TPSA)
- `data/jcim_novelty_v0/tables/*` (except `MASTER_RESULTS_TABLE.csv`) — canonical/figure/adjudication inputs and close-publication outputs
- `data/jcim_strengthen_t0t1_v0/tables/` — unified threshold, PM110 vs PM48, wrong-pocket, ECFP4
- `data/pik3ca_mtor_panel110_rdkit_v0/` — frozen PM110 scores behind Figure S2
- `data/jcim_supply_crossdb_v0/tables/` — BindingDB/PubChem count provenance
- `data/public_pair_selection/` — Table S9 pair-eligibility provenance
- `data/manuscript_lock/SI_TABLE_MERGE_MAP_v1.csv` — cited by SI
- `figures/jcim_article/scripts/inject_s2_boxes_from_json.py` — scanned by `check_current_chain.py`
- Production `*.pdbqt` under current score trees — deposited docking inputs, not abandoned panels

## 3. Current manuscript source

`docs/MANUSCRIPT_JCIM_EN.md`  
`docs/MANUSCRIPT_JCIM_ZH.md`

These are the only merged manuscripts. Section-level drafts were deleted so there is one publication text chain. `patch_publication_text.py` writes these files (extras are skipped when absent).

SHA-256 of both manuscripts is **unchanged** vs pre-cleanup `d824eb6b`.

## 4. Current SI source

`docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md`  
Chinese working copy: `docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md`

SHA-256 unchanged vs pre-cleanup.

## 5. Current canonical source

`results/canonical/`

`primary_summary_min.csv` (scheme B, θ = 6.0):

| Pair | summary_min | 95% CI |
|------|-------------|--------|
| EGFR/HER2 | 0.3237 | [0.1953, 0.4744] |
| JAK1/JAK2 | 0.5884 | [0.4482, 0.7158] |
| JAK1/TYK2 | 0.3649 | [0.2329, 0.5051] |
| PIK3CA/mTOR | 0.6921 | [0.4801, 0.8016] |
| AChE/BChE | 0.6058 | [0.4392, 0.7354] |
| F2/F10 | 0.3448 | [0.2157, 0.4819] |
| PPARG/PPARA | 0.6492 | [0.5111, 0.7460] |
| PPARA/PPARD | 0.4463 | [0.3008, 0.5898] |

All `results/canonical/*.csv` SHA-256 values match the pre-cleanup freeze.

## 6. Current figure source

Generator: `figures/jcim_article/scripts/update_figures_pr32.py`  
Shared loader/style: `plot_jcim_article_figures_v3.py`, `jcim_figure_style.py`  
Plotted values: `figures/jcim_article/plotted_values_postfix.json`

Nested plotted scientific values are identical to pre-cleanup except removal of the unused `fig1C_standalone` copy. Official Figure 1–5 / S1–S5 PDF binaries were restored from HEAD after a metadata-only regenerate.

## 7. Current submission pack

`submission_pack/` (100 files after rebuild):

- manuscripts EN/ZH + figure/table lock
- SI EN + ZH draft
- **Main:** Figure 1–5 (pdf/png/tif) + aliases `Figure1`–`Figure5`
- **SI:** Figure S1–S5 + aliases `FigureS1`–`FigureS5`
- TOC graphic
- captions + `plotted_values_postfix.json`
- `tables/canonical/*.csv`
- EGFR corrected boxes `3POZ_box_corrected.json`, `3RCD_box_corrected.json`

No Fig1_C, Fig8, FigS6+, or historical drafts.

## 8. `check_current_chain.py`

**PASS**

Zero-dock chain after cleanup (no Vina/GNINA):

`build_current_score_master.py` → `compute_canonical_results.py` → `compute_leave_one_document.py` → `compute_class_chemistry.py` → `compute_detectable_effect.py` → `fit_ecfp4_models.py` → `compute_descriptor_baselines.py` → `close_publication_from_canonical.py` → `patch_publication_text.py` → `update_figures_pr32.py` → `rebuild_submission_pack.py` → `check_current_chain.py`

Figure generator: `PASS: 11 figures` (Figures 1–5, S1–S5, TOC).

## 9. Fig6+ / FigS6+ / TableS11+ publication-facing leftovers

No leftover artwork files named Fig6+, Fig8, FigS6+, Fig1_C, or TableS11.

Remaining text mentions are not typeset figures:

- `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`: “There is no main-text Figure 6”
- captions: historical tokens “are not current results”
- `docs/FINAL_PRE_SUBMISSION_AUDIT.md` (historical / unavailable) and this report (cleanup notes)
- internal JSON keys `fig6A` inside `plotted_values_postfix.json` (function names; stems are `FigS4_*`)

---

Scientific results: unchanged  
Canonical results: unchanged  
Manuscript numbers: unchanged  
Zero-dock chain: PASS  
Submission pack: current official set only  

REPOSITORY CLEANUP STATUS:  
SUBMISSION-ONLY TREE READY
