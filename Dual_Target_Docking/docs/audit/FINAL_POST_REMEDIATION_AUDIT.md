# FINAL POST-REMEDIATION AUDIT

Date: 2026-09-19  
Baseline audit commit: `a87c846b006b6f83e84ddf2a3c5ea58b34982c36`  
Tree: Dual_Target_Docking on PR #39 (`cursor/scientific-freeze-c7cc`)  
Prior scan: `docs/audit/FINAL_FULL_PROJECT_AUDIT.md` (NOT_READY)  
This scan: fresh clean rebuild + `scripts/qa/run_final_forensic_audit.py` (not a reused verdict) + `scripts/qa/check_table23_and_stale_tokens.py`

No experiment was redesigned. Table 2 / Table 3 scientific values were not edited. Historical EGFR `0.3237` was not restored.

---

## 1. Executive verdict

**READY_FOR_FINAL_MANUSCRIPT_DRAFTING**

P0 = **0**  
actionable P1 = **0**  
Table 2 = **MATCH** (64/64 EN+ZH cells)  
Table 3 = **MATCH** (64/64 EN+ZH cells)  
EN/ZH/SI/pack numeric STALE_ERROR = **0**  
fresh clean rebuild = **PASS** (31/31 canonical scientific CSV identical)  
current path unresolved = **0**  
submission pack = **SYNC**  
PR #39 description = **CURRENT**  
primary replay = **MATCH** (192/192)

Current EGFR/HER2 truth (unchanged by this writing/provenance pass):

- uniform RDKit / ETKDGv3 / MMFF / Meeko
- production seed `20260727`
- primary n = 28 / 37 / 31
- D-vs-A = 0.6631
- D-vs-B = 0.3341
- `summary_min` = 0.3341 [0.1970, 0.4713]

---

## 2. Remediation closures

| Item | Status | Evidence |
|------|--------|----------|
| SI S6 prose/table contradiction | **CLOSED** | EN/ZH SI and submission-pack SI: AChE/BChE Δ = 0.177 [0.053, 0.291]; EGFR/HER2 Δ = 0.107 [0.006, 0.220]; both main-panel intervals exclude 0; other six main pairs and all scored holdouts include 0; no “matched pockets generally outperform” claim |
| Manuscript cluster CI stale | **CLOSED** | EN/ZH: scaffold [0.232, 0.637], document [0.117, 0.636]; both still exclude 0; source `results/canonical/cluster_bootstrap_sensitivity.csv` |
| SI TPSA stale footnote | **CLOSED** | `0.4389 → 0.439`, labeled as a descriptor AUROC, not Table 2 docking AUROC |
| Fig 5C fuzzy `five_seed_comparable` | **CLOSED** | plotted metadata uses `same_protocol_as_primary` / `same_membership_as_primary`; AChE = 1 / 0; y-values unchanged (AChE range 0.5527–0.6058, primary overlay 0.6058); caption cites n=95 / 94 / 88 and `results/qa/five_seed_fixed_membership_sensitivity.csv` |
| GNINA default historical path | **CLOSED** | `GNINA_SOURCES["EGFR/HER2"]` default is `data/egfr_her2_uniform_rdkit_v1/tables/gnina_dock_scores_EGFR_HER2.csv`; `current_gnina_sources()` returns that dict; PIK3CA/mTOR and JAK1/TYK2 unchanged |
| Stale 0.3237 current docs | **CLOSED** | `docs/PR39_FINAL_SCIENTIFIC_FREEZE_CHECK.md`, `docs/PR39_SCIENTIFIC_DATA_AUDIT.md`, `docs/SUBMISSION_ONLY_CLEANUP_REPORT.md` moved to `docs/archive/` with SUPERSEDED banners; current indexes no longer treat them as truth |
| Legacy EGFR box dual-authority | **CLOSED** | current files are `*_box_corrected.json`; `3POZ_box.json` / `3RCD_box.json` live only in `data/egfr_her2_panel120_v0/boxes/archive/` with superseded metadata; no redock |
| PR #39 description | **CURRENT** | https://github.com/1280602962-debug/gwj260531/pull/39 body now states uniform EGFR rebuild, n=28/37/31, 0.3341 [0.1970, 0.4713], and removes `different_box_realization` / `comparable_to_current_primary=0` / “Table 2 AUROCs are unchanged” |
| Duplicate processed master | **CLOSED (P2)** | authority is `results/canonical/current_score_master.csv`; pointer is `data/processed/CURRENT_SCORE_MASTER.md` |
| Table 3 surrounding prose (D-vs-neither CI / EF / Top-10) | **CLOSED** | EN/ZH + lock: 0.759 [0.551, 0.926]; EF 0.351; composition 1 / 4 / 6 / 0. Table 3 cells were already canonical and were not edited |
| Token scan 0.3237 / 0.554 / 0.354 / 0.4275 / old cluster CIs / 1 / 5 / 5 / 0 | **STALE_ERROR = 0** | `results/qa/stale_token_scan.csv`; remaining hits only in `docs/archive/` or prior-audit notes with explicit superseded context |

---

## 3. P0

None.

Independent replay (re-run 2026-09-19 21:06, not a reused verdict): **192/192 MATCH**. `new_issues` = 0.  
EGFR `summary_min` independent 0.3341013825 vs canonical 0.3341.  
Score-sample mismatches = 0.  
θ=6.0 primary class mismatches = 0.  
Second-scan `new_issues` = 0.

Table 2 / Table 3 cells were not rewritten. Clean rebuild vs `results/canonical`: **31/31 CSV files identical**.

---

## 4. Remaining true reproducibility boundaries (not P0, not actionable P1)

These stay in Methods / SI / protocol docs. They were not invented, completed, or upgraded:

- 14 receptor slots: original receptor-preparation method is **NOT_RECOVERABLE**
- PIK3CA / AChE / Track-B exact ligand PDBQT incomplete
- seven non-EGFR pairs are not a complete SMILES→dock-input archive
- Track-B activity kind = `chembl37_dump_panel`
- eight pairs are not independent biological replicates
- pointwise CIs have no multiplicity correction
- independent GNINA covers 3 pairs only
- no external docking validation (eligibility only)
- rigid single-receptor primary setup
- PIK3CA/mTOR neither n=4
- no systematic tautomer / protonation enumeration
- AChE five-seed available-case membership ≠ primary membership (documented; fixed-membership n=88 qualitative conclusion unchanged)
- scaffold parse fallback to `unparsed` / `acyclic` remains; folds were not rebuilt

---

## 5. Fresh clean rebuild

Directory: `/tmp/dual_target_final_remediation_rebuild`  
Env: `/tmp/dual_target_pr39_env` (Python 3.12.3; numpy 2.5.2; pandas 3.0.5; scipy 1.18.1; sklearn 1.9.0; RDKit 2026.03.5)

```
python3 scripts/check_analysis_env.py
python3 scripts/analysis/rebuild_freeze.py --outdir /tmp/dual_target_final_remediation_rebuild
python3 scripts/qa/verify_freeze_rebuild.py --outdir /tmp/dual_target_final_remediation_rebuild
python3 scripts/qa/compare_rebuild_to_canonical.py --rebuilt /tmp/dual_target_final_remediation_rebuild --canonical results/canonical
```

Results:

- `check_analysis_env.py`: PIN_MATCH
- `verify_freeze_rebuild.py`: PASS
- compare: **31/31 CSV identical** (header, row count, cell text)
- EGFR n_scored 28/37/31; `summary_min` 0.3341 [0.1970, 0.4713]
- EGFR MM CI excludes 0: [0.0058, 0.2201]
- AChE MM CI excludes 0: [0.0528, 0.2906]
- ECFP max |Δ| = 0.0234

Figure 5 PNG was regenerated for metadata/caption consistency. Numerical y-values (seed min/max/median/primary overlay) were not changed.

---

## 6. Submission pack

`scripts/analysis/rebuild_submission_pack.py` rewrote `submission_pack/` (74 files) from current docs/figures/canonical.

Confirmed in the pack:

- SI S6 corrected prose (EN + ZH)
- cluster CIs [0.232, 0.637] / [0.117, 0.636]
- TPSA footnote 0.4389 → 0.439
- Figure 5 caption / `plotted_values_postfix.json` flags
- 31/31 canonical CSV byte-identical to `results/canonical`

No hand-copied second version remains.

---

## 7. Second forensic scan gates

| Gate | Result |
|------|--------|
| P0 | 0 |
| remaining actionable P1 | 0 |
| Table 2 MATCH | 64/64 EN+ZH (`results/qa/table23_lock.csv`) |
| Table 3 MATCH | 64/64 EN+ZH |
| EN/ZH numeric stale errors | 0 (`results/qa/stale_token_scan.csv`) |
| SI S6 prose/table contradiction | CLOSED |
| cluster CI stale | CLOSED |
| TPSA stale footnote | CLOSED |
| Fig 5C fuzzy comparable flag | CLOSED |
| GNINA default historical path | CLOSED |
| stale 0.3237 current docs | CLOSED |
| PR #39 description | CURRENT |
| current broken paths | 0 |
| primary replay | MATCH (192/192) |
| canonical ↔ submission_pack | SYNC (31/31) |
| clean rebuild | PASS |

Non-blocking leftovers (same class as the first audit’s P2):

- Figure 3B displayed 0.023 vs raw 0.0234 (rounding)
- The token `9.505` still appears in SI/lock as a **forbidden historical RMSD**, with surrounding “not current” / “不得将历史” wording. The scanner’s English `historical` heuristic missed two of those contexts; they are not current RMSD claims.

---

## 8. What this phase did not do

- no new target pairs
- no AChE redock
- no 14-receptor re-prep
- no eight-pair GNINA
- no new receptor swaps
- no MD
- no new external docking
- no change of θ or PRIMARY_PAIRS
- no rewrite of Table 2 to chase a better number
- no restoration of EGFR 0.3237

---

## 9. Current authority

Write numbers from:

- `results/canonical/*.csv`
- `figures/jcim_article/plotted_values_postfix.json`
- this file

Do not write current EGFR numbers from:

- `docs/archive/PR39_SCIENTIFIC_DATA_AUDIT.md`
- `docs/archive/PR39_FINAL_SCIENTIFIC_FREEZE_CHECK.md`
- `docs/archive/SUBMISSION_ONLY_CLEANUP_REPORT.md`
- `data/jcim_independent_dock_v0/tables/gnina_dock_scores_EGFR_HER2.csv`
- `data/egfr_her2_panel120_v0/boxes/archive/*_box.json`

---

## 10. Verdict line

READY_FOR_FINAL_MANUSCRIPT_DRAFTING
