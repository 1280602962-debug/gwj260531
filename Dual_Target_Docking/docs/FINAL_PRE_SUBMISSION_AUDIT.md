# Final pre-submission consistency audit

Branch: `cursor/methods-sentence-audit-c7cc`  
Directory: `Dual_Target_Docking/`  
Date: 2026-09-17  
Authoritative science: `results/canonical/`  
Writing freeze log: `docs/WRITING_FREEZE_REPORT.md` (previous-run snapshot; not a second Table 2 source)

This audit did not add analyses, docking, statistics, models, figures, or SI items. Independent replay of Table 2 from `current_score_master.csv` (activity-eligible complete-case) matched the writing-freeze Table 2 points to four decimals.

---

# 1. Overall verdict

**READY_TO_FREEZE**

After the deterministic publication-facing syncs listed in §10, `check_current_chain.py` PASS, primary Table 2 matches canonical, EGFR corrected scores are used throughout current analyses, AChE AB_056 is not mixed into dump-gated max/median, JAK1/TYK2 document-cluster is not presented as a current recomputation, independent GNINA is distinguished from CNN rescoring, 4JT6 is 3.60 Å, PIK3CA/PIK3CB is not a current primary pair, and old EGFR tokens 0.430 / 0.808 / 0.378 are not used as current results.

**SCIENTIFIC RESULTS CAN BE FROZEN.**  
**NO NEW DOCKING OR NEW SCIENTIFIC ANALYSIS IS REQUIRED BEFORE MANUSCRIPT LANGUAGE FINALIZATION.**

---

# 2. BLOCKERS

None

---

# 3. MAJOR ISSUES

None remaining after the canonical-sourced syncs in §10.

Found and fixed in this audit (were MAJOR before the sync):

1. **Manuscript GNINA weaker-arm CI ≠ SI / canonical.**  
   `docs/MANUSCRIPT_JCIM_EN.md` / ZH Results 3.5 had EGFR independent-GNINA dual-versus-B-only `0.265 [0.144, 0.401]`. SI Table S7a and `computational_robustness.csv` `summary_min` CI are `0.2645 [0.1417, 0.4018]` → `0.265 [0.142, 0.402]`. Now synced.

2. **Manuscript detectable-effect prose ≠ SI Table S4.**  
   Main text still said “at most 0.068 … at least 0.817 … PIK3CA/mTOR (0.453)”. Current `detectable_effect_simulation.csv` and Table S4 are max \(P(\mathrm{CI}\setminus 0.5\mid 0.60)=0.071\) (F2/F10), min non-PIK3CA at 0.75 = 0.835 (AChE/BChE), PIK3CA/mTOR 0.460. `patch_publication_text.py` only inserted this sentence if absent, so a stale sentence was never replaced. Now replaced from canonical.

`docs/WRITING_FREEZE_REPORT.md` still records the previous-run 0.068 / 0.817 / 0.453 snapshot. That file is a freeze log, not a typeset source. Do not copy those three tokens back into the manuscript.

---

# 4. MINOR ISSUES

Only items that remain worth a language/packaging pass. None change Table 2 or the scientific freeze.

1. **Abstract first use of “single-target selectives”.**  
   `docs/MANUSCRIPT_JCIM_EN.md` line 5 uses the phrase before §2.1 (line 31) states that A-only / B-only are pair-and-threshold experimental states, not proteome-wide selectivity. Optional one-clause qualifier in the abstract; do not rewrite the paper.

2. **Methods §2.5.3 and §2.6.3 both mention the binormal detectable-effect simulation.**  
   After collapsing the five-fold duplicate, one sentence remains in each subsection. Redundant, not contradictory.

3. **Five-seed archive still uses complete-case EGFR 28/38/32 and D vs A 0.6607.**  
   File: `data/jcim_multiseed_v0/tables/multiseed_auroc_by_seed_v2.csv` production seed 20260727. Figure 5C plots `summary_min` (EGFR 0.3237, weaker arm unchanged). Do not quote that CSV’s D vs A 0.6607 as current Table 2.

4. **Leftover artwork in `figures/jcim_article/` is not typeset.**  
   `Fig8_diagnostic_workflow`, `FigS4_pocket_matched_forest`, `FigS7_posthoc_diagnostics`, `FigS8_bindingdb_native_slice`, `FigS_detectable_effect` remain on disk. `rebuild_submission_pack.py` copies only Fig1–5 and S1–S5. `ARCHIVED_FIGURES.md` still refers to historical “Figure 6” panel roles.

5. **Internal JSON keys in `plotted_values_postfix.json` retain pre-postfix function names.**  
   `fig4A` is manuscript Figure 5 (GNINA); `fig5A` is manuscript Figure 4 (matched–mismatched); `fig6A` is Figure S4. The written PNG/PDF stems are correct (`Fig4_mismatched_pocket`, `Fig5_computational_realization`, `FigS4_label_source_robustness`).

6. **`jcim_figure_style.py` `PAIR_ORDER` still lists historical PIK3CA/PIK3CB.**  
   Comment says not to use it for the eight-row primary set. Current generator uses `PRIMARY_PAIRS`. `plot_jcim_article_figures_v3.py` still contains withdrawn-pair strings in unused/historical helpers.

7. **GNINA and RTMScore have no numbered software citations.**  
   Methods name GNINA 1.3.2 and RTMScore. References (1)–(22) are consecutive and otherwise intact. Optional citation add only; do not restyle the bibliography.

8. **Receptor protonation / tautomer / missing-loop reconstruction remain unrecorded.**  
   Methods already state this as a limitation. Do not invent a uniform protocol.

---

# 5. Stale-value scan

Publication-facing and current generators only. Historical archives are listed as HISTORICAL ONLY.

| Token | File | Line / locus | Context | Class | Fix? |
|---|---|---|---|---|---|
| `28 / 38 / 32 / 12` | `docs/MANUSCRIPT_JCIM_EN.md` Table 1 | 51 | EGFR construction quota | LEGITIMATE CONTEXT | No |
| `28 / 37 / 32` | same Table 1 / Table 2 | 51, 164 | primary n_scored | CURRENT | No |
| `0.430` / `0.808` / `0.378` | `docs/MANUSCRIPT_JCIM_EN.md`, SI, tables | — | not present as current AUROCs | CURRENT (absent) | No |
| `0.430` / `0.808` / `0.378` | `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md` | 99 | “pre-fix tokens are not current” | LEGITIMATE CONTEXT | No |
| `0.8088` / SI `0.809` | `fixed_score_negative_class_delta.csv`; SI Table S4 | JAK1/TYK2 dual vs neither | current fixed-score neither AUROC | CURRENT | No |
| `-0.2199` / SI `−0.220` | same; SI Table S4 | PIK3CA/mTOR Δ | current negative delta, not stale GNINA 0.2199 | CURRENT | No |
| `0.2199` as GNINA AUROC | publication-facing | — | absent | CURRENT (absent) | No |
| `0.6607` | `data/jcim_multiseed_v0/tables/multiseed_auroc_by_seed_v2.csv` | EGFR seed 20260727 D vs A | five-seed archive complete-case | LEGITIMATE CONTEXT (archive; not Table 2) | No (do not cite as Table 2) |
| `0.6607` | `docs/WRITING_FREEZE_REPORT.md` | 11 | “do not reintroduce” | LEGITIMATE CONTEXT | No |
| `9.505` | SI Figure S3 caption | 378 | “historical 9.505 Å is not a current result”; EGFR 3POZ top-1 1.019 Å | LEGITIMATE CONTEXT | No |
| `9.505` | `current_score_master.csv` AB_038 `score_B` | AChE ligand score | CURRENT score, not RMSD | CURRENT | No |
| `1/96` | `scripts/analysis/patch_publication_text.py` | old-string needles | generator search text | LEGITIMATE CONTEXT | No |
| `1/94` | manuscript Results 3.6; SI Table S3 prose | AChE dump-gated max/median | CURRENT | No |
| `5/109` | same | EGFR max/median | CURRENT | No |
| `6/110` | publication-facing | — | absent as current | CURRENT (absent) | No |
| `PIK3CA/PIK3CB` | SI figure-list prose; FIGURE_TABLE_LOCK | “historical … figures are not typeset” | LEGITIMATE CONTEXT | No |
| `PIK3CA/PIK3CB` | `figures/jcim_article/scripts/jcim_figure_style.py` `PAIR_ORDER`; `plot_jcim_article_figures_v3.py` | historical four-pair / guard strings | HISTORICAL ONLY in unused helpers | No for typeset files |
| `PIK3CA/PIK3CB` | `remediation_outputs/`, `data/jcim_bench_v0/`, `FIGURE_AUDIT_PR32.md`, `MASTER_RESULTS_TABLE.csv` | superseded | HISTORICAL ONLY | Do not edit |
| `0.068` / `0.817` / `0.453` | `docs/WRITING_FREEZE_REPORT.md` | detectable snapshot | freeze log | HISTORICAL ONLY | Do not copy into MS |
| `Figure 6` | `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md` | “There is no main-text Figure 6” | LEGITIMATE CONTEXT | No |
| `1/96`, `0.6607`, `28/38/32` as n_scored | `remediation_outputs/freeze_audit/*.md` | superseded banners | HISTORICAL ONLY | Do not edit |

Independent ECFP4 OOF probabilities of 0.430 / 0.378 / 0.808 in `ecfp4_oof_predictions.csv` are ligand-level model outputs, not publication AUROCs.

**STALE AND MUST FIX after this audit:** none remaining in typeset manuscript / SI / current tables / current figures.

---

# 6. Numerical consistency matrix

Display rounding is three decimals (half-up) from canonical. Status is after §10 syncs.

| Result | Canonical | Manuscript | SI | Figure / Table | Status |
|---|---|---|---|---|---|
| 1. Table 1 panel composition | quota EGFR 28/38/32/12; n_scored 28/37/32 (`class_counts.csv`) | Table 1 same | Table S3 n_scored 28/37/32 | Fig 1 supply, not n_scored | PASS |
| 2. Table 2 directional AUROC | EGFR 0.6564 / 0.3237; JAK1/JAK2 0.5884 / 0.7275; JAK1/TYK2 0.5751 / 0.3649; PIK3CA/mTOR 0.7143 / 0.6921; AChE 0.6524 / 0.6058; F2/F10 0.4133 / 0.3448; PPARG 0.6492 / 0.7061; PPARA 0.6465 / 0.4463; n as freeze | Table 2 three-decimal match | Table S3 `summary_min` match | Fig 2B from canonical | PASS (independent replay) |
| 3. Table 3 D-vs-neither + Top10 + EF | EGFR 0.7589 [0.5535, 0.9256], n_neither 12, 1/5/5/0, EF 0.3539; AChE 0.6587, 5/3/1/1, EF 1.7593; n_ranked 109 / 95 | Table 3 0.759 / 0.354 / 1.759 | Table S10 same operating points | Fig 2C/D | PASS |
| 4. EGFR fixed-score Δ | 0.4621 [0.2596, 0.6406] | 0.462 [0.260, 0.641] | Table S4 0.462 [0.260, 0.641] | Fig 2A | PASS |
| 5. JAK1/TYK2 fixed-score Δ | 0.4438 [0.2605, 0.6299] | 0.444 [0.261, 0.630] | Table S4 0.444 [0.261, 0.630] | Fig 2A | PASS |
| 6. Descriptor baselines | AChE TPSA 0.7422 / 0.8009; PIK3CA heavy `summary_min` 0.4630; 5/8 Δ CIs include 0 | 0.742 / 0.801; 0.463 | Table S5 match | Fig S1 | PASS |
| 7. ECFP4 AUROC | EGFR D/A 0.7992; max \|Δ\| 0.0112 JAK1/JAK2 D/A | “at most 0.011” | Table S5 0.011; EGFR 0.799 / 0.793 | Fig 3A | PASS |
| 8. ECFP4+docking incremental Δ | max \|Δ\| 0.0112 | 0.011 | 0.011 | Fig 3B `fig3B_max_abs` 0.0112 | PASS |
| 9. Matched-minus-mismatched | EGFR 0.0558 [−0.0368, 0.1551] includes 0; AChE 0.1770 [0.0528, 0.2906] excludes 0; only AChE main-panel excludes 0 | 0.056 [−0.037, 0.155]; 0.177 [0.053, 0.291] | Table S6 match | Fig 4A `excl` true only for AChE | PASS |
| 10. Unused-pool holdout | 7 pairs; all `mm_ci_excludes_zero=0`; EGFR missing | all seven include 0; EGFR † | Table S6 | Fig 4B | PASS |
| 11. GNINA independent | EGFR `summary_min` 0.2645 [0.1417, 0.4018], n 28/37/32/11; not 0.2199 | 0.265 [0.142, 0.402] | Table S7a 0.265 [0.142, 0.402]; 28/37/32/11 | Fig 5A | PASS (after sync) |
| 12. Receptor substitution | 4JPS 0.4861 [0.2639, 0.6944]; 5DXT 0.5046; 4JSX 0.6389; PIK3CA/mTOR only | 0.486 / 0.505 / 0.639 | Table S7 | Fig 5B | PASS |
| 13. Multiseed | Fig 5C `summary_min` ranges; EGFR primary 0.3237 | “limited numerical fluctuation”; not independent validation | S7 prose | Fig 5C | PASS for plotted weaker-arm; archive D/A 0.6607 not Table 2 |
| 14. Exhaustiveness | PIK3CA/mTOR E=16; others 8 | Methods / Table 1 | Table S1 / S2 | Fig S2 | PASS |
| 15. Cognate redocking | EGFR 3POZ CalcRMS top-1 1.019 Å; not 9.505 Å | 1.019 / 1.947 | Table S2; Fig S3 caption | Fig S3 | PASS |
| 16. θ sensitivity | EGFR θ=6.0 28/37/32, 0.3237 | Results 3.6 | Table S3 | Fig S4A | PASS |
| 17. Max vs median | EGFR 5/109; AChE 1/94 CHEMBL659; PPARA 1/110; AB_056 not in dump-gated 94 | 5/109; 1/94 | Table S3 prose | — | PASS |
| 18. Cluster bootstrap | EGFR document 0.4621 [0.1246, 0.6441]; JAK document `unresolved_mapping_unavailable` | “not recomputed”; do not copy old CI | Table S4 “not recomputed” | Fig S4B skips empty JAK document CI | PASS |
| 19. Detectable-effect | n_mc=1000; EGFR n=28/37/32; AChE n_neg 26/28; P(0.60) max 0.071; P(0.75) min non-PIK 0.835; PIK 0.460 | 0.071 / 0.835 / 0.460 | Table S4 cells | not a typeset SI figure | PASS (after sync) |
| 20. BindingDB eligibility | 0/8 pairs | eligibility screen, not external validation | Table S8; Fig S5 | Fig S5 | PASS |

Conflict items A–G:

- **A EGFR corrected scores:** ablation 3POZ/3RCD affinities match master; Table 2, fixed-score, matched–mismatched, θ, max/median, cluster, two-pocket, top-10%, GNINA Vina reference, manuscript, SI, figures, tables all use 0.6564 / 0.3237, not 0.430 / 0.808 / 0.378. EH120_059 `activity_eligible=0`.
- **B EGFR class count:** Table 1 quota 28/38/32/12 vs n_scored 28/37/32. No publication-facing use of 38 as primary A-only n_scored.
- **C AChE AB_056:** primary 27/26/28 eligible; AB_056 `panel_pchembl_no_audit_rows` kept; dump-gated max/median n=94, A-only 25; AB_087 unresolved, not in primary neither=14.
- **D JAK document cluster:** status `unresolved_mapping_unavailable`; SI “not recomputed”; manuscript states map/sqlite unavailable. Point and scaffold cluster are current.
- **E GNINA:** Methods 2.4.5 distinguishes CNN rescoring of Vina poses from independent pose generation. Independent EGFR is 0.2645, not 0.2199.
- **F 4JT6:** SI Table S2 3.60 Å; Methods: ≤3.5 Å is pair-supply screen, 4JT6 is the docked receptor and is not rewritten as 3.5 Å.
- **G PIK3CA/PIK3CB:** not in canonical eight pairs, Table 1–3, or typeset figures. Mentions are historical-not-typeset only.

Figures: main 1–5 only (no typeset Figure 6). SI S1–S5 only (no generated S11/S12 files; `P['figS11']` is a JSON dump key). Figure 4: only AChE main-panel CI excludes 0; EGFR includes 0; seven holdouts include 0. Figure 5: independent GNINA; receptor swap PIK3CA/mTOR only; five-seed is robustness.

---

# 7. Reproducibility check

Zero-dock chain from `WRITING_FREEZE_REPORT.md`, this audit, no docking / no ChEMBL download / no fabricated substitutes.

| Step | Status |
|---|---|
| `python3 scripts/analysis/adjudicate_activity_records.py` | PASS |
| `python3 scripts/analysis/build_current_score_master.py` | PASS (1226 rows; EGFR eligible 109 = 28/37/32/12; AChE eligible 95 = 27/26/28/14) |
| `python3 scripts/analysis/compute_canonical_results.py` | PASS |
| `python3 scripts/analysis/compute_leave_one_document.py` | PASS |
| `python3 scripts/analysis/compute_class_chemistry.py` | PASS |
| `python3 scripts/analysis/compute_detectable_effect.py` | PASS (~19 min; class sizes 28/37/32/12 etc.) |
| `python3 scripts/analysis/fit_ecfp4_models.py` | PASS (max \|Δ\| = 0.0112 JAK1/JAK2 D_vs_A) |
| `python3 scripts/analysis/compute_descriptor_baselines.py` | PASS |
| `python3 scripts/analysis/close_publication_from_canonical.py` | PASS |
| `python3 scripts/analysis/patch_publication_text.py` | PASS (after idempotent detectable insert + current GNINA / detectable-prose needles) |
| `python3 figures/jcim_article/scripts/update_figures_pr32.py --source-root .` | PASS (12 figures; Arial) |
| `python3 scripts/analysis/rebuild_submission_pack.py` | PASS (103 files) |
| `python3 scripts/qa/check_current_chain.py` | PASS |

Independent Table 2 replay from activity-eligible master: **PASS** for all eight pairs.

JAK1/TYK2 document-cluster mapping / ChEMBL 37 sqlite: **unavailable**; not faked.

---

# 8. Evidence-boundary audit

No sentence in the typeset manuscript requires a mandatory rewrite. Reviewed claims stay inside the allowed boundary.

| File | Line | Sentence | Assessment |
|---|---|---|---|
| `MANUSCRIPT_JCIM_EN.md` | 198 | “In this retrospective setting, the docking score did not add a stable ranking increment beyond a ligand-based model of the same experimental classes.” | In-bound: ECFP4 / scaffold-grouped setting; not “docking contains no structural information”. |
| | 206 | “An interval that includes zero does not establish that the receptor contributes no information” | In-bound. |
| | 232 | BindingDB search “is an eligibility screen, not external validation” | In-bound. |
| | 239 | Dual-vs-neither vs dual-vs-selectives “answer different questions” | In-bound. |
| | 251 | High retrospective AUROC “cannot, by its numerical value alone, be attributed to three-dimensional complementarity” | In-bound. |
| | 257 | Holdouts “did not consistently reproduce a matched-pocket advantage”; AChE is not generalized | In-bound. |
| | 263 | Independent GNINA: task difference can reappear; not a matched-pocket proof | In-bound. |
| | 269 | Retrospective AUROC ≠ prospective hit rate; eight pairs not a leaderboard | In-bound. |
| | 275 | Eight pairs “should not be treated as independent statistical replicates” | In-bound. |
| | 283 | Dual-vs-neither “does not establish discrimination against single-target selectives” | In-bound. |
| Abstract L5 | “single-target selectives that remain active at one target” | Borderline first-use; §2.1 later restricts the definition. Optional qualifier only. |

No “docking generally does not work”, “chemistry causes all performance”, “AChE proves specificity generally”, “eight pairs are independent replicates”, or “BindingDB 0 eligible pairs = external validation failure of docking”.

Methods vs code: ChEMBL 37, pChEMBL max, median sensitivity, θ=6.0, strict 6.5/5.5 as supply, panel sampling, scaffold caps, PDB IDs, cognate AABB+5 Å / min 20 Å box, Meeko, Vina 1.2.7, n_modes 9 / energy_range 3, E=8 except PIK3CA/mTOR 16, rank-1 `REMARK VINA RESULT` (MODEL 1), score sign \(S=-E\), `summary_min`, two-pocket mean, B=2000 seed 20260729, matched/mismatched, GroupKFold ECFP4, GNINA independent vs CNN, BindingDB eligibility — all stated. Receptor protonation is explicitly unrecorded. **No Methods/code BLOCKER.**

References (1)–(22) are consecutive. Vina 1.2.0 paper cited for Vina 1.2.7 software is acceptable. ChEMBL 37, BindingDB [16], PDB IDs in Table 1/S2 are present. No obvious DOI/year/page errors inspected.

---

# 9. Final freeze recommendation

**A. Re-run docking?** No.  
**B. Add scientific analysis?** No.  
**C. Freeze scientific results?** Yes.  
**D. Remaining work?** Language, abstract qualifier, optional software citations, leftover untyped figure files, and submission formatting only.

---

# 10. Modifications made in this audit

All from canonical / current definitions. No new analysis, no docking, no historical-archive edits.

| File | Before | After | Source |
|---|---|---|---|
| `scripts/analysis/patch_publication_text.py` | Non-idempotent replace appended another “A binormal detectable-effect simulation used…” on every run (EN/ZH) | `_ensure_once`: keep a single copy | generator bug; manuscript §2.6.3 had five copies |
| same | GNINA needle only matched `0.265 [0.148, 0.394]` | also matches current `0.265 [0.144, 0.401]` and writes canonical CI | `computational_robustness.csv` 0.2645 [0.1417, 0.4018] |
| same | detectable prose inserted only if absent | regex-replace existing 0.068/0.817/0.453 sentence | `detectable_effect_simulation.csv` |
| `docs/MANUSCRIPT_JCIM_EN.md` L132 | five copies of the §2.6.3 detectable sentence | one copy | same |
| `docs/MANUSCRIPT_JCIM_EN.md` L214 | `0.265 [0.144, 0.401]` | `0.265 [0.142, 0.402]` | GNINA `summary_min` CI |
| `docs/MANUSCRIPT_JCIM_EN.md` L158 | at most 0.068; at least 0.817; PIK 0.453 | 0.071; 0.835; 0.460 | detectable-effect CSV / Table S4 |
| `docs/MANUSCRIPT_JCIM_ZH.md` | same three issues | same three corrections | same sources |
| `docs/METHODS_SECTION_JCIM_EN_V1.md`, `METHODS_DRAFT_ZH_JCIM_V1.md` | duplicate §2.6.3 sentence | collapsed | same |
| `docs/RESULTS_SECTION_JCIM_EN_V1.md` | directional range 0.345–0.728; Figure 6 image (missing file) | 0.324–0.728; Figure S4 / S5; no Fig6 include | Table 2; postfix numbering |
| `docs/RESULTS_DRAFT_ZH_JCIM_V1.md` | same Figure 6 / 0.345 range | S4/S5; 0.324 range | same |
| `submission_pack/` | rebuilt | 103 files from current docs/figures/canonical | `rebuild_submission_pack.py` |
| `figures/jcim_article/` PNG/PDF/TIFF | regenerated from current CSVs | PASS 12 figures | `update_figures_pr32.py` |

Not modified: `remediation_outputs/`, `FIGURE_AUDIT_PR32.md`, `ARTICLE_ASSET_INDEX_v1.csv`, `MASTER_RESULTS_TABLE.csv`, docking outputs, ChEMBL dumps, five-seed archive CSV (left as robustness archive).

Independent confirmation (do not overwrite canonical with the user-supplied list; the list matched):

| Pair | n_scored | D vs A (B) | D vs B (A) | summary_min [95% CI] |
|---|---|---|---|---|
| EGFR/HER2 | 28/37/32 | 0.6564 | 0.3237 | 0.3237 [0.1953, 0.4744] |
| JAK1/JAK2 | 32/32/32 | 0.5884 | 0.7275 | 0.5884 [0.4482, 0.7158] |
| JAK1/TYK2 | 31/32/32 | 0.5751 | 0.3649 | 0.3649 [0.2329, 0.5051] |
| PIK3CA/mTOR | 18/14/12 | 0.7143 | 0.6921 | 0.6921 [0.4801, 0.8016] |
| AChE/BChE | 27/26/28 | 0.6524 | 0.6058 | 0.6058 [0.4392, 0.7354] |
| F2/F10 | 31/32/32 | 0.4133 | 0.3448 | 0.3448 [0.2157, 0.4819] |
| PPARG/PPARA | 32/31/32 | 0.6492 | 0.7061 | 0.6492 [0.5111, 0.7460] |
| PPARA/PPARD | 32/32/32 | 0.6465 | 0.4463 | 0.4463 [0.3008, 0.5898] |
