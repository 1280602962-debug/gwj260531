# FINAL REVIEWER RED-TEAM AUDIT

Date: 2026-09-18  
Branch: `cursor/scientific-freeze-c7cc` (worktree `/tmp/pr39_fiveseed`)  
Companion tables: `docs/SCIENTIFIC_PROVENANCE_MATRIX.csv`, `results/qa/*`

## Verdict

**NOT_READY_FOR_WRITING_FREEZE**

P0 count: **1 open** (EGFR deposited ligand-prep chain).  
P1 items have explicit evidence boundaries below.  
Item 15 fresh rebuild + promotion was **not** run as an acceptance test: science state is not determined until the CASE 2 EGFR campaign is finished and wired.

This is not a manuscript-quality objection. It is a traceability objection. Internal AUROC replay of the deposited score master can still match `results/canonical/`. That consistency does not prove the EGFR ligand 3D inputs.

---

## P0 MUST FIX BEFORE WRITING

### P0-1 EGFR deposited ligand inputs are not recoverable

- **question:** Can Table 2 / Fig 5A / Fig 5C EGFR numbers be traced to a proven ligand-prep method and ligand PDBQT?
- **evidence:** `docs/EGFR_LIGAND_PREP_PROVENANCE_AUDIT.md` CASE 2. Panel120 has no `ligands_pdbqt/`. LigPrep `panel_v0_40_ligprep-out.maegz` was never in git. `EGFR_SCORE_SOURCE` is still `data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv`. Uniform rebuild ligands exist (`data/egfr_her2_uniform_rdkit_v1/ligands_pdbqt`, 110/110) but Vina 1100 jobs are not finished and are not promoted.
- **risk:** Methods text that says “all pairs used RDKit ETKDGv3 + Meeko” is not true of the deposited EGFR scores. Keeping 0.3237 to avoid a writing pass would be protocol shopping.
- **action:** Finish uniform Vina five-seed + EGFR independent GNINA (timeout_skipped allowed at 600 s). Write `EGFR_UNIFORM_PREP_IMPACT.md` with real deltas. Point `EGFR_SCORE_SOURCE`, `FIVE_SEED_EGFR_CORRECTED_BOX_SCORES`, and `GNINA_SOURCES['EGFR/HER2']` at the new tables. Rebuild canonical. New numbers win.
- **status:** OPEN

---

## P1 MUST DOCUMENT BEFORE WRITING

### P1-1 Receptor prep method is not_recoverable for all 14 slots

- **question:** How were receptors prepared?
- **evidence:** `docs/RECEPTOR_PREP_PROVENANCE.md`. PDBQT files exist. Waters absent. 9V8H chain B peptide present as claimed. No PPW/pH/loop script.
- **risk:** Reviewer asks for a common protonation protocol; inventing one would be worse than stating the gap.
- **action:** Keep the manuscript limitation already present (Table S2 records PDB/cognate/box, not a reconstructed PPW protocol). Do not re-prepare receptors to fill text.
- **status:** DOCUMENTED

### P1-2 Track-B labels are dump-gated, not assay-adjudicated

- **question:** Are activity labels harmonized?
- **evidence:** `build_current_score_master.py` tags Track-B `chembl37_dump_panel`. Track-A EGFR/PIK3CA/AChE are `chembl_assay_adjudicated` except five AChE `panel_pchembl_no_audit_rows` (only AB_056 in primary). EH120_059 and AB_087 `unresolved_missing_arm` are out of primary.
- **risk:** Calling all eight pairs assay-adjudicated.
- **action:** Keep the existing freeze-check wording. Do not change θ=6.0.
- **status:** DOCUMENTED

### P1-3 AChE five-seed membership ≠ Table 2 membership

- **question:** Does changing sample membership affect seed sensitivity?
- **evidence:** `results/qa/five_seed_fixed_membership_sensitivity.csv`. AChE primary n=95; production-seed five-seed n=94; five-seed intersection n=88 (drop AB_027, AB_028, AB_040, AB_046, AB_055, AB_056, AB_071). Max |Δ summary_min| vs current per-seed = 0.0161. Other seven pairs: intersection = primary n, Δ=0.
- **risk:** Over-claiming that five-seed uses the identical Table 2 sample. `comparable_to_current_primary=1` is hardcoded even for AChE n=94 vs 95.
- **action:** No redock. Report Fig 5C AChE as per-seed complete-case. Fixed-membership does not change the qualitative seed-range reading.
- **status:** DOCUMENTED (QA-only)

### P1-4 Ligand PDBQT missing for PIK3CA and AChE (method in yaml only)

- **question:** Can every pair be redocked from this tree?
- **evidence:** PIK3CA pack has 2 receptor PDBQT only. AChE pack has receptors + cognate QC, not panel ligand PDBQT. Track-B has production poses, not `ligands_pdbqt/`.
- **risk:** Claiming bit-for-bit docking reproduction rather than zero-dock score reproduction.
- **action:** Writing must say statistics rebuild from deposited scores. PIK3CA LigPrep predecessor is not current primary (`protocol.yaml` note neutralized this round).
- **status:** DOCUMENTED

### P1-5 Table S9 cites missing ladder/roster markdown

- **question:** Why these eight pairs?
- **evidence:** `pair_eligibility_audit_s14_v1.csv` exists and has no outcome metrics (`docs/TARGET_SELECTION_NO_OUTCOME_LEAKAGE_CHECK.md`). `TIER1_DOCKING_ROSTER_V1.md` and `FEASIBLE_PAIR_LADDER_V1.md` are cited and absent.
- **risk:** Footnote rot, not outcome leakage.
- **action:** Cite the CSV as the surviving audit; mark the md files historical/missing.
- **status:** DOCUMENTED

### P1-6 SI typeset Python 3.12.13 vs freeze 3.12.3

- **question:** Are results dependent on undocumented environment choices?
- **evidence:** `docs/PR39_SCIENTIFIC_DATA_AUDIT.md`: freeze env Python 3.12.3, pin-match for numpy/pandas/scipy/sklearn/rdkit. SI Table S1 typesets 3.12.13.
- **risk:** Reviewer treats 3.12.13 as the analysis pin.
- **action:** Writing pass: SI S1 should follow `check_analysis_env.py`, not a typeset leftover. ECFP max |Δ| is the pinned 0.0234, not 0.0112.
- **status:** DOCUMENTED (manuscript/SI body not edited this round)

### P1-7 Methods still claim uniform ligand prep; numbers still 0.3237

- **question:** Will the current manuscript survive CASE 2 promotion?
- **evidence:** `docs/MANUSCRIPT_JCIM_EN.md` methods; WRITING_INDEX EGFR 0.324 / 0.3237. This audit does not rewrite manuscript body.
- **risk:** After promotion, EN/ZH and WRITING_INDEX will be stale.
- **action:** Required writing pass after EGFR scores are wired. Out of this audit’s P0 docking work.
- **status:** DOCUMENTED

### P1-8 Descriptor baseline selects the best descriptor on the full panel

- **question:** Could simple ligand chemistry explain the docking signal?
- **evidence:** `compute_descriptor_baselines.py` picks the best descriptor on the full panel; manuscript already says this is not a selection-adjusted test. ECFP uses scaffold GroupKFold with shared `fold_id` across ECFP4 / docking / ECFP4+docking. Docking column is the same primary members.
- **risk:** Overstating descriptor “significance”.
- **action:** Keep the existing caveat. No new model fitting this round.
- **status:** DOCUMENTED

### P1-9 No external docking validation (0/8 gate)

- **question:** Why is there no external validation?
- **evidence:** `results/canonical/external_eligibility.csv`; Fig S5 is an eligibility screen.
- **risk:** Calling BindingDB counts external validation.
- **action:** Keep current caption language.
- **status:** DOCUMENTED

### P1-10 Eight pairs are not independent replicates

- **question:** Are the eight pairs independent observations?
- **evidence:** JAK1 (`6N7A`) and PPARA (`6LXA`) each used twice. Manuscript already warns. CIs not multiplicity-adjusted. Detectable-effect is binormal simulation, not observed power.
- **action:** Keep those sentences. Do not add a meta-analysis across eight pairs.
- **status:** DOCUMENTED

---

## P2 OPTIONAL / PRESENTATION

- FIGURE_TABLE_LOCK Unique sources remapped to canonical (done). Historical `data/jcim_*/*_v1.csv` files still exist as superseded packs.
- Track-B README pose wording corrected.
- EGFR box JSON status field corrected.
- Hardcoded `comparable_to_current_primary=1` could later encode AChE n mismatch; not required for this freeze if P1-3 is written.

---

## CLOSED this round

| Item | Closure |
|------|---------|
| PIK3CA protocol “Same receptors/boxes/E as LigPrep panel48” as fact | Rewritten to predecessor-panel + current RDKit/Meeko. Primary RDKit scores not changed. |
| EGFR protocol/MANIFEST presenting uniform RDKit as the deposited method | Deposited method = not_recoverable; rebuild isolated. |
| Two Unique sources for Fig 2–5 / S1–S5 | Lock + submission_pack copy now name `results/canonical/`. |
| `HISTORICAL_LEFTOVER_FILES.md` citation in submission_pack lock | Removed. |
| Pair-selection outcome leakage | CSV has no AUROC/summary_min/GNINA/ECFP/ranking. |
| Activity threshold / pair inclusion / ECFP 0.0112 chase | Unchanged. |
| Receptor re-prep to invent protonation text | Not done. |

---

## Reviewer questions (simulated)

**How were ligands prepared consistently across target pairs?**  
They were not proven to be. Track-B yaml + git prep script: RDKit ETKDGv3 seed 20260727, MMFF 200, Meeko 0.7.1, no LigPrep. PIK3CA/AChE: same method in protocol yaml, ligand PDBQT not deposited. EGFR deposited scores: method **not_recoverable**; CASE 2 rebuild is the intended unification and is not yet primary.

**How were receptors prepared?**  
Deposited PDBQT are the receptors (`docs/RECEPTOR_PREP_PROVENANCE.md`). Prep method not_recoverable. Peptide on 9V8H chain B is in the file. Waters stripped. Do not invent PPW.

**Why these eight target pairs?**  
`PRIMARY_PAIRS` plus `pair_eligibility_audit_s14_v1.csv`. Supply, holo, site, protocol compatibility. EGFR is a supply-limited exception. Supporting roster markdown is missing.

**Could pair selection be outcome-driven?**  
Not in the eligibility CSV. Downstream AUROC files are outcomes.

**Are activity labels harmonized?**  
θ=6.0 four-class everywhere in Table 2. Provenance is not: assay-adjudicated Track-A vs ChEMBL37 dump Track-B. Unresolved missing-arm ligands excluded. AB_056 is the only primary panel-fallback.

**How sensitive are results to activity thresholds?**  
`label_aggregation_sensitivity.csv` / Fig S4A. Primary remains θ=6.0. Not re-optimized this round.

**How sensitive are results to receptor choice?**  
PIK3CA substitution only (Fig 5B). Other pairs: one holo each. Limitation.

**How sensitive are results to random seed?**  
Five Vina seeds, Fig 5C, deposited scores. EGFR five-seed currently uses the same unproven ligand inputs as production (P0-1).

**Does changing sample membership affect seed sensitivity?**  
Only AChE: n 95 → 88 on the five-seed intersection; |Δ| ≤ 0.0161; qualitative range unchanged. QA-only; no redock.

**Could simple ligand chemistry explain the apparent docking signal?**  
Descriptor baselines and ECFP4 OOF are the control. Best-descriptor selection is full-panel. ECFP incremental max |Δ| = 0.0234 (pinned env), not 0.0112.

**Does docking add information beyond ECFP4?**  
Same ligands, shared `fold_id`, `ECFP4+docking` vs `ECFP4` in `ecfp4_incremental_information.csv`. Small increments. Do not claim a large docking increment.

**Is GNINA independent docking or rescoring?**  
Independent search: EGFR, PIK3CA, JAK1/TYK2 (`GNINA_SOURCES`). CNN best-of-9 and RTMScore: Vina-pose rescores (Table S7). EGFR independent GNINA still inherits P0-1 ligand inputs.

**Why is there no external validation?**  
0/8 passed the independent-source docking gate. Fig S5 is eligibility, not validation.

**Are the eight pairs independent observations?**  
No. Shared JAK1 and PPARA receptors. No multiplicity-adjusted CI. No eight-pair random-effects claim.

**Are any results dependent on undocumented software/environment choices?**  
Analysis pin is `requirements-analysis.txt` + Python 3.12.3 freeze env. ECFP 0.0112 was an unpinned-environment leftover. Receptor protonation is undocumented (P1-1). EGFR ligand prep for deposited scores is undocumented (P0-1).

**Can every table and figure be regenerated from supplied data?**  
Zero-dock statistics and figures: yes, from `results/canonical/` + plotter, once EGFR remains on deposited scores. Docking regeneration from ligand PDBQT: **no** for deposited EGFR, PIK3CA, AChE; Track-B poses exist. Item 15 acceptance rebuild is deferred until P0-1 is closed.

---

## Statistics / ML checks (no P0 findings)

- AUROC: higher `score_S` better; dual vs A-only uses `score_B`; dual vs B-only uses `score_A`; `summary_min` = min of those arms; shared-dual draw inside bootstrap; B=2000; seed 20260729.
- Pointwise percentile CIs are not multiplicity-adjusted (stated).
- Detectable-effect: n_mc=1000, not observed power (stated).
- Scaffold GroupKFold; three models share `fold_id`; no test-label feature selection on ECFP; docking scores concatenated, not used to pick folds.
- Primary membership: `is_primary_row` = eight pairs ∩ main ∩ complete_case ∩ activity_eligible.

---

## Item 15 reproducibility run

Not executed as the freeze-acceptance run.

Reason: promoting or comparing canonical while EGFR still sits on unproven ligand coordinates would freeze the wrong science state.

When P0-1 closes, run:

```
python3 scripts/check_analysis_env.py
rm -rf /tmp/dual_target_final_rebuild
python3 scripts/analysis/rebuild_freeze.py --outdir /tmp/dual_target_final_rebuild
python3 scripts/qa/verify_freeze_rebuild.py --outdir /tmp/dual_target_final_rebuild
python3 scripts/qa/compare_rebuild_to_canonical.py \
  --rebuilt /tmp/dual_target_final_rebuild \
  --canonical results/canonical
python3 scripts/qa/check_current_chain.py
git diff --exit-code -- data/processed/activity_adjudication
```

and confirm `submission_pack/tables/canonical` matches `results/canonical`.

---

## READY / NOT_READY rule (applied)

P0 = 0? **No.**  
All P1 have evidence boundaries? **Yes.**  
Fresh rebuild PASS on the final science state? **No (not run).**

**NOT_READY_FOR_WRITING_FREEZE**
