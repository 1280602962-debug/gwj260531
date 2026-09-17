# SECOND PASS FREEZE REPORT

**Not authoritative.** Superseded by `docs/WRITING_FREEZE_REPORT.md`. Do not copy residual 1/96 wording from this file as current.

Branch: `cursor/methods-sentence-audit-c7cc`  
Directory: `Dual_Target_Docking/`  
Date: 2026-09-17  
Scope: residual inconsistencies between `FINAL_FREEZE_REPORT.md` and the current repository. No full freeze re-audit. No new scientific analyses. No submission_pack regeneration. No redocking.

---

## Item 1. AChE/BChE max-vs-median and AB_056 provenance

### Current ligand-level aggregation input

Canonical generator:

- `scripts/analysis/compute_canonical_results.py` → `compute_max_median()`
- Join key: `(pair, ligand_id)` against
  `data/jcim_chembl_universe_v0/local_track_b_v0/tables/eight_pair_dump_gated_v1/max_vs_median_ligand_v1.csv`
- Labels used: dump `class_max` / `class_median` at θ=6.0
- Scores used: current `results/canonical/current_score_master.csv` (`analysis_set=main`, complete-case)
- Ligands absent from the dump-gated table are skipped. Panel pChEMBL is **not** copied.

Dump-gated AChE table: **95** ligands. Class_max = 27 dual / 25 A_only / 28 B_only / 15 neither.

### AB_056 / CHEMBL3960861

| source | status |
|---|---|
| `current_score_master` main complete-case | **present**. `construction_class=A_only`, `primary_class_theta6=A_only`, pA=6.94, pB=4.3, scores 12.777 / 9.05 |
| `panel_v0_strict.csv` / `ablation_ligand_scores.csv` | present; construction rule `strict_6.5_5.5` |
| `eight_pair_dump_gated_v1/max_vs_median_ligand_v1.csv` | **absent** (no n_act_A / n_act_B) |
| `high_confidence_activity_audit_v1.csv` | **0 rows** for CHEMBL3960861 |
| `mols_ACHE.json` / `mols_BCHE.json` | present as REST aggregates 6.94 / 4.3 (same numbers as the panel; not assay-level max/median) |
| `assay_max_vs_median_ligand_v1.csv` | present with n_act=1 copied from panel pChEMBL — **not used**; not a legal dump provenance |

The dump-gated table was harvested from the scored ablation set of that run (the five panel IDs missing from dump are exactly the then-failures: AB_001, AB_053, AB_054, AB_056, AB_097). AB_056 later gained both-end Vina scores and entered the primary complete-case (n=96). The dump harvest was not re-run, so there are still no ChEMBL 37 dump assay rows for this ligand.

**Decision: EXCLUDE AB_056 from the max/median aggregation population.**  
Do not copy panel or REST pChEMBL into max/median.

### Recomputed AChE/BChE (generator + dump-gated labels + current scores)

Bootstrap: class-stratified shared-dual, B=2000, seed=20260729.

| aggregation | n_dual | n_A_only | n_B_only | D_vs_A | D_vs_B | summary_min | 95% CI |
|---|---:|---:|---:|---:|---:|---:|---|
| max | 27 | 25 | 28 | 0.6504 | 0.6058 | 0.6058 | [0.4429, 0.7395] |
| median | 26 | 25 | 28 | 0.6754 | 0.6291 | 0.6291 | [0.4670, 0.7459] |

The one max→median flip among dump-gated ligands is AB_018 / CHEMBL659 (n_act_A=92, n_act_B=55; class_max=dual, class_median=neither). That flip is legal.

### Four-way consistency (this pass)

| piece | AChE max | D_vs_A |
|---|---|---|
| generator `compute_max_median()` | 27/25/28 | 0.6504 |
| ligand-level source dump-gated table | 27/25/28 (95 ligands; AB_056 absent) | — |
| `results/canonical/max_vs_median_sensitivity.csv` | 27/25/28 | 0.6504 |
| `FINAL_FREEZE_REPORT.md` (corrected header) | 27/25/28 | 0.6504 |

The first-pass claim 27/26/28 and D_vs_A=0.6524 mixed **primary** complete-case (which includes AB_056) with **dump-gated** max membership (which does not).

---

## Item 2. AChE primary label definition

`load_main()` sets `cls = primary_class_theta6`. That field is written by `build_current_score_master.row_of()` as `assign_fourclass(pA, pB)` at `THETA_PRIMARY=6.0`. `construction_class` is the panel-construction label from `panel_v0_strict.csv` and is **not** the AUROC class.

Independent recount on AChE main complete-case (n=96), from current `pA`/`pB`:

| definition | dual | A_only | B_only | neither | gray | disagreements vs `primary_class_theta6` |
|---|---:|---:|---:|---:|---:|---:|
| production / current (`primary_class_theta6`) | 27 | 26 | 28 | 15 | — | 0 |
| `construction_class` (strict panel) | 27 | 26 | 28 | 15 | — | 0 |
| true θ=6.0 from pA/pB | 27 | 26 | 28 | 15 | — | 0 |
| true θ=5.5 from pA/pB | 27 | 26 | 28 | 15 | — | 0 |
| true θ=6.5 from pA/pB | 27 | 26 | 28 | 15 | — | 0 |
| strict 6.5/5.5 from pA/pB | 27 | 26 | 28 | 15 | 0 | 0 |

Ligand-level disagreement list: **empty**.

### Answers

1. **Current AChE primary uses θ=6.0.** The class consumed by directional AUROC / `summary_min` is `primary_class_theta6`. It happens to equal strict 6.5/5.5 on this 96-ligand panel (0 flips), because the panel was sampled from the strict pool and contains no 5.5–6.5 gray-zone ligands.

2. **Why `FINAL_FREEZE_REPORT` wrote strict:** it treated `construction_class` / `panel_v0_strict.csv` as the production analysis label. That is the **sampling** rule, not the field `load_main()` uses. The first-pass Methods suggestion “θ=6.0 except AChE strict” mixed those two layers.

3. **Why the master field is still `primary_class_theta6`:** every pair’s analysis label is computed by the same `assign_fourclass(..., cut=6.0)` in `row_of()`. Construction may be θ=6.0 (EGFR, PIK3CA/mTOR) or strict 6.5/5.5 (the other six, including AChE). The analysis field name records the analysis cut, not the sampling pool.

4. **If all eight pairs are unified to θ=6.0:** AChE directional AUROC and `summary_min` **do not change**. They are already θ=6.0:

   - D_vs_A = 0.6524
   - D_vs_B = 0.6058
   - `summary_min` = 0.6058 [0.4392, 0.7354]
   - n = 27 / 26 / 28 (neither 15)

5. **If the manuscript instead kept strict as the *analysis* rule:** the current Methods already say primary analysis used θ=6.0 and that AChE was drawn from the strict 6.5/5.5 pool. Those sentences are correct for the code. What would have to change is any sentence that claims the AChE AUROC class is strict rather than θ=6.0 — specifically the first-pass freeze suggestion “eight pairs; θ=6.0 except AChE strict 6.5/5.5”. Table 1 “label_rule = strict 6.5/5.5” is the **construction** column and can stay. Primary labels were not modified in this pass.

---

## Item 3. AChE threshold sensitivity (real relabel)

Source: current master `pA`/`pB` for AChE main complete-case. Generator: `compute_label_sensitivity()` (now also writes `n_neither`, `n_gray`, `n_missing_pchembl`, `status`). Not copied production rows.

All 96 ligands have pA and pB, so every rule is **evaluable as a relabel of the current scored set**. This is not a reconstruction of a θ=5.5 or θ=6.5 *sampling pool* from ChEMBL (gray-zone ligands were never admitted to this panel).

| label_rule | n_dual | n_A_only | n_B_only | n_neither | n_gray | D_vs_A | D_vs_B | summary_min | 95% CI | status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| theta=5.5 | 27 | 26 | 28 | 15 | 0 | 0.6524 | 0.6058 | 0.6058 | [0.4392, 0.7354] | OK |
| theta=6.0 | 27 | 26 | 28 | 15 | 0 | 0.6524 | 0.6058 | 0.6058 | [0.4392, 0.7354] | OK |
| theta=6.5 | 27 | 26 | 28 | 15 | 0 | 0.6524 | 0.6058 | 0.6058 | [0.4392, 0.7354] | OK |
| strict 6.5/5.5 | 27 | 26 | 28 | 15 | 0 | 0.6524 | 0.6058 | 0.6058 | [0.4392, 0.7354] | OK |

Identity across the four rules is a property of this strict-constructed panel, not a copy of the production row. EGFR and PIK3CA/mTOR in the same file still change with θ (sanity check that the generator actually relabels).

Written to: `results/canonical/label_aggregation_sensitivity.csv`.

---

## Item 4. Stale auxiliary files

### `docking_failure_census_v1.csv`

Recomputed from current `ablation_ligand_scores.csv` via `claim_hardening_v1.docking_census()` (census only; aggregation/descriptor tables not rewritten).

AChE main panel, n_attempted=100:

| field | before | after |
|---|---:|---:|
| n_success_pocket_A | 96 | 97 |
| n_success_pocket_B | 95 | 96 |
| n_success_both_ends | 95 | **96** |
| n_fail_either_end | 5 | **4** |
| failed_ligands | AB_001, AB_053, AB_054, **AB_056**, AB_097 | AB_001, AB_053, AB_054, AB_097 |

AB_056 is a both-end success (`vina_ACHE=-12.777`, `vina_BCHE=-9.05`).

Also synced the Table S27 census row in `data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` (n_success_both_ends 95→96; AB_056 removed from the fail list).

### Canonical / figure-source files tied to items 1–2

- `results/canonical/max_vs_median_sensitivity.csv` — regenerated by the generator; AChE max remains 27/25/28, D_vs_A=0.6504.
- `results/canonical/label_aggregation_sensitivity.csv` — regenerated from pA/pB; AChE AUROCs unchanged, `n_neither` now explicit.
- `figures/jcim_article/plotted_values.json` AChE n=27/26/28 is the **primary** complete-case, not max-vs-median. Not modified.
- `data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv` AChE rows already equal the real relabel. Not rewritten (figure SHA would change without a numerical change).

Not regenerated (out of scope): `submission_pack/`, docking, other analyses.

Left stale on purpose (not the canonical max-vs-median source):

- `data/jcim_novelty_v0/tables/assay_max_vs_median_*.csv` still contain copied panel pChEMBL for AB_056 (first-pass overlay).
- `data/jcim_novelty_v0/tables/docking_failed_ligand_properties_v1.csv` still lists AB_056 as a torsdof skip.

---

## Item 5. Freeze checks

### A. max-vs-median canonical vs report

**No unresolved conflict.** Generator, dump-gated ligand table, `results/canonical/max_vs_median_sensitivity.csv`, and the corrected `FINAL_FREEZE_REPORT` header all use 27/25/28 and D_vs_A=0.6504. Primary 27/26/28 / 0.6524 is a different population (includes AB_056).

### B. AChE primary label definition

**Fully explicit.** Analysis = `primary_class_theta6` = θ=6.0. Construction = strict 6.5/5.5. On this panel they coincide (0 disagreements). Primary labels were not changed.

### C. AChE threshold sensitivity

**Real computation** from current pA/pB. All four rules evaluable on the current 96. Composition and AUROC are identical; that is a panel-construction result, not a copied row.

### D. Publication figures / tables

**Primary Table 2 / Figure plotted values are not affected.** AChE primary remains 27/26/28, D_vs_A=0.6524, `summary_min`=0.6058 [0.4392, 0.7354].

Affected only as SI/source hygiene:

- Table S27 docking census source now 96 both-end successes (pack copy not rebuilt).
- SI prose that says AChE max-vs-median “1/96” and “0 missing ends” still describes n=96 as if AB_056 had dump assay rows. The dump-gated comparison is **1/95** (CHEMBL659). That wording was not edited in this pass.

---

## CHANGED NUMBERS

| item | before (conflict / stale) | after |
|---|---|---|
| AChE max-vs-median membership (canonical + report) | report claimed 27/26/28, D_vs_A=0.6524; CSV was 27/25/28, 0.6504 | **27/25/28**, D_vs_A=**0.6504**, D_vs_B=0.6058, summary_min=0.6058 [0.4429, 0.7395] |
| AChE median-vs-median (recomputed, unchanged vs CSV) | 26/25/28, D_vs_A=0.6754 | **26/25/28**, D_vs_A=**0.6754**, D_vs_B=0.6291, summary_min=0.6291 [0.4670, 0.7459] |
| AChE docking census both-end success | 95; AB_056 listed as fail | **96**; remaining fails AB_001, AB_053, AB_054, AB_097 |
| AChE θ-grid schema | no n_neither; first pass called the identical rows “not a real reassignment” | real relabel; n_neither=15, n_gray=0, status=OK at all four rules |

## UNCHANGED PRIMARY NUMBERS

| item | value |
|---|---|
| AChE primary complete-case | 27 dual / 26 A_only / 28 B_only / 15 neither; n=96 |
| AChE D_vs_A | 0.6524 |
| AChE D_vs_B | 0.6058 |
| AChE summary_min | 0.6058 [0.4392, 0.7354] |
| Other seven pairs’ primary AUROC / class counts | unchanged |
| EGFR / PIK3CA/mTOR θ-grid numbers | unchanged (still actually move with θ) |
| Primary labels | not modified |

## REMAINING LIMITATIONS

1. Dump-gated max/median cannot be rebuilt here: ChEMBL 37 sqlite is absent. AB_056 exclusion rests on the committed dump-gated ligand table, not a live harvest.
2. Novelty `assay_max_vs_median_*.csv` still has copied panel pChEMBL for AB_056. It is **not** the canonical source.
3. `submission_pack/` was not regenerated; pack copies of census / max-vs-median remain older.
4. `docking_failed_ligand_properties_v1.csv` still lists AB_056 as a fail.
5. Manuscript/SI “1/96” max→median flip (CHEMBL659) should be 1/95 dump-gated ligands plus one complete-case ligand (AB_056) with no dump assay provenance. Not edited in this pass.
6. AChE θ=5.5/6.5/strict sensitivity cannot recover gray-zone ligands that were never sampled into `panel_v0_strict.csv`. Relabel of the current 96 is evaluable; a universe-level panel redraw is not.

## FINAL VERDICT

**READY_TO_FREEZE**

A–D have no unresolved conflict in the current generator, dump-gated ligand source, `results/canonical` tables, or the corrected freeze report. Primary AChE numbers are unchanged. Max-vs-median is explicitly the dump-gated 95-ligand population. AChE primary is explicitly θ=6.0. Threshold sensitivity is a real pA/pB relabel.
