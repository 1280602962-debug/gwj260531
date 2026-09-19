# Stale / current-truth token audit

Date: 2026-09-18  
Scope: non-archive current documents. Not a forbidden-word CI gate.

| Token | Classification | Where it remains | Action this round |
|-------|----------------|------------------|-------------------|
| 0.430 / 0.808 / 0.378 as EGFR AUROC | historical_with_context when listed as pre-fix; stale_error if claimed current | WRITING_INDEX / FIGURE_TABLE_LOCK “must not claim” lists | left as forbidden-current lists |
| 0.808 in SI Table S5 AChE ECFP+docking | current (different quantity) | SI EN/ZH | do not treat as EGFR leftover |
| 0.0112 | historical_with_context (superseded ECFP env) | audit docs; WRITING_INDEX | must not be chased; canonical is 0.0234 |
| 9.505 Å | historical_with_context (pre-fix EGFR RMSD) | WRITING_INDEX forbidden list; SI Fig S3 context | not current 3POZ RMSD (1.019 Å) |
| 0.760 | historical_with_context (PR35 reconstructed RMSD) | WRITING_INDEX; PR35 decisions | not current |
| PIK3CA/PIK3CB | historical_with_context / withdrawn primary | WRITING_INDEX forbidden; some data READMEs | not a current pair |
| different_box_realization | superseded helper string | `plot_jcim_article_figures_v3.py` only | helper, not current five-seed status (`current_or_compatible`) |
| LigPrep as current EGFR/PIK3CA primary | stale_error before this round | protocol notes / MANIFEST | neutralized; remaining LigPrep text is predecessor/sensitivity |
| plotted_values.json | must_be_absent | QA forbids; plotter deletes | absent |
| results/freeze_rebuild | must_be_absent in repo | rebuild writes `/tmp/...` | absent |
| HISTORICAL_LEFTOVER_FILES.md | stale_error citation | submission_pack lock | citation removed; file still absent (correct) |
| Python 3.12.13 | typeset SI S1 vs freeze 3.12.3 | SI | P1 document: typeset string is not the analysis pin (3.12.3) |
| EGFR box `corrected_not_yet_canonical` | stale_error | box JSON | changed to `current_primary_corrected_heavy_atom` |
| Track-B README “poses not in slim tree” | stale_error | README | corrected; poses exist |
| FIGURE_TABLE_LOCK historical Unique sources | stale_error (two current-looking sources) | lock + submission_pack copy | remapped to `results/canonical/` |
| EGFR protocol claiming uniform RDKit as deposited method | stale_error / two current truths | panel120 protocol/MANIFEST | deposited method = not_recoverable; rebuild separate |

Manuscript EN/ZH still describe a uniform RDKit/Meeko ligand pipeline and still typeset EGFR `summary_min` 0.3237. Those sentences are **out of date relative to the ligand-prep audit**. This round does not rewrite manuscript body; they are P1 writing-pass items after CASE 2 promotion.
