# DualFourClass scientific freeze inventory

Baseline: PR #38 `origin/cursor/methods-sentence-audit-c7cc`  
SHA (re-fetched 2026-09-18): `17435413410290960da3437de640509705f00b51`  
Working branch: `cursor/scientific-freeze-c7cc`  
Scope: previous-round scientific issues only. No new pairs, docking, models, or unbounded audit.

Config used by code: `scripts/analysis/analysis_config.py`.

---

## A. Confirmed errors (must close this round)

| ID | Issue | Evidence | Scope | Decision | Done when |
|---|---|---|---|---|---|
| C1 | Chinese Discussion still says ECFP4+docking max \|Δ\| = **0.023**; English Results/Discussion and ZH Abstract/Results say **0.011**. | `docs/MANUSCRIPT_JCIM_ZH.md` L251 vs EN L198/L251; `ecfp4_incremental_information.csv` max \|Δ\| = 0.0112 at JAK1/JAK2 D_vs_A. Patcher needle was `AUROC 最大绝对变化仅为 0.023` (text has no 仅). `check_current_chain.py` did not test this sentence. | Language only. Table S5 / Figure 3B unchanged. | Fix ZH Discussion to 0.011; add the actual needle to `patch_publication_text.py`; QA must fail if 0.023 remains as an ECFP incremental claim. | ZH/EN Discussion and Results all use the three-decimal max \|Δ\| from the rebuilt incremental table; QA PASS. |
| C2 | `compute_canonical_results.py` did not emit the Table 3 two-pocket-mean **CI**, JAK1/TYK2 GNINA row, or PIK3CA receptor-substitution rows. Those lived in `close_publication_from_canonical.py`, which also writes figure/pack files. | Canonical `two_pocket_mean_ranking.csv` has `ci_lo/ci_hi`; `write_two_pocket_ranking()` omitted them. `gnina_rows()` listed only EGFR and PIK3CA; canonical `computational_robustness.csv` has JAK1/TYK2 + 4JPS/5DXT/4JSX. | Rebuild completeness. Primary Table 2 points were already computed in the science script. | Move retained robustness into the science rebuild; do not call pack/manuscript generators during rebuild. | Freeze directory contains two-pocket CIs, three GNINA pairs, and three receptor substitutions, all recomputed from deposited scores + eligible classes. |
| C3 | `receptor_sub()` hard-coded PIK3CA primary `summary_min` as `0.6921`. | `close_publication_from_canonical.py` L324. | Receptor-substitution table note column. | Read the rebuilt PIK3CA `summary_min`. | Note matches rebuilt Table 2. |

---

## B. Risks checked this round (not primary-number bugs)

| ID | Issue | Evidence | Scope | Decision | Done when |
|---|---|---|---|---|---|
| R1 | AChE max/median n=94 vs primary eligible n=95. | Master main: eligible complete-case 95 (27/26/28/14); `activity_status=both_arms_present` 94; AB_087 ineligible (`unresolved_missing_arm`); AB_056 `panel_pchembl_no_audit_rows` kept in Table 2, dropped from max/median. Manuscript 1/94 matches the sensitivity denominator. | Sensitivity Table S3, not Table 2. | Keep. Writing index must say 94 is the qualified-record intersection. Do not force 95 into max/median without both-arm records. | Independent replay of max/median n=94; primary n_scored 27/26/28. |
| R2 | PR #35 RMSD numbers (EGFR 3POZ best 0.760 / top-1 9.505) vs current 14-receptor table (best/top-1 1.019). | `bb3f134f` `primary_cognate_rmsd_calcrrms_v1.csv` vs `all14_cognate_rmsd_calcrrms_v1.csv` (`canonical_heavy_atom_redock`). | SI Table S2 / Fig S3. | Do not replace corrected-box RMSD with PR #35 reconstructed-QC numbers. Adopt unification only. See `PR35_ITEM_DECISIONS.md`. | Current 14-row CalcRMS table remains the RMSD source; PR #35 CSV kept as evidence, not copied over. |
| R3 | Holdout loader did not apply `activity_eligible`. | `compute_holdout()` filtered complete-case only. Live holdout ineligible count = 0. | Holdout matched−mismatched. | Apply the same eligibility gate. Numbers should not move. | Rebuild holdout n unchanged; code path uses `is_primary_row(..., analysis_set="holdout")`. |
| R4 | Five-seed archive still has complete-case EGFR 28/38/32 and D vs A 0.6607. | `data/jcim_multiseed_v0/tables/multiseed_auroc_by_seed_v2.csv`. Figure 5C plots `summary_min` (EGFR 0.3237). | Archive only. | Keep archive. Do not cite 0.6607 as Table 2. | Writing index lists the archive as historical. |
| R5 | `check_current_chain.py` treated file presence + Table 2 point replay as enough. | Script did not check incremental 0.023, two-pocket CI provenance, or freeze-dir independence. | QA. | Independent verifier checks labels, n, score channel, folds, points, and CIs from the freeze directory. | `verify_freeze_rebuild.py` PASS without reading old canonical as the truth. |

---

## C. Research limitations (do not “fix” by invention)

| ID | Limitation | Effect on claims | This round |
|---|---|---|---|
| L1 | JAK1/TYK2 document-cluster map not deposited; ChEMBL 37 sqlite not used. | Document-cluster CI is `unresolved_mapping_unavailable`. Scaffold cluster and ligand-level Δ remain. | Keep status; do not copy an old CI. |
| L2 | Receptor protonation / tautomer / missing-loop protocol unrecorded. | Methods limitation. | No invented protocol. |
| L3 | Matched−mismatched intervals are not multiplicity-adjusted. | Only AChE main-panel CI excludes 0; that is descriptive. | Keep wording. |
| L4 | Detectable-effect is a binormal simulation (N_MC=1000), not observed power. | Caps how strongly we can talk about “underpowered”. | Rebuild from current class sizes; do not fill from old CSV. |
| L5 | GNINA intersection is smaller than Vina eligible n (EGFR 108 vs 109). | Independent pose-generation robustness, not a Vina replacement. | Recompute from deposited poses ∩ eligible classes. |
| L6 | AB_056 has panel pChEMBL but no dump audit rows. | In Table 2; out of max/median. | Document, do not drop from primary. |

---

## D. Explicitly out of scope

- New target pairs, new docking, new models, new power-simulation designs.
- Merging PR #35 or any earlier Dual_Target PR as a whole.
- Cleaning evidence directories this round.
- Editing `JNK1_Selectivity_Project/`.
- Merging to main, closing PRs, deleting branches.

---

## E. Primary locked claims (writing)

Eight pairs; θ=6.0; activity-eligible complete-case; Vina mode-1 higher-better; D vs A = pocket B; D vs B = pocket A; `summary_min` = weaker arm; two-pocket mean only for dual-vs-neither and ranking; EGFR corrected-box ablation; AChE corrected panel; scheme B B=2000 seed=20260729; no pre-fix 0.430 / 0.808 / 0.378 as current AUROCs.
