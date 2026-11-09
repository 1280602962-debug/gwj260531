# FINAL FULL-PROJECT SCIENTIFIC FORENSIC AUDIT

Date: 2026-09-19  
Subject: Dual_Target_Docking on PR #39 (`cursor/scientific-freeze-c7cc`)  
Mode: READ-ONLY except `docs/audit/`, `data/provenance/`, `results/qa/`  
Independent statistics: `scripts/qa/run_final_forensic_audit.py` (does not import `bootstrap_metrics` or `compute_canonical_results`)

---

## 1. Executive verdict

**NOT_READY_FOR_FINAL_MANUSCRIPT_DRAFTING**

Primary **numbers** independently replay. Score extraction samples match. Pair selection as recorded is not outcome-driven. EGFR uniform RDKit/Meeko is the only current EGFR production truth. Analysis pins are internally consistent (Python 3.12.3, not 3.12.13).

The freeze is **not** ready for a final manuscript draft because **current-facing SI and main-text claims still contradict the deposited tables** they cite. That is a P1 writing-integrity failure under a JCIM standard, not a reason to change Table 2.

P0 (data / primary arithmetic / score sign / hidden exclusion / outcome-driven roster) = **0**.  
P1 (must fix or explicitly archive before submission) remain.

---

## 2. P0 issues

None.

Independent replay of directional AUROC, summary_min, CIs, D-vs-neither, EF, top-k, AND dual counts, fixed-score Δ, and paired matched−mismatched Δ: **192/192 MATCH** (tolerance 1e-4).  
40/40 source→master score traces MATCH.  
θ=6.0 class reconstruction on primary rows: 0 mismatches.  
EH120_059 / AB_087 do not enter primary.

A first audit-only matched/mismatched CI used two independent draws and produced 16 false P0s. Production uses one paired draw. After the audit matched that protocol, the CIs matched. That error is documented in `NEWLY_DISCOVERED_ISSUES.md` so it is not recycled.

---

## 3. P1 issues

1. **SI Table S6 prose contradicts SI Table S6 cells and Figure 4.**  
   Prose: only AChE excludes 0; EGFR includes 0.  
   Table: EGFR 0.107 [0.006, 0.220], excludes 0 = yes.  
   Copied into `submission_pack/SI/`.

2. **Main-text EGFR cluster CIs are stale.**  
   Typeset [0.235, 0.665] and [0.125, 0.644].  
   Canonical [0.232, 0.637] and [0.117, 0.636]. SI S4 is already correct. Qualitative “excludes 0” unchanged.

3. **SI header footnote still cites EGFR TPSA 0.4275 → 0.428.**  
   Current TPSA D-vs-B / summary_min = 0.4389.

4. **Figure 5C lock still plots fuzzy `five_seed_comparable=1` for AChE** while `same_membership_as_primary=0`. Production CSV has the two flags; the plotted-value schema does not use them.

5. **`GNINA_SOURCES` default remains the historical EGFR GNINA path.** Runtime override exists.

6. **Root `docs/` freeze memos still treat 0.3237 as current Table 2** (superseded historical snapshots now under `docs/archive/`: `PR39_FINAL_SCIENTIFIC_FREEZE_CHECK.md`, `PR39_SCIENTIFIC_DATA_AUDIT.md`, `SUBMISSION_ONLY_CLEANUP_REPORT.md`).

7. **Receptor preparation method is NOT_RECOVERABLE** for 14 slots. Manuscript currently states this; any leftover “Meeko/PPW prepared all receptors” would be P1. None found in the current EN Methods paragraph.

8. **Ligand PDBQT absent for PIK3CA, AChE, Track-B.** Zero-dock analysis is reproducible; SMILES→dock is not. Must stay in Methods/SI.

9. **Track-B activity is `chembl37_dump_panel`, not assay-adjudicated.** Already distinguished in Methods; must not regress.

10. **Legacy EGFR box JSON files remain** next to corrected files (coordinates currently identical). Confusion risk.

---

## 4. P2 issues

- Duplicate identical `data/processed/current_score_master.csv` (superseded; now pointer-only).
- Nested descriptor fold AUROCs of 0 or 1 on tiny test folds (expected; do not over-interpret fold table).
- Figure 3B displayed 0.023 vs raw 0.0234 (rounding).
- Scaffold `except Exception` fallbacks to `"unparsed"` / `"acyclic"`.
- Recovered PIK3CA script and ENV_PIN contain personal absolute paths (historical snapshots).
- 9.505 appears as an AChE **score** (AB_038 pocket B) and as a forbidden historical RMSD token. Context-dependent; not a current RMSD.

---

## 5. Closed issues

| Item | Evidence |
|------|----------|
| EGFR uniform is current Table 2 | Master `score_source` = uniform five-seed CSV; n=28/37/31; summary_min 0.3341 [0.1970, 0.4713] |
| Old 0.3237 not in manuscript Table 2 | EN Table 2 typeset 0.334 |
| Primary estimand mapping | Config + independent replay: D-vs-A uses score_B; D-vs-B uses score_A; two-pocket mean only in Table 3 |
| θ=6.0 Table 2 labels | 0 class mismatches |
| EH120_059 / AB_087 excluded | Master flags; not in primary membership |
| AChE five-seed n=95 / 94 / 88 | Independent fixed-membership n=88; max \|Δ\| vs available-case ≈ 0.016; qualitative change = no |
| ECFP max \|Δ\| = 0.0234 | Canonical and independent max |
| Cognate EGFR/HER2 1.019 / 1.947 Å | `cognate_rmsd.csv`; 9.505 / 0.760 not current |
| External = eligibility only, 0/8 | `external_eligibility.csv`; no external docking |
| JAK1/TYK2 document cluster unresolved | status `unresolved_mapping_unavailable` |
| Path audit unresolved = 0 | `docs/CURRENT_PATH_AUDIT.md` (re-used as a pointer; not as scientific truth) |
| No outcome leakage in eligibility CSV | Token scan + structural PIK3CB fail_H3 |
| Analysis Python 3.12.3 | SI S1 and `requirements-analysis.txt`; 3.12.13 absent from SI S1 |

---

## 6. Newly discovered issues

See `docs/audit/NEWLY_DISCOVERED_ISSUES.md`.  
Non-empty. Highest new finding: SI S6 prose vs table.

---

## 7. Primary-result independent replay

Source: `results/canonical/current_score_master.csv` (treated as input, not as a computed truth).  
Implementation: Mann–Whitney mid-rank AUROC; class-stratified shared-dual bootstrap B=2000 seed=20260729; paired MM draw.

| Pair | n D/A/B | D-vs-A (B) | D-vs-B (A) | summary_min | MM Δ | Replay |
|------|---------|------------|------------|-------------|------|--------|
| EGFR/HER2 | 28/37/31 | 0.6631 | 0.3341 | 0.3341 | 0.1071 | MATCH |
| JAK1/JAK2 | 32/32/32 | 0.5884 | 0.7275 | 0.5884 | −0.0190 | MATCH |
| JAK1/TYK2 | 31/32/32 | 0.5751 | 0.3649 | 0.3649 | −0.0645 | MATCH |
| PIK3CA/mTOR | 18/14/12 | 0.7143 | 0.6921 | 0.6921 | 0.0903 | MATCH |
| AChE/BChE | 27/26/28 | 0.6524 | 0.6058 | 0.6058 | 0.1770 | MATCH |
| F2/F10 | 31/32/32 | 0.4133 | 0.3448 | 0.3448 | −0.0307 | MATCH |
| PPARG/PPARA | 32/31/32 | 0.6492 | 0.7061 | 0.6492 | 0.0296 | MATCH |
| PPARA/PPARD | 32/32/32 | 0.6465 | 0.4463 | 0.4463 | 0.0117 | MATCH |

Machine file: `results/qa/independent_statistics_replay.csv`.

---

## 8. Experimental-design assessment

The live experiments answer the stated question: directional discrimination of dual vs single-target-active ligands, plus negative-class dependence.

HARKing risk is **present as writing drift**, not as a changed primary estimand. The primary mapping is still D-vs-A / pocket B and D-vs-B / pocket A. Two-pocket mean is not used as a directional primary.

Post-hoc metrics (AND filter, ranking, descriptors, five-seed, GNINA) are labeled supporting/diagnostic when the lock documents are followed. SI S6 prose is the place where a **secondary narrative** (only AChE) was not updated after EGFR MM flipped.

EGFR supply exception and PIK3CA E=16 have design-period reasons in files (supply; 4JT6 coverage). They are not AUROC shopping in the eligibility CSV.

---

## 9. Data-provenance assessment

Activity: three kinds, correctly separated in the master. Track-B is dump-gated.  
Ligands: EGFR CONFIRMED; others METADATA_ONLY / PDBQT not in tree.  
Receptors: PDBQT CONFIRMED; method NOT_RECOVERABLE.  
Boxes: 14 current JSONs exist; EGFR corrected declared heavy-atom.  
Scores: deposited tables; 40/40 sample traces MATCH.

Cannot reconstruct protonation, pH, loops, tautomer enumeration, or a common PPW/Meeko receptor campaign. Write NOT_RECOVERABLE, not a guessed protocol.

---

## 10. Docking-protocol assessment

PRIMARY E=16 PIK3CA only; others E=8.  
E8 / PM110 / 4JPS / 5DXT / 4JSX are sensitivity.  
Independent GNINA: three pairs only.  
Same-pose RTM/CNN must not be called independent docking.  
EGFR current seed 20260727; five-seed 20260811–14 deposited.

No unexplained PRIMARY/SENSITIVITY collapse was found in `analysis_config.py` or protocol docs.

---

## 11. Statistical assessment

Bootstrap: class-stratified, shared dual, B=2000, seed 20260729, percentile CI, point ≠ bootstrap mean. Independent CIs MATCH.

Reference: AUROC 0.5; Δ 0. Pointwise intervals. Eight pairs are not independent replicates.

JAK1/TYK2 document cluster remains `unresolved_mapping_unavailable`.

Detectable-effect simulation is named correctly in SI/main. Must not be renamed observed power.

---

## 12. ML / descriptor assessment

ECFP: radius 2, 2048 bits, GroupKFold, unscaled primary, same folds across models. Max \|Δ\| 0.0234.

Full-panel best descriptor is a **descriptive univariate screen**. Nested scaffold-CV is the predictive chemistry-control number. SI already says this. Do not put full-panel best in the abstract as a predictive estimate.

TPSA on AChE uses the primary AChE membership (including AB_056). Max-vs-median AChE n=94 is a different set.

---

## 13. Reproducibility assessment

A (zero-dock): defendable from frozen scores + analysis pins. Independent arithmetic done. Full temp-dir rebuild not re-run in this pass.

B (full dock): only EGFR uniform is a deposited SMILES→PDBQT→Vina/GNINA campaign. Seven other pairs are score-archive reproducible, not input-archive reproducible.

Do not claim eight-pair full reproducibility.

---

## 14. Script / code assessment

See `ACTIVE_SCRIPT_AUDIT.md`.  
Risks: historical GNINA default path; fuzzy comparable flag; two master copies; promote vs live writers.  
No silent AUROC invention.

---

## 15. Figure assessment

`plotted_values_postfix.json`: 319 MATCH, 18 SOURCE_EXISTS, 1 rounding miss on 0.023 vs 0.0234, 1 FUZZY_COMPARABLE_FLAG (Fig 5C AChE).

Scientific-expression risks: Fig 5C can look like AChE five-seed equals primary membership; Fig S5 must remain eligibility, not validation; sensitivity panels must not be visually equal to Table 2.

Main Figure 4 caption already states AChE **and** EGFR exclude 0. That is correct.

---

## 16. Table assessment

Table 1–3 cells traced: MATCH or ROUNDING_MATCH for numeric AUROC/CI/EF; Table 1 quotas correctly flagged as construction quotas, not n_scored.

SI S1–S10 source files exist. S2 boxes/RMSD parse to deposited JSON/CSV. S4 cluster cells match canonical. S6 **cells** match; S6 **prose** does not.

---

## 17. Claim–evidence assessment

See `CLAIM_EVIDENCE_MATRIX.md`.

Core scientific claims about negative-class dependence, pair-specific AUROCs, chemistry baselines, and incremental ECFP are supported with limits.

Two current-facing sentences are **OVERSTATED / STALE**: “only AChE excludes 0” (SI) and the old EGFR cluster CI bounds (main text).

---

## 18. Reviewer red-team assessment

See `REVIEWER_RED_TEAM_REPORT.md`.

No POTENTIAL_FATAL_ISSUE in the deposited numbers.  
Writing-quality fatal risk: SI S6 self-contradiction. A reviewer can find it in one sitting.

True limitations (receptor prep, ligand PDBQT archive, Track-B activity kind, single rigid receptor, GNINA subset, shared targets, small n, no external docking) are already mostly in the paper and must stay.

---

## 19. Exact remediation plan

Remediation is a **separate writing/docs branch**. Do not change canonical numbers.

1. Rewrite SI S6 (EN + ZH + submission_pack) to match the table: AChE **and** EGFR main-panel CIs exclude 0.
2. Replace manuscript EN/ZH cluster CIs with SI S4 / `cluster_bootstrap_sensitivity.csv` values.
3. Delete or update the 0.4275 TPSA footnote.
4. Change Figure 5C plotted-value metrics from `five_seed_comparable` to `same_protocol_as_primary` / `same_membership_as_primary`; set AChE membership = 0.
5. Point `GNINA_SOURCES` default at the uniform EGFR file, or delete the stale key.
6. Move or banner the three root docs that still treat 0.3237 as current.
7. Keep or delete legacy EGFR box JSON, but do not let current scripts read the uncorrected names.
8. Optionally drop `data/processed/current_score_master.csv` or generate it only as a copy with a pointer.

No new docking. No receptor re-prep. No pair added or dropped. No AUROC target.

---

## 20. What must NOT be changed

- Table 2 / Table 3 canonical numbers
- EGFR uniform as the only current EGFR truth (do not restore 0.3237)
- θ=6.0 primary labels
- Directional estimand mapping
- PRIMARY_PAIRS roster
- PIK3CA E=16 as primary; E=8/PM110/alts as sensitivity
- Track-B yaml `do_not_now`
- Receptor PDBQT files
- AChE production scores (no redock)
- ECFP pipeline (except documenting 0.0234)
- Calling detectable-effect “observed power”
- Calling 0/8 eligibility an external validation

---

## 21. What does NOT require additional experiments

- Fixing SI/main-text contradictions
- Figure 5C metadata flags
- Archiving stale 0.3237 memos
- Documenting NOT_RECOVERABLE receptor/ligand boundaries
- Nested vs full-panel descriptor wording
- External-eligibility wording

New experiments that are **out of scope** and not required to close this audit: MD, extra receptor swaps, eight-pair GNINA, receptor re-prep, AChE redock, new target pairs.

---

## Limitations that belong in the paper (not automatic P0)

Must stay in main text or SI: receptor-prep provenance; ligand-input archive completeness; Track-A vs Track-B activity kinds; rigid single receptor; limited substitution; GNINA on three pairs only; shared targets / multiplicity; small n and PM neither n=4; no external docking; descriptor selection optimism if full-panel best is shown; ChEMBL assay heterogeneity; max-pChEMBL aggregation; no systematic tautomer enumeration; AChE seed-dependent membership; construction labels vs primary θ=6.0.

---

## Gate check (section 43)

| Gate | Status |
|------|--------|
| P0 = 0 | PASS |
| Primary numbers independently replayed | PASS |
| Current scientific paths exist | PASS (prior path audit unresolved=0; not re-scanned as truth) |
| Primary ligand/receptor/box/score traceable or bounded | PASS with NOT_RECOVERABLE boundaries |
| No unexplained protocol contradiction | PASS in code/protocol; FAIL in SI S6 prose |
| No outcome-driven pair selection evidence | PASS (NO_OUTCOME_LEAKAGE) |
| No hidden sample exclusions | PASS (IDs in membership registry) |
| No unresolved score-sign / pocket-map error | PASS |
| Tables/figures numeric cells traceable | PASS; SI S6 prose FAIL |
| AChE five-seed membership handled | PASS in CSV flags; FAIL in Fig 5C plotted comparable flag |
| EGFR uniform is only current EGFR truth | PASS in canonical/master |
| Analysis env clean-rebuildable | PARTIAL (pins exist; temp rebuild not re-run this pass) |
| Canonical ↔ submission_pack sync | FAIL (shared stale SI prose) |
| Core claims evidence-supported | FAIL until SI S6 + cluster CI sentences are fixed |

**NOT_READY_FOR_FINAL_MANUSCRIPT_DRAFTING**

Remaining P0: none.  
Remaining P1: SI S6 prose; manuscript cluster CIs; TPSA 0.4275 footnote; Fig 5C comparable flag; historical GNINA default; stale 0.3237 root docs; archive boundaries already mostly written but must not regress.
