# WRITING FREEZE REPORT

Branch: `cursor/methods-sentence-audit-c7cc`  
Directory: `Dual_Target_Docking/`  
Date: 2026-09-17  
Review commit at start of this freeze: `a91324bf`  
Verdict: **READY_FOR_WRITING**

`remediation_outputs/freeze_audit/FINAL_FREEZE_REPORT.md` and `SECOND_PASS_FREEZE_REPORT.md` are not authoritative. `figures/jcim_article/FIGURE_AUDIT_PR32.md` and `data/manuscript_lock/ARTICLE_ASSET_INDEX_v1.csv` are historical indexes, not current Table 2 sources.

Writing may proceed from `docs/MANUSCRIPT_JCIM_EN.md`, `docs/MANUSCRIPT_JCIM_ZH.md`, SI, `figures/jcim_article/`, and `results/canonical/`. Do not reintroduce 0.6607, 28/38/32 as n_scored, 6/110, 1/96, or a copied JAK document-cluster CI.

---

## Chain run (zero-dock)

```
python3 scripts/analysis/adjudicate_activity_records.py
python3 scripts/analysis/build_current_score_master.py
python3 scripts/analysis/compute_canonical_results.py
python3 scripts/analysis/compute_leave_one_document.py
python3 scripts/analysis/compute_class_chemistry.py
python3 scripts/analysis/compute_detectable_effect.py
python3 scripts/analysis/fit_ecfp4_models.py
python3 scripts/analysis/compute_descriptor_baselines.py
python3 scripts/analysis/close_publication_from_canonical.py
python3 scripts/analysis/patch_publication_text.py
python3 figures/jcim_article/scripts/update_figures_pr32.py --source-root .
python3 scripts/analysis/rebuild_submission_pack.py
python3 scripts/qa/check_current_chain.py
```

`check_current_chain.py`: **PASS**. Figures: **PASS** (12 figures; Arial). Submission pack rebuilt (103 files).

Primary uncertainty: class-stratified nonparametric percentile bootstrap, B = 2000, seed 20260729. Dual vs A-only uses pocket B; dual vs B-only uses pocket A.

---

## Primary Table 2 (activity-eligible n_scored)

| pair | n (D / A / B) | D vs A pocket B | D vs B pocket A | summary_min [95% CI] |
|---|---:|---:|---:|---|
| EGFR/HER2 | 28 / 37 / 32 | 0.6564 | 0.3237 | 0.3237 [0.1953, 0.4744] |
| JAK1/JAK2 | 32 / 32 / 32 | 0.5884 | 0.7275 | 0.5884 [0.4482, 0.7158] |
| JAK1/TYK2 | 31 / 32 / 32 | 0.5751 | 0.3649 | 0.3649 [0.2329, 0.5051] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.7143 | 0.6921 | 0.6921 [0.4801, 0.8016] |
| AChE/BChE | 27 / 26 / 28 | 0.6524 | 0.6058 | 0.6058 [0.4392, 0.7354] |
| F2/F10 | 31 / 32 / 32 | 0.4133 | 0.3448 | 0.3448 [0.2157, 0.4819] |
| PPARG/PPARA | 32 / 31 / 32 | 0.6492 | 0.7061 | 0.6492 [0.5111, 0.7460] |
| PPARA/PPARD | 32 / 32 / 32 | 0.6465 | 0.4463 | 0.4463 [0.3008, 0.5898] |

EGFR quota in Table 1 remains 28 / 38 / 32 / 12 (construction). n_scored is 28 / 37 / 32 after EH120_059 became unresolved.

---

## Confirmed downstream numbers (current scores)

- EGFR fixed-score Δ (pocket A, neither − B-only): 0.4621 [0.2596, 0.6406]
- EGFR matched − mismatched: 0.0558 [−0.0368, 0.1551] (includes 0)
- AChE matched − mismatched: 0.1770 [0.0528, 0.2906] (excludes 0; only main-panel interval excluding 0)
- EGFR two-pocket mean D-vs-neither: 0.7589 [0.5535, 0.9256]
- EGFR Top 10% 1 / 5 / 5 / 0; EF_dual,10% = 0.3539 on n_ranked = 109
- AChE Top 10% 5 / 3 / 1 / 1; EF_dual,10% = 1.7593 on n_ranked = 95
- Max vs median flips: EGFR 5/109; AChE 1/94 (CHEMBL659); PPARA/PPARD 1/110
- EGFR document cluster Δ: 0.4621 [0.1246, 0.6441], n_valid_boot = 1994
- JAK document cluster: `unresolved_mapping_unavailable` (not recomputed)
- EGFR leave-one-document: 38 groups, none sign-flip the directional claim; G001 n_dropped = 43, Δ change = −0.1766, still positive
- ECFP4 incremental max |Δ| = 0.0112 (JAK1/JAK2 D_vs_A)
- Detectable-effect: P(CI excludes 0.5 | true AUROC = 0.60) ≤ 0.068; at 0.75, min among non-PIK3CA pairs = 0.817; PIK3CA/mTOR = 0.453 (binormal; not observed power)
- Independent GNINA EGFR: 28 / 37 / 32 / 11
- 4JT6 resolution 3.60 Å vs census ≤ 3.5 Å supply gate: stated in Methods

---

## What writing may say

Retained conclusions that still have evidence:

- Dual-versus-neither two-pocket mean AUROCs can look favorable while a directional arm is near or below 0.5 (EGFR/HER2, JAK1/TYK2).
- Control class, not only score aggregation, moves the readout (Figure 2A).
- Adding the matched-pocket docking score to ECFP4 changes OOF AUROC by at most 0.011 on these panels.
- Matched-pocket advantage is not supported on EGFR/HER2; AChE/BChE is the only main-panel matched−mismatched interval excluding 0.
- No pair met the independent external docking gate.

Do not say JAK1/TYK2 document-cluster currently includes 0 from a recomputed interval. Do not say AChE max/median is 1/96. Do not say EGFR n_scored A-only is 38.

---

## Non-authoritative files

- `figures/jcim_article/FIGURE_AUDIT_PR32.md`
- `data/manuscript_lock/ARTICLE_ASSET_INDEX_v1.csv`
- `remediation_outputs/freeze_audit/FINAL_FREEZE_REPORT.md`
- `remediation_outputs/freeze_audit/SECOND_PASS_FREEZE_REPORT.md`
- `data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv` (pre-adjudication index)
