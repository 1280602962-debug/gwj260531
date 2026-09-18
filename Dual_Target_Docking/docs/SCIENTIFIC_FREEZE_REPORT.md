# DualFourClass freeze report

Baseline SHA (re-fetched): `17435413410290960da3437de640509705f00b51`  
Working branch: `cursor/scientific-freeze-c7cc` (PR #39, draft)  
Not merged to main. Other Dual_Target PRs not closed. Branches not deleted.

## Closed

| ID | Outcome |
|---|---|
| C1 | EN/ZH/SI incremental |Δ| = **0.023** from freeze rebuild (0.0234, PPARA/PPARD D vs B). Deposited 0.0112 is not reproducible on this env with the original fitter and identical master. |
| C2 | `compute_canonical_results.py --outdir` now writes two-pocket CIs, three GNINA pairs, and receptor substitution. Freeze rebuild did not call pack generators. |
| C3 | Receptor-substitution `primary_summary_min` read from rebuilt Table 2 (still 0.6921). |

## Checked, not a primary-number bug

| ID | Outcome |
|---|---|
| R1 | AChE max/median n=94 vs Table 2 n=95 kept. Qualified-record intersection vs eligible complete-case. |
| R2 | PR #35 EGFR 0.760/9.505 Å not adopted. Current 3POZ CalcRMS 1.019 Å. |
| R3 | Holdout eligibility gate added; holdout n unchanged (0 ineligible). |
| R4 | Five-seed archive 0.6607 remains archive-only. |
| R5 | `verify_freeze_rebuild.py` PASS on freeze dir (labels, n, score channel, folds, points, CIs). |

## Limitations unchanged

JAK document-cluster unresolved; no invented protonation protocol; no multiplicity adjustment; detectable-effect is simulation not observed power; GNINA n can be smaller than Vina n; AB_056 in Table 2 only.

## Result changes vs deposited canonical

Allowed because evidenced. Table 2 / ranking / GNINA / detectable-effect **did not move**.

- ECFP incremental max |Δ|: 0.0112 → **0.0234** (rdkit 2026.03.6, sklearn 1.9.1; 492/930 fold IDs differ). Claim: docking still adds no stable increment; the three-decimal cap is 0.023.
- JAK1/TYK2 cLogP `summary_min`: 0.5806 → 0.5796 (RDKit). SI three-decimal 0.581 → 0.580.

Figure 3A/B and SI Table S5 now use the freeze CSVs. Leftover 0.0112 files are listed in `docs/HISTORICAL_LEFTOVER_FILES.md` and are not current sources.

## Rebuild

```
python3 Dual_Target_Docking/scripts/check_analysis_env.py
python3 Dual_Target_Docking/scripts/analysis/rebuild_freeze.py
python3 Dual_Target_Docking/scripts/qa/verify_freeze_rebuild.py
python3 Dual_Target_Docking/scripts/analysis/promote_freeze_to_canonical.py
python3 Dual_Target_Docking/scripts/qa/check_current_chain.py
```

Env: `results/freeze_rebuild/ENV.txt` (Python 3.12.3, numpy 2.4.4, sklearn 1.9.1, rdkit 2026.03.6).  
No old CSV was copied into the freeze directory. Detectable-effect 192 rows matched the deposited file exactly.

## Writing

Use `docs/WRITING_INDEX_FREEZE.md`. Core claims (eight pairs, θ=6.0, eligible complete-case, directional pockets, scheme B, EGFR 0.324/0.656, control-class Δ, AChE-only matched−mismatched) stand. Incremental docking claim is 0.023 in this environment.

## Integration

`docs/MAIN_INTEGRATION_PLAN.md`. `docs/PR35_ITEM_DECISIONS.md`. Confirm before merge/close/delete.
