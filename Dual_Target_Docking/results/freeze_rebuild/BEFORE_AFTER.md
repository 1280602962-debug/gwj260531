# Freeze rebuild vs deposited canonical (informational)

Independent verification is PASS on the freeze directory. This file is not an acceptance test.

## Unchanged (identical to deposited canonical)

- `current_score_master.csv` (MD5 match)
- Table 2 directional AUROC / summary_min / scheme-B CIs
- Fixed-score Δ, ranking/EF, matched−mismatched, holdout, max/median, cluster, labels, GNINA, receptor substitution, detectable-effect (192/192 rows)

## Changed with evidence

| Item | Deposited | Freeze rebuild | Why |
|---|---|---|---|
| ECFP4 incremental max \|Δ\| | 0.0112 (JAK1/JAK2 D vs A) | 0.0234 (PPARA/PPARD D vs B) | Same master + original fitter in freeze env (rdkit 2026.03.6, sklearn 1.9.1). Fold IDs differ. Not an analysis-choice change. |
| JAK1/TYK2 best-descriptor summary_min | 0.5806 | 0.5796 | RDKit MolLogP on this rdkit version |
| two-pocket CI / JAK GNINA / receptor-sub | present in canonical only after pack closer | emitted by `compute_canonical_results.py` | C2 |

Manuscript incremental claim is now 0.023 to match the freeze environment.
Figure 3B JSON still has 0.0112 until a later figure pass.
