# plan_v2_redteam_v0 — Stage M measurement audit (v0)

Exploration-pool re-analysis of already-frozen score tables. **No new docking.**

- Narrative and redesign recommendations: [`../../docs/PLAN_V2_REDTEAM_AND_REDESIGN.md`](../../docs/PLAN_V2_REDTEAM_AND_REDESIGN.md)
- Script: `scripts/metric_decomposition_v0.py`
- Inputs (read-only): `data/egfr_her2_panel120_v0/tables/`, `data/pik3ca_mtor_panel48_v0/tables/`

## Outputs

| File | Content |
|------|---------|
| `tables/directional_auroc_v0.csv` | AUROC(dual vs A_only), AUROC(dual vs B_only) and the pooled value, per arm, per subset, for both pairs; includes trivial baselines (heavy atoms, MW, cLogP) |
| `tables/label_margin_v0.csv` | Per class, how many ligands have an endpoint activity within ±0.5 log of the pChEMBL 6.0 cutoff |

## Headline numbers

- EGFR/HER2 (N=110): `vina_mean` gives dual-vs-A_only **0.689** but dual-vs-B_only **0.311**; the pooled 0.516 is the average of a working and an inverted discrimination.
- EGFR/HER2: heavy-atom count alone reaches pooled **0.549**, above `vina_mean` (0.516) and `rtm_min_z` (0.464).
- PIK3CA/mTOR (N=48): no inversion (D/A 0.698, D/B 0.597 for `vina_mean`) and volume baselines are uninformative (~0.46) while docking arms reach 0.65–0.69.
- Hard negatives cluster at the activity cutoff: EGFR/HER2 A_only 60%, B_only 78% within ±0.5 log of 6.0, versus 7% of duals.

## Reproduce

```bash
cd Dual_Target_Docking/data/plan_v2_redteam_v0/scripts
python3 metric_decomposition_v0.py
```

Requires RDKit (used only for MW / cLogP / heavy-atom counts from panel SMILES).
