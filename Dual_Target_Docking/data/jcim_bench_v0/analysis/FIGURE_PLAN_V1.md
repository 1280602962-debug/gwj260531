# DualFourClass-Bench — figure / table plan (JCIM Article)

Maps analysis outputs to a lean evaluation/benchmark manuscript.  
Claim ceiling: [`../CLAIM_CEILING.md`](../CLAIM_CEILING.md).

## Main figures

| Fig | Content | Source | Notes |
|-----|---------|--------|-------|
| **Fig 1** | Task schematic: four-class labels + directional AUROC definition | draw new | dual / A_only / B_only / neither; `summary_min` |
| **Fig 2** | Forest: `summary_min` ± 95% CI for 3 dock arms + 4 baselines × K=4 | `figures/forest_summary_min_ci_v1.png` | primary result figure |
| **Fig 3** | Baseline gate: Δ(dock − best trivial) ± CI | `figures/baseline_gate_delta_ci_v1.png` | green/gray/red = CI>0 / spans0 / CI<0 |
| **Fig 4** | Failure anatomy (2-panel): (A) D/A vs D/B asymmetry; (B) AChE TPSA class means | `asymmetry_*.csv`, `ache_*` | optional combine with Top10 bars |

## Main tables

| Table | Content | Source |
|-------|---------|--------|
| **Table 1** | Pair inventory, receptors, prep, N by class, label rule | `inventory_v1.csv` + YAML roles |
| **Table 2** | Directional AUROC point + CI (vina/rtm/gnina + best baseline) | `forest_summary_min_ci_v1.csv` |
| **Table 3** | Baseline gate Δ + CI + pass/fail flags | `baseline_gate_bootstrap_v1.csv` |
| **Table 4** | Supply audit headline (49 pairs → strict Y rare) | `../jcim_j0j1_v0/` |

## Supporting / SI

| SI | Content | Source |
|----|---------|--------|
| SI Fig S1 | Prep LigPrep vs RDKit (PM48) | `pm48_directional_by_prep_v1.csv` |
| SI Fig S2 | Threshold sensitivity curves | `threshold_sensitivity_v1.csv` |
| SI Table S1 | Full bootstrap all arms | `bootstrap_directional_ci_v1.csv` |
| SI Table S2 | Top10 hardneg composition | `top10_hardneg_bootstrap_v1.csv` |
| SI Table S3 | Continuous Spearman | `continuous_spearman_v1.csv` |
| SI Table S4 | Descriptor-by-class all pairs | `descriptor_by_class_v1.csv` |
| Zenodo | `assembled_all_pairs_long.csv` + scripts + CLAIM_CEILING | this pack |

## Suggested Results section order

1. Task + metric + baselines (Fig 1)  
2. K=4 forest with CI (Fig 2, Table 2)  
3. Baseline gate honesty (Fig 3, Table 3) — **state PM CI spans 0**  
4. Asymmetry / pooled deception (EGFR + PIK3CB)  
5. AChE TPSA confound as negative control  
6. Prep protocol sensitivity  
7. Supply limits (Table 4) → why K is small  

## Forbidden figure framing

- Do not title any panel “RTM/GNINA wins”  
- Do not show only pooled AUROC as primary  
- Do not omit trivial baselines from the main forest
