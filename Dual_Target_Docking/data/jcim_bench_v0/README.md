# jcim_bench_v0 — DualFourClass-Bench pack

JCIM evaluation/benchmark aggregation root.

Authorization: [`../protocols/PAIR_ROLES_APPROVED_JCIM.yaml`](../protocols/PAIR_ROLES_APPROVED_JCIM.yaml)  
Claim ceiling: [`CLAIM_CEILING.md`](CLAIM_CEILING.md)  
GNINA: [`analysis/GNINA_STATUS.md`](analysis/GNINA_STATUS.md) (**SKIPPED**)

| Pair | Pack | Dock status |
|------|------|-------------|
| EGFR/HER2 | `../jcim_feasibility_v0/` + EH110 unified | ready (no new dock) |
| PIK3CA/mTOR | `../pik3ca_mtor_panel48_rdkit_v0/` | **Vina+RTM done**; prep Δ written |
| AChE/BChE | `../ache_bche_panel_v0/` | receptors frozen; **panel Vina running** |
| PIK3CA/PIK3CB | `../pik3ca_pik3cb_panel_v0/` | receptors frozen; **panel Vina running** |

Primary result workspaces (large poses):
- `/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_rdkit_v0/`
- `/mnt/d/CADD paper exercise/dual target docking/results/ache_bche_panel_v0/`
- `/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_pik3cb_panel_v0/`

Unified directional forest tables will be added when AChE + PIK3CA/B RTM finish.
