# Track A — Figure & claim pack (Mol. Inf. / JCAMD ready; JCIM = C1/C2/C4)

**Headline is measurement + public-data supply ceiling — not “same-series extrapolation failed.”**

## Claim ceiling

### Allowed

| Claim | Evidence pointer |
|-------|------------------|
| Pooled Dual-vs-A∪B AUROC can cancel opposing directions | M1; EH unified: D/A 0.680 vs D/B 0.282 |
| On EGFR/HER2, docking ≤ trivial baselines under **unified RDKit prep** | `eh110_unified_prep_directional.csv` — vina `min` 0.282 < cLogP 0.482 |
| On PIK3CA/mTOR, docking beats volume baselines (directional) | Stage M M3 — vina summary_min 0.597 > heavy 0.463 |
| Conclusions are pair-dependent; do not average pairs | M3 gates differ |
| RTM absolute scores are prep-sensitive | M4-min LigPrep→RDKit |
| Strict four-class panels are supply-limited in public ChEMBL | J0: 4/49 Y; EGFR B_only_strict=7 |
| Labels are in principle distinguishable (oracle ceiling) | Stage M M2 noise ceiling |

### Forbidden

- `rtm_min_z` / RTM as validated universal primary arm  
- Mixed-prep panel120 RTM split as method conclusion  
- “General dual-target decision ruler already validated”  
- Wet-lab / passenger-moiety as main claim  
- Track B Full-Go / mass new-pair docking without approval  

## Figure plan

| Fig | Title | Data file(s) | Axes / content | One-line caption |
|-----|-------|--------------|----------------|------------------|
| **1** | Task definition | schematic | Four classes dual / A_only / B_only / neither; θ vs strict | DualFourClass asks for ranking duals above single-end hard negatives on *both* ends. |
| **2** | Directional decomposition | `stage_m_v0/tables/m1_directional_auroc.csv`; `eh110_unified_prep_directional.csv` | Bars: D/A vs D/B for EH (unified) & PM; pooled as ghost | Pooling hides an EGFR end-inversion that survives unified ligand prep. |
| **3** | Baselines vs docking | `m3_baselines_vs_arms.csv`; EH110 directional | `summary_min` docking vs heavy/MW/cLogP/TPSA by pair | EGFR docking loses to physicochemical baselines; PIK3CA/mTOR does not. |
| **4** | Label gray / strict supply | `j0_strict_label_supply.csv` (+ feasibility 12-pair subset) | Forest of min_strict_hardneg; mark Y/T/−; highlight EGFR=7 | Across 49 audited pairs only three non-metal pairs support thick strict panels. |
| **5** (opt.) | Prep sensitivity | `stage_m_v0/tables/m4_directional_by_prep.csv` | EH40 LigPrep vs RDKit directional | Unifying prep changes RTM absolute discrimination; protocol must freeze prep. |
| **6** (opt.) | Failure typology pointer | PIK3CA/mTOR typology docs | T1/T2/T5 case counts | Hard negatives that still score well have recurring pose/score failure modes. |

## Venue mapping

| Venue | Use this pack as |
|-------|------------------|
| Mol. Inf. / JCAMD | Full diagnostic article skeleton |
| JCIM Article | C1 (task/pooling) + C2 (supply audit) + C4 (prep); C3 needs approved K=4 docking |

## Reproduce pointers

- Stage M: `data/stage_m_v0/`  
- EH110 unified: `data/jcim_feasibility_v0/` (copied into this pack’s tables)  
- J0/J1: `data/jcim_j0j1_v0/`
