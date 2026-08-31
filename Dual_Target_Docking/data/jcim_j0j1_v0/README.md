# jcim_j0j1_v0 — J0 supply audit + J1 pair draft + Track A pack

**Zero new docking.** Stage M Track B remains Weak.

## Outputs

| Path | Content |
|------|---------|
| `analysis/J0_SUPPLY_AUDIT.md` | How many pairs support strict panels |
| `analysis/J1_PAIR_SELECTION_DRAFT.md` | K=4 draft + budget (not executed) |
| `analysis/TRACK_A_FIGURE_CLAIM_PACK.md` | Claims + figure plan |
| `tables/j0_strict_label_supply.csv` | 49-pair strict audit |
| `tables/j0_candidate_pairs.csv` | 53 candidates + auditable flags |
| `tables/j0_fetch_queue.csv` | Targets to fetch when ChEMBL API recovers |
| `tables/eh110_unified_prep_directional.csv` | Copied from feasibility pack |

Related: `../protocols/PAIR_ROLES_DRAFT_J1.yaml`, `../track_a_starter_v0/`.

## Reproduce

```bash
python3 Dual_Target_Docking/data/jcim_j0j1_v0/scripts/run_j0_supply_audit.py
```

## Not done here

- J2/J3/J4/J5 docking or GNINA  
- Fetching new ChEMBL dictionaries (API HTTP 500 at run time)
