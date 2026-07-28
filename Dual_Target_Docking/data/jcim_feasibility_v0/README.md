# jcim_feasibility_v0 — JCIM route feasibility checks

Exploration-pool recomputation from frozen score tables and already-fetched ChEMBL
dictionaries. **No new docking.** Narrative:
[`../../docs/JCIM_ROUTE_ASSESSMENT_V1.md`](../../docs/JCIM_ROUTE_ASSESSMENT_V1.md).

## Scripts

| Script | Question |
|--------|----------|
| `scripts/audit_strict_label_supply.py` | Per candidate target pair, how many molecules exist in each class under θ=6 labels versus strict (6.5 / 5.5) margin labels? Does public data supply enough hard negatives for a strict four-class panel? |
| `scripts/unified_prep_eh110.py` | Assemble an all-RDKit EGFR/HER2 panel120 (M4-min re-dock for the legacy 40 + the 70 already RDKit) and recompute directional AUROCs and trivial baselines with RTM z recomputed on the assembled panel. |

## Outputs

| Table | Content |
|-------|---------|
| `tables/strict_label_supply.csv` | Class counts, gray fraction, and a `supports_strict_panel` flag (both hard-negative classes ≥ 50) for 12 candidate pairs |
| `tables/eh110_unified_prep_scores.csv` | Per-ligand unified-prep scores and physicochemical baselines for the 110 ligands |
| `tables/eh110_unified_prep_directional.csv` | Directional AUROCs, `min(D/A, D/B)` summary, and `fail_baseline` flags |

## Headline numbers

- Only **3 of 12** audited pairs have ≥50 strict hard negatives on both sides (PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB).
- **EGFR/HER2 has only 7 strict B_only molecules in all of ChEMBL**, which caps any strict four-class analysis on that pair regardless of docking budget.
- Under fully unified RDKit prep, every EGFR/HER2 docking arm still falls below the best non-docking baseline (`cLogP`, `min(D/A,D/B)` = 0.482): `vina_mean` 0.282, `rtm_min` 0.256.

## Reproduce

```bash
cd Dual_Target_Docking/data/jcim_feasibility_v0/scripts
python3 audit_strict_label_supply.py
python3 unified_prep_eh110.py
```

Requires RDKit (descriptors from panel SMILES only).
