# jcim_supply_crossdb_v0 — count-level BindingDB / PubChem check

Zero docking. Frozen K=4 pairs only (PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB, EGFR/HER2).

Recompute strict 6.5/5.5 class counts on BindingDB REST and PubChem PUG REST, then compare with the frozen ChEMBL `mols_*.json` cache (J0). Ligand identity is per database (BindingDB `monomerid`; PubChem CID). No InChIKey union and no panel rebuild.

```bash
python3 Dual_Target_Docking/data/jcim_supply_crossdb_v0/scripts/bindingdb_pubchem_strict_count_v1.py
```

Outputs:

- `tables/crossdb_strict_supply_v1.csv`
- `tables/fetch_log_v1.json`
- `analysis/SUPPLY_CROSSDB_VERDICT_V1.md`
- compact pmax JSON under `tables/` (as_is and equal_only)
- raw API dumps under `cache/` (gitignored)

Primary manuscript comparison is **equal_only** (closer to ChEMBL pChEMBL). `as_is` includes `>` censored values and can inflate hard-negatives.
