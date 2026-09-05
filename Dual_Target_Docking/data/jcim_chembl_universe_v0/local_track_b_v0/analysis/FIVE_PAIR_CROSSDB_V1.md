# Five-pair BindingDB / PubChem count-only

Same fetch and classify rules as the frozen K=4 script
`jcim_supply_crossdb_v0/scripts/bindingdb_pubchem_strict_count_v1.py`
(IC50/Ki/Kd/EC50; max converted p; `as_is` vs `equal_only`; per-database
IDs; no InChIKey merge).

ChEMBL column is a fresh ChEMBL 37 dump harvest (STANDARD_OK + max
pChEMBL), not the K=4 `mols_*.json` cache. Dump both-measured counts match
the frozen panel-summary θ-class sums (F2/F10 1868, JAK1/TYK2 3503,
JAK1/JAK2 8513, PPARG/PPARA 2058, PPARA/PPARD 1226).

All eight UniProt fetches returned HTTP 200. Raw dumps stay in
`jcim_supply_crossdb_v0/cache/` (gitignored). This file does **not**
overwrite `jcim_supply_crossdb_v0/tables/`.

**Zero docking. Not Table 2. Do not hard-dock BindingDB as external
validation.**

## equal_only (primary comparison)

| pair | ChEMBL both / min HN | BindingDB both / min HN | PubChem both / min HN | ≥50 gate flip? |
|---|---:|---:|---:|---|
| F2/F10 | 1868 / 117 | 1985 / 129 | 2163 / 147 | no |
| JAK1/TYK2 | 3503 / 94 | 4184 / 95 | 4117 / 93 | no |
| JAK1/JAK2 | 8513 / 53 | 9761 / 54 | 9700 / 54 | no |
| PPARG/PPARA | 2058 / 85 | 2026 / 84 | 2134 / 84 | no |
| PPARA/PPARD | 1226 / 84 | 1155 / 71 | 1231 / 81 | no |

Under the matched equal-relation rule, every pair stays above the ≥50
both-side hard-neg gate. Absolute paired counts are not identical to
ChEMBL (BindingDB JAK pairs are thicker; PPARA/PPARD BindingDB min HN is
lower, 71 vs 84). That is a count-level identity check, not a new panel.

`as_is` inflates hard-negatives (F2/F10 BindingDB min HN 195 vs equal_only
129) the same way it did for EGFR/HER2 in K=4. Do not treat `as_is` as a
new thick-panel justification.

PubChem concise and BindingDB overlap by deposition; similar counts are
expected and are not two independent censuses.

## What this is not

- Not a merged unique-structure census
- Not a leftover-holdout or year-split substitute
- Not a reason to dock BindingDB-only ligands
