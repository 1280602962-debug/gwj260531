# ChEMBL 37 pair census — reproducibility

This directory documents how the manuscript count **2,164,618** (unordered pairs with ≥1 ligand measured on both targets) was produced.

## Status (post-audit remediation)

The first-level census **currently relies on the archived summary output**. A full unordered pair list and the ChEMBL 37 SQLite dump are **not deposited in this repository**, so `2,164,618` cannot be independently recomputed from files in git alone.

| Quantity | Archived value | Independently recomputable from this repo? |
|----------|----------------:|--------------------------------------------|
| Human SINGLE PROTEIN (1 component) | 5,869 | No (needs dump) |
| With ≥1 qualifying molecule | 4,672 | No (needs dump) |
| Unordered pairs possible | 10,911,456 | No (needs dump) |
| Pairs with n_both ≥ 1 | **2,164,618** | **No — archived summary only** |
| Pairs with n_both ≥ 10 | 63,790 | Checked against archived summary row |
| θ = 6.0 directional (dual/A/B all ≥ 10) | 5,253 | Checked against archived summary row |
| Strict 6.5/5.5 thick | 86 | Checked against archived summary + annotated thick-pair table |

Do **not** describe 2,164,618 as “repository independently reproduced”.

## Required input (not in git)

- ChEMBL 37 SQLite dump (`chembl_37.db`)
- Upstream tarball SHA-256 recorded in `analysis/CHEMBL_UNIVERSE_PAIR_CENSUS_V1.md`:
  `33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`
- Release: ChEMBL_37 (2026-05-01)

Download the official ChEMBL 37 SQLite release from the ChEMBL downloads page. Do not commit the database.

## Script

`Dual_Target_Docking/data/jcim_chembl_universe_v0/scripts/chembl_exhaustive_pair_census_v1.py`

Example:

```bash
python3 data/jcim_chembl_universe_v0/scripts/chembl_exhaustive_pair_census_v1.py \
  --sqlite /path/to/chembl_37.db \
  --out data/jcim_chembl_universe_v0 \
  --archive-sha256 33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281
```

## Construction (from the deposited script and protocol)

1. Human `SINGLE PROTEIN` targets with exactly one component.
2. Qualifying activity maps: max pChEMBL; endpoints IC50 / Ki / Kd / EC50 / Potency / IC50app / Ki app.
3. Untested ≠ inactive.
4. Unordered pair enumeration among targets with ≥1 qualifying molecule.
5. `n_pairs_n_both_ge_1` counts pairs with ≥1 ligand measured on both members.

## Archived outputs in this repository

- `tables/universe_census_summary_v1.csv` — headline counts including 2,164,618
- `analysis/CHEMBL_UNIVERSE_PAIR_CENSUS_V1.md`
- `analysis/PROTOCOL_CHEMBL_UNIVERSE_CENSUS_V1.md`

The complete pair-level table that would hash-verify 2,164,618 row-by-row is not stored here.
