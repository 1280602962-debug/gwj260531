# Eight-pair ChEMBL 37 dump-gated check

Same STANDARD_OK harvest, θ=6.0 max-versus-median, match tolerance 0.015, leftover unused-pool counts, and dump-level pool counts on all eight Table 2 pairs. Frozen Vina scores are unchanged. This run does **not** restock Table 2 and does **not** draw holdouts.

sqlite: `/tmp/chembl/chembl_37/chembl_37_sqlite/chembl_37.db` (ChEMBL 37; tarball SHA-256 `33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`).

## Label-source sameness

All eight pairs used the same dump join. Production harvest history is **not** identical: EGFR/HER2, AChE/BChE, and PIK3CA/mTOR remain the 2026-07-23 REST panels; the other five remain the ChEMBL 37 extract. Dump agreement does not rewrite that history.

- REST-three scored ligands joined: **253**
- REST-three dump-missing ends: **0**
- REST-three dump-max vs panel mismatches (tol 0.015): **0**
- REST-three θ=6.0 max→median class flips: **7**
- dump-five missing / mismatch: **0** / **0**
- dump-five leftover vs frozen `track_b_panel_summary_v1.csv`: **MATCH**
- all eight dump-max `summary_min` vs locked Table 2 (3 d.p. half-up): **MATCH**
- 2018 test pairs with dual/A/B each n≥10: **2** / 8 (archived; not formal SI)

## Pair parity

| pair | source | scored | missing | mismatch | flips | dump max min | Table 2 | leftover D/A/B |
|---|---|---:|---:|---:|---:|---:|---:|---|
| EGFR/HER2 | 2026-07-23 ChEMBL REST | 110 | 0 | 0 | 6 | 0.430 | 0.430 | 865/19/0 |
| AChE/BChE | 2026-07-23 ChEMBL REST | 95 | 0 | 0 | 1 | 0.606 | 0.606 | 511/149/39 |
| PIK3CA/mTOR | 2026-07-23 ChEMBL REST | 48 | 0 | 0 | 0 | 0.692 | 0.692 | 1414/64/74 |
| F2/F10 | ChEMBL 37 dump extract | 107 | 0 | 0 | 0 | 0.345 | 0.345 | 312/76/245 |
| JAK1/TYK2 | ChEMBL 37 dump extract | 109 | 0 | 0 | 0 | 0.365 | 0.365 | 1874/59/80 |
| JAK1/JAK2 | ChEMBL 37 dump extract | 110 | 0 | 0 | 0 | 0.588 | 0.588 | 5953/76/21 |
| PPARG/PPARA | ChEMBL 37 dump extract | 109 | 0 | 0 | 0 | 0.649 | 0.649 | 408/50/59 |
| PPARA/PPARD | ChEMBL 37 dump extract | 110 | 0 | 0 | 1 | 0.446 | 0.446 | 187/50/68 |

Dump-level pool (same harvest, not panel membership):

- θ=6.0 D/A/B: EGFR/HER2 1182/207/46 / AChE/BChE 986/483/225 / PIK3CA/mTOR 2002/266/236 / F2/F10 645/258/607 / JAK1/TYK2 2808/293/314 / JAK1/JAK2 7399/522/286 / PPARG/PPARA 802/322/413 / PPARA/PPARD 494/236/286
- strict 6.5/5.5 D/A/B: EGFR/HER2 951/39/7 / AChE/BChE 687/189/78 / PIK3CA/mTOR 1552/80/81 / F2/F10 389/117/281 / JAK1/TYK2 2097/94/117 / JAK1/JAK2 6316/113/53 / PPARG/PPARA 445/85/102 / PPARA/PPARD 276/84/110

The 2026-08-26 API snapshot on the REST three remains a separate sensitivity (Table S3). It is not this dump join.
