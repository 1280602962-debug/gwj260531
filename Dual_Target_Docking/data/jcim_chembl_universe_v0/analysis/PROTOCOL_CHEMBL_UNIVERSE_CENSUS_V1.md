# Protocol — exhaustive ChEMBL human SINGLE PROTEIN pair census

**Status:** frozen before inspecting universe pair counts.  
**Does not dock. Does not replace Table 2, Table S44, or K = 4.**  
**Source dump:** ChEMBL 37 SQLite (`chembl_37_sqlite.tar.gz`; SHA-256 `33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`).  
**Script:** `data/jcim_chembl_universe_v0/scripts/chembl_exhaustive_pair_census_v1.py`

This census answers a different question from J0. J0 counted four-state supply on a **hand-listed literature/pathway candidate set** (49 auditable pairs among cached maps). This run enumerates **all unordered pairs** among human SINGLE PROTEIN targets in one frozen dump, using the same activity and label rules as `run_j0_supply_audit.py`.

## Universe

| Gate | Rule |
|------|------|
| Target type | `target_dictionary.target_type = 'SINGLE PROTEIN'` |
| Organism | `target_dictionary.organism = 'Homo sapiens'` |
| Components | exactly one `target_components` row (true single-chain mapping) |
| Activities | `pchembl_value IS NOT NULL` |
| Endpoints | `standard_type ∈ {IC50, Ki, Kd, EC50, Potency, IC50app, Ki app}` |
| Aggregation | **max** pChEMBL per (target, molecule) |
| Untested | not inactive |
| Primary slice | no assay `confidence_score` cut (matches J0 REST harvest) |
| Sensitivity | same maps restricted to `confidence_score ≥ 8` |

pChEMBL in ChEMBL is a curated −log10 scale on standardised concentration–response values. Assay systems are not harmonised. This is a supply census, not assay-equivalent ground truth.

## Pair labels (identical to J0)

θ = 6.0 (primary four-state, same as Table S44):

- dual: both ≥ 6.0
- A-only: A ≥ 6.0 and B < 6.0
- B-only: B ≥ 6.0 and A < 6.0
- neither: both < 6.0

Strict 6.5/5.5 (same as J0 thick-panel gate):

- strict dual: both ≥ 6.5
- strict A-only: A ≥ 6.5 and B ≤ 5.5
- strict B-only: symmetric
- strict neither: both ≤ 5.5
- gray: both measured, otherwise

## Reported gates

| Name | Definition |
|------|------------|
| `n_both ≥ 10` | at least 10 ligands measured on both targets |
| `directional_n10` | dual, A-only, and B-only all ≥ 10 at θ = 6.0 |
| `formulation_n10` | directional_n10 and neither ≥ 10 |
| `supports_strict_panel` | min(strict A-only, strict B-only) ≥ 50 |
| `supports_thin_panel` | min(strict A-only, strict B-only) ≥ 20 |

Metal-enzyme **flag** (name/class contains HDAC, carbonic anhydrase, or matrix metalloproteinase) is diagnostic. It does not drop pairs from the universe tables.

## Crosswalk

J0 P01–P52 and `j0_fetch_queue.csv` extra targets are located in this dump by ChEMBL ID or UniProt and reported even if they fail a gate. Version drift versus the 2026-07-23 REST caches is expected and must be named.

## Forbidden uses

- Do not write “ChEMBL contains only 49 dual-target pairs.”
- Do not write that K=4 is the unique best dual-target set in ChEMBL.
- Do not replace Table S44’s 49-pair J0 recount with this universe count without saying the sampling frame changed.
- Do not dock every pair that passes `directional_n10` or `supports_strict_panel`.
- Do not expand K = 4 or unfreeze Table 2 from these counts.
- Do not treat kinase-selectivity panels as independently designed dual-target datasets.
