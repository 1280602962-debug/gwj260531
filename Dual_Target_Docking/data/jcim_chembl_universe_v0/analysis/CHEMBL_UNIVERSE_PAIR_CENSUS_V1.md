# ChEMBL 37 universe pair census — verdict

**Dump:** ChEMBL 37 SQLite (`chembl_37.db`; tarball SHA-256 `33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`).  
**Release:** ChEMBL_37 (2026-05-01).  
**Script:** `scripts/chembl_exhaustive_pair_census_v1.py`  
**Protocol:** `PROTOCOL_CHEMBL_UNIVERSE_CENSUS_V1.md`  
**Docking:** none. **Table 2 / K = 4:** unchanged.

This census is the ChEMBL-wide counterpart of J0. J0 counted four-state supply on a **hand-listed** literature/pathway set (49 auditable pairs). This run enumerates **all unordered pairs** among human `SINGLE PROTEIN` targets with exactly one component.

## Headline counts (primary slice)

Human SINGLE PROTEIN, one component, max pChEMBL, endpoints IC50/Ki/Kd/EC50/Potency/IC50app/Ki app, untested ≠ inactive.

| Quantity | N |
|----------|--:|
| Human SINGLE PROTEIN (1 component) | 5,869 |
| With ≥1 qualifying molecule | 4,672 |
| Unordered pairs possible | 10,911,456 |
| Pairs with ≥1 ligand measured on both | 2,164,618 |
| Pairs with n_both ≥ 10 | 63,790 |
| θ = 6.0 directional (dual/A/B all ≥ 10) | 5,253 |
| Plus neither ≥ 10 (`formulation_n10`) | 4,564 |
| Strict 6.5/5.5 thick (min hard-neg ≥ 50 both sides) | **86** |
| Thick and not metal-flagged | 78 |
| Thin or thick (min hard-neg ≥ 20) | 275 |

A `confidence_score ≥ 8` sensitivity is **identical**. In this dump every human SINGLE PROTEIN pChEMBL row already has confidence 8 or 9 (2,529,783 / 2,529,783).

## What the 86 thick pairs are

They are **not** 86 independent dual-target drug-design problems. Annotated in `universe_pairs_strict_thick_annotated_v1.csv`:

| Bucket | N | Meaning |
|--------|--:|---------|
| qHTS / common counter-screen | 46 | Huge PubChem-style maps (MAPT, LMNA, SMN1, ALDH1A1, POLB, …) sharing screened libraries |
| Same protein class | 22 | Isoforms / close homologs (CNR1/2, SLC6, adenosine, opioid, JAK, PI3K, AChE/BChE, …) |
| Metal enzyme | 8 | HDAC and carbonic anhydrase pairs |
| CYP ADME panel | 6 | CYP2D6/3A4/2C9/2C19 |
| Cross-class | 4 | Includes **PIK3CA/mTOR** (min HN = 80) and **CREBBP/BRD4** (min HN = 270) |

J0’s four thick pairs reappear with **the same counts as the 2026-07-29 REST caches**:

| Pair | n_both | min strict HN |
|------|-------:|--------------:|
| HDAC1/HDAC6 | 3,987 | 93 |
| PIK3CA/mTOR | 2,713 | 80 |
| AChE/BChE | 2,537 | 78 |
| PIK3CA/PIK3CB | 1,990 | 56 |

EGFR/HER2 is still supply-limited: n_both = 1,751, min HN = **7**, directional_n10 = yes, thick = no.

**K=4 is not a ChEMBL-wide top-4 ranking.** See `K4_UNIVERSE_SUITABILITY_V1.md`. Three of four paper pairs are thick under this dump; EGFR/HER2 is not. The thickest conventional cross-class pair in the dump is CREBBP/BRD4 (min HN = 270), which was never on the J0 list.

## Fetch queue (the ~20 targets that never entered `mols_*.json`)

All 20 extra fetch-queue targets resolve in ChEMBL 37. The API outage did **not** mean the proteins are absent; it meant the REST harvest never wrote maps.

Selected intended literature pairs (`fetch_queue_intended_pairs_v1.csv`):

| Pair | n_both | dual/A/B | min HN | directional_n10 | thick |
|------|-------:|----------|-------:|----------------:|------:|
| AXL/MERTK | 1,384 | 943/71/151 | 2 | yes | no |
| FLT3/KDR | 1,003 | 588/213/108 | 28 | yes | no |
| FGFR1/KDR | 1,749 | 1344/77/178 | 14 | yes | no |
| SRC/ABL1 | 814 | 453/125/93 | 13 | yes | no |
| BTK/EGFR | 641 | 382/21/215 | 6 | yes | no |
| BACE1/AChE | 274 | 91/66/23 | 5 | yes | no |
| BRAF/MEK1 | 37 | 22/6/3 | 2 | no | no |
| KDR/HDAC1 | 66 | 35/6/21 | 3 | no | no |
| SERT/ESR1 | 18 | 5/10/0 | 0 | no | no |
| SYK/HDAC1, WEE1/HDAC1, TOP1/HDAC1, ROCK1/HDAC1, PIM1/HDAC1 | <10 | — | — | no | no |

MERTK resolved via UniProt Q12866 to `CHEMBL5331` (the fetch-queue file listed `CHEMBL3983`).

## How to say this in the paper

**Allowed**

- Frozen J0 is a literature/pathway candidate list of 49 pairs, not an exhaustive ChEMBL search.
- In ChEMBL 37, among human SINGLE PROTEIN pairs, 86 unordered pairs meet the same strict thick hard-neg gate as J0; most are qHTS counter-screens, CYP panels, metal enzymes, or close homologs.
- After that filter, J0’s three non-metal development pairs (PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB) remain among the scarce **cross-family / conventional dual** thick pairs. They are **not** a universe top-4 ranking. EGFR/HER2 remains supply-limited (min HN = 7).
- CREBBP/BRD4 (min HN = 270) is a conventional thick pair that the J0 candidate list never included.
- Fetch-queue targets exist in the dump; several kinase pairs (AXL/MERTK, FLT3/KDR, FGFR1/KDR, SRC/ABL1) are directional at θ = 6.0 but **none of the queued HDAC hybrids or BRAF/MEK are thick**.

**Forbidden**

- “ChEMBL contains only 49 dual-target pairs.”
- “Public data have only four thick pairs” as a **universe** claim. The universe number is 86 under the same numeric gate, before scientific filters.
- “These four are the unique best dual-target pairs in ChEMBL.”
- Replacing Table S44 (49-pair J0 recount) with 5,253 directional pairs as if K were expanded.
- Docking the 86, the 5,253, or CREBBP/BRD4. This census does not unfreeze Table 2 or K = 4.

## Files

| File | Role |
|------|------|
| `universe_census_summary_v1.csv` | Primary + conf8 headline counts |
| `universe_targets_all.csv` | 4,672 mapped human SINGLE PROTEIN targets |
| `universe_pairs_n_both_ge10_all.csv` | All 63,790 pairs with n_both ≥ 10 |
| `universe_pairs_directional_n10_all.csv` | 5,253 directional pairs |
| `universe_pairs_strict_thick_all.csv` | 86 thick pairs |
| `universe_pairs_strict_thick_annotated_v1.csv` | Thick pairs + supply bucket |
| `j0_universe_crosswalk_v1.csv` | P01–P52 in this dump |
| `fetch_queue_universe_targets_v1.csv` | 20 queued targets resolved |
| `fetch_queue_intended_pairs_v1.csv` | Literature pairs the queue was meant to unlock |
| `k4_vs_universe_suitability_v1.csv` | Frozen K=4 vs universe ranks |
| `K4_UNIVERSE_SUITABILITY_V1.md` | Whether K=4 is a ChEMBL-wide optimum |

The SQLite dump is not in git.
