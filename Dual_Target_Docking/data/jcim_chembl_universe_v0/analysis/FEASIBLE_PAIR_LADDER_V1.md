# How many dual-target pairs are actually benchmarkable?

The answer is a **gate ladder with a number at every rung**, not a single figure. Every rung is reproducible from a committed script; no pair is dropped without a written reason.

| Rung | Gate | Pairs | Script |
|------|------|------:|--------|
| 0 | Human `SINGLE PROTEIN`, one component, ≥1 qualifying pChEMBL molecule | 4,672 targets → 10,911,456 possible pairs | `chembl_exhaustive_pair_census_v1.py` |
| 1 | ≥1 ligand measured on both ends | 2,164,618 | same |
| 2 | `n_both ≥ 10` | 63,790 | same |
| 3 | Four-state directional at θ = 6.0 (dual / A-only / B-only each ≥ 10) | 5,253 | same |
| 4 | + neither ≥ 10 (`formulation_n10`) | 4,564 | same |
| 5 | Strict 6.5/5.5 thick: min hard-negative ≥ 50 both sides | **86** | same |
| 6 | Drop qHTS hub proteins, CYP ADME panels, Zn-metal enzymes | 26 | `universe_structure_feasibility_v1.py` |
| 7 | Human holo PDB ≥ 5 per end (≤ 3.5 Å, non-polymer ligand) | **19** | same |
| 8 | Hard negatives still ≥ 50 after keeping only drug-like small molecules | **17** | `pair_ligand_identity_qc_v1.py` |
| 9 | Independent target systems among those 17 (shared targets collapsed) | **12** | same |

**Headline answer: 17 pairs, spanning 12 independent target systems.** Everything above rung 5 is arithmetic; rungs 6–8 are declared rules, not taste.

## Rung 8 — ligand identity check (new)

Molecules were kept only if ChEMBL `structure_type = MOL`, `molecule_type = Small molecule`, MW 150–750, heavy atoms 10–60, and no metal in the formula. Two pairs fall below the thick gate once non-small-molecules are removed:

| Pair | min HN raw → small-molecule only |
|------|----------------------------------|
| OPRM1/OPRK1 | 56 → **46** |
| JAK3/TYK2 | 51 → **48** |

Scaffold diversity of the dual sets was also checked (RDKit Bemis–Murcko). No surviving pair is a single congeneric series: the largest single-scaffold share is S1PR3/S1PR1 at 0.195, and 15 of 17 are below 0.08. Thick supply here is chemically diverse, not one patent family.

## The 17

| Pair | min HN (small-molecule) | dual scaffolds | Pocket class |
|------|------------------------:|---------------:|--------------|
| CREBBP/BRD4 | 251 | 36 | HAT+bromodomain / bromodomain |
| CNR1/CNR2 | 225 | 479 | GPCR |
| HCRTR1/HCRTR2 | 208 | 814 | GPCR |
| SLC6A4/SLC6A3 | 130 | 429 | SLC6 |
| F2/F10 | 108 | 192 | serine protease |
| SLC6A2/SLC6A4 | 92 | 481 | SLC6 |
| JAK1/TYK2 | 91 | 766 | kinase ATP |
| PPARG/PPARA | 82 | 205 | NR LBD |
| PPARA/PPARD | 82 | 84 | NR LBD |
| **PIK3CA/mTOR** | 71 | 536 | kinase ATP |
| **AChE/BChE** | 65 | 327 | hydrolase gorge |
| OPRM1/OPRD1 | 57 | 744 | GPCR |
| CTSK/CTSS | 57 | 255 | cysteine protease |
| F2/PRSS1 | 56 | 87 | serine protease (trypsin is usually an antitarget) |
| OPRD1/OPRK1 | 55 | 598 | GPCR |
| S1PR3/S1PR1 | 55 | 92 | GPCR |
| JAK1/JAK2 | 53 | 2,275 | kinase ATP |

Bold = already docked in this project.

## Why 17 is not "n = 17"

The 17 collapse into **12** independent systems: SLC6A4 appears twice, F2 twice, JAK1 twice, PPARA twice, and the opioid pairs share receptors. Any cross-pair statistic must be computed on the 12 systems (or with a system-level cluster bootstrap), not on 19 or 17 as if independent.

Composition matters too: only **PIK3CA/mTOR** and **CREBBP/BRD4** are genuinely cross-family. The rest are paralog/isoform selectivity pairs. A suite built from these measures **within-family selectivity**, which is a narrower question than cross-pathway polypharmacology. Say that explicitly rather than calling it a general dual-target benchmark.

## What still cannot be automated

Receptor construct choice (active vs inactive GPCR state, truncation, chimera risk), cognate redocking RMSD, and CREBBP domain selection remain human decisions. Rung 7 counts structures; it does not prove any single structure is dockable. See `RECEPTOR_IDENTITY_AUDIT_V1.md` for what happens when that step is skipped.

## Forbidden

- Do not present 17 (or 19, or 86) as a docked benchmark. Nothing here was docked.
- Do not call the 17 independent replicates; use the 12 systems.
- Do not claim the 17 measure cross-pathway dual-target design.
