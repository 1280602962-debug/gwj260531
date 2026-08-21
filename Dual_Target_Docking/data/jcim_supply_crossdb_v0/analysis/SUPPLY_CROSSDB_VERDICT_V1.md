# SUPPLY_CROSSDB_VERDICT_V1 — BindingDB / PubChem count-level check

**Date (UTC):** 2026-08-21T05:11:40Z
**Script:** `scripts/bindingdb_pubchem_strict_count_v1.py`
**Scope:** frozen K=4 pairs only; **zero docking**; per-database ligand IDs (no cross-DB structure merge).

## Rule (identical to J0 construction gate)

- Types: IC50 / Ki / Kd / EC50.
- Representative potency = **max** converted p-activity per ligand–target.
- p from nM: `9 − log10(nM)`; PubChem concise values are µM → `6 − log10(µM)`.
- dual: both ≥ 6.5; A_only: A ≥ 6.5 and B ≤ 5.5; B_only: converse.
- Gate: `min(strict_A_only, strict_B_only) ≥ 50` (thin: ≥ 20).
- **as_is:** strip `>`/`<` and use the numeric value (can inflate hard-negs; `>10000 nM` → p = 5.0).
- **equal_only:** keep only `=` (or unqualified) records. This is the apples-to-apples comparison with ChEMBL pChEMBL, which usually requires a standard `=` relation.
- ChEMBL column is the frozen `mols_*.json` max-pChEMBL cache (J0); no inequality mode.

## Fetch status (not fabricated)

- BindingDB REST `getLigandsByUniprots` (cutoff = 1000000 nM): **OK**
- PubChem PUG REST `protein/accession/…/concise/CSV`: **OK**

Per-target fetch metadata is in `tables/fetch_log_v1.json`. Raw dumps stay in `cache/` (gitignored); compact pmax maps are in `tables/`.

## Counts

| pair | source | rule | both | strict dual | strict A/B | min HN | ≥50 both-side | ≥20 thin |
|------|--------|------|-----:|------------:|------------|-------:|:-------------:|:--------:|
| PIK3CA/MTOR | ChEMBL_cache | pChEMBL | 2713 | 1552 | 80/81 | 80 | Y | Y |
| PIK3CA/MTOR | BindingDB | as_is | 3673 | 1610 | 389/151 | 151 | Y | Y |
| PIK3CA/MTOR | BindingDB | equal_only | 2739 | 1579 | 76/96 | 76 | Y | Y |
| PIK3CA/MTOR | PubChem | as_is | 3929 | 1637 | 405/153 | 153 | Y | Y |
| PIK3CA/MTOR | PubChem | equal_only | 2955 | 1602 | 86/93 | 86 | Y | Y |
| ACHE/BCHE | ChEMBL_cache | pChEMBL | 2537 | 687 | 189/78 | 78 | Y | Y |
| ACHE/BCHE | BindingDB | as_is | 3610 | 710 | 228/141 | 141 | Y | Y |
| ACHE/BCHE | BindingDB | equal_only | 2711 | 698 | 181/92 | 92 | Y | Y |
| ACHE/BCHE | PubChem | as_is | 3885 | 756 | 275/153 | 153 | Y | Y |
| ACHE/BCHE | PubChem | equal_only | 2916 | 742 | 214/97 | 97 | Y | Y |
| PIK3CA/PIK3CB | ChEMBL_cache | pChEMBL | 1990 | 602 | 56/67 | 56 | Y | Y |
| PIK3CA/PIK3CB | BindingDB | as_is | 3858 | 924 | 208/129 | 129 | Y | Y |
| PIK3CA/PIK3CB | BindingDB | equal_only | 2545 | 855 | 58/75 | 58 | Y | Y |
| PIK3CA/PIK3CB | PubChem | as_is | 4525 | 976 | 212/144 | 144 | Y | Y |
| PIK3CA/PIK3CB | PubChem | equal_only | 2860 | 908 | 61/74 | 61 | Y | Y |
| EGFR/HER2 | ChEMBL_cache | pChEMBL | 1751 | 951 | 39/7 | 7 | N | N |
| EGFR/HER2 | BindingDB | as_is | 3032 | 1519 | 85/92 | 85 | Y | Y |
| EGFR/HER2 | BindingDB | equal_only | 2269 | 1336 | 34/31 | 31 | N | Y |
| EGFR/HER2 | PubChem | as_is | 3213 | 1545 | 88/92 | 88 | Y | Y |
| EGFR/HER2 | PubChem | equal_only | 2068 | 1121 | 43/30 | 30 | N | Y |

## Does the public-data ceiling move?

**Primary comparison = ChEMBL pChEMBL vs BindingDB/PubChem `equal_only`.** The `as_is` rows are a sensitivity to censored `>` values, not the matched-rule headline.

- **PIK3CA/MTOR** ChEMBL min HN = 80 (A/B 80/81); BindingDB equal_only min HN = 76 (A/B 76/96; both=2739); as_is min HN = 151 (A/B 389/151); PubChem equal_only min HN = 86 (A/B 86/93; both=2955); as_is min HN = 153 (A/B 405/153).
- **ACHE/BCHE** ChEMBL min HN = 78 (A/B 189/78); BindingDB equal_only min HN = 92 (A/B 181/92; both=2711); as_is min HN = 141 (A/B 228/141); PubChem equal_only min HN = 97 (A/B 214/97; both=2916); as_is min HN = 153 (A/B 275/153).
- **PIK3CA/PIK3CB** ChEMBL min HN = 56 (A/B 56/67); BindingDB equal_only min HN = 58 (A/B 58/75; both=2545); as_is min HN = 129 (A/B 208/129); PubChem equal_only min HN = 61 (A/B 61/74; both=2860); as_is min HN = 144 (A/B 212/144).
- **EGFR/HER2** ChEMBL min HN = 7 (A/B 39/7); BindingDB equal_only min HN = 31 (A/B 34/31; both=2269); as_is min HN = 85 (A/B 85/92); PubChem equal_only min HN = 30 (A/B 43/30; both=2068); as_is min HN = 88 (A/B 88/92).

## EGFR/HER2 BindingDB qualifier diagnostic (as_is B_only)

- as_is strict B_only (HER2-selective): **92**
- of which EGFR records are **only** `>`: **49**
- of which EGFR has at least one `=` record: **43**

as-is B_only ligands whose EGFR records are exclusively '>' (typical IC50 > 10 µM panel values). These inflate hard-neg counts relative to ChEMBL pChEMBL, which usually needs '='.

### One-line verdict

**Under the matched equal-relation rule, BindingDB/PubChem do not flip the ≥50 both-side thick-panel membership of this K=4 set.** The three frozen thick pairs remain above the gate; EGFR/HER2 remains a supply-limited case rather than a newly thick panel. EGFR/HER2 remains below the ≥50 thick-panel gate under equal_only (BindingDB min HN = 31; PubChem = 30) versus ChEMBL min HN = 7; it does reach the thin ≥20 pool. The as_is counts that would pass ≥50 (BindingDB min HN = 85) are driven largely by `>` censored EGFR measurements, not by a new pool of equality-bounded HER2-selective ligands. Absolute paired counts are higher than the ChEMBL cache, so this is not a claim that the databases are identical — only that the construction gate used for the docked panels does not change. Count-only; **no new docking** and no panel rebuild.

## What this is not

- Not a merged unique-structure census (no RDKit/InChIKey union).
- Not a new docking panel and not a change to frozen K=4 ligands.
- Not species / assay-confidence filtering (matches J0 max-p aggregation).
- PubChem concise and BindingDB overlap (deposition); similar counts are expected and are not two independent censuses.

## Reproduce

```bash
python3 Dual_Target_Docking/data/jcim_supply_crossdb_v0/scripts/bindingdb_pubchem_strict_count_v1.py
```

Re-runs reuse `cache/` if present. Delete `cache/` to force a live refetch.

