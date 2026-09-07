# Structure feasibility of every ChEMBL 37 strict-thick pair

**Does not dock. Does not unfreeze Table 2 or K = 4.**  
**H3** is the same operational gate as `docs/PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`: RCSB experimental entries, resolution ≤ 3.5 Å, ≥1 non-polymer ligand instance; suggested pass **≥ 5 holos on each end**.  
**Script:** `scripts/universe_structure_feasibility_v1.py`  
**Tables:** `tables/universe_target_holo_v1.csv` (71 UniProts) · `tables/universe_pairs_structure_feasibility_v1.csv` (**all 86** thick pairs; none omitted)

Starting set = every pair with `supports_strict_panel` in the dump (86). No pair is dropped before a written reason.

## Layer 0 — all 86, with a reason

| Decision | N | What it is |
|----------|--:|------------|
| `exclude_qhts_counterscreen` | 46 | MAPT / LMNA / SMN1 / ALDH1A1 / TP53 / POLB / GAA / HTT / … shared PubChem libraries |
| `exclude_cyp_adme` | 6 | CYP2D6 / 3A4 / 2C9 / 2C19 panels |
| `exclude_metal` | 8 | HDAC or carbonic anhydrase (original H4) |
| `fail_H3_holo` | 7 | Remaining scientific pairs that fail ≥5 holos / end |
| `include_candidate_*` | **19** | H3 pass and not qHTS / CYP / metal |

The 46+6+8 = 60 exclusions are **identity of the proteins**, not missing PDB lookups. Every UniProt among the 86 was queried (71 unique accessions).

## The 7 H3 failures (among the 26 non-qHTS / non-CYP / non-metal)

These are the only scientific thick pairs that **cannot** enter a human-holo DualFourClass freeze under the written H3 gate.

| Pair | min HN | holo A/B (≤3.5 Å) | Human UniProt PDB xrefs | Why it fails |
|------|-------:|-------------------:|-------------------------|--------------|
| MAOA/MAOB | 128 | **4** / 57 | MAOA: 4 (all ≤3.15 Å) | Gate is ≥5. Borderline: all four MAOA entries are inhibitor holos. Would pass at ≥4. |
| ADORA1/ADORA2A | 119 | **4** / 183 | ADORA1: 5 xrefs, one EM at 3.60 Å dropped | A1 end < 5 at ≤3.5 Å |
| PIK3CG/PIK3CB | 96 | 104 / **0** | **P42338 has 0 PDB xrefs** | No human PI3Kβ crystal. Mouse 2Y3A (3.30 Å) / 4BFR (2.80 Å) only. |
| ADORA1/ADORA3 | 77 | 4 / **3** | ADORA3: 3 ≤3.5 Å + two 3.60 Å EM | Both ends short |
| ADORA2A/ADORA3 | 65 | 183 / **3** | same A3 bottleneck | A3 end < 5 |
| HTR7/HTR6 | 62 | **1** / **3** | 1 and 3 | Not enough to freeze two receptors |
| PIK3CA/PIK3CB | 56 | 106 / **0** | **P42338 has 0 PDB xrefs** | Same human-PI3Kβ hole. The current paper’s 2Y3A is a **mouse** surrogate. |

**PIK3CA/PIK3CB cannot be in a from-scratch human-H3 roster** unless you explicitly allow an ortholog receptor. That is an exception, not a H3 pass.

## The 19 that can be included (H3 pass)

Grouped by pocket class. This is the complete include set. Nothing else from the 86 remains.

### A. Conventional soluble small-molecule pockets (closest to DualFourClass H4)

| Pair | min HN | holo A/B | Shared non-solvent CCD (≤3.5 Å, upper bound) | Role if you rebuilt K from scratch |
|------|-------:|----------|----------------------------------------------|-------------------------------------|
| **F2/F10** | 117 | 429/182 | 11 | Thrombin / fXa duals; serine-protease homolog |
| **JAK1/TYK2** | 94 | 47/51 | 6 | JAK isoform dual |
| **PPARG/PPARA** | 85 | 363/78 | 20 | PPAR dual agonists (same NR1C family; not cross-class) |
| **PPARA/PPARD** | 84 | 78/54 | 8 | PPAR dual agonists |
| **PIK3CA/mTOR** | 80 | 106/41 | 0 at ≤3.5* | Pathway dual; **already in K=4** |
| **AChE/BChE** | 78 | 74/110 | 7 | Cholinesterase dual; **already in K=4** |
| **F2/PRSS1** | 62 | 429/10 | 0 | Structurally dockable; trypsin is usually an **antitarget**, not a dual-drug pair |
| **CTSK/CTSS** | 58 | 65/60 | 0 | Cathepsin isoform dual |
| **JAK1/JAK2** | 53 | 47/167 | 9 | JAK isoform dual |
| **JAK3/TYK2** | 51 | 43/51 | 4 | JAK isoform dual |

\*PI-103 (CCD **X6K**) is on 4L23 (PIK3CA, 2.50 Å) and **4JT6 (mTOR, 3.60 Å)**. The gold mTOR structure is 0.10 Å outside H3, so the automated ligand-overlap column is 0. 4JT5 (3.45 Å) is a different mTOR holo (P2X). A from-scratch freeze can still use 4JT6 if you document a 3.6 Å exception, as the current paper already did.

### B. Class A GPCR homologs (membrane; dockable, not kinase-ATP)

| Pair | min HN | holo A/B | Note |
|------|-------:|----------|------|
| CNR1/CNR2 | 246 | 26/13 | CB1/CB2 duals exist |
| HCRTR1/HCRTR2 | 210 | 14/11 | Dual orexin antagonists (e.g. suvorexant chemotypes); shared CCD includes SUV |
| OPRM1/OPRD1 | 64 | 31/16 | Opioid subtype duals |
| OPRD1/OPRK1 | 62 | 16/28 | Opioid subtype duals |
| OPRM1/OPRK1 | 56 | 31/28 | Opioid subtype duals |
| S1PR3/S1PR1 | 55 | 6/19 | S1P receptor duals; S1PR3 just meets ≥5 |

### C. SLC6 transporters (membrane, detergent/lipid holos)

| Pair | min HN | holo A/B |
|------|-------:|----------|
| SLC6A4/SLC6A3 | 131 | 14/16 |
| SLC6A2/SLC6A4 | 93 | 36/14 |

### D. Epigenetic, domain-choice required

| Pair | min HN | holo A/B | Note |
|------|-------:|----------|------|
| **CREBBP/BRD4** | 270 | 113/593 | Thickest conventional-looking cross-family pair. CREBBP = HAT **and** bromodomain; BRD4 = bromodomain. Must freeze which CREBBP domain. Not equivalent to PI3K/mTOR ATP docking. |

## Cross-check against the original 12-pair H3 table

Where the same UniProt appears, RCSB counts match the 2026-07-23 table within archive drift:

| Target | 2026-07-23 | This run |
|--------|-----------:|---------:|
| PIK3CA | 106 | 106 |
| mTOR | 41 | 41 |
| AChE | 74 | 74 |
| BChE | 108 | 110 |
| BRD4 | 584 | 593 |
| HDAC1 | 6 | 6 |

PIK3CB was **never** in `pdb_holo_counts.csv`. This run is the first time that end was counted: human = 0.

## What a from-scratch DualFourClass roster would actually be

If you keep the original scientific intent (conventional soluble pockets, human holo ≥5, no metal, no qHTS, no CYP):

**Must-include soluble set (8 pairs, dropping F2/PRSS1 as antitarget):**  
F2/F10, JAK1/TYK2, PPARG/PPARA, PPARA/PPARD, PIK3CA/mTOR, AChE/BChE, CTSK/CTSS, JAK1/JAK2  
(JAK3/TYK2 is a ninth JAK isoform pair if you want another control.)

If you also allow membrane proteins: add the 6 GPCR + 2 SLC6 above.

If you also allow epigenetic domain-matched docking: add CREBBP/BRD4 after a domain freeze.

**Cannot include without a written exception:** PIK3CA/PIK3CB (no human β structure), EGFR/HER2 (not even thick), MAOA/MAOB (MAOA has 4 holos), adenosine A1/A3 pairs, HTR6/HTR7.

This census does **not** authorize docking any of the 19, expanding K, or replacing Table 2. It closes the “never looked at PDB for the 26” hole.

## Reproduce

```bash
python3 Dual_Target_Docking/data/jcim_chembl_universe_v0/scripts/universe_structure_feasibility_v1.py
```

Queries `https://search.rcsb.org/rcsbsearch/v2/query`. Caches under `cache/rcsb_holo_v1/` (gitignored).
