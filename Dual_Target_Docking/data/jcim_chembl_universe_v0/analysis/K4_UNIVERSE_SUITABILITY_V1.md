# Are the four paper pairs the most suitable dual-target pairs in ChEMBL?

**Short answer: no.** The four DualFourClass pairs were **not** chosen as the top four pairs from an exhaustive ChEMBL search. They are the pairs that survived a **pre-specified literature/pathway candidate list** (J0: 49 auditable pairs among P01–P52) plus metal / isoform / pose-gold filters. A ChEMBL 37 dump-level census with the **same numeric gates** shows both that three of the four sit in a scarce conventional thick set **and** that the candidate list was incomplete.

This note does not dock anything, does not unfreeze Table 2, and does not expand K = 4.  
Table: `tables/k4_vs_universe_suitability_v1.csv`  
Script: `scripts/k4_universe_suitability_v1.py` (reads existing universe CSVs; no SQLite).

## What was actually screened before K=4

| Step | Sampling frame | Result |
|------|----------------|--------|
| 2026-07-23 public-pair audit | ~12 literature / pathway / pose-gold pairs with REST caches | Hard-gate table for those pairs only |
| 2026-07-29 J0 | Hand-listed P01–P52 among already-cached `mols_*.json` | **49** auditable pairs; **4** thick (min hard-neg ≥ 50); HDAC1/HDAC6 dropped as metal |
| Frozen K=4 | J0 survivors + one intentional thin case | EGFR/HER2, AChE/BChE, PIK3CA/PIK3CB, PIK3CA/mTOR |

J0 is a candidate list. It is **not** “ChEMBL only has 49 dual-target pairs.”

## Universe with the same numeric gates (ChEMBL 37)

Human `SINGLE PROTEIN`, one component, max pChEMBL, same endpoints and 6.5/5.5 strict labels as J0:

| Gate | N |
|------|--:|
| Mapped targets | 4,672 |
| Unordered pairs possible | 10,911,456 |
| n_both ≥ 10 | 63,790 |
| Directional (dual/A/B ≥ 10 at θ = 6.0) | 5,253 |
| Strict thick (min hard-neg ≥ 50) | **86** (78 non-metal) |

Of the 86, 46 are qHTS / common counter-screens, 6 CYP ADME, 8 metal enzymes, 22 same protein class, and **4** annotated `cross_class`. Two of those four are PPAR subtype pairs (same nuclear-receptor family, different ChEMBL subclass strings). The remaining conventional cross-class thick pairs in this dump are **CREBBP/BRD4** (min HN = 270) and **PIK3CA/mTOR** (min HN = 80).

## Pair-by-pair versus the universe

| Pair | In J0? | K=4 role | Thick? | Rank / 86 | min HN | Universe reading |
|------|:------:|----------|:------:|----------:|-------:|------------------|
| **PIK3CA/mTOR** | yes | development | yes | 34 | 80 | Among the scarce conventional cross-class thick pairs, and the only such pair that was **both** on the J0 list **and** pose-gold (PI-103). Not unique in the dump. |
| **AChE/BChE** | yes | development (homolog) | yes | 37 | 78 | Valid literature dual with thick supply. Several GPCR/transporter homologs are thicker (e.g. CNR1/CNR2 min HN = 246). J0 did not rank-1 homologs. |
| **PIK3CA/PIK3CB** | yes | isoform **control** | yes | 74 | 56 | Thick, but labelled `too_close_for_primary`. PIK3CG/PIK3CB is thicker (min HN = 96) and was never a J0 candidate. |
| **EGFR/HER2** | yes | **case / supply-limited** | **no** | — | **7** | Directional at θ = 6.0 (757th of 5,253 by min HN) but fails the thick gate. Kept for pose-gold (TAK-285), not as a ChEMBL-wide thick optimum. |

J0’s four thick pairs (HDAC1/HDAC6, PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB) **reproduce exactly** in the dump. HDAC1/HDAC6 (min HN = 93) was correctly excluded from K=4 as Zn-metal + isozyme.

## What is a real hole (incomplete screening)

These are limitations of the **candidate list and universe definition**, not evidence that Table 2 was fabricated.

1. **Candidate-list incompleteness.** CREBBP/BRD4 is a conventional thick pair (min HN = 270, 3rd of 86) and was **never in J0**. BRD4 was cached; CREBBP was not in `EXTRA_TARGETS_FOR_QUEUE`. This census did **not** score holo/cognate/metal docking feasibility for CREBBP/BRD4. Do not dock it in this paper.
2. **Homologs were not ranked.** CNR, HCRTR, SLC6, MAO, adenosine, opioid, and JAK isoform pairs can match or beat AChE/BChE on min HN. K=4 did not pick the thickest homolog.
3. **Isoform controls were not ranked.** PIK3CG/PIK3CB is thicker than PIK3CA/PIK3CB. The isoform seat is by design, not by universe rank.
4. **EGFR/HER2 would lose a thick-supply ranking** to dozens of GPCR/isoform pairs. Its seat is a pose-gold case, not a supply optimum.
5. **Universe definition.** Complexes, protein-family targets, and non-human proteins are out of frame. Kinase selectivity panels inflate the 5,253 directional pairs; those are not designed dual-target sets.
6. **No PDB / cognate / metal prefilter on the 86.** A “most suitable for DualFourClass docking” ranking was never computed. Fetch-queue literature pairs (AXL/MERTK, FLT3/KDR, FGFR1/KDR, SRC/ABL1, BRAF/MEK1, HDAC hybrids) exist in the dump; **none are thick**. HDAC hybrids would be metal-excluded anyway.

## What is not a hole

- J0’s 49-pair recount (Table S44) is internally consistent with this dump for the pairs it listed.
- Metal HDAC pairs were excluded on purpose.
- K=4 was never advertised in the English methods as a ChEMBL-wide top-4 (3.1 already says “49 candidate pairs”).
- This census is post-hoc. It does not rebuild panels.

## How to say this in the paper

**Allowed**

- Frozen J0 is a literature/pathway candidate list of 49 pairs, not an exhaustive ChEMBL search.
- In ChEMBL 37, 86 unordered human SINGLE PROTEIN pairs meet the same strict thick gate; most are qHTS, CYP, metal, or close homologs.
- After those scientific filters, conventional thick duals remain scarce. Three of the four K=4 pairs sit in that scarce set. EGFR/HER2 is an intentional supply-limited case (min HN = 7).
- CREBBP/BRD4 is a conventional thick pair that the candidate list missed.

**Forbidden**

- “These four are the unique best dual-target pairs in ChEMBL.”
- “ChEMBL contains only 49 dual-target pairs.”
- “Public data have only four thick pairs” as a **universe** claim.
- Docking CREBBP/BRD4, the 86, or the 5,253. This note does not unfreeze Table 2 or K = 4.
