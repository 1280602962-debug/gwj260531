# Tier-1 dockable roster — definitive

One number, one rule, one list. Supersedes the loose "8 / 12 / 17 / 19" language in earlier notes.

## The rule

A pair enters Tier 1 only if it passes **all five** pre-declared gates:

| Gate | Rule | Script |
|------|------|--------|
| G1 supply | strict 6.5/5.5 min hard-negative ≥ 50 on **both** ends | `chembl_exhaustive_pair_census_v1.py` |
| G2 target class | not a qHTS hub protein, not a CYP ADME panel, not a Zn-metal enzyme | `universe_structure_feasibility_v1.py` |
| G3 structure | ≥ 5 human holo entries per end (≤ 3.5 Å, ≥1 non-polymer ligand) | same |
| G4 ligand identity | hard negatives still ≥ 50 after keeping only drug-like small molecules (`structure_type=MOL`, `molecule_type=Small molecule`, MW 150–750, heavy 10–60, no metal) | `pair_ligand_identity_qc_v1.py` |
| G5 pocket regime | conventional soluble small-molecule pocket, rigid-receptor docking defensible; pharmacologically a designed dual | `receptor_shortlist_v1.py` + declared |

## The answer

| Set | Pairs | Independent target systems |
|-----|------:|---------------------------:|
| Passes G1–G4 (full feasible universe) | 17 | 12 |
| Census G5 (conventional soluble pocket) | **8** | **6** |
| **Track B production Vina (ordinary noncovalent)** | **5 new** + 2 already docked | **3 new** + 2 already docked |
| Tier 2: deferred frontier (6 GPCR + 2 SLC6 + CREBBP/BRD4) | 9 | 6 |

Census G5 is still 8 pairs / 6 systems. **Production docking is not that list.** CTSK/CTSS passes G1–G5 on labels and human holos, but both crystals are reversible-covalent, so it is **out of the unified AutoDock Vina campaign** (`DOCKING_PLAN_V1.md`). All 14 receptors stay frozen (`RECEPTOR_FREEZE_V1.md`); 4X6H/9GJ2 remain structure records only.

| # | Pair | min HN (small-mol) | System | Status | Receptor A (verified) | Receptor B (verified) |
|---|------|-------------------:|--------|--------|----------------------|----------------------|
| 1 | PIK3CA/mTOR | 71 | PI3K–mTOR | **already docked** (ordinary thick) | 4L23 (frozen) | 4JT6 (frozen, 3.60 Å ATP site) |
| 2 | AChE/BChE | 65 | cholinesterase | **already docked** (ordinary thick) | 4EY7 (frozen) | 4BDS (frozen) |
| 3 | F2/F10 | 108 | coagulation protease | **Track B Vina** | **4UDW** 1.16 Å / N6L (not auto 5AFY fragment) | 2JKH 1.25 Å / BI7 |
| 4 | JAK1/TYK2 | 91 | JAK family | **Track B Vina** | 6N7A 1.33 Å / KEV | 3LXP 1.65 Å / IZA (**JH1**, not JH2) |
| 5 | JAK1/JAK2 | 53 | JAK family | **Track B Vina** | 6N7A 1.33 Å / KEV | 8BXH 1.30 Å / C87 |
| 6 | PPARG/PPARA | 82 | PPAR family | **Track B Vina** | **9V8H** 1.39 Å / BRL + PG08-NL peptide (not auto 9F7W BPA) | 6LXA 1.23 Å / EPA |
| 7 | PPARA/PPARD | 82 | PPAR family | **Track B Vina** | 6LXA 1.23 Å / EPA | 5U3Q 1.50 Å / 7UJ |
| 8 | CTSK/CTSS | 57 | cathepsin | Layer-2 PASS; **not Track B Vina** (covalent) | 4X6H 1.00 Å / I37 (structure record) | 9GJ2 1.15 Å / KH0 (structure record) |

Historical K=4 also docked **EGFR/HER2** (supply-limited special case) and **PIK3CA/PIK3CB** (receptor-identity special case: 2WXF = mouse p110δ). Those two are not on this G5 roster. Destination identity (`PROJECT_IDENTITY_LOCK_V1.md`): EGFR/HER2 remains a main-table row; PIK3CA/PIK3CB is withdrawn to SI as a documented receptor-identity failure. Do not re-dock either. The five Track B pairs are not a side chapter — they join the same article's analysis stack (8 main-table rows after PIK3CB withdrawal). CTSK/CTSS is still not ordinary Vina.

Receptor candidates: `tables/tier1_receptor_shortlist_v1.csv` (12 ranked candidates per target, with accession, organism, entity length, mutation count, cognate). Pair-level plan: `tables/tier1_pair_receptor_plan_v1.csv`. Locked execution: `DOCKING_PLAN_V1.md`.

Track B needs **eight** new receptors: F2, F10, JAK1, TYK2, JAK2, PPARG, PPARA, PPARD. JAK1 and PPARA are each reused. CTSK and CTSS stay frozen but are not production receptors.

## Excluded, with reasons (no silent drops)

| Pair | Fails | Why |
|------|-------|-----|
| PIK3CA/PIK3CB | **G3 and G4** (new selection) | Human PIK3CB (P42338) has zero PDB entries; hard negatives fall 56 → 40 under G4. The protocol would never have selected it for a *new* ordinary thick pair. **Existing docking is kept** as a declared special case, like EGFR/HER2: human-PIK3CB labels scored in mouse p110δ (2WXF / O35904). Not an isoform-control result. Do not re-dock. |
| EGFR/HER2 | G1 | min hard-negative = 7. Retained only as a declared supply-limited case, never as a thick panel. |
| F2/PRSS1 | G5 | Passes G1–G4, but trypsin is an antitarget, not a designed dual partner. |
| OPRM1/OPRK1, JAK3/TYK2 | G4 | 56 → 46 and 51 → 48 once non-small-molecules are removed. |
| 6 GPCR + 2 SLC6 pairs | G5 | Membrane proteins; construct/conformational-state choice is a separate methodological problem. |
| CREBBP/BRD4 | G5 | CREBBP has both a HAT site and a bromodomain; the target domain must be declared before it is dockable. |
| MAOA/MAOB, adenosine A1/A3 pairs, HTR6/HTR7, PIK3CG/PIK3CB | G3 | Too few human holo structures per end (MAOA has 4, A3 has 3, HTR7 has 1, PIK3CB has 0). |

## Receptor selection cannot be fully automated

The shortlist verifies accession, organism, entity length, mutation count, and cognate drug-likeness. It **cannot** verify that the cognate sits in the intended pocket. Two Tier-1 targets prove it:

- **mTOR**: every entry at ≤ 2.5 Å is an FKBP–rapamycin/FRB complex. The automated top pick, 8PPZ (1.85 Å, cognate 0AN), is the allosteric FRB site, **not** the ATP pocket. The correct ATP-site receptor is 4JT6 at 3.60 Å — worse resolution, right site.
- **PIK3CA**: ranks 4–10 at ≤ 2.5 Å are HLA class I + β2-microglobulin complexes in which "PIK3CA" is a presented **9-mer peptide**, and the top pick 9CMK is a molecular-glue structure already blacklisted in `FROZEN_PUBLIC_PAIRS.yaml`. The correct receptor remains the frozen 4L23.

Two entity-handling traps were also found and fixed: presented peptides pass an accession match (rejected by a ≥ 80-residue protein cut), and zymogen-derived proteases split catalytic and light chains across entities sharing one accession, so the **longest** matching entity must be used (thrombin 4UDW: **258 aa** PDB heavy 364–621, not the 28 aa light chain; UniProt mature heavy is 364–622).

Therefore the receptor protocol is **three layers, in this order**:

1. **Identity** — accession == the accession the labels came from; human; protein entity ≥ 80 aa (longest matching entity); record mutation count and complex partners.
2. **Site** — declared `intended_site` per target, verified by a human against the cognate's location. Metadata pack: `SITE_VERIFICATION_EVIDENCE_V1.md`. Human + literature PASS: `LAYER2_LITERATURE_SIGN_OFF_V1.md`. Log: `tables/site_verification_log_v1.csv` (all ten new ends `PASS` as of 2026-09-04). Numbering: `tables/receptor_span_registry_v1.csv`. Covalent ligand rule: `COVALENT_LIGAND_PREP_V1.md`.
3. **Cognate redocking** — best-of-9 heavy-atom RMSD gate, run only after 1 and 2 pass.

Layer 3 alone passed the mouse p110δ receptor at 0.405 Å. Layers 1 and 2 are what make layer 3 meaningful.

2024–2026 literature URLs (DualDiff, PLINDER, LIT-PCBA audit, TopU-LBVS, DTDL review): `LITERATURE_2024_2026_DUAL_BENCHMARKS_V1.md`.

## Not authorised by this document

Production Vina of the five ordinary pairs (local Vina + Layer-3 cognate RMSD). Expanding K or editing Table 2. Covalent docking of CTSK/CTSS. Re-docking PIK3CA/PIK3CB. Swapping to 7AWC, 4P6G, 8EXL, or 8PPZ without a new written amendment. See `DOCKING_PLAN_V1.md`.
