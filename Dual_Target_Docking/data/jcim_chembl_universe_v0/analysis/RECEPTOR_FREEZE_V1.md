# Receptor freeze (2026-09-04) — keep all 14; do not swap

Decision recorded after the human RCSB/UniProt/paper check, the independent mmCIF geometry pass, and a written survey of higher-resolution or chemically cleaner alternatives. **The 14 receptors stay as they are.** Changing any of them now requires a new written amendment, not a silent re-rank.

Locked table: `tables/receptor_freeze_v1.csv`.

## The 14 frozen receptors

| Protein | UniProt | PDB | Cognate for Vina | Pair(s) | Status |
|---------|---------|-----|------------------|---------|--------|
| PIK3CA | P42336 | **4L23** | X6K (PI-103) | PIK3CA/mTOR | already docked |
| mTOR | P42345 | **4JT6** | X6K (PI-103) | PIK3CA/mTOR | already docked; ATP not FRB |
| AChE | P22303 | **4EY7** | E20 (donepezil) | AChE/BChE | already docked |
| BChE | P06276 | **4BDS** | THA (tacrine) | AChE/BChE | already docked |
| F2 | P00734 | **4UDW** | N6L | F2/F10 | new |
| F10 | P00742 | **2JKH** | BI7 | F2/F10 | new |
| JAK1 | P23458 | **6N7A** | KEV | JAK1/TYK2, JAK1/JAK2 | new; reused |
| TYK2 | P29597 | **3LXP** | IZA | JAK1/TYK2 | new; JH1 not JH2 |
| JAK2 | O60674 | **8BXH** | C87 (momelotinib) | JAK1/JAK2 | new |
| PPARG | P37231 | **9V8H** | BRL (rosiglitazone) | PPARG/PPARA | new; ternary + PG08-NL |
| PPARA | Q07869 | **6LXA** | EPA | PPARG/PPARA, PPARA/PPARD | new; reused |
| PPARD | Q03181 | **5U3Q** | 7UJ | PPARA/PPARD | new |
| CTSK | P43235 | **4X6H** | **I37** (not 3XT) | CTSK/CTSS | frozen **structure record only**; not Track B Vina |
| CTSS | P25774 | **9GJ2** | reconstructed α-ketoamide 13b (not KH0) | CTSK/CTSS | frozen **structure record only**; not Track B Vina |

JAK1 and PPARA are each used on two pairs. Census G5 is 8 pairs / 6 systems. **Track B production Vina uses the eight ordinary new receptors** (F2–PPARD). CTSK/CTSS stay in this freeze as structure records (`DOCKING_PLAN_V1.md`).

## What “most suitable” meant here

Not “highest resolution in PDB” and not “most famous crystal.” The declared objective, in order:

1. Accession = the UniProt the ChEMBL labels came from; *Homo sapiens*; longest matching polymer ≥ 80 aa.
2. Cognate sits in the **declared** site (ATP / S1 / LBD / papain active site), verified by literature **and** by mmCIF geometry (`GEOMETRIC_POCKET_VERIFICATION_V1.md`).
3. Among remaining OK candidates, take the best resolution (`scripts/receptor_shortlist_v1.py`, top 12 per target in `tables/tier1_receptor_shortlist_v1.csv`).
4. Already-docked receptors stay frozen unless identity or site is wrong.
5. Written amendments only when auto-rank fails identity, site, or cognate size (5AFY fragment; 9F7W BPA; 8PPZ FRB; 9CMK glue).

Under **that** objective the 14 are the correct freeze. A different objective (binary PPAR LBD only; strictly noncovalent cathepsin holos) would pick different crystals; that path was considered and **declined**.

## Alternatives considered and declined

| Current | Alternative | Why it looks better | Why it was not taken |
|---------|-------------|---------------------|----------------------|
| 9V8H PPARG | **7AWC** 1.74 Å binary LBD + BRL (also 4EMA / 4XLD) | No PG08-NL peptide; fairer PPAR-subtype comparison | 9V8H is still the highest-resolution rosiglitazone LBD (1.39 Å); pocket occupancy confirmed; peptide is a Methods note, not a fail |
| 9GJ2 CTSS | **4P6G** 1.58 Å or **2HHN** 1.55 Å (title: non-covalent) | Ordinary noncovalent Vina; published papers | 9GJ2 is 1.15 Å at the catalytic Cys; identity and pocket are correct; ligand-prep rule already handles the adduct |
| 4X6H CTSK | none in the current top-12 | — | 6ASH is a **non-active-site** ligand; 5JA7 is **C25S + allosteric**. No equal-rule noncovalent orthosteric swap exists. I37 is already deposited non-bonded (Cys25 Sγ 2.93 Å) in 4X6H |
| 4L23 PIK3CA | **8EXL** 1.99 Å taselisib (human p110α ATP) | Better resolution than 2.50 Å | Would unfreeze a docked receptor and drop the shared PI-103 pose-gold with 4JT6 |
| 4JT6 mTOR | 4JSP 3.3 Å ATPγS; 4JSX 3.5 Å Torin2; **not** 8PPZ | Slightly better resolution | Still ~3 Å kinase holos; 8PPZ is FRB. Shared PI-103 with 4L23 is kept |
| 4EY7 AChE | 4M0E 2.0 Å dihydrotanshinone | Better resolution | Frozen clinical donepezil gorge template; 4M0E is a natural-product soak |
| 4BDS BChE | 6ZWI 1.85 Å | Better resolution | Frozen tacrine gorge template; 6ZWI cognate is a large atypical ligand |
| 4UDW F2 | 5AFY / 4UD9 | Better resolution | S1 **fragments** (MW 156 / 162). 4UDW is the first inhibitor-sized cognate |

Do **not** revive: 8PPZ, 9CMK, 5AFY, 9F7W, 6ASH, 5JA7, 2WXF, HLA-presented 9-mers.

## Verification already on this branch (do not re-do to “confirm keep”)

| Record | What it proves |
|--------|----------------|
| `tables/tier1_receptor_shortlist_v1.csv` | How auto-rank listed 12 candidates per target |
| `tables/tier1_pair_receptor_plan_v1.csv` | Pair-level top pick before human amendment |
| `analysis/SITE_VERIFICATION_EVIDENCE_V1.md` | RCSB/PDBe metadata pack |
| `tables/site_verification_log_v1.csv` | Layer-2 `PASS` on all ten new ends |
| `analysis/LAYER2_LITERATURE_SIGN_OFF_V1.md` | Human RCSB + UniProt + primary-paper check |
| `tables/receptor_span_registry_v1.csv` | Construct vs resolved vs UniProt domain |
| `analysis/COVALENT_LIGAND_PREP_V1.md` | I37 / reconstructed 13b rule |
| `analysis/GEOMETRIC_POCKET_VERIFICATION_V1.md` | mmCIF distances and pocket contacts for all 14 |
| `scripts/verify_geometric_pockets_v1.py` | Reproducible geometry script |
| this file | No-swap freeze after the alternative survey |

## Still not authorised by this freeze

Layer-3 cognate RMSD and production docking of the **five** ordinary new pairs (local Vina; `DOCKING_PLAN_V1.md`). Expanding K or editing Table 2. Covalent docking of CTSK/CTSS. Re-docking PIK3CA/PIK3CB. Docking CREBBP/BRD4, GPCRs, or SLC6.
