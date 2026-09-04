# Site verification checklist (Layer 2) — before any Vina job

Mandatory gate between Layer 1 (identity) and Layer 3 (cognate best-of-9 RMSD). Metadata ranking, resolution, and `site_verified_by_human=0` in `tables/tier1_receptor_shortlist_v1.csv` are **not** a pass. DualDiff-style P2Rank / AlphaFold / one-reference-ligand pockets are **not** a pass.

Log every receptor end in `tables/site_verification_log_v1.csv`. A pair is dockable only when **both** ends are `PASS`. Do not start Track B docking, do not freeze a new PDBQT, and do not run cognate redocking until that log row exists.

Declared sites live in `scripts/receptor_shortlist_v1.py` (`INTENDED_SITE`). Changing a declared site requires a written protocol amendment, not a silent swap.

Box construction (same as Methods): axis-aligned box on **that cognate’s heavy atoms**, expand 5 Å per axis, any edge < 20 Å raised to 20 Å.

## Stop rule

Open the structure in RCSB 3D / PyMOL / ChimeraX. If you cannot point to the cognate inside the declared pocket, the receptor **fails**, even if it is rank-1 by resolution.

## Per receptor (one PDB, one end)

Work top to bottom. Record the evidence in the log; do not tick from memory.

### A. Identity already passed (Layer 1 — re-confirm, do not skip)

These can be assisted by SIFTS / RCSB metadata. A human still signs the log.

1. **Accession.** Polymer entity UniProt = the accession the ChEMBL labels came from (pair plan columns `uniprot_A` / `uniprot_B`). Reject orthologs, paralogs, and “close isoform.”
2. **Species.** Source organism *Homo sapiens*. Mouse / rat / chimeric surrogates fail unless a written species-surrogate exception exists (none is authorised for Tier 1).
3. **SIFTS coverage.** Record coverage of the matched entity vs UniProt. Coverage ≈ 1 of the **wrong** protein is how 2WXF looked fine.
4. **Longest matching protein entity ≥ 80 aa.** Reject presented peptides / epitopes. For zymogens (thrombin, fXa), use the **catalytic (heavy) chain**, not the light chain that shares the same accession.
5. **Mutation count.** Record engineered substitutions; flag catalytic-site or gatekeeper mutants.
6. **Complex partners.** List every other polymer (p85, mLST8, FKBP12, HLA/β2m, antibodies, nanobodies). Partners are allowed only if they do not redefine the pocket.
7. **Not a chimera / graft.** Reject human/mouse kinase chimeras and “PI3K” entries that are another isoform (worked fail: 3T8M = human PIK3CG chimera used as if it were PIK3CB).

### B. Cognate ligand is a real inhibitor at the declared site

8. **Name + CCD.** Write the inhibitor name and three-letter (or CCD) ID. Example gold: 4L23 / X6K (PI-103); 4JT6 / X6K; 4EY7 / E20; 4BDS / THA.
9. **Drug-like small molecule.** MW 150–750, non-polymer, not a metal. Matches G4 spirit.
10. **Not a dummy occupant.** Fail if the “cognate” is solvent, detergent, buffer, crystallization additive, a nucleotide cofactor used only as a placeholder, a peptide antigen, a molecular glue at an interface, rapamycin (or analogue) on FKBP–FRB, or an HLA-presented oligopeptide.
11. **Pharmacology matches the labels.** The cognate should be an orthosteric (or declared) inhibitor of **this** target’s assay endpoint class (kinase ATP-site inhibitor, cholinesterase gorge inhibitor, S1 protease inhibitor, PPAR LBD agonist/antagonist, cathepsin active-site inhibitor). A fragment in a crystal soaking site that is not the labelled pocket fails.

### C. Site — open the 3D (this is the step metadata cannot do)

12. **Declared `intended_site`.** Copy the string from `INTENDED_SITE` into the log **before** looking at the structure, then confirm or fail.
13. **Cognate sits in that pocket.** Visually: ligand heavy atoms in the declared cavity, not an adjacent allosteric site, not a crystal-contact niche, not a peptide-binding groove.
14. **Orthosteric vs allosteric vs interface.** State which of the three it is. For a dual **isoform / paralog** pair, both ends must be the **same class of site** (ATP vs ATP, S1 vs S1, LBD vs LBD). ATP vs FRB on mTOR is an automatic fail for PIK3CA/mTOR.
15. **Construct / domain.** Record what is actually in the asymmetric unit:
    - kinases: JH1 catalytic kinase vs JH2 pseudokinase (JAK/TYK2: declare JH1; JH2 is a different experiment);
    - PI3K: p110 catalytic subunit vs p85-only or peptide;
    - mTOR: kinase domain ± mLST8 vs FRB-only / FKBP–rapamycin;
    - PPAR: LBD vs DNA-binding domain;
    - thrombin / fXa: catalytic domain heavy chain vs Gla / EGF / light chain;
    - cholinesterase: catalytic subunit gorge, not a peptide-only construct.
16. **Coverage vs UniProt length.** Domain-only constructs are allowed if the domain **is** the declared site. A 9-mer with accession match is not a domain construct.
17. **Blacklist traps already seen in this project (fail immediately):**
    - **8PPZ** and any ≤ 2.5 Å mTOR entry that is FKBP–rapamycin / FRB (cognate 0AN etc.) — not the ATP site; keep **4JT6**.
    - **9CMK** and other PIK3CA molecular-glue / SMARCA degrader complexes already blacklisted in `FROZEN_PUBLIC_PAIRS.yaml`.
    - HLA class I + β2m entries where “PIK3CA” is a **presented 9-mer**.
    - **2WXF** (murine p110δ / Pik3cd O35904) and mouse Pik3cb 2Y3A / 4BFR.
18. **Resolution is not a vote.** Worse resolution at the right site beats a 1.8 Å structure at the wrong site (4JT6 3.60 Å ATP vs 8PPZ 1.85 Å FRB).

### D. Docking box (only after A–C pass)

19. **Box from this cognate.** Heavy-atom AABB + 5 Å/axis, min edge 20 Å. Do not copy a box from a different PDB, a P2Rank centroid, or an AlphaFold pocket.
20. **Same protocol as frozen Methods.** Water and cognate removed; Meeko PDBQT; no LigPrep. Record box center and sizes in the log.
21. **Both-end consistency.** For the pair, write one sentence: why these two pockets are a fair dual-target comparison (same site class, comparable construct).

### E. Sign-off (blocks Layer 3 and production docking)

22. **Pass / fail** for this PDB + CCD.
23. **PDB ID + cognate CCD** written explicitly (the identifiers that will go into the receptor freeze).
24. **Operator + date.** Initials; ISO date.
25. **If fail:** either pick the next shortlist row and restart A–E, or **drop the pair** from Track B. Do not “dock anyway and check RMSD.” Layer 3 passed 2WXF at 0.405 Å.

## Per pair (after both ends PASS)

26. **Shared-ligand pose-gold (optional, record if used).** Same CCD on both ends (PI-103 / X6K on 4L23 and 4JT6; TAK-285 / 03P on EGFR/HER2) is a plus, not a requirement.
27. **No silent receptor reuse across the wrong pair.** JAK1 6N7A and PPARA 6LXA may be reused only on the pairs listed in the Tier-1 roster.
28. **Independent systems.** JAK1/TYK2 and JAK1/JAK2 share JAK1: they are two pairs, **one** JAK system. Same for PPARA across the two PPAR pairs. Statistics later count 6 systems, not 8 replicates.
29. **Written freeze.** PDB, CCD, box, and log SHA go into the receptor manifest **before** the first production ligand is docked.

## Worked failures (do not rediscover these in Track B)

| PDB | Looked like | Actually | Which items catch it |
|-----|-------------|----------|----------------------|
| 2WXF | “PIK3CB” holo, cognate RMSD 0.405 Å | Mouse p110δ (O35904), wrong gene and species | A1, A2, A7 |
| 8PPZ | Highest-resolution mTOR | FKBP–rapamycin/FRB, not kinase ATP | B10, C13, C14, C17 |
| 9CMK | Rank-1 PIK3CA | Molecular glue; blacklisted | B10, C17 |
| HLA–PIK3CA 9-mers | Accession match, high resolution | Presented peptide, not p110α | A4, C16, C17 |
| Thrombin light chain | Same UniProt, short entity | Not the S1 catalytic domain | A4, C15 |

## Already-frozen receptors (retrospective log, still required)

Fill the log even though docking already happened. Identity audit (`RECEPTOR_IDENTITY_AUDIT_V1.md`) is Layer 1 only.

| Pair | End | PDB | CCD | Declared site | Layer-1 identity | Layer-2 site (this checklist) |
|------|-----|-----|-----|---------------|------------------|-------------------------------|
| PIK3CA/mTOR | PIK3CA | 4L23 | X6K | class I PI3K ATP (p110α) | PASS (P42336, human) | **must still be signed** in the log |
| PIK3CA/mTOR | mTOR | 4JT6 | X6K | kinase ATP, not FRB | PASS (P42345, human) | **must still be signed**; 4JT6 is the ATP-site exception to auto-rank |
| AChE/BChE | AChE | 4EY7 | E20 | catalytic gorge | PASS (P22303) | **must still be signed** |
| AChE/BChE | BChE | 4BDS | THA | catalytic gorge | PASS (P06276) | **must still be signed** |

Do **not** replace 4L23 with auto-rank 9CMK, or 4JT6 with auto-rank 8PPZ, during this log. Auto-rank 4M0E / 6ZWI for AChE/BChE are not the frozen receptors.

## New Tier-1 receptors (no Vina until both ends PASS)

Shortlist proposals only — **not** verified:

| Pair | End | Proposed PDB / CCD | Declared site | Extra visual checks |
|------|-----|--------------------|---------------|---------------------|
| F2/F10 | F2 | 5AFY / WCE | thrombin S1, catalytic domain | Heavy chain ≥ 80 aa; S1 not exosite; not the 28 aa light chain |
| F2/F10 | F10 | 2JKH / BI7 | fXa S1/S4, catalytic domain | Same |
| JAK1/TYK2 | JAK1 | 6N7A / KEV | JH1 ATP | Confirm JH1 not JH2; Type I/II ATP geometry is fine, JH2 is not |
| JAK1/TYK2 | TYK2 | 3LXP / IZA | JH1 ATP (declare if JH2) | **Explicit JH1 vs JH2 call** |
| JAK1/JAK2 | JAK2 | 8BXH / C87 | JH1 ATP | Same JH1 rule; JAK1 reuses 6N7A |
| PPARG/PPARA | PPARG | 9F7W / 2OH | LBD | LBD agonist/antagonist pocket, not a coactivator-peptide-only interface |
| PPARG/PPARA | PPARA | 6LXA / EPA | LBD | Fatty-acid / fibrate LBD; EPA is a legitimate LBD occupant if it sits in LBD |
| PPARA/PPARD | PPARD | 5U3Q / 7UJ | LBD | Same; PPARA reuses 6LXA |
| CTSK/CTSS | CTSK | 4X6H / 3XT | S2 papain-fold active site | Catalytic Cys/His pair facing the ligand; not an occluding-loop-only crystal |
| CTSK/CTSS | CTSS | 9GJ2 / KH0 | S2 papain-fold active site | Same |

## What this checklist is not

- Not a substitute for G1–G5 pair gates.
- Not authorisation to dock CREBBP/BRD4, GPCRs, SLC6, or PIK3CA/PIK3CB.
- Not a DualDiff pocket protocol (no P2Rank, no AF-only receptors, no synergy-derived pairs).
- Not Layer 3. Cognate RMSD runs **after** PASS, and a pass RMSD never overrides a site fail.
