# Layer-2 metadata evidence pack (not a visual PASS)

Collected 2026-09-04 from RCSB entry records + PDBe UniProt mappings + CCD definitions. Script: `scripts/site_verification_evidence_v1.py`. Table: `tables/site_verification_evidence_v1.csv`.

This pack fills Layer 1 and the non-visual half of Layer 2. **`pass_fail` stays `AWAITING_VISUAL` on new receptors until a human opens the 3D** (`HUMAN_VISUAL_SIGN_OFF_V1.md`).

## Identity (all proposed receptors)

Every chosen PDB maps to the ChEMBL-label UniProt and *Homo sapiens*. Longest matching span ≥ 80 aa. No 2WXF-class species swap.

| PDB | Acc | Human | UniProt span | Domain window | Partners |
|-----|-----|-------|--------------|---------------|----------|
| 4L23 | P42336 | yes | 1–1068 | p110 | p85α P27986 (niSH2; allowed) |
| 4JT6 | P42345 | yes | 1376–2549 | FRB **and** kinase ATP | mLST8 Q9BVC4 (allowed) |
| 4EY7 | P22303 | yes | 33–574 | catalytic | none |
| 4BDS | P06276 | yes | 29–557 | catalytic | none |
| 4UDW | P00734 | yes | H 364–621; L 333–360 | **heavy catalytic** + light | hirudin P09945 (standard soak) |
| 2JKH | P00742 | yes | A 235–475; L 126–180 | **heavy catalytic** + light | none other |
| 6N7A | P23458 | yes | 854–1154 | **JH1** (289 aa); JH2 overlap 2 aa | none |
| 3LXP | P29597 | yes | 888–1182 | **JH1 only** | none |
| 8BXH | O60674 | yes | 840–1132 | **JH1** | none |
| 9V8H | P37231 | yes | 231–505 | **LBD** | synthetic PG08-NL peptide |
| 6LXA | Q07869 | yes | 200–468 | **LBD** | none |
| 5U3Q | Q03181 | yes | 170–441 | **LBD** | none |
| 4X6H | P43235 | yes | 115–329 | papain | none |
| 9GJ2 | P25774 | yes | 114–331 | papain | none |

## Cognates

| PDB | CCD | Name | MW | In entry | Pharmacology note |
|-----|-----|------|---:|:--------:|-------------------|
| 4L23 | X6K | PI-103 | 348 | yes | class I PI3K/mTOR ATP inhibitor |
| 4JT6 | X6K | PI-103 | 348 | yes | same ligand, mTOR ATP (not FRB) |
| 4EY7 | E20 | donepezil | 379 | yes | gorge inhibitor |
| 4BDS | THA | tacrine | 198 | yes | gorge inhibitor |
| 4UDW | N6L | D-Phe-N-(2,5-dichlorobenzyl)-Pro | 420 | yes | thrombin S1 peptidomimetic |
| 2JKH | BI7 | fXa cation inhibitor | 464 | yes | S1/S4-class fXa inhibitor |
| 6N7A | KEV | JAK1 compound 39 | 397 | yes | JH1 ATP-site ligand |
| 3LXP | IZA | CMP-6 | 309 | yes | TYK2 JH1 ATP-site ligand |
| 8BXH | C87 | momelotinib | 414 | yes | JAK2 JH1 clinical inhibitor |
| 9V8H | BRL | rosiglitazone | 357 | yes | PPARG LBD agonist |
| 6LXA | EPA | eicosapentaenoic acid | 302 | yes | PPARA LBD fatty-acid agonist |
| 5U3Q | 7UJ | PPARD agonist 1 | 460 | yes | PPARD LBD agonist |
| 4X6H | 3XT | homocycloleucyl-glycinonitrile | 320 | yes | CTSK active-site; **possible covalent nitrile** |
| 9GJ2 | KH0 | ketoamide 13b | 596 | yes | CTSS active-site; **possible covalent ketoamide** |

## Written receptor amendments (auto-rank rejected)

| Target | Auto-rank | Why rejected | Chosen instead |
|--------|-----------|--------------|----------------|
| F2 | 5AFY / WCE | 3-chlorobenzamide, MW 156, fragment soak | **4UDW / N6L** (shortlist rank 3) |
| PPARG | 9F7W / 2OH | bisphenol A + PGC1α peptide | **9V8H / BRL** rosiglitazone (rank 2) |
| mTOR | 8PPZ / 0AN | FKBP–rapamycin/FRB | keep frozen **4JT6** |
| PIK3CA | 9CMK | molecular glue | keep frozen **4L23** |

## JH1 vs JH2 (the call the shortlist could not make)

- **3LXP is JH1.** PDBe maps P29597 888–1182 (JH1 window 888–1176). Not the JH2 structures (3NYX / 4GVJ class).
- **6N7A is JH1.** 854–1154; title “JAK1 kinase domain.”
- **8BXH is JH1.** Title “JAK2 JH1 in complex with momelotinib”; map 840–1132.

## What metadata still cannot do

Cognate **occupancy of the cavity** (vs a crystal-contact niche), covalent Cys–warhead geometry, and the mTOR ATP-vs-FRB visual on 4JT6. Those are the visual sign-off items.
