# Layer-2 metadata evidence pack

Collected 2026-09-04 from RCSB entry records + PDBe UniProt mappings + CCD definitions. Script: `scripts/site_verification_evidence_v1.py`. Table: `tables/site_verification_evidence_v1.csv`.

This pack filled Layer 1 and the non-visual half of Layer 2. **Human + literature PASS is now recorded** in `LAYER2_LITERATURE_SIGN_OFF_V1.md` and `tables/site_verification_log_v1.csv`. Numbering columns: `tables/receptor_span_registry_v1.csv`. Covalent ligand rule: `COVALENT_LIGAND_PREP_V1.md`.

## Identity (all proposed receptors)

Every chosen PDB maps to the ChEMBL-label UniProt and *Homo sapiens*. Longest matching span ≥ 80 aa. No 2WXF-class species swap.

| PDB | Acc | Human | UniProt span | Domain window | Partners |
|-----|-----|-------|--------------|---------------|----------|
| 4L23 | P42336 | yes | 1–1068 | p110 | p85α P27986 (niSH2; allowed) |
| 4JT6 | P42345 | yes | 1376–2549 | FRB **and** kinase ATP | mLST8 Q9BVC4 (allowed) |
| 4EY7 | P22303 | yes | 33–574 | catalytic | none |
| 4BDS | P06276 | yes | 29–557 | catalytic | none |
| 4UDW | P00734 | yes | H **364–621** (UniProt mature heavy **364–622**) | **heavy catalytic** + light | hirudin P09945 (standard soak) |
| 2JKH | P00742 | yes | A **235–475** (UniProt activated heavy **235–488**) | **heavy catalytic** + light | none other |
| 6N7A | P23458 | yes | construct ~854–1154 | **JH1** (UniProt 875–1153); JH2 overlap 2 aa | none |
| 3LXP | P29597 | yes | construct 888–1182; resolved ~888–1178 | **JH1 only** (UniProt 897–1176) | none |
| 8BXH | O60674 | yes | 840–1132 | **JH1** | none |
| 9V8H | P37231 | yes | 231–505 | **LBD** | synthetic PG08-NL peptide (**ternary**, not binary LBD) |
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
| 4X6H | I37 (dock) / 3XT (adduct) | glycinonitrile / 2-iminoethyl | 318 / 320 | yes | CTSK active-site; **reversible covalent nitrile — dock I37, not 3XT** |
| 9GJ2 | KH0 (bound) | ketoamide 13b | 596 | yes | CTSS active-site; **reversible covalent α-ketoamide — reconstruct 13b, not KH0**; journal still To Be Published |

## Written receptor amendments (auto-rank rejected)

| Target | Auto-rank | Why rejected | Chosen instead |
|--------|-----------|--------------|----------------|
| F2 | 5AFY / WCE | 3-chlorobenzamide, MW 156, fragment soak | **4UDW / N6L** (shortlist rank 3) |
| PPARG | 9F7W / 2OH | bisphenol A + PGC1α peptide | **9V8H / BRL** rosiglitazone (rank 2) |
| mTOR | 8PPZ / 0AN | FKBP–rapamycin/FRB | keep frozen **4JT6** |
| PIK3CA | 9CMK | molecular glue | keep frozen **4L23** |

## JH1 vs JH2 (the call the shortlist could not make)

- **3LXP is JH1.** Construct 888–1182; resolved ≈888–1178; UniProt catalytic JH1 897–1176. Not the JH2 structures (3NYX / 4GVJ class).
- **6N7A is JH1.** Construct ~854–1154; UniProt catalytic JH1 875–1153; title “JAK1 kinase domain.”
- **8BXH is JH1.** Title “JAK2 JH1 in complex with momelotinib”; PDB 840–1132; UniProt JH1 849–1124.

## What this pack no longer blocks

Cavity occupancy and covalent chemistry were signed from RCSB + UniProt + primary papers on 2026-09-04. Remaining gate is Layer-3 cognate RMSD (local Vina), with I37 as the 4X6H gold and **no** KH0 RMSD gold on 9GJ2.
