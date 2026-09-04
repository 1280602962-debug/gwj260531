# What you must check by eye (Layer 2 visual)

Metadata for every proposed receptor is already filled (`SITE_VERIFICATION_EVIDENCE_V1.md`, `tables/site_verification_evidence_v1.csv`). **That is not a PASS.** Open each structure below in RCSB 3D (or PyMOL). If the cognate is in the declared pocket, write `PASS` and your initials in `tables/site_verification_log_v1.csv`. If not, write `FAIL` and stop that pair.

Do **not** start Vina / cognate RMSD until every **new** row is `PASS`. Frozen pairs are already docked; a 5-minute confirm is enough.

## The three leftover items (plain language)

| Leftover | What it is | This turn |
|----------|------------|-----------|
| Write conclusions into the JCIM main manuscript | Track A: rewrite Methods/Results around the census ladder, 2WXF, Tier-1 | **Not done — you asked to leave the main text alone** |
| Sign the site-verification log | Layer 2: a person looks at each holo and confirms the cognate sits in the declared pocket | Metadata packed; **your eyes still required on the 10 new PDBs** |
| Optional: dock 6 new pairs | Track B: local Vina/GNINA after Layer 3 cognate RMSD | **Not started.** Cloud has no Vina. This is the next compute step **after** you sign the log |

So: next scientific work after you sign is **not** immediately production docking. Order is: **you sign 3D → local cognate best-of-9 RMSD (Layer 3) → then dock the 6 new pairs.** Frozen PIK3CA/mTOR and AChE/BChE stay as they are.

## Unique proteins / PDBs you must open

Ten new receptor files (JAK1 and PPARA are reused, so 10 freezes cover 6 pairs). Click the 3D link; confirm the one sentence in the last column.

| # | Protein (UniProt) | Pair(s) | PDB / CCD | Declared site | Open 3D | You confirm |
|---|-------------------|---------|-----------|---------------|---------|-------------|
| 1 | Thrombin F2 **P00734** | F2/F10 | **4UDW / N6L** | S1 catalytic heavy chain | https://www.rcsb.org/3d-view/4UDW | N6L in **S1**, chain H (heavy, UniProt 364–621), not the 28-residue light chain. Hirudin peptide may stay. |
| 2 | Factor Xa F10 **P00742** | F2/F10 | **2JKH / BI7** | S1/S4 catalytic | https://www.rcsb.org/3d-view/2JKH | BI7 spans **S1/S4** on the heavy chain (235–475), not light-chain-only. |
| 3 | JAK1 **P23458** | JAK1/TYK2 and JAK1/JAK2 | **6N7A / KEV** | JH1 ATP | https://www.rcsb.org/3d-view/6N7A | KEV in the **JH1 ATP cleft**. Construct is UniProt 854–1154 (JH1). Two residues touch the JH2 window — ignore that; this is not a JH2 structure. |
| 4 | TYK2 **P29597** | JAK1/TYK2 | **3LXP / IZA** | JH1 ATP | https://www.rcsb.org/3d-view/3LXP | **JH1, not JH2.** Mapping is 888–1182 = catalytic kinase. IZA (CMP-6) in the ATP site. Title also mentions JAK3; this entry’s polymer is TYK2. |
| 5 | JAK2 **O60674** | JAK1/JAK2 | **8BXH / C87** | JH1 ATP | https://www.rcsb.org/3d-view/8BXH | Title already says JH1 + momelotinib. Confirm C87 in the ATP cleft. |
| 6 | PPARG **P37231** | PPARG/PPARA | **9V8H / BRL** | LBD | https://www.rcsb.org/3d-view/9V8H | **Rosiglitazone (BRL) in the LBD** (231–505). Peptide PG08-NL may be present; ignore it if BRL is in the LBD cavity. |
| 7 | PPARA **Q07869** | PPARG/PPARA and PPARA/PPARD | **6LXA / EPA** | LBD | https://www.rcsb.org/3d-view/6LXA | EPA (fatty acid) in the **LBD** (200–468), not a crystal-contact groove. |
| 8 | PPARD **Q03181** | PPARA/PPARD | **5U3Q / 7UJ** | LBD | https://www.rcsb.org/3d-view/5U3Q | Specific agonist 7UJ in the **LBD** (170–441). Extra PEG/ions are not the cognate. |
| 9 | Cathepsin K **P43235** | CTSK/CTSS | **4X6H / 3XT** | papain S2 / active site | https://www.rcsb.org/3d-view/4X6H | 3XT facing **Cys25/His** (papain 115–329). Note: glycinonitrile series can be **covalent**; if the warhead is bonded to Cys, record `covalent=yes` — Vina will still treat it noncovalently. |
| 10 | Cathepsin S **P25774** | CTSK/CTSS | **9GJ2 / KH0** | papain S2 / active site | https://www.rcsb.org/3d-view/9GJ2 | Ketoamide 13b (KH0) in the **active site** (115–331). Same covalent flag as CTSK. |

## Optional 5-minute confirms (already docked)

| Protein | PDB / CCD | 3D | Already known |
|---------|-----------|-----|----------------|
| PIK3CA P42336 | 4L23 / X6K (PI-103) | https://www.rcsb.org/3d-view/4L23 | Human p110α + p85 niSH2; PI-103 in class I ATP site. Do not switch to 9CMK. |
| mTOR P42345 | 4JT6 / X6K (PI-103) | https://www.rcsb.org/3d-view/4JT6 | ΔN-mTOR + mLST8 **includes FRB+kinase** (1376–2549); PI-103 is in the **kinase ATP** site, not FRB. Do not switch to 8PPZ. |
| AChE P22303 | 4EY7 / E20 (donepezil) | https://www.rcsb.org/3d-view/4EY7 | Catalytic gorge. Keep 4EY7, not 4M0E. |
| BChE P06276 | 4BDS / THA (tacrine) | https://www.rcsb.org/3d-view/4BDS | Catalytic gorge. Keep 4BDS, not 6ZWI. |

## Why two shortlist PDBs were swapped (written amendment)

Auto-rank is not a freeze.

- **F2:** 5AFY / WCE (3-chlorobenzamide, MW 156) is an S1 **fragment**. Replaced by **4UDW / N6L** (D-Phe-N-(2,5-dichlorobenzyl)-Pro, MW 420), first shortlist row with an inhibitor-sized cognate.
- **PPARG:** 9F7W / 2OH is **bisphenol A** plus a PGC1α peptide. Replaced by **9V8H / BRL (rosiglitazone)**, a canonical LBD agonist.

## How to sign

In `tables/site_verification_log_v1.csv`, for each new row set:

- `site_matches_declared` = 1 or 0 (after 3D)
- `operator` = your initials
- `date_iso` = YYYY-MM-DD
- `pass_fail` = `PASS` or `FAIL`
- `notes` = one sentence (e.g. “N6L in S1 on chain H”)

A pair is dockable only when **both** ends are `PASS`. JAK1/TYK2 and JAK1/JAK2 share the JAK1 row; PPARG/PPARA and PPARA/PPARD share the PPARA row.

## After you sign — is the next step docking?

**Almost, but not quite.** On a machine that has AutoDock Vina (this cloud environment does not):

1. Freeze the signed PDB + CCD (do not re-pick from the shortlist).
2. Build the box from **that** cognate (AABB + 5 Å/axis, min edge 20 Å) — same as Methods.
3. **Layer 3:** cognate best-of-9 heavy-atom RMSD. Fail the receptor if it cannot generate a near-crystal pose (do not treat this as identity proof; 2WXF passed RMSD).
4. **Then** production-dock the six new pairs (F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, PPARA/PPARD, CTSK/CTSS), ~1,300–1,500 Vina jobs at the current 110-ligand depth, plus five-seed and failure typology if you want the same audit depth as the frozen pairs.

Do **not** dock CREBBP/BRD4, GPCRs, SLC6, or PIK3CA/PIK3CB. Do not expand K or edit Table 2 without a separate written amendment.
