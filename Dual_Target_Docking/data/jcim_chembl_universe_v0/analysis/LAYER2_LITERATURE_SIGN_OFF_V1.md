# Layer-2 sign-off (2026-09-04) — what can be used

Human check against RCSB PDB, UniProt, and the primary structure papers. Independent RCSB REST confirmation of the two covalent CCD names and of 9GJ2’s citation status is in this file. This **is** the Layer-2 PASS for the ten new receptors. It does **not** authorise production Vina (Layer 3 cognate RMSD is still required).

Log: `tables/site_verification_log_v1.csv` (`PASS` on all ten new ends).  
Numbering registry: `tables/receptor_span_registry_v1.csv`.  
Covalent ligand protocol: `COVALENT_LIGAND_PREP_V1.md`.

## Verdict in one sentence

**Keep all ten PDBs as receptors.** Pairs 1–8 can go forward with numbering and Methods notes only. Pairs 9–10 keep the same crystals; the **free-ligand chemistry** used for Vina must be the pre-reaction form, not the crystallographic adduct.

## Usable as written (receptor + extracted cognate)

| # | Protein / PDB / CCD | Use? | Locked wording |
|---|---------------------|------|----------------|
| 1 | F2 P00734 / **4UDW / N6L** | **Yes** | N6L in thrombin S1 on heavy chain H. PDB mapping **364–621**; UniProt mature heavy chain **364–622**. Rühmann *J Med Chem* 2015, 10.1021/acs.jmedchem.5b00812. |
| 2 | F10 P00742 / **2JKH / BI7** | **Yes** | BI7 in Factor Xa active site spanning S1/S4. PDB heavy construct **235–475** (not UniProt activated heavy 235–488). Salonen *Angew Chem Int Ed* 2009, 10.1002/anie.200804695. |
| 3 | JAK1 P23458 / **6N7A / KEV** | **Yes** | JAK1 **JH1** ATP pocket; KEV = compound 39. Construct ~**854–1154**; UniProt catalytic JH1 **875–1153**. Do not write 854–1154 = JH1 domain. Zak *Bioorg Med Chem Lett* 2019, 10.1016/j.bmcl.2019.04.008. |
| 4 | TYK2 P29597 / **3LXP / IZA** | **Yes** | **JH1, not JH2.** IZA = CMP-6, ATP-site type-I. Construct **888–1182**; resolved ≈ **888–1178**; UniProt catalytic JH1 **897–1176**. Chrencik *J Mol Biol* 2010, 10.1016/j.jmb.2010.05.020. |
| 5 | JAK2 O60674 / **8BXH / C87** | **Yes** | JAK2 JH1, C87 = momelotinib, ATP site. PDB **840–1132**; UniProt JH1 **849–1124**. Miao *J Med Chem* 2024, 10.1021/acs.jmedchem.4c00197. |
| 6 | PPARG P37231 / **9V8H / BRL** | **Yes, with a Methods sentence** | BRL = rosiglitazone in PPARγ LBD **231–505**. Structure is **ternary**: LBD + rosiglitazone + **PG08-NL peptide**. Usable as the rosiglitazone LBD template; not a binary LBD. Sigal *JACS* 2025, 10.1021/jacs.5c13803. |
| 7 | PPARA Q07869 / **6LXA / EPA** | **Yes** | EPA in PPARα LBD, residues **200–468**. Kamata *iScience* 2020, 10.1016/j.isci.2020.101727. |
| 8 | PPARD Q03181 / **5U3Q / 7UJ** | **Yes** | 7UJ = specific agonist 1 in PPARδ LBD **170–441**; not PEG. Wu *PNAS* 2017, 10.1073/pnas.1621513114. |

## Usable as receptors only after ligand-prep change

| # | Protein / PDB / CCD | Receptor? | Ligand for Vina |
|---|---------------------|-----------|-----------------|
| 9 | CTSK P43235 / **4X6H / 3XT** | **Keep 4X6H** | **Do not Meeko-cut 3XT.** 3XT is the 2-iminoethyl reacted form (C16H21FN4O2, MW 320). The pre-reaction nitrile is already in the same entry as **I37** (cyanomethyl, C16H19FN4O2, MW 318). Dock **I37**. Boríšek *J Med Chem* 2015, 10.1021/acs.jmedchem.5b00746. |
| 10 | CTSS P25774 / **9GJ2 / KH0** | **Keep 9GJ2** | **Do not Meeko-cut KH0.** KH0 is ketoamide 13b in the bound (thiohemiketal-like) state. Reconstruct **pre-reaction α-ketoamide 13b**. Primary journal paper is still **To Be Published** (RCSB 2026-09-04: Falke/Turk et al.; no DOI). Cite the PDB deposition, not a journal article. |

Do **not** swap 4X6H or 9GJ2 for a different PDB. The identity and pocket occupancy are correct; only the free-ligand representation was wrong for a uniform noncovalent Vina protocol. The 2026-09-04 freeze (`RECEPTOR_FREEZE_V1.md`) records that 7AWC (binary PPARG) and 4P6G/2HHN (noncovalent CTSS) were considered and declined.

## Three numbering columns (do not collapse)

The 364–621 / 235–475 / 888–1182 disputes were all the same error: mixing PDB construct, resolved coordinates, and UniProt domain. Locked table: `tables/receptor_span_registry_v1.csv`.

## Independent RCSB REST checks on this date

- **4X6H** CCD **I37** name contains cyanomethyl; formula C16 H19 F N4 O2. CCD **3XT** name contains (2Z)-2-iminoethyl; formula C16 H21 F N4 O2. Primary citation title matches Boríšek et al. 2015.
- **9GJ2** `rcsb_journal_abbrev` = **To be published**; authors Falke, Karnicar, Usenik, Lindic, Sekirnik, Reinke, Guenther, Turk, Meents; no DOI. CCD **KH0** formula C31 H41 N5 O7, MW 595.687.
- **9V8H** title: *PPARgamma ligand-binding domain in complex with PG08-NL and rosiglitazone*.

## Independent geometric re-check (same day)

The table above is a text/literature check. `GEOMETRIC_POCKET_VERIFICATION_V1.md` re-does it from the raw mmCIF coordinates: real Cys25 bond lengths (3XT 1.83 Å covalent, I37 2.93 Å non-bonded, KH0 1.78 Å covalent) and real ligand-residue contact lists for all 14 receptors. No pocket assignment failed; no accession or organism failed. Use that file when a reviewer wants numbers instead of a title.

## What this PASS does not authorise

- Layer 3 cognate best-of-9 RMSD on the eight ordinary new receptors (needs local Vina).
- Production docking of the five ordinary pairs in `DOCKING_PLAN_V1.md`.
- Covalent / ordinary Vina of CTSK/CTSS. Re-docking PIK3CA/PIK3CB. Expanding K, editing Table 2, or docking CREBBP/BRD4, GPCR, or SLC6.
- Citing a journal paper for 9GJ2.
- Treating 9V8H as a binary PPARγ–rosiglitazone complex when comparing PPAR subtype LBD conformations.
