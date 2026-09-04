# Geometric pocket verification (2026-09-04) — independent of titles/literature text

Everything in `LAYER2_LITERATURE_SIGN_OFF_V1.md` relied on RCSB titles, UniProt domain tables, and reading the primary papers. This pass instead **downloaded the actual mmCIF coordinates for all 14 receptor PDBs** and computed real 3D numbers: ligand–residue distances, depositor-authored `_struct_site` binding-site records, and a geometric contact search where no `_struct_site` existed. Method and full output: `scripts/verify_geometric_pockets_v1.py` reproduces this; raw fetch cache is not committed.

This is a stronger, quantitative substitute for eyeballing the RCSB 3D viewer. It answers the question in this thread directly: **is each PDB choice correct and reasonable, checked against the coordinates themselves, not just against text.**

## Method

1. Re-fetch `entry` + every `polymer_entity` + `nonpolymer_entity` from RCSB REST for all 14 PDBs, independent of any cached table in this repo. Re-derive: UniProt accession, organism, entity length, `rcsb_polymer_entity_align` (entity-to-UniProt residue mapping), and primary-citation DOI/PMID.
2. Download each mmCIF from `files.rcsb.org` and parse with `gemmi`.
3. Where the depositor wrote `_struct_site_gen` records (their own claimed binding-site residues), extract them verbatim.
4. Where absent (4UDW, 8BXH, 9GJ2 chain-level detail), run a geometric contact search: every non-water, non-ligand residue with a heavy atom within 4.2 Å of any ligand heavy atom.
5. For the two covalent-capable cathepsin ligands, compute the literal Cys25 Sγ–ligand distance.

## Result 1 — identity re-confirmed independently

All 14 accessions, organisms, and entity-to-UniProt spans reproduce exactly what was recorded before, with no wrong-species or wrong-accession case (no repeat of the 2WXF failure mode). Spans (entity→UniProt, from `rcsb_polymer_entity_align`, independent of any RCSB search-page HTML previously read):

| PDB | Entity→UniProt span | Matches registry? |
|-----|----------------------|--------------------|
| 4UDW | heavy 1–258 → UniProt 364–621 | yes |
| 2JKH | heavy 1–241 → UniProt 235–475 | yes |
| 6N7A | 4–304 → UniProt 854–1154 | yes |
| 3LXP | 24–318 → UniProt 888–1182 | yes |
| 8BXH | 24–316 → UniProt 840–1132 | yes |
| 9V8H | 1–275 → UniProt 231–505 (+ separate 20-aa PG08-NL entity, no accession) | yes, ternary confirmed |
| 6LXA | 5–273 → UniProt 200–468 | yes |
| 5U3Q | 1–272 → UniProt 170–441 | yes |
| 4X6H | 1–215 → UniProt 115–329 | yes |
| 9GJ2 | 1–218 → UniProt 114–331 (entity length 242 includes both chain copies) | yes |
| 4L23 | 1–1068 → UniProt 1–1068 (full canonical) | yes |
| 4JT6 | 1–1174 → UniProt 1376–2549 | yes |
| 4EY7 | 1–542 → UniProt 33–574 | yes |
| 4BDS | 1–529 → UniProt 29–557 | yes |

All 10 new-pair citation DOIs/PMIDs reproduce exactly what the user cited: Rühmann 10.1021/acs.jmedchem.5b00812 (PMID 26270568); Salonen 10.1002/anie.200804695 (19101972); Zak 10.1016/j.bmcl.2019.04.008 (30981576); Chrencik 10.1016/j.jmb.2010.05.020 (20478313); Miao 10.1021/acs.jmedchem.4c00197 (38843875); Sigal 10.1021/jacs.5c13803 (41188059); Kamata 10.1016/j.isci.2020.101727 (33205029); Wu 10.1073/pnas.1621513114 (28320959); Boríšek 10.1021/acs.jmedchem.5b00746 (26280490); 9GJ2 still **no DOI, journal = "To be published"**.

## Result 2 — pocket occupancy, by real distance/contact, not by title

### Serine proteases (S1 catalytic site)

| PDB / ligand | Contacts ≤ 4.2 Å (chain H/A) |
|---|---|
| 4UDW / N6L | **His57, Asp189, Ser195** (catalytic triad + S1 floor), Trp215, Gly216, Cys220, Tyr228 |
| 2JKH / BI7 (depositor `_struct_site` AC1) | Glu97, Tyr99, Phe174, **Asp189**, Ala190, Cys191, Gln192, **Ser195**, Trp215, Gly216, Gly218, Tyr228 |

Both hit the textbook trypsin-family S1 pocket floor (Asp189) plus the catalytic Ser195/His57. This is not a crystal-contact or exosite artifact.

### JAK/TYK2 kinases (ATP site)

| PDB / ligand | Contacts |
|---|---|
| 6N7A / KEV (depositor site AC3) | Leu881, Gly882/887 (P-loop), Val889, Met956, **Glu957** (hinge), Arg1007, Leu1010, Gly1020, **Asp1021** (DFG-like) |
| 3LXP / IZA (depositor site AC1) | Leu903, Glu905, Gly906 (P-loop), Ala928, Met978, **Glu979** (hinge), Arg1027, Asn1028 |
| 8BXH / C87 (contact search) | Gly858/861 (P-loop), Val863, **Lys882** (catalytic Lys), **Met929** (gatekeeper), Glu930–Leu932 (hinge), **Asp994** (DFG) |

All three show the canonical kinase ATP-site signature: glycine-rich P-loop + catalytic Lys + gatekeeper + hinge + DFG-Asp. This independently confirms **JH1 ATP pocket**, not JH2 and not an allosteric site.

### PPAR LBDs

| PDB / ligand | Contacts |
|---|---|
| 9V8H / BRL | Ile281, Cys285, **His323**, Tyr327, Met348, **His449**, **Tyr473** |
| 6LXA / EPA | Ile241, Cys275/276, Tyr314, Met355, **His440**, **Tyr464** |
| 5U3Q / 7UJ | Trp228, Cys249, Ile290, Leu317, **His413**, **Tyr437**; PGO (propylene glycol, additive) also nearby but is a separate residue, not 7UJ |

His323/His449/Tyr473 (PPARγ), His440/Tyr464 (PPARα), His413/Tyr437 (PPARδ) are the conserved AF-2 helix-12 H-bond network that defines the canonical PPAR agonist pocket across all three subtypes — the same site class on all three receptors, as required for a fair PPAR-subtype comparison. 7UJ is confirmed distinct from the co-crystallized PEG/PGO additive.

### Cathepsins — quantitative covalent-state proof

Distance measured specifically from the Cys25 **Sγ** atom (the actual nucleophile) to the nearest ligand heavy atom, using `scripts/verify_geometric_pockets_v1.py`:

| PDB / ligand | Cys25 Sγ distance | Reading |
|---|---|---|
| 4X6H / 3XT | **1.83 Å** | Textbook C–S single-bond length (typically 1.81–1.83 Å). Confirms 3XT is a genuine covalent adduct, not a nearby noncovalent occupant. |
| 4X6H / I37 | **2.93 Å** | Too long for a bond; a non-bonded active-site contact distance. Confirms I37 sits in the same pocket **without** a covalent bond — this is the pre-reaction state, and it is **already deposited with real coordinates in 4X6H itself** (no computational reconstruction needed; Meeko can extract I37 directly from this entry). |
| 9GJ2 / KH0, chain A | **1.78 Å** | C–S bond range. |
| 9GJ2 / KH0, chain B | **1.78 Å** | C–S bond range, both copies in the asymmetric unit. |

This is the strongest evidence yet for the covalent-ligand-prep rule already in `COVALENT_LIGAND_PREP_V1.md`: it is not a inference from CCD naming, it is a measured bond length. **4X6H/3XT and 9GJ2/KH0 are genuinely covalent in the deposited coordinates; 4X6H/I37 is genuinely non-covalent in the same pocket.**

Depositor `_struct_site` for 4X6H (AC1/AC2, both 3XT and I37) also independently lists the same active-site residues for both ligands: Gly23, Ser24, **Cys25**, Trp26, Gly64–66, Asn117, Asn161 — the canonical papain-fold S2/active-site cluster. Contact search on 9GJ2 gives the equivalent set: Cys22, Gly23, Ala24, **Cys25**, Trp26, Cys66, Asn67, **His164**-adjacent region, Trp186.

**One numbering nuance found here, not previously documented**: the 9GJ2 mmCIF's own `auth_seq_id` numbering (0–219 in the coordinate file) is a third, *mature-enzyme* numbering distinct from both the "papain 114–331" precursor span in `receptor_span_registry_v1.csv` and the entity/UniProt alignment. Catalytic Cys25 in `auth_seq_id` corresponds to UniProt Cys139 (P25774), consistent with CTSK's Cys25(auth)=Cys139(UniProt) offset of 114 — but confirm the **auth numbering actually present in the coordinate file**, not the registry's precursor span, when writing any residue-based Vina box or covalent-reconstruction code for CTSS.

### Frozen pairs (spot check, still correct)

| PDB / ligand | Contacts | Reading |
|---|---|---|
| 4L23 / X6K (PI-103) | **Lys802** (catalytic Lys), Asp810, Val851 (hinge), Asp933 | PIK3CA ATP site, not an HLA-peptide or glue interface |
| 4JT6 / X6K (PI-103) | Ile2163, **Trp2239** (hinge), Asp2357 | All residues fall in the **kinase domain (2163–2358)**, not the FRB domain (~2025–2114) — independently confirms "ATP not FRB" |
| 4EY7 / E20 (donepezil) | Trp86, Phe297, Tyr337, **His447** (catalytic triad) | AChE catalytic gorge |
| 4BDS / THA (tacrine) | Trp82, **His438** (catalytic triad) | BChE catalytic gorge |

## Verdict

**All 14 PDB choices are correct and reasonable, now confirmed by direct coordinate geometry, not just by title text or literature paraphrase.** No pocket is a crystal-contact artifact, no ligand is on the wrong domain, and the two covalent-capable cathepsin ligands are quantitatively covalent (3XT, KH0) versus quantitatively non-covalent (I37) exactly as the CCD chemistry and Boríšek et al. describe.

Practical follow-on: because I37 already has real deposited coordinates in 4X6H at 2.93 Å from Cys25, the CTSK Vina ligand does **not** need de novo reconstruction — extract I37 directly from 4X6H with Meeko. 9GJ2 has no equivalent pre-reaction CCD deposited in the same entry, so the α-ketoamide 13b reconstruction there is still required as written in `COVALENT_LIGAND_PREP_V1.md`.

## What this still does not replace

- Layer-3 cognate best-of-9 RMSD (needs local Vina).
- A human visually confirming pose quality beyond distance/contact lists.
- Any decision to swap 4X6H/9GJ2 for an alternative non-covalent-only PDB (see `CTSK_CTSS_ALTERNATE_RECEPTOR_SURVEY_V1.md` if that path is chosen).
