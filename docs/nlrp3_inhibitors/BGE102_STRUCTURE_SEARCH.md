# BGE-102 structure search (updated)

## Verdict

**As of this search, no free public source labels a specific chemical structure as “BGE-102.”**  
PubChem (CID/SID), Synapse free pages, ClinicalTrials.gov, and BioAge press releases still give only the code name.

What **is** newly public is BioAge’s issued composition-of-matter patent **US 12,509,459 B2** (granted **2025-12-30**), which BioAge highlighted in its **2026-01-12** Phase 1 update (“Patent issued covering additional composition of matter and novel NLRP3 binding site”). That patent discloses full example structures for the 5-azaindazole / biphenyl clinical chemotype, including cryo-EM of **Compound 007** and obesity (DIO) pharmacology for **Compound 040 / 096 / 211**.

Equating any example to clinical BGE-102 remains **inference**, not an official company identity statement.

## Sources checked

| Source | Structure labeled BGE-102? |
| --- | --- |
| PubChem name/substance | No (404) |
| Synapse / 智慧芽 free drug page | Synonym only: BGE-102 |
| ClinicalTrials NCT07656727 | No structure / CAS |
| BioAge IR press releases | Mechanism / PK / PD only |
| Corporate / R&D Day PDFs (`ir.bioagelabs.com`) | Download timed out in this environment (~20–23 MB) |
| PDB “BGE” / 9MIG / 9MGY | Unrelated ligands (not BioAge clinical) |
| **US 12,509,459 B2** | Full example structures; **no “BGE-102” string** |

## US 12,509,459 B2 — key disclosed structures

### Compound 007 (cryo-EM exemplar)

Used for NLRP3 epitope / cryo-EM (Figs. 1–6). From Table 1 + Example 3 text:

- East: **1H-pyrazolo[4,3-c]pyridin-7-yl** (5-azaindazole)
- Chiral benzylic carbon: **(R)** (cryo-EM); 2-hydroxyethyl side chain
- Central: **N-methyl amide**
- West: **2-ethoxy-[1,1′-biphenyl]-4-carboxamide** with **3′-fluoro** on the distal ring
- LCMS: **m/z = 449.1 [M+H]+**; THP-1 IC50 bin **A** (&lt;100 nM)

Approximate name:

`(R)-N-(2-hydroxy-1-(1H-pyrazolo[4,3-c]pyridin-7-yl)ethyl)-2-ethoxy-N-methyl-3'-fluoro-[1,1'-biphenyl]-4-carboxamide`

Compound **009** is the opposite enantiomer (same connectivity; weaker or matched depending on assay context in Table 4: both listed A in the extracted table).

### Compound 040 (DIO / weight-loss exemplar)

Prepared from Compound 007 alcohol → aldehyde → reductive amination with dimethylamine (Example synthesis). Table 1 / pharmacology:

- Same biphenyl–ethoxy–3′-F west and 5-azaindazole east
- Chiral side chain: **2-(dimethylamino)ethyl** instead of 2-hydroxyethyl
- LCMS: **m/z = 476.3 [M+H]+**
- Used in DIO mouse studies vs semaglutide / NT-0796 / VTX3232 (Figs. 8–13, 17)

This is the patent example whose **preclinical narrative** (oral NLRP3 inhibitor, DIO weight loss ± GLP-1) most closely tracks public BGE-102 messaging — still **not labeled** BGE-102.

### Related in vivo examples

- **Compound 096**, **Compound 211**: additional DIO agents in the same patent
- Binding-site residues claimed for the series include Y143, R147, F257, Y258, H260, E263, V264, L272, L275, I276, C279, F299, G328, L331, L332, L335, C514 (SEQ ID NO:1 numbering)

## What this does *not* settle

1. Official identity **BGE-102 ≡ Compound 040** (or 007 / 096 / 211 / other) — not stated.
2. Whether clinical BGE-102 is a salt, deuterated analog (e.g. Compound 211 uses bis(methyl-d₃)amine), or a later example outside the highlighted DIO set.
3. IR deck slides may show a structure drawing without naming CAS/SMILES; those large PDFs could not be retrieved here.

## Practical takeaway for SAR work

- Treat **US12509459** (and related WO/US filings in the same family) as the public structure corpus for BioAge’s clinical chemotype.
- Do **not** write “BGE-102 structure = …” in tables unless a primary source equates the code to a specific example.
- Closest working hypothesis for the clinical DC chemotype: **5-azaindazole east + N-methyl amide + chiral C with basic amine (often dimethylaminoethyl) + 2-ethoxy-3′-fluoro-biphenyl west**, i.e. Compound **040**-like — pending official confirmation.

## Primary links

- Patent: https://patents.google.com/patent/US12509459B2/en  
- USPTO PDF: https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12509459  
- Press note on issued patent: BioAge 2026-01-12 Phase 1 update / 8-K Exhibit 99.1  
