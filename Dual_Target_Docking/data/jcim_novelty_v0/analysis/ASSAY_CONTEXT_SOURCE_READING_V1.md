# Assay-context source reading v1 (7 uncertain ligands)

Date: 2026-08-27.
Scope: the 9 decisive ChEMBL documents plus the original Pictilisib paper
(Folkes 2008). Frozen DualFourClass labels are **not** changed in this
pass. Table 2 is **not** recomputed.

Public text used without local PDFs: US 9,181,263 B2 (Google Patents),
Hong et al. 2024 (PMC11284801), Folkes values via ChEMBL document
CHEMBL1140078 plus secondary citations of the 3 nM / 580 nM pair,
BindingDB assay CHEMBL5133798, and a Lookchem HTML extract of Sang 2020.
Paywalled PDFs still needed: Ma 2022, Yang 2023 SI, Elmenier 2019,
Bass 2021, Folkes 2008 body/SI, Sang 2020 Table 1, Cheng 2021 SI,
Shi 2024 SI. Drop them in `literature_sources/pdfs/` using the README
filenames.

Label rule used below: DualFourClass θ = 6.0 on the per-target maximum
pChEMBL. Class flips are recorded only as *conditional* outcomes.

## Frozen vs source-checked maxima

| ligand | pair | frozen class | frozen pA / pB | decisive docs | source-checked status |
|--------|------|--------------|----------------|---------------|------------------------|
| EH120_045 | EGFR/HER2 | dual | 9.30 / 9.12 | patent CHEMBL3886356; Ma 2022 CHEMBL5131445 | EGFR 0.5 nM located in patent Table 1. HER2 9.12 is a 0.76 vs 76 nM conflict. Class likely remains dual even if HER2 is 76 nM. |
| PM48_04 | PIK3CA/mTOR | dual | 8.85 / 9.35 | Hong 2024 CHEMBL5500428; Yang 2023 CHEMBL5620391 | 8.85 is a **cited** 1.4 nM, not a Hong measurement. mTOR 0.45 nM still needs ATP/n/error from Yang SI. Class remains dual if 8.85 is dropped. |
| **PM48_05** | PIK3CA/mTOR | dual | **10.00 / 8.52** | Bass 2021 review; Elmenier 2019 review; **must-read Folkes 2008** | Review maxima are not original Pictilisib values. Folkes 3 nM / 580 nM is already in ChEMBL as CHEMBL1140078 (8.52 / 6.24). At θ=6.0 this remains dual unless further mTOR rows are excluded. |
| AB_089 | AChE/BChE | neither | 4.78 / 5.12 | Sang 2020 CHEMBL4680246 | Paper used eeAChE + hAChE and eqBChE + hBChE. ChEMBL type A on AChE is a metadata error relative to Ellman potency. Class stays neither. |
| AB_091 | AChE/BChE | neither | 4.76 / 5.19 | same | same |
| AB_094 | AChE/BChE | neither | 4.71 / 5.29 | same | same |
| PM48_22 | PIK3CA/mTOR | A_only | 8.77 / 5.83 | Cheng 2021 CHEMBL4765307; Shi 2024 CHEMBL5579880 | 5.83 is cellular p70S6K ELISA (type A). Biochemical mTOR IC50 5.52 already exists. Class stays A_only. |

---

## 1. EH120_045 / Ibrutinib / CHEMBL1873475

### EGFR — US 9,181,263 B2 (CHEMBL3886356)

Located in patent Table 1 (and restated in Table 2 / Example 1c):

- Compound 1 EGFR IC50 = **0.5 nM** (table header is misspelled `EFGR`).
- Same row: BTK 0.5 nM, HER2 **9.4 nM**, HER4 0.1 nM.
- 0.5 nM → pChEMBL **9.30**, matching the frozen EGFR maximum.
- Assay: in vitro HotSpot kinase assay; **purified enzymes**; ³³P-ATP;
  “an appropriate substrate”; **1 μM ATP**; **1 hour** inhibitor
  incubation; 10-point curve 10 μM to 0.0005 μM; Prism IC50.
- EGFR mutation status is **unspecified**. The panel lists EGFR, not
  L858R/T790M/exon-19. Treat as unspecified WT commercial kinase, not
  a documented full-length cellular EGFR.
- Construct: not given beyond “purified enzymes”. Commercial HotSpot
  panels are typically kinase-domain / catalytic constructs, but the
  patent does not print residue ranges.

Compound 1 identity: the patent calls Compound 1 “our highly selective
BTK inhibitor” and uses the ibrutinib chemotype
`1-(3-(4-amino-3-(4-phenoxyphenyl)-1H-pyrazolo[3,4-d]pyrimidin-1-yl)piperidin-1-yl)prop-2-en-1-one`.
ChEMBL maps the document to CHEMBL1873475. Ibrutinib as a drug is the
**(R)** piperidine (PCI-32765). Table 1 does not print the R/S label on
Compound 1. This is a residual identity caveat, not a 100-fold potency
error.

Same-document HER2 IC50 9.4 nM → pChEMBL 8.03, matching audit row
CHEMBL3887534. That is **not** the frozen HER2 maximum.

### HER2 — Ma et al., 2022 BMCL (CHEMBL5131445)

Frozen HER2 maximum is Kd pChEMBL **9.12** (assay CHEMBL5133798), which
is **0.76 nM**. BindingDB for the same ChEMBL assay ID reports Kd
**76 nM** (pChEMBL 7.12) and attaches a **non-ibrutinib** ligand
(BDBM250082, amide-pyrazole SMILES, patent compound 27). So there is
both a 100-fold decimal conflict and a possible ligand mismatch.

Ibrutinib is expected to appear in Ma 2022 as a **kinase-selectivity
comparator**, not as a designed HER2 inhibitor. The local PDF must show:

1. the table cell for ibrutinib vs HER2 (0.76 nM vs 76 nM vs other);
2. whether the protein is full-length HER2 or a kinase-domain construct
   (KinomeScan-style assays are usually T7-tagged kinase domains).

Even if the 9.12 row is demoted to 7.12 or excluded, remaining HER2
rows include Kd 8.92 and IC50 8.82. EGFR remains 9.30. **Conditional
class at θ=6.0: still dual.**

---

## 2. PM48_04 / Dactolisib / CHEMBL1879463

### PIK3CA — Hong et al., 2024 JMC (PMC11284801, CHEMBL5500428)

The 1.4 nM PI3Kα value is **cited, not measured in this paper**:

> “Dactolisib (1) possesses functional features to potently inhibit
> DNA-PK kinase activity (biochemical IC50: 0.8 nM) and also potently
> inhibits other PIKK and PI3K enzymes (e.g., IC50: PI3Kα, 1.4 nM;
> mTOR, 4.3 nM).”

Table 7 reports measured IC50s for peposertib (2), AZD7648 (3), and
compounds 53/78. Dactolisib is **not** in that table. The 8.85
pChEMBL (1.4 nM) deposited on CHEMBL5500428 should not be treated as a
Hong experimental result.

Hong methods (Reaction Biology) for **their** compounds:

- PI3Kα: human **p110α/p85α**, PIP2 10 μM, ATP **10 μM**, ADP-Glo.
- mTOR: human mTOR, 4EBP1 1 μM, ATP **10 μM**, HotSpot ³³P.
- DNA-PK: human DNA-PK, 10 μM ATP, 20 μM peptide, 10 μg/mL DNA.
- n = 2 unless stated.

Primary literature for the 1.4 nM citation is likely Maira 2008
(*Mol. Cancer Ther.*, Hong ref 35). Optional PDF:
`Maira_2008_MCT_BEZ235.pdf`.

If the 8.85 row is excluded as a review-style citation, remaining
PIK3CA maxima are 8.40 and 8.22. **Class remains dual.**

### mTOR — Yang et al., 2023 EJMECH (CHEMBL5620391)

Public ChEMBL/GtoPdb description already matches the user brief:
human mTOR, N-terminal FLAG, residues 1362–end, 4E-BP1 Thr37/46,
LANCE Ultra, IC50 0.45 nM (pChEMBL 9.35). Local PDF/SI still required
for ATP concentration, replicate count, and error. This does not
control the class: even a much weaker retained mTOR value stays ≥6.

---

## 3. PM48_05 / Pictilisib / CHEMBL521851  (highest risk)

Frozen maxima **both come from reviews**:

| arm | frozen | ChEMBL document | what it is |
|-----|--------|-----------------|------------|
| PIK3CA | 10.00 (0.1 nM) | CHEMBL5214883 Bass 2021 HDAC-hybrid **review** | not the discovery paper |
| mTOR | 8.52 (3 nM) | CHEMBL4373732 Elmenier 2019 PI3K **review** | not the discovery paper |

### Original Folkes 2008 (already in ChEMBL as CHEMBL1140078)

Secondary citations (Raynaud/Workman 2009; Sarker 2015 CCR; APSB 2016)
and the ChEMBL Folkes document agree:

- PI3Kα IC50 **3 nM** → pChEMBL **8.52** (audit rows CHEMBL967868,
  CHEMBL967881, CHEMBL980850).
- mTOR **580 nM** class value; ChEMBL stores **Ki 6.24** on
  CHEMBL980857 (~575 nM).
- Sarker 2015: “193-fold less activity against mTOR compared to
  p110α”.
- DNA-PK ~1230 nM in the same selectivity profile.

Local Folkes PDF is still required to pin the exact table (compound 17
= GDC-0941), ATP, and whether 580 nM is IC50 or Ki.

### What happens if the two review maxima are dropped

Remaining PIK3CA (not from Bass): Kd 9.12, Kd 9.07, then Folkes IC50
8.52. Remaining mTOR (not from Elmenier 8.52): Kd 7.32, IC50 6.87, Kd
6.70, IC50 6.51, Folkes Ki 6.24.

| rule | pA | pB | class at θ=6.0 |
|------|----|----|----------------|
| keep frozen review maxima | 10.00 | 8.52 | dual (current, **not acceptable**) |
| drop both reviews, keep all other rows | 9.12 (Kd) | 7.32 (Kd) | dual |
| drop reviews; biochemical IC50/Ki only | 8.52 (Folkes) | 6.87 | dual |
| Folkes original pair only | 8.52 | 6.24 | **dual, barely** |
| Folkes + treat 580 nM as inactive because of 193-fold PI3K selectivity | 8.52 | <6 by policy | A_only (**only if a new rule is added**) |

**A silent dual → A_only flip is not justified at the frozen θ=6.0
max-pChEMBL rule.** The review 10.00 / 8.52 pair **must not** remain
the reported maxima. The realistic outcomes after source correction
are (i) still dual with pA ≈ 8.5–9.1 and pB ≈ 6.2–7.3, or (ii) A_only
only after an explicit extra incomparability/selectivity rule.

Elmenier depositing mTOR = 3 nM is the same number as Folkes **PI3Kα**,
which is the table-shift hypothesis. Bass depositing PI3Kα = 0.1 nM is
30-fold tighter than Folkes 3 nM and is treated as a review
mis-assignment until the Bass table is sighted.

---

## 4. AB_089 / AB_091 / AB_094 — Sang et al., 2020 EJMECH

One document, CHEMBL4680246. Frozen class **neither** for all three
(pA 4.71–4.78, pB 5.12–5.29).

From the Sang text (Lookchem extract; Table 1 still needs the PDF):

- Ellman assays on **eeAChE and hAChE** and **eqBChE and hBChE**.
- Compounds were first screened on electric eel / equine enzymes, then
  re-evaluated on human enzymes.
- Best compound is **4f** (eqBChE 0.92 μM; hBChE 0.97 μM), which is
  **not** these three neither ligands.
- ChEMBL assay type **A** on AChE (CHEMBL4686681) is inconsistent with
  a primary Ellman IC50 table. That is a ChEMBL metadata problem, not
  evidence that the measurement is ADME.

Local PDF Table 1 is required to map CHEMBL4792013 / 4787165 /
4761179 onto compound numbers and to confirm whether the deposited
~16 μM AChE / ~6–8 μM BChE values are **human** or **eel/horse**.
Species mixing would change construct notes, not the neither label:
all current pChEMBL values stay <6.

---

## 5. PM48_22 / Alpelisib / CHEMBL2396661

### mTOR — Cheng et al., 2021 JMC (CHEMBL4765307)

Audit row: assay type **A**, IC50 pChEMBL **5.83**. This is the
cellular MCF7 **TSC1-null p70S6K Thr389 ELISA**, not purified mTOR.
It is already marked `uncertain` / `incomparable_record=1`.

A separate biochemical-looking mTOR row already exists:
CHEMBL4706640 IC50 pChEMBL **5.52** (3000 nM). If Cheng 5.83 is
excluded, pB becomes 5.52. **Still A_only.**

### PIK3CA — Shi et al., 2024 JMC (CHEMBL5579880)

Alpelisib as PI3Kα control: IC50 **1.7 nM** → pChEMBL **8.77**,
ADP-Glo, 60 min ATP incubation. Local SI still needed for
p110α/p85α vs p110α-only and ATP concentration. Class is insensitive
to that detail.

---

## Actions taken / not taken

- Frozen `human_reviewed_class` left unchanged (7 uncertain, 0 exclude,
  0 class flips).
- Table 2 **not** recomputed.
- Next step after PDFs land: fill `protein_construct` /
  `wildtype_or_mutant` / `incomparable_record` on the nine decisive
  rows; recompute maxima; recompute Table 2 **only** if a frozen class
  actually changes.

## PDF still required (highest first)

1. Folkes 2008 body + SI
2. Elmenier 2019 (mTOR 3 nM table)
3. Bass 2021 (PI3Kα 0.1 nM table)
4. Ma 2022 body + SI (ibrutinib HER2 0.76 vs 76 nM)
5. Sang 2020 Table 1 (compound IDs + species)
6. Yang 2023 SI (ATP / n / SD)
7. Cheng 2021 SI (confirm cellular ELISA)
8. Shi 2024 SI (PI3Kα construct)
9. Optional: Maira 2008 for the dactolisib 1.4 nM citation
