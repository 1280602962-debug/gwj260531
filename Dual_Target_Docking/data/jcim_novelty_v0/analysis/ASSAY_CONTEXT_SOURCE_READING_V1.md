# Assay-context source reading v1 (7 uncertain ligands)

Date: 2026-08-27 (BindingDB article-dump + patent PDF + Hong PMC pass).
Scope: the 9 decisive ChEMBL documents plus the original Pictilisib
paper (Folkes 2008). Frozen DualFourClass labels are **not** changed
in this pass. Table 2 is **not** recomputed.

What was actually read:

| source | what landed | usable as |
|--------|-------------|-----------|
| US 9,181,263 B2 | local PDF `literature_sources/pdfs/US9181263B2.pdf` (93 pp) | primary patent Tables 1–2 and Example 4a |
| Hong et al. 2024 | PMC11284801 HTML (publisher PDF blocked) | primary paper intro + methods; dactolisib 1.4 nM is cited |
| Folkes, Elmenier, Bass, Ma, Sang, Yang, Cheng, Shi | BindingDB **article pages**, not publisher PDFs | deposited numbers + assay blurbs for the named document |
| Workman 2011 (*Chem. Biol.*) PMC3242038 | OA review quoting Folkes | independent restatement of the Folkes isoform / mTOR profile |

DOI downloads named `folkes.pdf` / `elmenier.pdf` / etc. resolved to
BindingDB HTML or a reCAPTCHA page, not ACS/Elsevier PDFs. Those dumps
are stored under `literature_sources/bindingdb_article_dumps/`. Publisher
page locators still need the real PDFs listed in
`literature_sources/README.md`.

Label rule used below: DualFourClass θ = 6.0 on the per-target maximum
pChEMBL. Class flips are recorded only as *conditional* outcomes.

## Frozen vs source-checked maxima

| ligand | pair | frozen class | frozen pA / pB | decisive docs | source-checked status |
|--------|------|--------------|----------------|---------------|------------------------|
| EH120_045 | EGFR/HER2 | dual | 9.30 / 9.12 | patent CHEMBL3886356; Ma 2022 CHEMBL5131445 | EGFR 0.5 nM confirmed in patent Table 1. HER2 Kd **0.75 nM** confirmed on Ma BindingDB (KINOMEscan). The 76 nM conflict was a **wrong-ligand** BindingDB page. Class remains dual. |
| PM48_04 | PIK3CA/mTOR | dual | 8.85 / 9.35 | Hong 2024 CHEMBL5500428; Yang 2023 CHEMBL5620391 | 8.85 is a **cited** 1.4 nM, not a Hong measurement. Yang mTOR 0.45 nM is a truncated FLAG construct (1362–end). Class remains dual if 8.85 is dropped. |
| **PM48_05** | PIK3CA/mTOR | dual | **10.00 / 8.52** | Bass 2021 review; Elmenier 2019 review; **must-read Folkes 2008** | Both review maxima are **wrong**. Folkes original is PI3Kα **3 nM** / mTOR **Ki 580 nM** (8.52 / 6.24). At θ=6.0 this remains dual. |
| AB_089 | AChE/BChE | neither | 4.78 / 5.12 | Sang 2020 CHEMBL4680246 | Frozen pB matches **human serum BuChE** (7600 nM), not equine. ChEMBL species=horse is a metadata error. Class stays neither. |
| AB_091 | AChE/BChE | neither | 4.76 / 5.19 | same | Human BuChE 6500 nM. Class stays neither. |
| AB_094 | AChE/BChE | neither | 4.71 / 5.29 | same | Human BuChE 5100 nM. Class stays neither. |
| PM48_22 | PIK3CA/mTOR | A_only | 8.77 / 5.83 | Cheng 2021 CHEMBL4765307; Shi 2024 CHEMBL5579880 | Shi PI3Kα 1.70 nM confirmed. Cheng cellular p70S6K row is not on the BindingDB Cheng dump. Class stays A_only. |

No frozen class flips at θ=6.0. Table 2 stays as published.

---

## 1. EH120_045 / Ibrutinib / CHEMBL1873475

### EGFR — US 9,181,263 B2 (CHEMBL3886356) — **verified on local PDF**

Patent pages 99–100 (OCR page 63) and restated as Example 3
(OCR page 86), Table 1:

| kinase | Compound 1 IC50 (nM) | Compound 9 IC50 (nM) |
|--------|----------------------|----------------------|
| BTK | 0.5 | 1.0 |
| **EFGR** (header misspelled) | **0.5** | 20.6 |
| **HER2** | **9.4** | 1536 |
| HER4 | 0.1 | 3.2 |

0.5 nM → pChEMBL **9.30**, matching the frozen EGFR maximum.

Assay protocol (same table caption + Example 4a, OCR pages 86–87):

- in vitro **HotSpot** kinase assay
- **purified enzymes**
- ³³P-ATP phosphorylation of an appropriate peptide substrate
- **1 μM ATP**
- **1 hour** inhibitor incubation
- 10-point curve 10 μM to 0.0005 μM; Prism IC50

EGFR mutation status is **unspecified**. The panel lists EGFR, not
L858R / T790M / exon-19. Treat as unspecified WT commercial kinase,
not a documented full-length cellular EGFR. Residue ranges are not
printed.

Compound 1 identity: “our highly selective BTK inhibitor” with the
ibrutinib chemotype. ChEMBL maps the document to CHEMBL1873475.
Ibrutinib as a drug is the **(R)** piperidine (PCI-32765). Table 1
does not print the R/S label on Compound 1. Residual identity caveat,
not a 100-fold potency error.

Same-document HER2 IC50 9.4 nM → pChEMBL 8.03, matching audit row
CHEMBL3887534. That is **not** the frozen HER2 maximum.

### HER2 — Ma et al., 2022 BMCL (CHEMBL5131445) — **BindingDB dump**

Frozen HER2 maximum is Kd pChEMBL **9.12** (assay CHEMBL5133798) =
**0.76 nM**. BindingDB for the **ibrutinib** ligand BDBM50357312
(`IBRUTINIB | US9181263, 1`) on this article deposits:

| target | Kd (nM) | assay |
|--------|---------|-------|
| BTK | 0.640 | KINOMEscan |
| **ERBB2 / HER2** | **0.75** | KINOMEscan |
| EGFR | 3.90 | KINOMEscan |
| ERBB4 | 0.910 | KINOMEscan |

0.75 nM ≈ frozen 9.12. **This is not a 100-fold error.**

An earlier BindingDB assay page showing 76 nM attached a
**different ligand** (BDBM250082, amide-pyrazole, patent compound 27),
not ibrutinib. That page must not be used for EH120_045.

Ibrutinib is a KINOMEscan **selectivity comparator**. KINOMEscan is a
kinase-domain competition binding assay (typically T7-tagged catalytic
domain), not full-length cellular HER2. Construct caveat remains;
the number itself matches.

Even if the 9.12 row were dropped, remaining HER2 rows include Kd 8.92
and IC50 8.82. EGFR remains 9.30. **Conditional class at θ=6.0: still
dual.**

---

## 2. PM48_04 / Dactolisib / CHEMBL1879463

### PIK3CA — Hong et al., 2024 JMC (PMC11284801, CHEMBL5500428) — **verified on PMC HTML**

The 1.4 nM PI3Kα value is **cited, not measured in this paper**:

> “Dactolisib (1) possesses functional features to potently inhibit
> DNA-PK kinase activity (biochemical IC50: 0.8 nM) and also potently
> inhibits other PIKK and PI3K enzymes (e.g., IC50: PI3Kα, 1.4 nM;
> mTOR, 4.3 nM).”

Hong Table 7 / SAR tables measure peposertib (2), AZD7648 (3), hit 53,
lead 78, and analogues. Dactolisib is the **scaffold-hopping starting
point**, not a compound in those IC50 tables. The 8.85 pChEMBL
(1.4 nM) deposited on CHEMBL5500428 should not be treated as a Hong
experimental result. Primary literature for the citation is likely
Maira 2008 (*Mol. Cancer Ther.*, Hong ref 35).

Hong methods for **their** compounds (Reaction Biology):

- PI3Kα: human **p110α/p85α**, PIP2 10 μM, ATP **10 μM**, ADP-Glo
- mTOR: human mTOR, 4EBP1 1 μM, ATP **10 μM**, HotSpot ³³P
- DNA-PK: human DNA-PK, 10 μM ATP, 20 μM peptide, 10 μg/mL DNA
- n = 2 unless stated

If the 8.85 row is excluded as a citation-not-measured, remaining
PIK3CA maxima are 8.40 and 8.22. **Class remains dual.**

### mTOR — Yang et al., 2023 EJMECH (CHEMBL5620391) — **BindingDB dump**

BEZ235 / dactolisib BDBM92862 on this article:

| target | IC50 (nM) | assay blurb |
|--------|-----------|-------------|
| **mTOR** | **0.450** | N-terminal FLAG-tagged recombinant human mTOR **(1362 to end)**; ULight-4E-BP1 Thr37/46 |
| PI3Kα | 19 | PIP2, 60 min ATP, ADP-Glo |
| PI3Kδ | 78 | same ADP-Glo |
| PI3Kγ | 267 | same ADP-Glo |
| PI3Kβ | 1000 | same ADP-Glo |

0.450 nM → pChEMBL **9.35**, matching the frozen mTOR maximum.

Construct: truncated catalytic fragment (residues 1362–end), not
full-length mTORC1/mTORC2. ATP concentration, n, and SD are still
truncated in the BindingDB blurb (“incubated for …”). Publisher SI
is still required for those three fields. Class is insensitive:
even a much weaker retained mTOR value stays ≥6.

Yang PI3Kα 19 nM (p 7.72) is **not** the frozen 8.85; that 8.85 is
the Hong citation.

---

## 3. PM48_05 / Pictilisib / CHEMBL521851  (highest risk)

Frozen maxima **both come from reviews** and are **not** the Folkes
discovery values.

| arm | frozen | ChEMBL document | BindingDB on that document |
|-----|--------|-----------------|----------------------------|
| PIK3CA | 10.00 (0.1 nM) | CHEMBL5214883 Bass 2021 HDAC-hybrid **review** | BDBM25028 PI3Kα IC50 **0.100 nM** |
| mTOR | 8.52 (3 nM) | CHEMBL4373732 Elmenier 2019 PI3K **review** | BDBM25028 mTOR IC50 **3 nM** |

### Original Folkes 2008 (BindingDB entry 2822 = JMC 51:5522–32)

Ligand BDBM25028 (GDC-0941 / pictilisib), Piramed Pharma deposition,
pH 7.5, T = 2°C:

| target | value | type |
|--------|-------|------|
| PI3Kα WT | **3 nM** | IC50 (Ysi polylysine SPA) |
| PI3Kα H1047R | 3 nM | IC50, same SPA |
| PI3Kα E545K | 3 nM | IC50, same SPA |
| PI3Kδ | 3 nM | IC50 |
| PI3Kβ | 33 nM | IC50 |
| PI3Kγ | 75 nM | IC50 |
| **mTOR** | **580 nM** | **Ki**, GFP-4EBP HTRF |
| PI3KC2β | 670 nM | IC50 |
| Vps34 | 10 000 nM | IC50 |

Already in ChEMBL as `CHEMBL1140078` (pA 8.52 / pB Ki 6.24).
Workman 2011 (PMC3242038) restates the same Folkes profile: p110α
3 nM, p110δ 3 nM, p110β 33 nM, p110γ 75 nM, mTOR 580 nM, DNA-PK
1230 nM, and “193-fold less activity against mTOR compared to
p110α” in later clinical citations.

### Elmenier 2019 — **confirmed column shift**

BindingDB article dump of Elmenier vs Folkes original for the same
ligand BDBM25028:

| isoform | Folkes 2008 | Elmenier 2019 deposit |
|---------|-------------|------------------------|
| PI3Kα | **3 nM** | 33 nM (= Folkes β) |
| PI3Kβ | 33 nM | **3 nM** (= Folkes α/δ) |
| PI3Kδ | 3 nM | 75 nM (= Folkes γ) |
| PI3Kγ | 75 nM | **580 nM** (= Folkes **mTOR**) |
| mTOR | **Ki 580 nM** | **IC50 3 nM** (= Folkes **PI3Kα**) |

Frozen mTOR 8.52 **is Folkes PI3Kα mis-assigned to mTOR**. It must
not remain the reported maximum. Assay descriptions on the Elmenier
dump are empty (“Inhibition of mTOR (unknown origin)”) because this
is a review table, not a measured assay.

### Bass 2021 — review 0.100 nM is not Folkes

Bass BindingDB dump, same ligand BDBM25028:

| isoform | Bass deposit | Folkes original |
|---------|--------------|-----------------|
| PI3Kα | **0.100 nM** | 3 nM |
| PI3Kδ | 4 nM | 3 nM |
| PI3Kβ | 31 nM | 33 nM |

δ and β approximately match Folkes; α is **30-fold too tight**.
Treat 0.100 nM as a review transcription error until the Bass table
is sighted on the publisher PDF. Frozen PIK3CA 10.00 must not remain
the reported maximum.

### What happens if the two review maxima are dropped

Remaining PIK3CA (not from Bass): Kd 9.12, Kd 9.07, then Folkes IC50
8.52. Remaining mTOR (not from Elmenier 8.52): Kd 7.32, IC50 6.87, Kd
6.70, IC50 6.51, Folkes Ki 6.24.

| rule | pA | pB | class at θ=6.0 |
|------|----|----|----------------|
| keep frozen review maxima | 10.00 | 8.52 | dual (current, **not acceptable as reported max**) |
| drop both reviews, keep all other rows | 9.12 (Kd) | 7.32 (Kd) | dual |
| drop reviews; biochemical IC50/Ki only | 8.52 (Folkes) | 6.87 | dual |
| Folkes original pair only | 8.52 | 6.24 | **dual, barely** |
| Folkes + treat 580 nM as inactive because of 193-fold PI3K selectivity | 8.52 | <6 by policy | A_only (**only if a new rule is added**) |

**A silent dual → A_only flip is not justified at the frozen θ=6.0
max-pChEMBL rule.** The review 10.00 / 8.52 pair **must not** remain
the reported maxima. After source correction the ligand is still dual
at θ=6.0 (pA ≈ 8.5–9.1, pB ≈ 6.2–7.3), unless an extra
incomparability / selectivity-ratio rule is added later.

---

## 4. AB_089 / AB_091 / AB_094 — Sang et al., 2020 EJMECH

One document, CHEMBL4680246. Frozen class **neither** for all three
(pA 4.71–4.78, pB 5.12–5.29).

BindingDB article dump (Ellman, 15 min). Paper used **eeAChE + hAChE**
and **eqBChE + hBChE**. Substrates: equine BuChE uses **BTC**; human
BuChE uses **acetylcholine iodide** (same substrate as AChE).

ChEMBL AChE assay type **A** is a metadata error relative to a primary
Ellman IC50 table.

Frozen pB matches **human serum BuChE**, not equine. ChEMBL
species=horse on the BChE assay is therefore a **metadata error**:

| ligand | ChEMBL | frozen pB | human BuChE (BindingDB) | equine BuChE (BindingDB) |
|--------|--------|-----------|-------------------------|--------------------------|
| AB_089 | CHEMBL4792013 | 5.12 | **7600 nM** | 1300 nM |
| AB_091 | CHEMBL4787165 | 5.19 | **6500 nM** | 3900 nM |
| AB_094 | CHEMBL4761179 | 5.29 | **5100 nM** | 2800 nM |

This BindingDB dump did **not** list AChE rows for these three
CHEMBL IDs (only BuChE + MAO-A/B). Frozen pA ~4.71–4.78 (~16–20 μM)
still needs Sang Table 1 from the publisher PDF for compound-number
mapping. Best compound 4f is **not** these three.

**Class stays neither** (all current pChEMBL values stay <6). Species
correction for BChE (horse → human) changes construct notes, not the
label.

---

## 5. PM48_22 / Alpelisib / CHEMBL2396661

### mTOR — Cheng et al., 2021 JMC (CHEMBL4765307)

Audit row: assay type **A**, IC50 pChEMBL **5.83**. This is the
cellular MCF7 **TSC1-null p70S6K Thr389 ELISA**, not purified mTOR.
It is already marked `uncertain` / `incomparable_record=1`.

The BindingDB Cheng 2021 article dump is a PI3Kα biochemical page
(full-length p110α/p85α, residues 322–600 of p85, Sf21) and does
**not** mention alpelisib, BYL719, p70S6K, ELISA, MCF7, or TSC1.
The cellular surrogate is therefore ChEMBL-only / likely SI, not in
this BindingDB extract.

A separate biochemical-looking mTOR row already exists:
CHEMBL4706640 IC50 pChEMBL **5.52** (3000 nM). If Cheng 5.83 is
excluded, pB becomes 5.52. **Still A_only.**

### PIK3CA — Shi et al., 2024 JMC (CHEMBL5579880) — **BindingDB dump**

Alpelisib CHEMBL2396661 / BDBM50436459:

- PI3Kα IC50 **1.70 nM** → pChEMBL **8.77**
- ADP-Glo, 60 min ATP incubation
- construct still “unknown origin” in BindingDB

Class is insensitive to the remaining p110α/p85α vs p110α-only /
ATP-concentration detail.

---

## Actions taken / not taken

- Frozen `human_reviewed_class` left unchanged (7 uncertain, 0 exclude,
  0 class flips).
- Table 2 **not** recomputed.
- BindingDB article dumps archived under
  `literature_sources/bindingdb_article_dumps/`.
- Verdict table: `tables/uncertain_ligand_source_verdicts_v1.csv`.
- Next step after publisher PDFs land: fill `protein_construct` /
  `wildtype_or_mutant` / `incomparable_record` on the nine decisive
  rows; recompute maxima; recompute Table 2 **only** if a frozen class
  actually changes.

## PDF still required (highest first)

1. Folkes 2008 body + SI (table locator for compound 17 = GDC-0941)
2. Elmenier 2019 publisher table (column-shift already proven from
   BindingDB vs Folkes, page locator still missing)
3. Bass 2021 publisher table (0.100 nM vs Folkes 3 nM)
4. Ma 2022 body + SI (KINOMEscan 0.75 nM already matches; kinase-panel
   SI would pin the construct)
5. Sang 2020 Table 1 (AChE mapping for the three neither ligands)
6. Yang 2023 SI (ATP / n / SD for the 0.45 nM mTOR row)
7. Cheng 2021 SI (confirm cellular ELISA)
8. Shi 2024 SI (PI3Kα construct / ATP)
9. Optional: Maira 2008 for the dactolisib 1.4 nM citation
