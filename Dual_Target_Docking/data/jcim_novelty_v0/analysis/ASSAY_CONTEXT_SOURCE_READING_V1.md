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

Local PDFs supplied on 2026-08-27 were identity/integrity screened and the
decisive tables in Ma 2022, Hong 2024, Yang 2023, Elmenier 2019, Bass 2021,
Sang 2020, Cheng 2021, Shi 2024, and Folkes 2008 were read. Supporting
information for Ma, Yang, Cheng, and Shi is still not local. The Folkes SI
is also absent, but the paper states that it contains compound-purity tables
and HPLC traces rather than the decisive assay protocol.
The main papers are sufficient to decide the recorded values and evidence
tier, but not to invent unreported residue ranges, ATP concentrations, or
mutation states.

Label rule used below: DualFourClass θ = 6.0 on the per-target maximum
pChEMBL. Class flips are recorded only as *conditional* outcomes.

## Frozen vs source-checked maxima

| ligand | pair | frozen class | frozen pA / pB | decisive docs | source-checked status |
|--------|------|--------------|----------------|---------------|------------------------|
| EH120_045 | EGFR/HER2 | dual | 9.30 / 9.12 | patent CHEMBL3886356; Ma 2022 CHEMBL5131445 | EGFR 0.5 nM is in patent Table 1. Ma Table 6 prints EGFR Kd 3.9 nM and HER2 Kd **0.75 nM**; the BindingDB 76 nM page belongs to a different ligand and is not evidence against this value. |
| PM48_04 | PIK3CA/mTOR | dual | 8.85 / 9.35 | Hong 2024 CHEMBL5500428; Yang 2023 CHEMBL5620391 | 8.85 is a **cited** 1.4 nM, not a Hong measurement. Yang Table 1 directly measures dactolisib mTOR IC50 **0.45 nM**, duplicate dose response; ATP/construct remain unreported in the main PDF. |
| **PM48_05** | PIK3CA/mTOR | dual | **10.00 / 8.52** | Bass 2021 review; Elmenier 2019 review; Folkes 2008 CHEMBL1140078 | Review maxima are not original Pictilisib values. Folkes Table 4 directly gives p110α IC50 3 nM and mTOR Kiapp 580 nM (8.52 / 6.24). At θ=6.0 this remains dual. |
| AB_089 | AChE/BChE | neither | 4.78 / 5.12 | Sang 2020 CHEMBL4680246 | Table 1 maps it to 4c: human erythrocyte AChE 16.6 μM and human serum BChE 7.6 μM. ChEMBL type A on AChE is a metadata error. |
| AB_091 | AChE/BChE | neither | 4.76 / 5.19 | same | Table 1 compound 4b: hAChE 17.3 μM, hBChE 6.5 μM. |
| AB_094 | AChE/BChE | neither | 4.71 / 5.29 | same | Table 1 compound 4a: hAChE 19.7 μM, hBChE 5.1 μM. |
| PM48_22 | PIK3CA/mTOR | A_only | 8.77 / 5.83 | Cheng 2021 CHEMBL4765307; Shi 2024 CHEMBL5579880 | 5.83 is cellular TSC1-null p-p70S6(Thr389) ELISA. Shi directly measures alpelisib PI3Kα IC50 1.7 nM by ADP-Glo using Carna 11-101, 10 μM ATP, 60 min. |

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

Ma Table 6 directly prints ibrutinib Kd **3.9 nM for EGFR** and
**0.75 nM for ERBB2/HER2**. These convert to pKd 8.41 and 9.12 and
exactly support the ChEMBL rows. The paper describes the DiscoverX
KINOMEscan platform and reports the experiment as a comparator panel,
not as a designed HER2 program. The main PDF does not report a residue
range or explicit WT/mutant state. Record `commercial KINOMEscan;
construct not reported` and `mutation unspecified`; do not infer a
full-length receptor or a T7-tagged kinase-domain construct from
platform convention alone.

The BindingDB 76 nM/non-ibrutinib record is therefore treated as a
BindingDB mapping/decimal inconsistency, not as evidence to demote the
ChEMBL value. **Class remains dual.**

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

Yang Table 1 directly reports dactolisib as the reference control:
mTOR IC50 **0.45 nM** (pChEMBL 9.35). The table footnote states that
8- or 10-point dose responses were run in duplicate, and the methods
state that the mTOR enzyme screen was performed by Wuxi Bioduro
Biologics. The main PDF does not print the ATP concentration or protein
residue range; those fields remain `not reported in main paper` pending
SI. This does not control the class: even a much weaker retained mTOR
value stays ≥6.

The archived BindingDB assay blurb additionally describes an N-terminal
FLAG-tagged recombinant human mTOR fragment (residues 1362–end). Because
that construct detail is not printed in the main PDF, it is retained as
database metadata rather than attributed directly to Yang's article.

---

## 3. PM48_05 / Pictilisib / CHEMBL521851  (highest risk)

Frozen maxima **both come from reviews** and are **not** the Folkes
discovery values.

| arm | frozen | ChEMBL document | BindingDB on that document |
|-----|--------|-----------------|----------------------------|
| PIK3CA | 10.00 (0.1 nM) | CHEMBL5214883 Bass 2021 HDAC-hybrid **review** | BDBM25028 PI3Kα IC50 **0.100 nM** |
| mTOR | 8.52 (3 nM) | CHEMBL4373732 Elmenier 2019 PI3K **review** | BDBM25028 mTOR IC50 **3 nM** |

### Original Folkes 2008 (BindingDB entry 2822 = JMC 51:5522–32)

The local Folkes paper directly establishes:

- PI3Kα IC50 **3 nM** → pChEMBL **8.52** (audit rows CHEMBL967868,
  CHEMBL967881, CHEMBL980850).
- mTOR **Kiapp 0.58 μM**; ChEMBL stores **Ki 6.24** on CHEMBL980857
  (~575 nM). It must not be described as an IC50.
- Sarker 2015: “193-fold less activity against mTOR compared to
  p110α”.
- DNA-PK ~1230 nM in the same selectivity profile.

Folkes Table 4 further shows wild-type p110α, E545K and H1047R all at
0.003 μM. Values are averages of at least two determinations with typical
variation below ±30%. The class-I PI3K SPA used recombinant human
p110α/β/δ coexpressed with p85α as purified GST fusions, 1 μM ATP,
33P-ATP, 1 h at room temperature. The mTOR assay used 8 μM ATP,
GFP-4E-BP1, 30 min at 25 °C, duplicate dose-response curves, and a
competitive tight-binding fit. The exact mTOR construct is not printed.

Elmenier section 4.1.6 explicitly prints PI3Kα/β/δ/γ IC50 values
3/33/3/75 nM and mTOR IC50 580 nM. The deposited mTOR 3 nM is therefore
the PI3Kα table value shifted onto the wrong target. Bass is a broad
HDAC-hybrid review; it discusses pictilisib as compound 114 but does
not provide an original assay protocol supporting a 0.1 nM PI3Kα
measurement. Both review maxima are excluded from a primary-source
sensitivity analysis. Folkes now resolves the endpoint ambiguity: the
580 nM mTOR value is **Kiapp**, not IC50. Missing Folkes SI affects the
purity/HPLC audit, not the reported Table 4 assay values or methods.

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

Sang Table 1 and its footnotes resolve the mapping:

- CHEMBL4792013 = 4c: hAChE 16.6 ± 0.17 μM; hBChE 7.6 ± 0.26 μM.
- CHEMBL4787165 = 4b: hAChE 17.3 ± 0.24 μM; hBChE 6.5 ± 0.28 μM.
- CHEMBL4761179 = 4a: hAChE 19.7 ± 0.32 μM; hBChE 5.1 ± 0.12 μM.
- Values are means ± SD of three experiments. hAChE is from human
  erythrocytes and hBChE from human serum; eel AChE and equine BChE are
  separate columns.

Both target arms are biochemical Ellman enzyme IC50 measurements.
ChEMBL assay type **A** on the AChE rows is a metadata error and should
not make those rows incomparable. All values remain below θ=6.0, so
the three labels remain `neither`.

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

Alpelisib as PI3Kα control: IC50 **1.7 nM** → pChEMBL **8.77**.
The main methods specify Promega ADP-Glo, PI3Kα Carna #11-101 at
1.25 nM, **10 μM ATP**, 25 μM PIP2:3PS, 60 min at room temperature,
and a six-point curve. The commercial catalogue number is the most
specific construct identifier printed; the paper does not print a
residue range or mutant state. Class is insensitive to that detail.

---

## Actions taken / not taken

- Frozen `human_reviewed_class` was left unchanged; paper-level findings are
  stored separately in `decisive_source_human_review_v1.csv` rather than
  silently overwriting the frozen benchmark.
- Source-verified construct, mutation, assay-context, and comparability fields
  were recorded for the decisive rows where the papers provide them;
  genuinely unreported fields remain `unspecified`.
- A decision-targeted primary-only label and model sensitivity analysis was
  run and deposited in `PRIMARY_ONLY_LABEL_MODEL_SENSITIVITY_V1.md` and its
  companion CSV files. The frozen Table 2 remains the prespecified primary
  analysis; the new result is an SI sensitivity analysis.
- One ligand, EH120_059, becomes `unresolved_missing_arm` after its HER2
  cellular surrogate is removed. Vina `summary_min` is unchanged for all four
  pairs; this does not constitute a claim that every ChEMBL source record has
  been paper-level adjudicated.
- BindingDB article dumps remain archived under
  `literature_sources/bindingdb_article_dumps/`; the earlier verdict table is
  retained as `tables/uncertain_ligand_source_verdicts_v1.csv`.

## Material still required (highest first)

1. Yang 2023 SI (ATP / construct / error for the 0.45 nM mTOR assay)
2. Ma 2022 SI, if available (KINOMEscan construct identifiers)
3. Cheng 2021 SI (replicate statistics for the cellular mTOR surrogate)
4. Shi 2024 SI (replicate statistics and any construct detail beyond catalogue ID)
5. Optional: Folkes 2008 SI for compound purity/HPLC only
6. Optional: Maira 2008 for the primary dactolisib 1.4/4.3 nM citation
