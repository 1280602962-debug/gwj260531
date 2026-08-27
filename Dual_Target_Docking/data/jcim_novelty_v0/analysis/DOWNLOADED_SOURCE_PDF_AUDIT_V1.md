# Downloaded source PDF audit v1

Date: 2026-08-27
Input folder: `D:\CADD paper exercise\dual target docking\check`
Machine-readable inventory: `../tables/downloaded_source_pdf_inventory_v1.csv`

## Scope and interpretation

This pass performs two different checks that must not be conflated:

1. **File identity/integrity screening** for every mandatory document.
2. **Assay-level human reading** for the decisive documents controlling
   the seven previously uncertain ligand labels.

Passing the first check means only that a readable PDF appears to be the
intended publication. It does not mean that every potency value, construct,
mutation, ATP concentration, replicate count, or supporting-information
method has been verified.

## Inventory result

- Mandatory document records: **135**
- Downloaded files mapped to those records after replacements: **130**
- Exact DOI matches: **105**
- Exact patent-number matches: **6**
- Probable title matches: **9**
- Manual identity checks: **10**
- Invalid/non-PDF files in the active folder: **0**
- Missing database/deposited-dataset records: **5**
- Exact SHA-256 duplicate groups: **0**

All conventional article and patent records in the mandatory queue now have
an active local file. The five unmatched records are database/deposited-
dataset sources rather than conventional papers.

## Previously detected files that were replaced or renamed

| supplied filename | original finding | resolution |
|---|---|---|
| `CHEMBL1134862` | HTML security-verification page, not a PDF | correct DOI/PII article subsequently supplied and visually confirmed |
| `CHEMBL1131301.pdf` | contains the CHEMBL1129291 “Tyrosine kinase inhibitors 8” article | replaced with the article corresponding to DOI `10.1021/jm970641d` |
| `CHEMBL1145585.pdf` | unrelated 2025 Dimroth-rearrangement review | replaced with DOI `10.1016/j.bmcl.2008.05.024` |
| `CHEMBL3638541.pdf` | contains US8623883B2, which belongs to CHEMBL3639119 | provenance corrected and US8772480B2 separately supplied for CHEMBL3638541 |
| `CHEMBL1125266].pdf` | content appears correct; filename contains an extra `]` | identity accepted; filename normalization remains recommended |

On 2026-08-27 the first three files were moved to
`check\_rejected_wrong_sources\` rather than permanently erased, and the
misnamed `CHEMBL3638541.pdf` was renamed to `CHEMBL3639119.pdf`. The active
folder therefore no longer exposes those wrong files to the audit script.
Correct replacements were subsequently supplied for CHEMBL1131301,
CHEMBL1145585, and CHEMBL3638541 and passed title/DOI/patent inspection.
CHEMBL1134862 was subsequently supplied and visually confirmed against the
expected title and PII `S0960-894X(02)00364-5`.

## Repository bibliography corrections

- **CHEMBL1141539**: the supplied PDF title matches the expected article,
  but its DOI is `10.1016/j.bmcl.2007.08.073`. The repository's current
  DOI `10.1016/j.bmcl.2008.02.035` is inconsistent with the article and
  should be treated as a cached mapping error pending regeneration.
- **CHEMBL5532641** is resolved to “MNK, mTOR or eIF4E—selecting the best
  anti-tumor target for blocking translation initiation”, DOI
  `10.1016/j.ejmech.2023.115781`.
- **CHEMBL5584230** is resolved to “Exploring fluorine-substituted
  piperidines as potential therapeutics for diabetes mellitus and
  Alzheimer's diseases”, DOI `10.1016/j.ejmech.2024.116523`.

These are bibliography/provenance repairs. They do not by themselves alter
the frozen model labels.

## Missing records: what to obtain

| record | interpretation | action |
|---|---|---|
| CHEMBL1201862 | ChEMBL kinase-profiling dataset/source | export the ChEMBL document/assay records; do not search indefinitely for a journal PDF |
| CHEMBL1909046 | DrugMatrix kinase-inhibition dataset | preserve the database record and assay metadata as the source artifact |
| CHEMBL5210307 | EUbOPEN wave-3 deposited dataset, DOI `10.6019/CHEMBL5210307` | save the repository landing page/data attachment and version/date |
| CHEMBL5446079 | EUbOPEN deposited dataset | save landing page/data attachment and version/date |
| CHEMBL5465560 | EUbOPEN selectivity dataset | save landing page/data attachment and version/date |

## Decisive assay findings

### Sun EGFR/HER2 compound 26f (EH120_059)

The replacement CHEMBL1134862 article is correct. Table 1 reports an
EGF-R biochemical kinase IC50 of 0.029 μM (pIC50 7.54). The deposited
HER2 value of 49.3 μM (pIC50 4.31) comes from the Table 2 `3T3
Proliferation / Her2` column and is a ligand-driven BrdU cellular response,
not direct HER2 enzyme inhibition. It is therefore marked as a cellular
surrogate/incomparable with biochemical potency. The frozen `A_only` class
does not change because the HER2 value was already below threshold.

### Ibrutinib (EH120_045)

Ma 2022 Table 6 supports EGFR Kd 3.9 nM and ERBB2/HER2 Kd 0.75 nM.
The ChEMBL pKd values 8.41/9.12 are supported. The BindingDB 76 nM
record is not supported by this table and is coupled to a ligand-mapping
discrepancy. The main PDF identifies DiscoverX KINOMEscan but does not
print the exact construct or mutation state. Label remains `dual`.

### Dactolisib (PM48_04)

Hong 2024 merely cites PI3Kα 1.4 nM and mTOR 4.3 nM for dactolisib; it
does not newly measure those dactolisib values. Those rows are secondary
evidence. Yang 2023 Table 1 directly measures dactolisib mTOR IC50
0.45 nM in duplicate dose responses. ATP concentration and construct are
not stated in the main PDF. Label remains `dual`.

### Pictilisib (PM48_05)

Elmenier 2019 states PI3Kα/β/δ/γ IC50 = 3/33/3/75 nM and mTOR IC50 =
580 nM. The ChEMBL mTOR 3 nM review maximum is a target/table-shift error.
Bass 2021 is also a review and does not establish an original 0.1 nM
PI3Kα assay. Excluding both review maxima leaves the class `dual` at the
frozen θ=6 rule, but with much less impressive and more defensible maxima.
The Folkes 2008 primary paper is now present as `CHEMBL1140078.pdf`.
Table 4 directly reports p110α IC50 0.003 μM and mTOR Kiapp 0.58 μM.
The assay methods specify recombinant p110α/p85α GST-fusion SPA with
1 μM ATP and an mTOR GFP-4E-BP1 assay with 8 μM ATP. The missing SI is
described by the paper as compound-purity/HPLC material and does not block
the potency or assay-context conclusion.

### Sang AChE/BChE compounds (AB_089, AB_091, AB_094)

The deposited values map exactly to Sang Table 1 compounds 4c, 4b, and
4a, respectively. They are human erythrocyte AChE and human serum BChE
Ellman assays with mean ± SD from three experiments. The AChE assay-type
`A` flag is a ChEMBL metadata error; both arms are biochemical potency
measurements. All three labels remain `neither`.

### Alpelisib (PM48_22)

Cheng 2021 mTOR pChEMBL 5.83 is a cellular TSC1-null p-p70S6(Thr389)
ELISA surrogate and must not be merged with purified-mTOR potency. Shi
2024 directly measures alpelisib PI3Kα IC50 1.7 nM using ADP-Glo, Carna
11-101 enzyme, 10 μM ATP, 25 μM PIP2:3PS, 60 min, and a six-point curve.
Label remains `A_only`.

## What this improves for the manuscript

This audit materially strengthens source traceability and prevents three
reviewer-visible errors: using a review-derived table-shift as a primary
measurement, mixing a cellular pathway surrogate with biochemical potency,
and claiming protein constructs that the article never reports. It improves
the credibility of the dataset and limitations section, but it does not
replace independent validation, robust scaffold/time splits, negative
controls, uncertainty analysis, or external docking validation.

## Remaining evidence gaps

1. SI for Yang 2023, Ma 2022, Cheng 2021, and Shi 2024 where available.
2. Optional Folkes 2008 SI for purity/HPLC provenance.
3. Primary Maira 2008 source for dactolisib 1.4/4.3 nM if those values are
   retained in a primary-source sensitivity analysis.
4. Explicit snapshot/version files for ChEMBL, DrugMatrix, and EUbOPEN
   dataset documents.
5. Complete paper-level source-tier adjudication for every record capable of
   determining a ligand-target maximum. The scripted decision-targeted
   primary-only sensitivity analysis is now deposited as
   `PRIMARY_ONLY_LABEL_MODEL_SENSITIVITY_V1.md`; it does not silently
   overwrite frozen labels.
