# Local PDFs for the 7 uncertain DualFourClass ligands

Drop publisher PDFs (and supporting information) into `pdfs/`.
Exact filenames below are matched by the source-reading notes in
`analysis/ASSAY_CONTEXT_SOURCE_READING_V1.md`. Extra files are fine;
wrong names just slow matching.

Do **not** commit PDFs. The directory is gitignored.

## Filename checklist (priority order)

| file | ligand / role |
|------|----------------|
| `US9181263B2.pdf` | EH120_045 EGFR (patent) |
| `Ma_2022_BMCL_128549.pdf` | EH120_045 HER2 Kd |
| `Ma_2022_BMCL_128549_SI.pdf` | same, kinase-panel SI |
| `Hong_2024_JMC_DNA-PK.pdf` | PM48_04 PIK3CA |
| `Hong_2024_JMC_DNA-PK_SI.pdf` | same |
| `Yang_2023_EJMECH_115543.pdf` | PM48_04 mTOR |
| `Yang_2023_EJMECH_115543_SI.pdf` | ATP / n / error |
| `Elmenier_2019_EJMECH_111718.pdf` | PM48_05 mTOR review table |
| `Bass_2021_EJMECH_112904.pdf` | PM48_05 PIK3CA review table |
| `Folkes_2008_JMC_jm800295d.pdf` | **PM48_05 original Pictilisib** |
| `Folkes_2008_JMC_SI.pdf` | same |
| `Sang_2020_EJMECH_112265.pdf` | AB_089 / AB_091 / AB_094 |
| `Sang_2020_EJMECH_112265_SI.pdf` | Table 1 species columns |
| `Cheng_2021_JMC_0c01652.pdf` | PM48_22 mTOR cellular ELISA |
| `Cheng_2021_JMC_0c01652_SI.pdf` | same |
| `Shi_2024_JMC_4c00992.pdf` | PM48_22 PIK3CA |
| `Shi_2024_JMC_4c00992_SI.pdf` | protein construct / ATP |

Optional primary for the dactolisib PI3Kα citation in Hong 2024:

| `Maira_2008_MCT_BEZ235.pdf` | PM48_04 original BEZ235 / dactolisib |

ACS/Elsevier HTML exports (`.html`) and page-range screenshots of the
activity tables are also usable if a PDF is not available.

DOI downloads on 2026-08-27 did **not** return publisher PDFs (BindingDB
HTML or reCAPTCHA). The BindingDB article pages that were actually read
are archived in `bindingdb_article_dumps/`. Hong 2024 was read from
PMC11284801. The ibrutinib patent PDF is present locally and gitignored.
