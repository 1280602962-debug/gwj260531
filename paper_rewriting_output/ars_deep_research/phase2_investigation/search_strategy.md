# Phase 2 — Search Strategy（可复现）

**Search date:** 2026-07-17  
**Primary index:** OpenAlex API (`api.openalex.org/works`)  
**Filters:** `publication_year:2019-2026`, `type:article|review`  
**Polite pool:** mailto header set  

## Boolean / queries

| ID | Query focus | OpenAlex hits (meta.count) | Retrieved |
|----|-------------|----------------------------|-----------|
| Q1 | covalent + VS/docking + acrylamide/warhead + kinase | 1406 | 25 |
| Q2 | scaffold hopping + VS/fingerprint/pharmacophore + kinase/covalent | 930 | 25 |
| Q3 | AlphaFold3 + covalent + ligand/screening | 1072 | 25 |
| Q4 | JNK2/JNK + covalent + inhibitor/acrylamide | 3183 | 25 |
| Q5 | active learning + VS + docking/enrichment | 523 | 25 |
| Q6 | acrylamide/warhead + GSH/reactivity + covalent | 6763 | 25 |
| Q7 | CovDock/DOCKovalent + VS/enrichment | 622 | 25 |
| Q8 | enumeration/building-block acrylamide / REAL covalent | — | **429 rate-limit** |

**Dedup:** 160 unique works across Q1–Q7.  
**Supplementary channels:** ACS full-text fetches (RosettaAMRLD, CovalentLab, extended warheads JCIM); prior project corpus (Lu 2023, Wydra 2025, COValid/JACS AF3, Phase0 notes).

## Screening
- Pass 1: title/abstract relevance to pre-docking triage / covalent library / JNK / enrichment.  
- Pass 2: keep methodologically transferable sources even if not JNK-specific.  
- Grade: peer-reviewed non-OA preferred for “field practice”; OA reviews for warhead landscape.
