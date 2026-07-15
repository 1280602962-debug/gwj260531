# Citation Bank (Option A — JCIM-oriented)

Local-first; map each claim to a citable source. Prefer ACS/high-quality venues already curated in `docs/REFERENCES.md` and `LITERATURE_COMPARISON_zh.md`.

## A. Field motivation (JNK1 / fibrosis / isoform problem)

| ID | Claim in manuscript | Source |
|----|---------------------|--------|
| CIT-A1 | JNK pathway therapeutically relevant | Manning & Davis, *Nat Rev Drug Discov* 2003 (REFERENCES #6) |
| CIT-A2 | CC-90001 JNK1-biased clinical candidate | Bennett et al., *J Med Chem* 2021 (REFERENCES #7; doi:10.1021/acs.jmedchem.1c01716) |
| CIT-A3 | E1 potent JNK1 inhibitor / IPF | Pan et al., *J Med Chem* 2024 (doi:10.1021/acs.jmedchem.4c01764) |
| CIT-A4 | SP600125 pan-JNK tool | Bennett et al., *PNAS* 2001 (REFERENCES #8) |
| CIT-A5 | ATP-site conservation / selectivity hard | Duong et al., *Comput Struct Biotechnol J* 2020 (REFERENCES #11); Kinase-Bench context |

## B. Method peers (selectivity computation)

| ID | Claim | Source |
|----|-------|--------|
| CIT-B1 | Docking/IFP can guide subtype selectivity in some kinases | Bajusz et al., *JCIM* 2016 (doi:10.1021/acs.jcim.5b00634) |
| CIT-B2 | Need tailored VS benchmarks for kinase selectivity | Kinase-Bench, *JCIM* 2024 (doi:10.1021/acs.jcim.4c01830) |
| CIT-B3 | FEP selectivity readiness / error cancellation | Albanese et al., *JCIM* 2020 (doi:10.1021/acs.jcim.0c00815) |
| CIT-B4 | ML+structure VS → buy → enzyme assay pipeline exemplar | *JCIM* 2020 successive modeling (doi:10.1021/acs.jcim.9b01204) |
| CIT-B5 | Honest µM / limited selectivity reporting style | CCR2/CCR5, *JCIM* 2025 (doi:10.1021/acs.jcim.5c01596) |

## C. Tools & data

| ID | Claim | Source |
|----|-------|--------|
| CIT-C1 | ChEMBL bioactivity | Zdrazil et al., *NAR* 2024 (REFERENCES #1) |
| CIT-C2 | ECFP / Morgan fingerprints | Rogers & Hahn, *JCIM* 2010 (REFERENCES #21) |
| CIT-C3 | Murcko scaffolds | Bemis & Murcko, *J Med Chem* 1996 (REFERENCES #22) |
| CIT-C4 | PAINS filters | Baell & Holloway, *J Med Chem* 2010 (REFERENCES #23) |
| CIT-C5 | RDKit | Landrum (REFERENCES #20) |
| CIT-C6 | XGBoost | Chen & Guestrin, KDD 2016 (REFERENCES #28) |
| CIT-C7 | Glide docking (if licensed Methods) | Friesner et al. Glide papers (add from REFERENCES docking section) |
| CIT-C8 | AutoDock Vina (C2 fallback) | Eberhardt et al., *JCIM* 2021 / Trott & Olson 2010 |

## D. JNK structural panel PDBs (Methods)

| ID | PDB | Role |
|----|-----|------|
| CIT-D1 | 3ELJ | JNK1 primary docking/MD |
| CIT-D2 | 3E7O | JNK2 sole drug co-crystal used |
| CIT-D3 | 3TTI | JNK3 / CC-930 |
| CIT-D4 | 4L7F / 4WHZ | Ensemble secondary (benchmark Δsel) |

## E. Project-internal evidence (not bibliography — cite as data/SI)

| ID | Artifact | Supports |
|----|----------|----------|
| DAT-E1 | `results/selectivity_autopsy/` | RQ-C negative result |
| DAT-E2 | `results/chemotype_novelty/` | Chemotype distance |
| DAT-E3 | `results/assay_analysis/` | Pre-registered IC50 rules |
| DAT-E4 | `results/pose_consensus/` | C2 Vina multi-seed |
| DAT-E5 | `results/c11_2231_comparison/` | Unbought 2231 opportunity cost |
| DAT-E6 | `results/purchase_risk/` | PAINS/physchem |
| DAT-E7 | `docs/JNK1_PROJECT_REPORT.md` | Funnel numbers |

## Do-not-cite-as-peer-anchors

Low-tier MDPI OA hit papers (*Molecules* Tricin; *Pharmaceuticals* JNK3 DL) — optional one-line background only; not contribution templates.
