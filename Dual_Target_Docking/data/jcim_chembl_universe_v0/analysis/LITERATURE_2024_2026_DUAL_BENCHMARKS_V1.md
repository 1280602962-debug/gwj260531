# 2024–2026 dual-target / docking-benchmark literature (verified URLs)

Positioning for this project. Nothing here docks, unfreezes Table 2, or expands K. Cite these URLs in Track A Methods/Discussion; do not paraphrase DualDiff as having four-state labels or verified holo receptors.

## DualDiff — current de-facto “dual-target benchmark”

Zhou, Guan, Zhang, Peng, Wang, Ma. *Reprogramming Pretrained Target-Specific Diffusion Models for Dual-Target Drug Design.* NeurIPS 2024.

| What | URL |
|------|-----|
| NeurIPS abstract | https://proceedings.neurips.cc/paper_files/paper/2024/hash/9ebf3213e396975cce47f2762e87e166-Abstract.html |
| NeurIPS PDF | https://proceedings.neurips.cc/paper_files/paper/2024/file/9ebf3213e396975cce47f2762e87e166-Paper-Conference.pdf |
| OpenReview | https://openreview.net/forum?id=Y79L45D5ts |
| arXiv abs | https://arxiv.org/abs/2410.20688 |
| arXiv PDF | https://arxiv.org/pdf/2410.20688 |
| DOI | https://doi.org/10.48550/arXiv.2410.20688 |
| Code | https://github.com/zhouxiangxin1998/dualdiff |
| Dataset | https://huggingface.co/datasets/zhouxiangxin/DualDiff |

**What DualDiff actually built (§3.1 and §4.1, NeurIPS PDF):** 12,917 synergistic drug combinations → treated as target pairs; §4.1 also writes “438 unique targets”; GitHub maps indices 0–437 to **reference-ligand SMILES** (one reference ligand per target). Pairing source is **DrugCombDB** (ZIP, Bliss, Loewe and HSA all positive in at least one cell line), not paired experimental four-state activity. Structures: PDBBind if present, else PDB, else **AlphaFold DB** with pLDDT < 70 dropped. Pockets: **P2Rank** “most possible pocket”, then **AutoDock Vina** to place the reference ligand. Evaluation: Vina Dock scores (P-1 / P-2 / Dual High Affinity). No paired four-class labels, no per-end verified holo, no hard-negative class.

Synergy / structure sources DualDiff itself cites:

| Source | URL |
|--------|-----|
| DrugCombDB (Liu et al., *NAR* 2020) | https://doi.org/10.1093/nar/gkz1007 — http://drugcombdb.denglab.org/ |
| P2Rank (Krivák & Hoksza, 2018) | https://doi.org/10.1186/s13321-018-0285-8 — https://github.com/rdk/p2rank |
| AlphaFold DB | https://alphafold.ebi.ac.uk/ |
| PDBBind | http://www.pdbbind.org.cn/ |

## Tightening single-target / docking standards (same window)

### PLINDER — leakage-aware protein–ligand splits + PoseBusters

Durairaj et al. *PLINDER: The protein-ligand interactions dataset and evaluation resource.* bioRxiv 2024.

| What | URL |
|------|-----|
| DOI | https://doi.org/10.1101/2024.07.17.603955 |
| bioRxiv v3 | https://www.biorxiv.org/content/10.1101/2024.07.17.603955v3 |
| Docs / citation | https://plinder-org.github.io/plinder/citation.html |
| Site | https://www.plinder.sh/ |
| Code | https://github.com/plinder-org/plinder |

### PoseBusters — physical validity of poses (PLINDER uses this bar)

Buttenschoen, Morris, Deane. *PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences.* *Chem. Sci.* 2024, 15, 3130–3139.

| What | URL |
|------|-----|
| Publisher | https://doi.org/10.1039/D3SC04185A |
| PMC | https://pmc.ncbi.nlm.nih.gov/articles/PMC10901501/ |
| arXiv | https://arxiv.org/abs/2308.05777 |
| Code | https://github.com/maabuu/posebusters |
| Docs | https://posebusters.readthedocs.io/ |
| Data | https://zenodo.org/records/8278563 |

### LIT-PCBA leakage audit — explicit call to rebuild benchmarks

Huang, Knight, Naprienko. *Data Leakage and Redundancy in the LIT-PCBA Benchmark.* arXiv:2507.21404 (2025).

| What | URL |
|------|-----|
| arXiv abs | https://arxiv.org/abs/2507.21404 |
| HTML | https://arxiv.org/html/2507.21404 |
| PDF | https://arxiv.org/pdf/2507.21404 |
| Code (canonical) | https://github.com/sievestack/LIT-PCBA-audit |
| Mirror | https://github.com/mireklzicar/lit-pcba-audit |

Original LIT-PCBA (Tran-Nguyen, Jacquemard, Rognan, *JCIM* 2020) — the set being audited, not a dual-target benchmark:

| What | URL |
|------|-----|
| DOI | https://doi.org/10.1021/acs.jcim.0c00155 |
| PubMed | https://pubmed.ncbi.nlm.nih.gov/32282202/ |
| Dataset | http://drugdesign.unistra.fr/LIT-PCBA |

The 2025 audit’s conclusion is that LIT-PCBA in its current form does not measure recovery of novel chemotypes and **should not be treated as evidence of methodological progress**; they do not release a “cleaned” LIT-PCBA and call for new benchmark development.

### TopU-LBVS — property-matched hard negatives (1:40); EF@1% drops ~4× vs random decoys

Kumar, Zhou, Shiralkar, Huang, Coskunuzer. *TopU-LBVS: A Realistic Multi-Target Benchmark for Ligand-Based Virtual Screening.* NeurIPS 2026 Datasets & Benchmarks, **under review** (no journal DOI yet).

| What | URL |
|------|-----|
| Code | https://github.com/topu-benchmark/topu-lbvs |
| Dataset | https://huggingface.co/datasets/topu-benchmark/topu-lbvs |
| Org | https://huggingface.co/topu-benchmark |

Cite as unpublished / under review. Do not invent a DOI.

## Independent medicinal-chemistry evidence: public dual space is paralog-biased

Lembo & Bottegoni. *Systematic Investigation of Dual-Target-Directed Ligands.* *J. Med. Chem.* 2024, 67, 10374–10385. (DTDLs, not Azure DTDL.)

| What | URL |
|------|-----|
| DOI | https://doi.org/10.1021/acs.jmedchem.4c00838 |
| Publisher | https://pubs.acs.org/doi/10.1021/acs.jmedchem.4c00838 |
| PMC (OA) | https://pmc.ncbi.nlm.nih.gov/articles/PMC11215722/ |
| PubMed | https://pubmed.ncbi.nlm.nih.gov/38843874/ |

158 rationally conceived DTDLs from ChEMBL. Finding used here: target-pair choice and chemistry are driven by known pathology associations, **existing scaffolds**, and **binding-pocket similarity** — independent confirmation that public dual-target space is biased toward paralogs / shared chemotypes, not novel cross-pathway pairs.

## How this project differs (one paragraph, claim-ceiling safe)

DualDiff asserts 12,917 synergy-derived pairs with predicted or docked pockets and Vina evaluation. Under paired experimental four-state labels, human holo receptors, and drug-like hard negatives, the ChEMBL 37 census yields **17** feasible pairs, of which **8** are conventionally dockable over **6** independent systems (`TIER1_DOCKING_ROSTER_V1.md`). That comparison is a **feasibility audit**, not a DualDiff bake-off (CLAIM_CEILING item 44 still forbids re-scoring DualDiff molecules on DualFourClass-Bench). Site verification (`SITE_VERIFICATION_CHECKLIST_V1.md`) is the layer DualDiff does not have: P2Rank + AlphaFold + one reference ligand cannot substitute for a human-confirmed cognate in a declared pocket.
