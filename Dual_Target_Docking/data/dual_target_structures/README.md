# Dual-target structural benchmark (starter curation)

Curated seed set of **published dual-/multi-target ligands with PDB evidence**, for evaluating dual-target docking protocols.

## Files

| File | Description |
|------|-------------|
| `dual_target_cocrystal_catalog.csv` | Master catalog (one row per ligand–target-pair case) |
| `README.md` | This file: classification scheme, use cases, expansion protocol |

## Classification scheme

### 1. Structure completeness (`structure_tier`)

| Tier | Meaning | Use for docking eval |
|------|---------|----------------------|
| **A_both_ends** | Same dual ligand has co-crystal (or clear complex) structures **for both targets** | Gold standard for dual pose recovery (RMSD on A and B) |
| **B_single_end** | Molecule is a published dual inhibitor, but PDB currently covers **only one** target | Strict pose eval on crystallized end only; other end = docking-only |
| **C_series_related** | Dual-target series with structures of **related** ligands on each target (not necessarily identical chem_comp on both) | Weaker pose transfer / cross-docking stress test |
| **D_claim_only** | Dual activity claimed; structural work uses single-target templates only | Do **not** use for pose RMSD of the dual ligand |

### 2. Design type (`design_type`, Morphy-style)

| Type | Definition used here |
|------|----------------------|
| **linked** | Two pharmacophores joined by an explicit linker / tethering |
| **fused** | Pharmacophores joined with short connection / partial atom sharing |
| **merged** | Highly overlapping / shared scaffold serving both pockets |
| **unclear** | Not enough public description to assign confidently |

### 3. Intended evaluation roles

| Dataset role | Recommended metrics |
|--------------|---------------------|
| Tier A | Self-docking RMSD ≤ 2 Å on **both** ends; interaction fingerprint recovery |
| Tier B | RMSD on crystallized end; optional cross-docking |
| Activity labels (separate) | Dual vs single enrichment (ChEMBL / paper IC50) |
| Your PROTAC set | Linked bifunctional / ternary protocol validation |

**Using this catalog + PROTAC activity data together is appropriate** for a methods paper: structures → pose quality; PROTAC labels → functional discrimination.

## Seed inventory summary (v0.1)

| Tier | n cases (approx.) | Examples |
|------|-------------------|----------|
| A_both_ends | 5 ligand systems | Mcl-1/Bcl-xL compound 10; LpxA/LpxD Q5M; PknA/PknB CJJ; EGFR/HER2 TAK-285 |
| B_single_end | 10+ | BET–HDAC (BRD4 only); ER dual modulators (ER only); MurD/MurE (MurD only); PD-L1/VISTA (PD-L1 only) |
| C_series_related | few | GyrB / ParE pyrrolopyrimidine series (different ligands per structure) |

Counts will grow as the catalog is expanded; see expansion protocol below.

## Expansion protocol (recommended)

1. Literature keywords: `dual inhibitor` / `dual-target` / `MTDL` / `framework combination` + `crystal` / `PDB`.
2. For each paper, extract ligand identifiers and **all** deposited PDBs from the same study.
3. Confirm **same ligand** (chem_comp ID or chemical identity) on both targets before assigning Tier A.
4. Assign `design_type` from paper text/figures (linker vs merge).
5. Record affinities if reported (nM / assay type).
6. Prefer resolution ≤ 3.0 Å for pose benchmarks.

RCSB full-text search hints used for this seed: `dual inhibitor`, `BET HDAC`, `Mcl-1 Bcl-xL`, `PknA PknB`, `LpxA LpxD`.

## Caveats

- PDB titles containing “dual” are **noisy** (dual-specificity phosphatases, dual binding modes of single-target ligands, etc.). Manual curation is mandatory.
- Many dual inhibitors **lack** both-end crystals; do not inflate Tier A by pairing unrelated ligands.
- BET–HDAC duals often crystallize only with BRD4; HDAC end remains a docking hypothesis unless a dual-ligand–HDAC structure exists.
- This seed set is for **protocol development**, not a claim of exhaustive PDB coverage.

## Version

- **v0.1** (2026-07-19): initial curated seed for dual-target docking method evaluation.
