# PDB structures for STAD-AIDD / MASFL ensemble docking

See `docking_ensemble_pdb.csv` and `config/docking_ensemble.yaml` for the authoritative mapping.

## URAT1 three-state docking (Teacher M-CPDL / S_trap)

**One PDB per conformational state** — Schrödinger/Maestro does **not** auto-split states from a single file.

| State | Primary PDB | Ligand in structure | EMDB | Use |
|-------|-------------|---------------------|------|-----|
| **Inward-open** | [9DKB](https://www.rcsb.org/structure/9DKB) | lesinurad | EMD-46950 | `grid_inward` — inhibitor docking, redock gate |
| **Occluded** | [9B1K](https://www.rcsb.org/structure/9B1K) | urate | EMD-44082 | `grid_occluded` |
| **Outward-open** | [9B1L](https://www.rcsb.org/structure/9B1L) | urate | EMD-44083 | `grid_outward` |

Source: Dai & Lee, *Cell Res* 2024 (transport cycle); inward inhibitor grid: Fedor/Suo *Nat Commun* 2025 (9DKB, highest resolution lesinurad).

### Common mistake

| PDB | RCSB title | **Not** |
|-----|------------|---------|
| 9JDZ | Human URAT1 bound to **lesinurad** | occluded or outward — it is inward-open only (Wu *Cell Discov* 2025) |

Wu et al. solved urate outward/occluded states cryo-EM, but those conformations are **not** deposited as separate PDB entries alongside 9JDZ. Use **9B1K / 9B1L** from Dai 2024.

### Preparation workflow (per PDB)

1. Download from RCSB; remove co-crystallized ligand for docking (retain for grid center)
2. Protein Preparation Wizard: pH 7.4, fill loops, optimize H-bond network
3. Receptor Grid Generation: center on co-crystal ligand; box 22×22×22 Å
4. Glide SP → XP on `distill_manifest.csv` (or screening library) against all three grids

### Structural sanity check (PyMOL)

Align 9B1L, 9B1K, 9DKB on NTD (Dai 2024 method). Outward → occluded → inward should show progressive CTD (TM7–TM12) closure of the extracellular cavity.

## URAT1 supplementary (benchmark / redock)

| PDB | Description |
|-----|-------------|
| [9B1H](https://www.rcsb.org/structure/9B1H) | lesinurad inward (Dai 2024) |
| [9DKA](https://www.rcsb.org/structure/9DKA) | benzbromarone inward — redock only |
| [9JDY](https://www.rcsb.org/structure/9JDY) | verinurad (Wu 2025) |
| [9JE1](https://www.rcsb.org/structure/9JE1) | dotinurad (Wu 2025) |

## NLRP3 NACHT domain

| PDB | Description | DOI |
|-----|-------------|-----|
| [7ALV](https://www.rcsb.org/structure/7ALV) | NACHT + MCC950-class analog NP3-146, 2.84 Å X-ray | 10.2210/pdb7ALV/pdb |
| [8ETR](https://www.rcsb.org/structure/8ETR) | NACHT + GDC-2394, cryo-EM | 10.2210/pdb8ETR/pdb |

## Key references

- Dai Y, Lee CH. *Cell Res* 2024 — URAT1 transport cycle (9B1H–9B1O, **9B1K** occluded, **9B1L** outward)
- Suo Y et al. *Nat Commun* 2025 — **9DKB** lesinurad inward
- Wu C et al. *Cell Discov* 2025 — native URAT1 drug structures (9JDZ inward lesinurad)
- Dekker A et al. *J Mol Biol* 2021 — NLRP3 7ALV
- McBride CJ et al. *J Med Chem* 2022 — GDC-2394 / 8ETR

Full methodology: `docs/URAT1_THREE_STATE_DOCKING.md`
