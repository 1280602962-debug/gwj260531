# PDB structures for STAD-AIDD ensemble docking

See `docking_ensemble_pdb.csv` for the full table.

## URAT1 (SLC22A12)

Download from RCSB PDB:

| PDB | Description | DOI |
|-----|-------------|-----|
| [9B1H](https://www.rcsb.org/structure/9B1H) | URAT1 + lesinurad, inward-open | 10.2210/pdb9B1H/pdb |
| [9DKB](https://www.rcsb.org/structure/9DKB) | URAT1 + lesinurad, 2.55 Å | 10.2210/pdb9DKB/pdb |
| [9JDZ](https://www.rcsb.org/structure/9JDZ) | Native URAT1, multiple states | 10.2210/pdb9JDZ/pdb |

**Preparation workflow**:
1. Remove co-crystallized ligand, waters, ions (retain structural ions if needed)
2. Protonation at pH 7.4 (PDB2PQR or Schrödinger Protein Preparation Wizard)
3. For MD: embed in POPC membrane via CHARMM-GUI

## NLRP3 NACHT domain

| PDB | Description | DOI |
|-----|-------------|-----|
| [7ALV](https://www.rcsb.org/structure/7ALV) | NACHT + MCC950-class analog NP3-146, 2.84 Å X-ray | 10.2210/pdb7ALV/pdb |
| [8ETR](https://www.rcsb.org/structure/8ETR) | NACHT + GDC-2394, cryo-EM | 10.2210/pdb8ETR/pdb |

## Key references

- Dai et al., *Cell Res* 2024 — URAT1 transport cycle (PDB 9B1H 等)
- Dekker A et al., *J Mol Biol* 2021 — NLRP3 NACHT + inhibitor analog（PDB 7ALV；配体为 NP3-146，非 MCC950）
- McBride CJ et al., *J Med Chem* 2022 — GDC-2394（PDB 8ETR）
