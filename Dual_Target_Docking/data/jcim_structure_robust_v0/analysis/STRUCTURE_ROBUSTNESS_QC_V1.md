# Structure robustness cognate QC v1

P0 alternative-receptor gate for the PIK3CA/mTOR PM study. Cognate redocking QC for 4JPS / 5DXT / 4JSX; panel-level redock results are in `STRUCTURE_ROBUSTNESS_VERDICT_V1.md`.

## Frozen protocol

- AutoDock Vina 1.2.7; Meeko 0.7.1; seed `20260727`; exhaustiveness `16`; `n_poses=9`; `cpu=1`.
- Box: cognate-ligand heavy-atom AABB center, `5` Å padding on each side, `20` Å minimum edge (same construction as frozen PM boxes).
- RMSD: heavy atoms only, Meeko SMILES-index mapping with template automorphisms, no superposition. Pass gate: best-of-9 < 2.0 Å.
- Receptor preparation: Meeko `mk_prepare_receptor.py --read_pdb -p -a --default_altloc A`; ligand PDBQT via Meeko.

## Results

| target | PDB | cognate ligand | RMSD mode 1 (Å) | best-of-9 RMSD (Å) | best mode | verdict |
|---|---|---|---:|---:|---:|---|
| PIK3CA | 4JPS | 1LT | 0.607 | 0.607 | 1 | PASS |
| PIK3CA | 5DXT | 5H5 | 0.624 | 0.624 | 1 | PASS |
| mTOR | 4JSX | 17G | 0.515 | 0.515 | 1 | PASS |

## Polymer-entity identity audit

### 4JPS — PIK3CA
- entity 1: Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha isoform
- entity 2: Phosphatidylinositol 3-kinase regulatory subunit alpha
- Selected cognate: `1LT` chain `A`.
- PASS: receptor PDBQT/box in `receptors/4JPS_*`.

### 5DXT — PIK3CA
- entity 1: Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha isoform
- Selected cognate: `5H5` chain `A`.
- PASS: receptor PDBQT/box in `receptors/5DXT_*`.

### 4JSX — mTOR
- entity 1: Serine/threonine-protein kinase mTOR
- entity 2: Target of rapamycin complex subunit LST8
- Selected cognate: `17G` chain `B` (Torin2).
- PASS: receptor PDBQT/box in `receptors/4JSX_*`.

## Interpretation

A PASS means the frozen docking prep + cognate-centered box recovers the deposited ligand. It does not by itself prove panel-level structural robustness; that requires the separate PM48 pocket-swap redocks summarized in `STRUCTURE_ROBUSTNESS_VERDICT_V1.md`.

Note: an interrupted `--pdb 4JSX` refresh briefly overwrote this report to show only 4JSX; this file was restored to the full three-structure table. Underlying cognate outputs remain under `cognate_qc/{4JPS,5DXT,4JSX}/`.
