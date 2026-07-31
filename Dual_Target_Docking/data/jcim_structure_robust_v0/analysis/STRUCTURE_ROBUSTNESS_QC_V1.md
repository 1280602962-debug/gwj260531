# Structure robustness cognate QC v1

P0 alternative-receptor gate for the PIK3CA/mTOR PM study. This is cognate redocking QC only; no frozen-panel redocking was started.

## Frozen protocol

- AutoDock Vina 1.2.7; Meeko 0.7.1; seed `20260727`; exhaustiveness `16`; `n_poses=9`; `cpu=1`.
- Box: cognate-ligand heavy-atom AABB center, `5` Å padding on each side, `20` Å minimum edge (the same construction used for the frozen PM boxes).
- RMSD: heavy atoms only, Meeko SMILES-index mapping with all template automorphisms, no superposition. Pass gate: best-of-9 < 2.0 Å.
- Receptor preparation: Meeko `mk_prepare_receptor.py --read_pdb -p -a --default_altloc A`; ligand PDBQT: Meeko. Open Babel only converts the crystal PDB residue to SDF to retain its bond graph.

## Results

| target | PDB | cognate ligand | RMSD mode 1 (Å) | best-of-9 RMSD (Å) | best mode | verdict |
|---|---|---|---:|---:|---:|---|
| mTOR | 4JSX | 17G | 0.515 | 0.515 | 1 | PASS |

## Polymer-entity identity audit

### 4JSX — mTOR
- entity 1: Serine/threonine-protein kinase mTOR
- entity 2: Target of rapamycin complex subunit LST8
- Selected cognate: `17G` chain `B`, residue `2601_`. Selection was the largest non-solvent organic HET residue; full candidate audit is in the JSON record.
- PASS: promoted receptor PDBQT and box to `data/jcim_structure_robust_v0/receptors/4JSX_*`.

## Interpretation

A PASS establishes that the frozen docking preparation and cognate-centered box can recover the deposited ligand for that alternative crystal structure. It does not itself establish panel-level structural robustness; only PASS receptors may be used in the later, separate frozen-panel redocking step.
