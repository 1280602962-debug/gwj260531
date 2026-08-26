#!/usr/bin/env python3
"""Re-audit ranked cognate RMSD for receptors with deposited SDF/PDBQT pairs.

Meeko reconstructs each multi-pose PDBQT as an RDKit molecule.  The reference
and pose heavy-atom graphs must have identical canonical isomeric SMILES before
RDKit CalcRMS is used.  CalcRMS enumerates symmetry-equivalent atom mappings and
does not superpose the docking coordinates.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"

SPECS = (
    (
        "AChE", "4EY7", 8,
        "data/ache_bche_panel_v0/cognate_qc/4EY7_E20_crystal.sdf",
        "data/ache_bche_panel_v0/cognate_qc/4EY7_cognate_out.pdbqt",
        None,
    ),
    (
        "BChE", "4BDS", 8,
        "data/ache_bche_panel_v0/cognate_qc/4BDS_THA_crystal.sdf",
        "data/ache_bche_panel_v0/cognate_qc/4BDS_cognate_out_E8.pdbqt",
        "data/ache_bche_panel_v0/cognate_qc/4BDS_THA_E8.pdbqt",
    ),
    (
        "PIK3CB", "2WXF", 8,
        "data/pik3ca_pik3cb_panel_v0/cognate_qc/2WXF_039_crystal.sdf",
        "data/pik3ca_pik3cb_panel_v0/cognate_qc/2WXF_cognate_out_E8.pdbqt",
        "data/pik3ca_pik3cb_panel_v0/cognate_qc/2WXF_039_E8.pdbqt",
    ),
)


def read_reference(path: Path):
    mol = Chem.SDMolSupplier(str(path), removeHs=False)[0]
    if mol is None:
        raise ValueError(f"could not parse reference SDF: {path}")
    return Chem.RemoveHs(mol)


def read_poses_meeko(path: Path):
    pdbqt = PDBQTMolecule.from_file(str(path), skip_typing=True)
    made = RDKitMolCreate.from_pdbqt_mol(pdbqt)
    if len(made) != 1:
        raise ValueError(f"expected one reconstructed molecule, found {len(made)}: {path}")
    return Chem.RemoveHs(made[0])


def pdbqt_element(line: str) -> str:
    atom_type = line.split()[-1]
    if atom_type.upper().startswith("H"):
        return "H"
    if atom_type in {"A", "C"}:
        return "C"
    if atom_type in {"N", "NA"}:
        return "N"
    if atom_type in {"O", "OA"}:
        return "O"
    if atom_type in {"S", "SA", "S6"}:
        return "S"
    if atom_type.lower() == "cl":
        return "Cl"
    if atom_type.lower() == "br":
        return "Br"
    return atom_type.title()


def pdbqt_models(path: Path):
    models, current = [], []
    in_model = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL"):
            current, in_model = [], True
        elif line.startswith("ENDMDL"):
            models.append(current)
            current, in_model = [], False
        elif line.startswith(("ATOM", "HETATM")):
            element = pdbqt_element(line)
            if element != "H":
                current.append(
                    (element, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
                )
    if current and not models:
        models = [current]
    return models


def read_poses_legacy(ref, input_path: Path, output_path: Path):
    """Recover legacy PDBQT poses on the reference SDF topology.

    The prepared cognate input is mapped to the reference by element-constrained
    coordinate assignment. Vina preserves PDBQT atom order across output poses.
    """
    prepared = pdbqt_models(input_path)[0]
    output = pdbqt_models(output_path)
    ref_conf = ref.GetConformer()
    ref_elements = [atom.GetSymbol() for atom in ref.GetAtoms()]
    ref_xyz = np.asarray([list(ref_conf.GetAtomPosition(i)) for i in range(ref.GetNumAtoms())])
    prep_elements = [x[0] for x in prepared]
    prep_xyz = np.asarray([x[1] for x in prepared])
    if len(ref_elements) != len(prep_elements):
        raise ValueError(f"heavy-atom count mismatch: reference={len(ref_elements)}, PDBQT={len(prep_elements)}")
    cost = ((ref_xyz[:, None, :] - prep_xyz[None, :, :]) ** 2).sum(axis=2)
    for i, ref_element in enumerate(ref_elements):
        for j, prep_element in enumerate(prep_elements):
            if ref_element != prep_element:
                cost[i, j] = 1e6
    ref_idx, prep_idx = linear_sum_assignment(cost)
    mapped_distances = np.sqrt(cost[ref_idx, prep_idx])
    if float(mapped_distances.max()) > 0.15:
        raise ValueError(f"unsafe reference-to-PDBQT atom map; max distance={mapped_distances.max():.3f} A")
    prep_for_ref = dict(zip(ref_idx, prep_idx))

    poses = Chem.Mol(ref)
    poses.RemoveAllConformers()
    for model in output:
        if len(model) != len(prepared):
            raise ValueError("output pose heavy-atom count differs from prepared input")
        if [x[0] for x in model] != prep_elements:
            raise ValueError("output PDBQT changed heavy-atom order or types")
        conf = Chem.Conformer(ref.GetNumAtoms())
        for ri in range(ref.GetNumAtoms()):
            xyz = model[prep_for_ref[ri]][1]
            conf.SetAtomPosition(ri, xyz)
        poses.AddConformer(conf, assignId=True)
    return poses, float(mapped_distances.max())


def main():
    rows = []
    for target, pdb, exhaustiveness, ref_rel, pose_rel, input_rel in SPECS:
        ref = read_reference(ROOT / ref_rel)
        if input_rel is None:
            poses = read_poses_meeko(ROOT / pose_rel)
            mapping_method = "Meeko topology reconstruction"
            mapping_max_distance = ""
        else:
            poses, mapping_max_distance = read_poses_legacy(
                ref, ROOT / input_rel, ROOT / pose_rel
            )
            mapping_method = "element-constrained crystal-coordinate map to reference SDF"
        ref_smiles = Chem.MolToSmiles(ref, isomericSmiles=True)
        pose_smiles = Chem.MolToSmiles(poses, isomericSmiles=True)
        if ref_smiles != pose_smiles:
            raise ValueError(f"topology mismatch for {pdb}: {ref_smiles} != {pose_smiles}")
        rmsd = [
            float(rdMolAlign.CalcRMS(poses, ref, prbId=i, refId=0))
            for i in range(poses.GetNumConformers())
        ]
        for rank, value in enumerate(rmsd, 1):
            rows.append(
                {
                    "target": target,
                    "pdb": pdb,
                    "exhaustiveness": exhaustiveness,
                    "n_modes_deposited": len(rmsd),
                    "pose_rank": rank,
                    "rmsd_A": round(value, 4),
                    "best_top1_A": round(min(rmsd[:1]), 4),
                    "best_top3_A": round(min(rmsd[:3]), 4),
                    "best_all_deposited_A": round(min(rmsd), 4),
                    "pass_top1_lt2": int(min(rmsd[:1]) < 2.0),
                    "pass_top3_lt2": int(min(rmsd[:3]) < 2.0),
                    "pass_all_deposited_lt2": int(min(rmsd) < 2.0),
                    "mapping_method": mapping_method,
                    "mapping_max_distance_A": (
                        "" if mapping_max_distance == "" else round(mapping_max_distance, 4)
                    ),
                    "rmsd_method": "RDKit symmetry-aware CalcRMS; no superposition",
                }
            )
        print(target, pdb, [round(x, 3) for x in rmsd])

    output = TAB / "cognate_rank_rmsd_reaudit_v1.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(output)


if __name__ == "__main__":
    main()
