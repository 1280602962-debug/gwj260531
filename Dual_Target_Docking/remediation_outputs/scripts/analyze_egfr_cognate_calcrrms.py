#!/usr/bin/env python3
"""Official CalcRMS cognate RMSD on corrected-box 03P redocks."""
from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
SRC = ROOT / "data/jcim_novelty_v0/scripts/cognate_rank_qc_v1.py"
OUT = ROOT / "remediation_outputs/phase1_cognate_redock/cognate_rank_rmsd_corrected_box.csv"


def load():
    spec = importlib.util.spec_from_file_location("cognate_rank_qc_v1", SRC)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    L = load()
    rows = []
    jobs = [
        (
            "EGFR",
            "3POZ",
            ROOT / "data/egfr_her2_panel40_v0/cognate_qc/3POZ_03P_crystal.sdf",
            ROOT / "remediation_outputs/phase1_cognate_redock/3POZ_cognate_corrected_out.pdbqt",
            ROOT / "remediation_outputs/phase1_cognate_redock/3POZ_03P_ligand.pdbqt",
        ),
        (
            "HER2",
            "3RCD",
            ROOT / "data/egfr_her2_panel40_v0/cognate_qc/3RCD_03P_crystal.sdf",
            ROOT / "remediation_outputs/phase1_cognate_redock/3RCD_cognate_corrected_out.pdbqt",
            ROOT / "remediation_outputs/phase1_cognate_redock/3RCD_03P_ligand.pdbqt",
        ),
    ]
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    for target, pdb, cry, pose, inp in jobs:
        ref = L.read_reference(cry)
        poses, map_max = L.read_poses_legacy(ref, inp, pose)
        ref_smiles = Chem.MolToSmiles(ref, isomericSmiles=True)
        pose_smiles = Chem.MolToSmiles(poses, isomericSmiles=True)
        if ref_smiles != pose_smiles:
            raise SystemExit(f"topology mismatch {pdb}: {ref_smiles} != {pose_smiles}")
        rmsd = [float(rdMolAlign.CalcRMS(poses, ref, prbId=i, refId=0)) for i in range(poses.GetNumConformers())]
        print(target, pdb, [round(x, 3) for x in rmsd], "best", round(min(rmsd), 3), "top1", round(rmsd[0], 3))
        for rank, value in enumerate(rmsd, 1):
            rows.append(
                {
                    "target": target,
                    "pdb": pdb,
                    "box": "canonical_heavy_atom",
                    "exhaustiveness": 8,
                    "n_modes_deposited": len(rmsd),
                    "pose_rank": rank,
                    "rmsd_A": round(value, 4),
                    "best_top1_A": round(min(rmsd[:1]), 4),
                    "best_top3_A": round(min(rmsd[:3]), 4),
                    "best_all_deposited_A": round(min(rmsd), 4),
                    "pass_top1_lt2": int(min(rmsd[:1]) < 2.0),
                    "pass_top3_lt2": int(min(rmsd[:3]) < 2.0),
                    "pass_all_deposited_lt2": int(min(rmsd) < 2.0),
                    "mapping_method": "element-constrained crystal-coordinate map to reference SDF",
                    "mapping_max_distance_A": round(map_max, 4),
                    "rmsd_method": "RDKit symmetry-aware CalcRMS; no superposition",
                }
            )
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
