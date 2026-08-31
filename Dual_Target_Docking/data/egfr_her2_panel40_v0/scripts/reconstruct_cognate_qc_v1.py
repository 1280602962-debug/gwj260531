#!/usr/bin/env python3
"""Reconstructed EGFR/HER2 cognate redock QC (local).

Original nine-mode cognate PDBQTs are absent from git and the local results
tree. This script rebuilds them under the frozen protocol and labels outputs
**reconstructed QC** (not historical production artifacts).

Protocol: Vina 1.2.7, seed 20260727, exhaustiveness 8, num_modes 9,
energy_range 3; frozen panel40 receptor PDBQT + boxes.
"""
from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from vina import Vina

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/egfr_her2_panel40_v0/cognate_qc"
LOCAL_RESULTS = Path(
    "/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0"
)
MAESTRO_03P = Path(
    "/mnt/d/CADD paper exercise/dual target docking/Maestro doc/vina_docking/rtmscore_3POZ/03P_crystal.sdf"
)
OBABEL = Path("/mnt/d/CADD paper exercise/gnina/conda_env/bin/obabel")
MEKO_PY = Path("/home/gwj/miniconda3/bin/python")
MK_PREP = Path("/home/gwj/miniconda3/bin/mk_prepare_ligand.py")
SEED = 20260727
EXHAUST = 8
N_MODES = 9

SPECS = (
    {
        "target": "3POZ",
        "ligand_name": "03P",
        "crystal_pdb": LOCAL_RESULTS
        / "analysis/exhaustiveness_sensitivity_v1/tables/3POZ_cocrystal_03P.pdb",
        "receptor": LOCAL_RESULTS / "receptors/3POZ_receptor.pdbqt",
        "box": LOCAL_RESULTS / "boxes/3POZ_box.json",
    },
    {
        "target": "3RCD",
        "ligand_name": "03P",  # TAK-285 residue code in both crystals
        "crystal_pdb": LOCAL_RESULTS
        / "analysis/exhaustiveness_sensitivity_v1/tables/3RCD_cocrystal_03P.pdb",
        "receptor": LOCAL_RESULTS / "receptors/3RCD_receptor.pdbqt",
        "box": LOCAL_RESULTS / "boxes/3RCD_box.json",
    },
)


def pdb_het_xyz(path: Path, resname: str = "03P") -> np.ndarray:
    xyz = []
    for line in path.read_text().splitlines():
        if not line.startswith(("HETATM", "ATOM")):
            continue
        if line[17:20].strip() != resname:
            continue
        name = line[12:16].strip()
        if name.startswith("H"):
            continue
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(xyz, dtype=float)


def template_to_crystal_sdf(crystal_pdb: Path, out_sdf: Path) -> Chem.Mol:
    """Apply Maestro 03P bond orders onto crystal heavy-atom coordinates."""
    tmpl = Chem.RemoveHs(Chem.SDMolSupplier(str(MAESTRO_03P), removeHs=False)[0])
    tmp_sdf = out_sdf.with_suffix(".obabel.sdf")
    subprocess.run(
        [str(OBABEL), str(crystal_pdb), "-O", str(tmp_sdf), "-d"],
        check=True,
        capture_output=True,
    )
    cry = Chem.RemoveHs(Chem.SDMolSupplier(str(tmp_sdf), removeHs=False)[0])
    if cry is None:
        raise RuntimeError(f"obabel failed to parse {crystal_pdb}")
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(tmpl, cry)
    except Exception:
        mol = cry
    Chem.SanitizeMol(mol)
    mol.SetProp("_Name", crystal_pdb.stem)
    mol_h = Chem.AddHs(mol, addCoords=True)
    w = Chem.SDWriter(str(out_sdf))
    w.write(mol_h)
    w.close()
    return Chem.RemoveHs(Chem.Mol(mol))


def prepare_pdbqt(sdf: Path, pdbqt: Path) -> None:
    proc = subprocess.run(
        [str(MEKO_PY), str(MK_PREP), "-i", str(sdf), "-o", str(pdbqt)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not pdbqt.exists():
        raise RuntimeError(f"meeko failed: {proc.stderr[-400:]}")


def dock(receptor: Path, ligand: Path, box: dict, out_pdbqt: Path) -> None:
    v = Vina(sf_name="vina", cpu=4, seed=SEED, verbosity=1)
    v.set_receptor(str(receptor))
    v.set_ligand_from_file(str(ligand))
    v.compute_vina_maps(
        center=[box["center_x"], box["center_y"], box["center_z"]],
        box_size=[box["size_x"], box["size_y"], box["size_z"]],
    )
    v.dock(exhaustiveness=EXHAUST, n_poses=N_MODES)
    v.write_poses(str(out_pdbqt), n_poses=N_MODES, overwrite=True, energy_range=3)


def pdbqt_models_xyz(path: Path) -> list[np.ndarray]:
    models, cur = [], []
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            models.append(np.asarray(cur, float))
            cur = []
        elif line.startswith(("ATOM", "HETATM")):
            atom_type = line.split()[-1]
            if atom_type.upper().startswith("H"):
                continue
            cur.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if cur and not models:
        models = [np.asarray(cur, float)]
    return models


def rmsd_hungary(ref_xyz: np.ndarray, pose_xyz: np.ndarray) -> float:
    if len(ref_xyz) != len(pose_xyz):
        return float("nan")
    # no superposition — docking frame RMSD
    # element-agnostic optimal assignment (symmetry)
    from scipy.optimize import linear_sum_assignment

    cost = ((ref_xyz[:, None, :] - pose_xyz[None, :, :]) ** 2).sum(axis=2)
    r, c = linear_sum_assignment(cost)
    return float(np.sqrt(cost[r, c].mean()))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in SPECS:
        target = spec["target"]
        print(f"=== {target} reconstructed cognate QC ===", flush=True)
        crystal_sdf = OUT / f"{target}_{spec['ligand_name']}_crystal.sdf"
        lig_pdbqt = OUT / f"{target}_{spec['ligand_name']}_E{EXHAUST}.pdbqt"
        out_pdbqt = OUT / f"{target}_cognate_out_E{EXHAUST}.pdbqt"
        log_path = OUT / f"{target}_cognate_vina_E{EXHAUST}.log"

        template_to_crystal_sdf(spec["crystal_pdb"], crystal_sdf)
        prepare_pdbqt(crystal_sdf, lig_pdbqt)
        box = json.loads(spec["box"].read_text())
        try:
            dock(spec["receptor"], lig_pdbqt, box, out_pdbqt)
            status = "success"
            reason = ""
        except Exception as exc:  # noqa: BLE001
            status = "fail"
            reason = str(exc)[:300]
            log_path.write_text(reason)
            rows.append(
                {
                    "pair": "EGFR/HER2",
                    "target": target,
                    "pdb": target,
                    "cognate_ligand": spec["ligand_name"],
                    "protocol": "reconstructed_qc",
                    "exhaustiveness": EXHAUST,
                    "seed": SEED,
                    "status": status,
                    "reason": reason,
                    "rmsd_top1": "",
                    "rmsd_top3_min": "",
                    "rmsd_best9_min": "",
                    "pass_top1_lt2": "",
                    "pass_top3_lt2": "",
                    "n_modes": 0,
                }
            )
            continue

        ref_xyz = pdb_het_xyz(spec["crystal_pdb"])
        poses = pdbqt_models_xyz(out_pdbqt)
        rmsds = [rmsd_hungary(ref_xyz, p) for p in poses]
        top1 = rmsds[0] if rmsds else float("nan")
        top3 = min(rmsds[:3]) if rmsds else float("nan")
        best = min(rmsds) if rmsds else float("nan")
        (OUT / f"{target}_cognate_rmsd_E{EXHAUST}.csv").write_text(
            "mode,rmsd\n"
            + "\n".join(f"{i+1},{r:.4f}" for i, r in enumerate(rmsds))
            + "\n"
        )
        rows.append(
            {
                "pair": "EGFR/HER2",
                "target": target,
                "pdb": target,
                "cognate_ligand": spec["ligand_name"],
                "protocol": "reconstructed_qc",
                "exhaustiveness": EXHAUST,
                "seed": SEED,
                "status": "success",
                "reason": "",
                "rmsd_top1": round(top1, 4),
                "rmsd_top3_min": round(top3, 4),
                "rmsd_best9_min": round(best, 4),
                "pass_top1_lt2": int(top1 < 2.0),
                "pass_top3_lt2": int(top3 < 2.0),
                "n_modes": len(rmsds),
            }
        )
        print(
            f"  {target}: top1={top1:.3f} top3min={top3:.3f} best9={best:.3f}",
            flush=True,
        )

    out_csv = OUT / "cognate_reconstructed_qc_summary_v1.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    verdict = OUT / "EGFR_HER2_COGNATE_RECONSTRUCTED_QC_VERDICT.md"
    lines = [
        "# EGFR/HER2 cognate reconstructed QC verdict",
        "",
        "Original nine-mode cognate PDBQTs were not found on the local calculation",
        "disk or in git. Poses below were **re-redocked** under the frozen protocol",
        f"(Vina, seed {SEED}, exhaustiveness {EXHAUST}, 9 modes) and must be cited",
        "as reconstructed QC, not as the original production gate artifact.",
        "",
        "| target | top-1 RMSD | top-3 min | best-of-9 | pass top-1 (<2Å) | pass top-3 |",
        "|--------|------------:|----------:|----------:|:----------------:|:----------:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['target']} | {r['rmsd_top1']} | {r['rmsd_top3_min']} | "
            f"{r['rmsd_best9_min']} | {r['pass_top1_lt2']} | {r['pass_top3_lt2']} |"
        )
    lines += [
        "",
        "Crystal ligand: residue `03P` (TAK-285) in both 3POZ and 3RCD.",
        f"Artifacts: `{OUT.relative_to(ROOT)}/`.",
        "",
    ]
    verdict.write_text("\n".join(lines))
    print(f"wrote {out_csv}", flush=True)
    print(f"wrote {verdict}", flush=True)


if __name__ == "__main__":
    main()
