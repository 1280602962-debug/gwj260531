#!/usr/bin/env python3
"""JCIM Step 2B (robust): freeze receptors + cognate QC via Open Babel for ligands."""
from __future__ import annotations

import json
import math
import shutil
import subprocess
from pathlib import Path

from scipy.optimize import linear_sum_assignment
import numpy as np

DUAL = Path("/home/gwj/repos/gwj260531/Dual_Target_Docking")
PDB_CACHE = Path("/tmp/jcim_pdb")
SEED = 20260727
EXHAUST = 8
N_MODES = 9
VINA = "/home/gwj/miniconda3/bin/vina"
PY = "/home/gwj/miniconda3/bin/python"
OBABEL = shutil.which("obabel") or "/home/gwj/miniconda3/envs/cadd_tools/bin/obabel"
PAD = 5.0
MIN_EDGE = 20.0

SPECS = [
    {"pair_pack": "ache_bche_panel_v0", "name": "ACHE", "pdb": "4EY7", "lig_resname": "E20", "note": "human AChE + donepezil"},
    {"pair_pack": "ache_bche_panel_v0", "name": "BCHE", "pdb": "6ZWI", "lig_resname": "QRH", "note": "human BChE holo"},
    {"pair_pack": "pik3ca_pik3cb_panel_v0", "name": "PIK3CB", "pdb": "2Y3A", "lig_resname": "GD9", "note": "PI3Kβ ATP-site holo"},
]


def extract_protein_pdb(src: Path, dst: Path):
    lines = [ln for ln in src.read_text().splitlines() if ln.startswith(("ATOM", "TER"))]
    dst.write_text("\n".join(lines) + "\nEND\n")


def extract_first_ligand_pdb(src: Path, resname: str, dst: Path):
    resid = None
    lines = []
    for line in src.read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        if line[17:20].strip() != resname:
            continue
        rid = line[22:26]
        if resid is None:
            resid = rid
        if rid == resid:
            lines.append(line)
    if not lines:
        raise RuntimeError(f"no {resname} in {src}")
    dst.write_text("\n".join(lines) + "\nEND\n")
    return resid


def ligand_heavy_xyz(path: Path):
    xyz = []
    for line in path.read_text().splitlines():
        if not line.startswith(("HETATM", "ATOM")):
            continue
        name = line[12:16].strip()
        if name.upper().startswith("H"):
            continue
        el = (line[76:78].strip() if len(line) >= 78 else "") or name[0]
        if el.upper() == "H":
            continue
        xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return xyz


def box_from_xyz(xyz):
    xs, ys, zs = zip(*xyz)
    return {
        "center_x": round((min(xs) + max(xs)) / 2, 3),
        "center_y": round((min(ys) + max(ys)) / 2, 3),
        "center_z": round((min(zs) + max(zs)) / 2, 3),
        "size_x": round(max(max(xs) - min(xs) + 2 * PAD, MIN_EDGE), 3),
        "size_y": round(max(max(ys) - min(ys) + 2 * PAD, MIN_EDGE), 3),
        "size_z": round(max(max(zs) - min(zs) + 2 * PAD, MIN_EDGE), 3),
        "n_heavy_atoms": len(xyz),
    }


def prepare_receptor(protein_pdb: Path, out_pdbqt: Path):
    base = out_pdbqt.with_suffix("")
    cmd = [
        PY,
        "/home/gwj/miniconda3/bin/mk_prepare_receptor.py",
        "--read_pdb",
        str(protein_pdb),
        "-o",
        str(base),
        "-p",
        "-a",
        "--default_altloc",
        "A",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    produced = base.with_suffix(".pdbqt")
    if not produced.exists():
        cands = list(base.parent.glob(base.name + "*.pdbqt"))
        if not cands:
            raise RuntimeError(proc.stderr[-500:] or proc.stdout[-500:])
        produced = cands[0]
    if produced.resolve() != out_pdbqt.resolve():
        out_pdbqt.write_bytes(produced.read_bytes())
    return out_pdbqt


def prepare_ligand_obabel(lig_pdb: Path, lig_pdbqt: Path):
    sdf = lig_pdb.with_suffix(".sdf")
    p1 = subprocess.run([OBABEL, str(lig_pdb), "-O", str(sdf), "-h"], capture_output=True, text=True)
    if p1.returncode != 0 or not sdf.exists():
        raise RuntimeError(f"obabel pdb->sdf failed: {p1.stderr[-300:]}")
    # Prefer meeko if single fragment; else obabel pdbqt
    p2 = subprocess.run(
        [PY, "/home/gwj/miniconda3/bin/mk_prepare_ligand.py", "-i", str(sdf), "-o", str(lig_pdbqt)],
        capture_output=True,
        text=True,
    )
    if p2.returncode == 0 and lig_pdbqt.exists():
        return "meeko"
    p3 = subprocess.run([OBABEL, str(sdf), "-O", str(lig_pdbqt)], capture_output=True, text=True)
    if p3.returncode != 0 or not lig_pdbqt.exists():
        raise RuntimeError(f"ligand prep failed meeko={p2.stderr[-200:]} obabel={p3.stderr[-200:]}")
    return "obabel"


def vina_dock(rec, lig, box, out_pdbqt, log):
    conf = out_pdbqt.with_suffix(".txt")
    conf.write_text(
        "\n".join(
            [
                f"receptor = {rec}",
                f"ligand = {lig}",
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                f"exhaustiveness = {EXHAUST}",
                f"num_modes = {N_MODES}",
                "energy_range = 3",
                "cpu = 1",
                f"seed = {SEED}",
                f"out = {out_pdbqt}",
            ]
        )
        + "\n"
    )
    proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True)
    log.write_text(proc.stdout + "\n" + proc.stderr)
    if proc.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(proc.stderr[-300:] or proc.stdout[-300:])


def pdbqt_models_xyz(path: Path):
    text = path.read_text().splitlines()
    models, cur = [], []
    for line in text:
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            models.append(cur)
            cur = []
        elif cur is not None and line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            if name.upper().startswith("H"):
                continue
            cur.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not models:
        xyz = []
        for line in text:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16].strip()
                if name.upper().startswith("H"):
                    continue
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        models = [xyz]
    return models


def rmsd_assign(ref, mob):
    R = np.asarray(ref, float)
    M = np.asarray(mob, float)
    n = min(len(R), len(M))
    R, M = R[:n], M[:n]
    cost = ((R[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
    ri, ci = linear_sum_assignment(cost)
    return float(math.sqrt(cost[ri, ci].sum() / len(ri)))


def process_one(spec):
    pack = DUAL / "data" / spec["pair_pack"]
    rec_dir = pack / "receptors"
    box_dir = pack / "boxes"
    qc_dir = pack / "cognate_qc"
    for d in (rec_dir, box_dir, qc_dir):
        d.mkdir(parents=True, exist_ok=True)
    pdb_id = spec["pdb"]
    src = PDB_CACHE / f"{pdb_id}.pdb"
    protein = rec_dir / f"{pdb_id}_protein.pdb"
    lig_pdb = qc_dir / f"{pdb_id}_{spec['lig_resname']}_crystal.pdb"
    extract_protein_pdb(src, protein)
    resid = extract_first_ligand_pdb(src, spec["lig_resname"], lig_pdb)
    xyz = ligand_heavy_xyz(lig_pdb)
    box = box_from_xyz(xyz)
    box.update({"pdb": pdb_id, "ligand": spec["lig_resname"], "resid": resid})
    (box_dir / f"{pdb_id}_box.json").write_text(json.dumps(box, indent=2) + "\n")

    rec_pdbqt = rec_dir / f"{pdb_id}_receptor.pdbqt"
    prepare_receptor(protein, rec_pdbqt)
    lig_pdbqt = qc_dir / f"{pdb_id}_{spec['lig_resname']}.pdbqt"
    how = prepare_ligand_obabel(lig_pdb, lig_pdbqt)
    out_dock = qc_dir / f"{pdb_id}_cognate_out.pdbqt"
    log = qc_dir / f"{pdb_id}_cognate_vina.log"
    vina_dock(rec_pdbqt, lig_pdbqt, box, out_dock, log)
    models = pdbqt_models_xyz(out_dock)
    best_rmsd, best_mode = 999.0, None
    for i, mob in enumerate(models, 1):
        if not mob:
            continue
        r = rmsd_assign(xyz, mob)
        if r < best_rmsd:
            best_rmsd, best_mode = r, i
    passed = best_rmsd < 2.0
    return {
        "name": spec["name"],
        "pdb": pdb_id,
        "ligand": spec["lig_resname"],
        "ligand_prep": how,
        "n_heavy_ref": len(xyz),
        "best_mode": best_mode,
        "best_of_9_rmsd": round(best_rmsd, 3),
        "pass_rmsd_lt_2": passed,
        "status": "PASS" if passed else "FAIL_RMSD",
        "note": spec["note"],
    }


def main():
    rows = []
    for spec in SPECS:
        print("processing", spec["name"], spec["pdb"], flush=True)
        try:
            rec = process_one(spec)
        except Exception as e:
            rec = {"name": spec["name"], "pdb": spec["pdb"], "status": "ERROR", "error": str(e)}
        rows.append(rec)
        print(rec, flush=True)
    lines = ["# Cognate QC — new receptor freeze (JCIM Step 2B)\n\n"]
    lines.append("| target | PDB | ligand | best_of_9 RMSD (Å) | pass | status |\n")
    lines.append("|--------|-----|--------|--------------------:|:----:|--------|\n")
    for r in rows:
        lines.append(
            f"| {r.get('name')} | {r.get('pdb')} | {r.get('ligand','')} | "
            f"{r.get('best_of_9_rmsd','')} | {r.get('pass_rmsd_lt_2','')} | {r.get('status')} |\n"
        )
    lines.append(
        "\nProtocol: Vina E=8 seed=20260727 n_modes=9; box=cognate AABB+5Å (min 20).\n"
        "Receptor: meeko mk_prepare_receptor (-a, default_altloc A).\n"
        "Cognate ligand: first residue copy → Open Babel (+meeko if possible).\n"
        "\n**Gate:** FAIL blocks Step-3 for that end until retuned.\n"
        "\nPIK3CA end reuses frozen **4L23** from PM48 (copied into pik3ca_pik3cb pack).\n"
    )
    md = "".join(lines)
    for pack in ("ache_bche_panel_v0", "pik3ca_pik3cb_panel_v0"):
        p = DUAL / "data" / pack / "cognate_qc" / "COGNATE_QC.md"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(md)
    (DUAL / "data/ache_bche_panel_v0/cognate_qc/cognate_qc_summary.json").write_text(
        json.dumps(rows, indent=2) + "\n"
    )
    print("done")


if __name__ == "__main__":
    main()
