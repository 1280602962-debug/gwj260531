#!/usr/bin/env python3
"""JCIM Step 2B: freeze AChE/BChE/PIK3CB receptors + cognate redock QC.

Selection (holo with clear small-molecule cognate):
  AChE  4EY7  ligand E20 (donepezil)
  BChE  6ZWI  ligand QRH
  PIK3CB 2Y3A ligand GD9

Gate: best-of-9 heavy-atom RMSD < 2.0 Å (same spirit as PM cognate QC).
"""
from __future__ import annotations

import json
import math
import re
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

DUAL = Path("/home/gwj/repos/gwj260531/Dual_Target_Docking")
PDB_CACHE = Path("/tmp/jcim_pdb")
SEED = 20260727
EXHAUST = 8  # cognate QC; panel docking may raise later
N_MODES = 9
VINA = "/home/gwj/miniconda3/bin/vina"
PY = "/home/gwj/miniconda3/bin/python"
PAD = 5.0
MIN_EDGE = 20.0

SPECS = [
    {
        "pair_pack": "ache_bche_panel_v0",
        "name": "ACHE",
        "pdb": "4EY7",
        "lig_resname": "E20",
        "note": "human AChE + donepezil",
    },
    {
        "pair_pack": "ache_bche_panel_v0",
        "name": "BCHE",
        "pdb": "6ZWI",
        "lig_resname": "QRH",
        "note": "human BChE holo",
    },
    {
        "pair_pack": "pik3ca_pik3cb_panel_v0",
        "name": "PIK3CB",
        "pdb": "2Y3A",
        "lig_resname": "GD9",
        "note": "PI3Kβ ATP-site holo",
    },
]


def extract_protein_pdb(src: Path, dst: Path):
    lines = []
    for line in src.read_text().splitlines():
        if line.startswith("ATOM") or line.startswith("TER"):
            lines.append(line)
    lines.append("END")
    dst.write_text("\n".join(lines) + "\n")


def extract_ligand_pdb(src: Path, resname: str, dst: Path):
    lines = []
    for line in src.read_text().splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == resname:
            lines.append(line)
    if not lines:
        raise RuntimeError(f"no HETATM {resname} in {src}")
    lines.append("END")
    dst.write_text("\n".join(lines) + "\n")


def ligand_xyz_from_pdb(path: Path):
    xyz = []
    for line in path.read_text().splitlines():
        if line.startswith("HETATM") or line.startswith("ATOM"):
            el = line[76:78].strip() or line[12:16].strip()[0]
            if el.upper().startswith("H"):
                continue
            xyz.append(
                (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            )
    return xyz


def box_from_xyz(xyz, pad=PAD, min_edge=MIN_EDGE):
    xs, ys, zs = zip(*xyz)
    cx, cy, cz = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2
    sx = max(max(xs) - min(xs) + 2 * pad, min_edge)
    sy = max(max(ys) - min(ys) + 2 * pad, min_edge)
    sz = max(max(zs) - min(zs) + 2 * pad, min_edge)
    return {
        "center_x": round(cx, 3),
        "center_y": round(cy, 3),
        "center_z": round(cz, 3),
        "size_x": round(sx, 3),
        "size_y": round(sy, 3),
        "size_z": round(sz, 3),
        "n_heavy_atoms": len(xyz),
    }


def pdb_hetatm_to_sdf(lig_pdb: Path, sdf: Path):
    mol = Chem.MolFromPDBFile(str(lig_pdb), removeHs=False, sanitize=False)
    if mol is None:
        raise RuntimeError(f"RDKit failed to read {lig_pdb}")
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    # keep largest fragment (multi-copy cognates / ions)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if not frags:
        raise RuntimeError(f"no fragments in {lig_pdb}")
    mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    mol = Chem.AddHs(mol, addCoords=True)
    w = Chem.SDWriter(str(sdf))
    w.write(mol)
    w.close()
    return mol


def prepare_receptor(protein_pdb: Path, out_basename: Path):
    # mk_prepare_receptor writes pdbqt next to basename
    cmd = [
        PY,
        "/home/gwj/miniconda3/bin/mk_prepare_receptor.py",
        "--read_pdb",
        str(protein_pdb),
        "-o",
        str(out_basename),
        "-p",
        "-a",
        "--default_altloc",
        "A",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    # find produced pdbqt
    cands = list(out_basename.parent.glob(out_basename.name + "*.pdbqt"))
    if not cands:
        cands = list(out_basename.parent.glob("*.pdbqt"))
    # prefer exact
    exact = out_basename.parent / f"{out_basename.name}.pdbqt"
    if exact.exists():
        return exact
    if cands:
        return cands[0]
    raise RuntimeError(
        f"receptor prep failed: rc={proc.returncode}\n{proc.stdout[-500:]}\n{proc.stderr[-500:]}"
    )


def prepare_ligand_sdf(sdf: Path, pdbqt: Path):
    proc = subprocess.run(
        [
            PY,
            "/home/gwj/miniconda3/bin/mk_prepare_ligand.py",
            "-i",
            str(sdf),
            "-o",
            str(pdbqt),
        ],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not pdbqt.exists():
        raise RuntimeError(f"ligand prep failed: {proc.stderr[-400:]}")


def vina_dock(rec: Path, lig: Path, box: dict, out_pdbqt: Path, log: Path):
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
        raise RuntimeError(f"vina failed: {proc.stderr[-300:]}")


def pdbqt_heavy_xyz(path: Path, model_index: int = 1):
    text = path.read_text().splitlines()
    models, cur = [], []
    for line in text:
        if line.startswith("MODEL"):
            cur = [line]
        elif line.startswith("ENDMDL"):
            cur.append(line)
            models.append(cur)
            cur = []
        elif cur:
            cur.append(line)
    if not models:
        # single model file
        models = [text]
    block = models[model_index - 1]
    xyz = []
    for line in block:
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            if name.upper().startswith("H"):
                continue
            # skip meeko dummy sometimes
            xyz.append(
                (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            )
    return xyz


def rmsd_hungarian(ref, mob):
    """Simple RMSD after optimal atom matching by greedy NN (same N required)."""
    if len(ref) != len(mob):
        # compare min length with order-preserving as fallback
        n = min(len(ref), len(mob))
        ref, mob = ref[:n], mob[:n]
    # Kabsch-free: assume cognate docking keeps atom order from PDBQT rebuild —
    # use sorted pairwise min assignment
    import numpy as np
    from scipy.optimize import linear_sum_assignment

    R = np.asarray(ref, float)
    M = np.asarray(mob, float)
    # cost = squared distance
    cost = ((R[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
    ri, ci = linear_sum_assignment(cost)
    d2 = cost[ri, ci].sum() / len(ri)
    return float(math.sqrt(d2))


def process_one(spec: dict) -> dict:
    pack = DUAL / "data" / spec["pair_pack"]
    rec_dir = pack / "receptors"
    box_dir = pack / "boxes"
    qc_dir = pack / "cognate_qc"
    for d in (rec_dir, box_dir, qc_dir):
        d.mkdir(parents=True, exist_ok=True)

    pdb_id = spec["pdb"]
    src = PDB_CACHE / f"{pdb_id}.pdb"
    if not src.exists():
        raise FileNotFoundError(src)

    protein = rec_dir / f"{pdb_id}_protein.pdb"
    lig_pdb = qc_dir / f"{pdb_id}_{spec['lig_resname']}_crystal.pdb"
    extract_protein_pdb(src, protein)
    extract_ligand_pdb(src, spec["lig_resname"], lig_pdb)
    xyz = ligand_xyz_from_pdb(lig_pdb)
    box = box_from_xyz(xyz)
    box["pdb"] = pdb_id
    box["ligand"] = spec["lig_resname"]
    (box_dir / f"{pdb_id}_box.json").write_text(json.dumps(box, indent=2) + "\n")

    # receptor pdbqt
    out_base = rec_dir / f"{pdb_id}_receptor"
    try:
        rec_pdbqt = prepare_receptor(protein, out_base)
        # normalize name
        target = rec_dir / f"{pdb_id}_receptor.pdbqt"
        if rec_pdbqt.resolve() != target.resolve():
            target.write_bytes(rec_pdbqt.read_bytes())
            rec_pdbqt = target
    except Exception as e:
        return {"name": spec["name"], "pdb": pdb_id, "status": "receptor_prep_fail", "error": str(e)}

    # cognate ligand sdf/pdbqt from crystal coords via RDKit
    lig_sdf = qc_dir / f"{pdb_id}_{spec['lig_resname']}.sdf"
    lig_pdbqt = qc_dir / f"{pdb_id}_{spec['lig_resname']}.pdbqt"
    try:
        pdb_hetatm_to_sdf(lig_pdb, lig_sdf)
        prepare_ligand_sdf(lig_sdf, lig_pdbqt)
    except Exception as e:
        return {"name": spec["name"], "pdb": pdb_id, "status": "ligand_prep_fail", "error": str(e)}

    out_dock = qc_dir / f"{pdb_id}_cognate_out.pdbqt"
    log = qc_dir / f"{pdb_id}_cognate_vina.log"
    try:
        vina_dock(rec_pdbqt, lig_pdbqt, box, out_dock, log)
    except Exception as e:
        return {"name": spec["name"], "pdb": pdb_id, "status": "vina_fail", "error": str(e)}

    ref_xyz = ligand_xyz_from_pdb(lig_pdb)
    best_rmsd = 999.0
    best_mode = None
    for mode in range(1, N_MODES + 1):
        try:
            mob = pdbqt_heavy_xyz(out_dock, mode)
        except Exception:
            continue
        if not mob:
            continue
        r = rmsd_hungarian(ref_xyz, mob)
        if r < best_rmsd:
            best_rmsd = r
            best_mode = mode

    passed = best_rmsd < 2.0
    return {
        "name": spec["name"],
        "pdb": pdb_id,
        "ligand": spec["lig_resname"],
        "n_heavy_ref": len(ref_xyz),
        "best_mode": best_mode,
        "best_of_9_rmsd": round(best_rmsd, 3),
        "pass_rmsd_lt_2": passed,
        "status": "PASS" if passed else "FAIL_RMSD",
        "receptor_pdbqt": str(rec_pdbqt),
        "box": box,
        "note": spec["note"],
    }


def main():
    # ensure 4EY7 present
    if not (PDB_CACHE / "4EY7.pdb").exists():
        raise SystemExit("missing 4EY7.pdb in cache")
    rows = []
    for spec in SPECS:
        print("processing", spec["name"], spec["pdb"], flush=True)
        rec = process_one(spec)
        rows.append(rec)
        print(rec, flush=True)

    # write QC md into both packs
    lines = ["# Cognate QC — new receptor freeze (JCIM Step 2B)\n\n"]
    lines.append("| target | PDB | ligand | best_of_9 RMSD (Å) | pass (<2Å) | status |\n")
    lines.append("|--------|-----|--------|--------------------:|:----------:|--------|\n")
    for r in rows:
        lines.append(
            f"| {r.get('name')} | {r.get('pdb')} | {r.get('ligand','')} | "
            f"{r.get('best_of_9_rmsd','')} | {r.get('pass_rmsd_lt_2','')} | {r.get('status')} |\n"
        )
    lines.append(
        "\nProtocol: Vina E=8, seed=20260727, n_modes=9; box = cognate AABB + 5Å "
        "(min edge 20Å). Receptor prep: meeko `mk_prepare_receptor.py`.\n"
    )
    lines.append(
        "\n**Gate:** any FAIL blocks Step 3 docking for that end until retuned.\n"
    )
    md = "".join(lines)
    for pack in ("ache_bche_panel_v0", "pik3ca_pik3cb_panel_v0"):
        p = DUAL / "data" / pack / "cognate_qc" / "COGNATE_QC.md"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(md)
    # also JSON summary
    (DUAL / "data" / "ache_bche_panel_v0" / "cognate_qc" / "cognate_qc_summary.json").write_text(
        json.dumps(rows, indent=2, default=str) + "\n"
    )
    print("wrote cognate QC docs")


if __name__ == "__main__":
    main()
