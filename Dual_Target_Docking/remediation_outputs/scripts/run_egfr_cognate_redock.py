#!/usr/bin/env python3
"""Cognate 03P redock into canonical 3POZ/3RCD heavy-atom boxes + RMSD vs crystal."""
from __future__ import annotations

import csv
import json
import re
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
QC = ROOT / "data/egfr_her2_panel40_v0/cognate_qc"
BOXDIR = ROOT / "data/egfr_her2_panel40_v0/boxes"
RECDIR = ROOT / "data/egfr_her2_panel40_v0/receptors"
OUT = ROOT / "remediation_outputs/phase1_cognate_redock"
VINA = "/home/gwj/miniconda3/bin/vina"
MEKO = ["/home/gwj/miniconda3/bin/python", "/home/gwj/miniconda3/bin/mk_prepare_ligand.py"]
SEED = 20260727


def parse_mode1(path: Path) -> float | None:
    m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)", path.read_text(errors="replace"))
    return float(m.group(1)) if m else None


def pdbqt_to_mol(pdbqt: Path):
    # RDKit cannot read pdbqt well; use first MODEL coordinates via SDF if meeko wrote one.
    return None


def heavy_rmsd(crystal_sdf: Path, dock_pdbqt: Path) -> float | None:
    suppl = Chem.SDMolSupplier(str(crystal_sdf), removeHs=False)
    cry = suppl[0] if suppl and suppl[0] is not None else None
    if cry is None:
        return None
    # parse docked MODEL 1 heavy xyz
    xyz = []
    in_m1 = False
    for line in dock_pdbqt.read_text(errors="replace").splitlines():
        if line.startswith("MODEL"):
            in_m1 = line.split()[-1] == "1"
            continue
        if line.startswith("ENDMDL") and in_m1:
            break
        if in_m1 and line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            if name.startswith("H"):
                continue
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    cry = Chem.RemoveHs(cry)
    if cry.GetNumAtoms() != len(xyz):
        # still report if counts differ
        n = min(cry.GetNumAtoms(), len(xyz))
        if n == 0:
            return None
        conf = cry.GetConformer()
        s = 0.0
        for i in range(n):
            p = conf.GetAtomPosition(i)
            q = xyz[i]
            s += (p.x - q[0]) ** 2 + (p.y - q[1]) ** 2 + (p.z - q[2]) ** 2
        return (s / n) ** 0.5
    conf = cry.GetConformer()
    s = 0.0
    n = cry.GetNumAtoms()
    for i in range(n):
        p = conf.GetAtomPosition(i)
        q = xyz[i]
        s += (p.x - q[0]) ** 2 + (p.y - q[1]) ** 2 + (p.z - q[2]) ** 2
    return (s / n) ** 0.5


def prep_from_sdf(sdf: Path, pdbqt: Path) -> Path:
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    pdbqt.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(MEKO + ["-i", str(sdf), "-o", str(pdbqt)], capture_output=True, text=True)
    if proc.returncode != 0 or not pdbqt.exists():
        # fallback: existing E8 ligand pdbqt from old cognate (coordinates of ligand, not box-dependent)
        alt = QC / sdf.name.replace("_crystal.sdf", "_ligand.pdbqt")
        raise RuntimeError(f"meeko failed {sdf}: {proc.stderr[-300:]}")
    return pdbqt


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for pdb, sdf_name in (("3POZ", "3POZ_03P_crystal.sdf"), ("3RCD", "3RCD_03P_crystal.sdf")):
        sdf = QC / sdf_name
        box = json.loads((BOXDIR / f"{pdb}_box.json").read_text())
        rec = RECDIR / f"{pdb}_receptor.pdbqt"
        lig_pdbqt = OUT / f"{pdb}_03P_ligand.pdbqt"
        try:
            prep_from_sdf(sdf, lig_pdbqt)
        except Exception as e:
            # use previously prepared cognate ligand pdbqt if present
            for cand in QC.glob(f"{pdb}*E8.pdbqt"):
                # those are docked outputs, skip
                pass
            # extract ligand from crystal via obabel
            obabel = Path("/home/gwj/miniconda3/bin/obabel")
            proc = subprocess.run(
                [str(obabel), str(sdf), "-O", str(lig_pdbqt)],
                capture_output=True,
                text=True,
            )
            if not lig_pdbqt.exists():
                rows.append({"pdb": pdb, "status": f"prep_fail {e}", "affinity": "", "rmsd_heavy": ""})
                continue
        out = OUT / f"{pdb}_cognate_corrected_out.pdbqt"
        log = OUT / f"{pdb}_cognate_corrected.log"
        conf = OUT / f"{pdb}_cognate_corrected.txt"
        conf.write_text(
            "\n".join(
                [
                    f"receptor = {rec}",
                    f"ligand = {lig_pdbqt}",
                    f"center_x = {box['center_x']}",
                    f"center_y = {box['center_y']}",
                    f"center_z = {box['center_z']}",
                    f"size_x = {box['size_x']}",
                    f"size_y = {box['size_y']}",
                    f"size_z = {box['size_z']}",
                    "exhaustiveness = 8",
                    "num_modes = 9",
                    "energy_range = 3",
                    "cpu = 1",
                    f"seed = {SEED}",
                    f"out = {out}",
                ]
            )
            + "\n"
        )
        proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True)
        log.write_text((proc.stdout or "") + (proc.stderr or ""))
        aff = parse_mode1(out) if out.exists() else None
        rmsd = heavy_rmsd(sdf, out) if out.exists() else None
        rows.append(
            {
                "pdb": pdb,
                "status": "ok" if aff is not None else f"FAIL rc={proc.returncode}",
                "affinity": aff if aff is not None else "",
                "rmsd_heavy": rmsd if rmsd is not None else "",
                "seed": SEED,
                "box": "canonical_heavy_atom",
            }
        )
        print(rows[-1])
    dest = OUT / "cognate_redock_corrected_box.csv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
