#!/usr/bin/env python3
"""Dock CHEMBL3960861 into AChE/BChE with the frozen official protocol.

Protocol (data/ache_bche_panel_v0/scripts/dock_panel.py + protocol.yaml):
  RDKit ETKDG seed 20260727 + MMFF + meeko
  Vina 1.2.7, exhaustiveness 8, num_modes 9, energy_range 3, cpu=1, seed 20260727
  receptors ACHE_receptor.pdbqt / BCHE_receptor.pdbqt
  boxes ACHE_box.json / BCHE_box.json (identical to 4EY7 / 4BDS heavy-atom boxes)

Box is unchanged for this pair. Only the newly selected A_only ligand is docked.
"""
from __future__ import annotations

import csv
import json
import re
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
SRC = ROOT / "data/ache_bche_panel_v0"
OUT = ROOT / "remediation_outputs/phase2_ache_bche/vina_new_ligand"
PANEL = ROOT / "remediation_outputs/phase2_ache_bche/panel_v0_strict_no_id_prefix_cap.csv"
VINA = "/home/gwj/miniconda3/bin/vina"
MEKO_PY = "/home/gwj/miniconda3/bin/python"
MKO = "/home/gwj/miniconda3/bin/mk_prepare_ligand.py"
SEED = 20260727
N_MODES = 9
EXHAUST = 8
NEW_ID = "CHEMBL3960861"


def prep_ligand(panel_id: str, smiles: str) -> Path:
    sdf = OUT / "ligands_sdf" / f"{panel_id}.sdf"
    pdbqt = OUT / "ligands_pdbqt" / f"{panel_id}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"bad smiles {panel_id}")
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        AllChem.EmbedMolecule(mol, randomSeed=SEED)
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    sdf.parent.mkdir(parents=True, exist_ok=True)
    pdbqt.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(sdf))
    mol.SetProp("_Name", panel_id)
    w.write(mol)
    w.close()
    proc = subprocess.run(
        [MEKO_PY, MKO, "-i", str(sdf), "-o", str(pdbqt)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not pdbqt.exists():
        raise RuntimeError(f"meeko {panel_id}: {proc.stderr[-400:]}")
    return pdbqt


def parse_mode1(out_pdbqt: Path) -> float | None:
    text = out_pdbqt.read_text(errors="replace")
    m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)", text)
    return float(m.group(1)) if m else None


def dock(target: str, lig: str, ligand: Path, box: dict) -> dict:
    rec = SRC / "receptors" / f"{target}_receptor.pdbqt"
    conf_dir = OUT / "confs"
    out_dir = OUT / "out"
    log_dir = OUT / "logs"
    for d in (conf_dir, out_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{target}_{lig}_out.pdbqt"
    log = log_dir / f"{target}_{lig}.log"
    conf = conf_dir / f"{target}_{lig}.txt"
    conf.write_text(
        "\n".join(
            [
                f"receptor = {rec}",
                f"ligand = {ligand}",
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
                f"out = {out}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    if out.exists() and out.stat().st_size > 0 and parse_mode1(out) is not None:
        aff = parse_mode1(out)
        return {"target": target, "ligand": lig, "status": "cached", "affinity": aff, "seed": SEED}
    proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True)
    log.write_text((proc.stdout or "") + (proc.stderr or ""), encoding="utf-8")
    if proc.returncode != 0 or not out.exists():
        return {
            "target": target,
            "ligand": lig,
            "status": f"FAIL rc={proc.returncode}",
            "affinity": "",
            "seed": SEED,
        }
    aff = parse_mode1(out)
    return {
        "target": target,
        "ligand": lig,
        "status": "ok" if aff is not None else "NO_MODE1",
        "affinity": aff if aff is not None else "",
        "seed": SEED,
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    row = None
    with PANEL.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["molecule_chembl_id"] == NEW_ID:
                row = r
                break
    if row is None:
        raise SystemExit(f"{NEW_ID} not in reconstructed panel")
    lig = row["panel_id"]
    smiles = row["smiles"]
    ligand = prep_ligand(lig, smiles)
    print(f"prepared {lig} {ligand} torsdof check")
    boxes = {
        "ACHE": json.loads((SRC / "boxes/ACHE_box.json").read_text()),
        "BCHE": json.loads((SRC / "boxes/BCHE_box.json").read_text()),
    }
    results = [dock(t, lig, ligand, boxes[t]) for t in ("ACHE", "BCHE")]
    dest = OUT / "scores_vina_mode1_new_ligand.csv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()) + ["molecule_chembl_id", "smiles"])
        w.writeheader()
        for r in results:
            r = dict(r)
            r["molecule_chembl_id"] = NEW_ID
            r["smiles"] = smiles
            w.writerow(r)
    print("wrote", dest)
    for r in results:
        print(r)
    return 0 if all(r["status"] in {"ok", "cached"} for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
