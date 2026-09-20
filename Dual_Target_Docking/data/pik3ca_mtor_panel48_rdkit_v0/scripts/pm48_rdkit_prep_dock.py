#!/usr/bin/env python3
"""JCIM Step 2A: re-prep PM48 with RDKit ETKDG+meeko and dock 4L23/4JT6 @ E=16."""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(
    "/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_rdkit_v0"
)
PANEL = ROOT / "tables" / "panel_v0_48.csv"
SEED = 20260727
EXHAUST = 16
N_MODES = 9
TARGETS = ["4L23", "4JT6"]
VINA = "/home/gwj/miniconda3/bin/vina"
MEKO_PY = "/home/gwj/miniconda3/bin/python"
MAX_WORKERS = 6


def prep_ligand(panel_id: str, smiles: str) -> Path:
    sdf = ROOT / "ligands_sdf" / f"{panel_id}.sdf"
    pdbqt = ROOT / "ligands_pdbqt" / f"{panel_id}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"bad smiles {panel_id}")
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
        [
            MEKO_PY,
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
        raise RuntimeError(f"meeko failed {panel_id}: {proc.stderr[-500:]}")
    return pdbqt


def load_box(target: str) -> dict:
    return json.loads((ROOT / "boxes" / f"{target}_box.json").read_text())


def write_conf(target: str, lig: str, ligand_pdbqt: Path) -> Path:
    box = load_box(target)
    conf_dir = ROOT / "logs" / "vina_confs"
    conf_dir.mkdir(parents=True, exist_ok=True)
    (ROOT / "logs" / "vina").mkdir(parents=True, exist_ok=True)
    conf = conf_dir / f"{target}_{lig}.txt"
    rec = ROOT / "receptors" / f"{target}_receptor.pdbqt"
    out = ROOT / "logs" / "vina" / f"{target}_{lig}_out.pdbqt"
    conf.write_text(
        "\n".join(
            [
                f"receptor = {rec}",
                f"ligand = {ligand_pdbqt}",
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
        + "\n"
    )
    return conf


def split_modes(out_pdbqt: Path, dest_dir: Path) -> int:
    dest_dir.mkdir(parents=True, exist_ok=True)
    text = out_pdbqt.read_text().splitlines()
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
    for i, m in enumerate(models, 1):
        (dest_dir / f"mode_{i:02d}.pdbqt").write_text("\n".join(m) + "\n")
    return len(models)


def run_one(target: str, lig: str, ligand_pdbqt: Path) -> dict:
    pose_dir = ROOT / "poses" / target / lig
    if (pose_dir / "mode_01.pdbqt").exists():
        n = len(list(pose_dir.glob("mode_*.pdbqt")))
        return {"target": target, "ligand": lig, "status": "exists", "n_modes": n}
    conf = write_conf(target, lig, ligand_pdbqt)
    log = ROOT / "logs" / "vina" / f"{target}_{lig}.log"
    out = ROOT / "logs" / "vina" / f"{target}_{lig}_out.pdbqt"
    proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True)
    log.write_text(proc.stdout + "\n" + proc.stderr)
    if proc.returncode != 0 or not out.exists():
        return {
            "target": target,
            "ligand": lig,
            "status": "fail",
            "reason": (proc.stderr or proc.stdout)[-300:],
            "n_modes": 0,
        }
    n = split_modes(out, pose_dir)
    return {"target": target, "ligand": lig, "status": "success", "n_modes": n}


def main():
    rows = list(csv.DictReader(PANEL.open()))
    print(f"PM48 ligands: {len(rows)}", flush=True)
    for r in rows:
        print("prep", r["panel_id"], flush=True)
        prep_ligand(r["panel_id"], r["smiles"])

    jobs = []
    for r in rows:
        lig = r["panel_id"]
        pdbqt = ROOT / "ligands_pdbqt" / f"{lig}.pdbqt"
        for target in TARGETS:
            jobs.append((target, lig, pdbqt))

    results = []
    print(f"dock jobs: {len(jobs)} workers={MAX_WORKERS} E={EXHAUST}", flush=True)
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {ex.submit(run_one, t, l, p): (t, l) for t, l, p in jobs}
        done = 0
        for fut in as_completed(futs):
            res = fut.result()
            results.append(res)
            done += 1
            print(
                f"[{done}/{len(jobs)}] {res['status']} {res['target']} {res['ligand']} n={res.get('n_modes')}",
                flush=True,
            )

    status_path = ROOT / "tables" / "job_status.csv"
    with status_path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["target", "ligand", "status", "n_modes", "reason"]
        )
        w.writeheader()
        for r in results:
            w.writerow(
                {
                    "target": r.get("target"),
                    "ligand": r.get("ligand"),
                    "status": r.get("status"),
                    "n_modes": r.get("n_modes"),
                    "reason": r.get("reason", ""),
                }
            )
    ok = sum(1 for r in results if r["status"] in ("success", "exists"))
    print(f"done ok={ok}/{len(results)}")
    return 0 if ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
