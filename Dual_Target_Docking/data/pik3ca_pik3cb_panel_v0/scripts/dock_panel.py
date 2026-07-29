#!/usr/bin/env python3
"""Generic Step-3 panel dock: RDKit+meeko → Vina for a frozen pair pack."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

SEED = 20260727
N_MODES = 9
VINA = "/home/gwj/miniconda3/bin/vina"
MEKO_PY = "/home/gwj/miniconda3/bin/python"
MAX_WORKERS = 4  # leave room for PM jobs


def prep_ligand(root: Path, panel_id: str, smiles: str) -> Path:
    sdf = root / "ligands_sdf" / f"{panel_id}.sdf"
    pdbqt = root / "ligands_pdbqt" / f"{panel_id}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"bad smiles {panel_id}")
    # drop salts / counterions — keep largest organic fragment
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
        raise RuntimeError(f"meeko {panel_id}: {proc.stderr[-400:]}")
    return pdbqt


def write_conf(root, target, lig, ligand_pdbqt, box, exhaust):
    conf_dir = root / "logs" / "vina_confs"
    conf_dir.mkdir(parents=True, exist_ok=True)
    (root / "logs" / "vina").mkdir(parents=True, exist_ok=True)
    conf = conf_dir / f"{target}_{lig}.txt"
    rec = root / "receptors" / f"{target}_receptor.pdbqt"
    out = root / "logs" / "vina" / f"{target}_{lig}_out.pdbqt"
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
                f"exhaustiveness = {exhaust}",
                f"num_modes = {N_MODES}",
                "energy_range = 3",
                "cpu = 1",
                f"seed = {SEED}",
                f"out = {out}",
            ]
        )
        + "\n"
    )
    return conf, out


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


def torsdof(ligand_pdbqt: Path) -> int:
    for line in ligand_pdbqt.read_text().splitlines():
        if line.startswith("TORSDOF"):
            return int(line.split()[1])
    return 0


def run_one(root, target, lig, ligand_pdbqt, box, exhaust, timeout_s=600):
    pose_dir = root / "poses" / target / lig
    if (pose_dir / "mode_01.pdbqt").exists():
        return {
            "target": target,
            "ligand": lig,
            "status": "exists",
            "n_modes": len(list(pose_dir.glob("mode_*.pdbqt"))),
        }
    td = torsdof(ligand_pdbqt)
    # Peptide-scale ligands thrash Vina; skip rather than hang workers
    if td >= 25:
        return {
            "target": target,
            "ligand": lig,
            "status": "fail",
            "reason": f"skip_torsdof={td}_ge_25",
            "n_modes": 0,
        }
    use_e = exhaust
    use_to = timeout_s
    if td >= 20:
        use_e = min(exhaust, 4)
        use_to = max(timeout_s, 900)
    conf, out = write_conf(root, target, lig, ligand_pdbqt, box, use_e)
    log = root / "logs" / "vina" / f"{target}_{lig}.log"
    try:
        proc = subprocess.run(
            [VINA, "--config", str(conf)],
            capture_output=True,
            text=True,
            timeout=use_to,
        )
        log.write_text(proc.stdout + "\n" + proc.stderr)
        rc, err = proc.returncode, (proc.stderr or proc.stdout)[-300:]
    except subprocess.TimeoutExpired as e:
        log.write_text(f"TIMEOUT after {use_to}s torsdof={td} E={use_e}\n{(e.stdout or b'').decode(errors='ignore')}")
        return {
            "target": target,
            "ligand": lig,
            "status": "fail",
            "reason": f"timeout_{use_to}s_torsdof={td}_E={use_e}",
            "n_modes": 0,
        }
    if rc != 0 or not out.exists():
        return {
            "target": target,
            "ligand": lig,
            "status": "fail",
            "reason": err,
            "n_modes": 0,
        }
    n = split_modes(out, pose_dir)
    return {"target": target, "ligand": lig, "status": "success", "n_modes": n}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--targets", nargs=2, required=True)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--workers", type=int, default=MAX_WORKERS)
    ap.add_argument("--timeout", type=int, default=600, help="per-job Vina timeout seconds")
    args = ap.parse_args()
    root = Path(args.root)
    panel = list(csv.DictReader(open(args.panel)))
    boxes = {
        t: json.loads((root / "boxes" / f"{t}_box.json").read_text())
        for t in args.targets
    }
    print(f"panel n={len(panel)} targets={args.targets} E={args.exhaustiveness}", flush=True)
    for r in panel:
        if not r.get("smiles"):
            raise SystemExit(f"missing smiles {r['panel_id']}")
        print("prep", r["panel_id"], flush=True)
        prep_ligand(root, r["panel_id"], r["smiles"])

    jobs = []
    for r in panel:
        lig = r["panel_id"]
        pdbqt = root / "ligands_pdbqt" / f"{lig}.pdbqt"
        for t in args.targets:
            jobs.append((t, lig, pdbqt, boxes[t]))
    # dock easy ligands first so high-TORSDOF cannot monopolize the pool
    jobs.sort(key=lambda j: torsdof(j[2]))

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {
            ex.submit(run_one, root, t, l, p, b, args.exhaustiveness, args.timeout): (t, l)
            for t, l, p, b in jobs
        }
        done = 0
        for fut in as_completed(futs):
            res = fut.result()
            results.append(res)
            done += 1
            print(
                f"[{done}/{len(jobs)}] {res['status']} {res['target']} {res['ligand']}",
                flush=True,
            )

    status = root / "tables" / "job_status.csv"
    status.parent.mkdir(parents=True, exist_ok=True)
    with status.open("w", newline="") as fh:
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
