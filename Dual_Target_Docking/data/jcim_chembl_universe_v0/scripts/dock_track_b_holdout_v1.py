#!/usr/bin/env python3
"""Dock Track B holdout panels with the same Vina protocol as primary production."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_track_b_v0"
HOLD = LOCAL / "allpairs_stack" / "holdout"
# Frozen IDs from dump-gated (do not re-draw).
FROZEN_HOLDOUT_DIR = LOCAL / "tables" / "five_pair_dump_gated_v1"
VINA = "/home/gwj/miniconda3/bin/vina"
PY = sys.executable
MK = "/home/gwj/miniconda3/bin/mk_prepare_ligand.py"
SEED = 20260727
N_MODES = 9
ENERGY = 3
EXHAUST = 8
TIMEOUT = 600  # timeout → skip (no retry)
SKIP_TORS = 25

PAIR_TARGETS = {
    "F2/F10": ["4UDW", "2JKH"],
    "JAK1/TYK2": ["6N7A", "3LXP"],
    "JAK1/JAK2": ["6N7A", "8BXH"],
    "PPARG/PPARA": ["9V8H", "6LXA"],
    "PPARA/PPARD": ["6LXA", "5U3Q"],
}


def prep_ligand(panel_id: str, smiles: str) -> Path:
    sdf = HOLD / "ligands_sdf" / f"{panel_id}.sdf"
    pdbqt = HOLD / "ligands_pdbqt" / f"{panel_id}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    sdf.parent.mkdir(parents=True, exist_ok=True)
    pdbqt.parent.mkdir(parents=True, exist_ok=True)
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
    w = Chem.SDWriter(str(sdf))
    mol.SetProp("_Name", panel_id)
    w.write(mol)
    w.close()
    proc = subprocess.run([PY, MK, "-i", str(sdf), "-o", str(pdbqt)], capture_output=True, text=True)
    if proc.returncode != 0 or not pdbqt.exists():
        raise RuntimeError(proc.stderr[-300:] or "meeko fail")
    return pdbqt


def torsdof(path: Path) -> int:
    for line in path.read_text().splitlines():
        if line.startswith("TORSDOF"):
            return int(line.split()[1])
    return 0


def parse_e(path: Path) -> float | None:
    for line in path.read_text().splitlines():
        if "VINA RESULT" in line:
            return float(line.split()[3])
    return None


def split_modes(out: Path, dest: Path) -> int:
    dest.mkdir(parents=True, exist_ok=True)
    models, cur = [], []
    for line in out.read_text().splitlines():
        if line.startswith("MODEL"):
            cur = [line]
        elif line.startswith("ENDMDL"):
            cur.append(line)
            models.append(cur)
            cur = []
        elif cur:
            cur.append(line)
    for i, m in enumerate(models, 1):
        (dest / f"mode_{i:02d}.pdbqt").write_text("\n".join(m) + "\n")
    return len(models)


def dock_one(pair, target, lig, lig_pdbqt, box):
    pose = HOLD / "poses" / target / lig
    if (pose / "mode_01.pdbqt").exists():
        e = parse_e(pose / "mode_01.pdbqt")
        return {"pair": pair, "target": target, "ligand": lig, "status": "exists",
                "mode1_energy": e, "score_S": None if e is None else -e, "reason": ""}
    td = torsdof(lig_pdbqt)
    if td >= SKIP_TORS:
        return {"pair": pair, "target": target, "ligand": lig, "status": "skip",
                "mode1_energy": None, "score_S": None, "reason": f"skip_torsdof={td}"}
    rec = LOCAL / "receptors" / f"{target}_receptor.pdbqt"
    conf_dir = HOLD / "logs" / "confs"
    conf_dir.mkdir(parents=True, exist_ok=True)
    out = HOLD / "logs" / f"{target}_{lig}_out.pdbqt"
    conf = conf_dir / f"{target}_{lig}.txt"
    conf.write_text(
        "\n".join([
            f"receptor = {rec}", f"ligand = {lig_pdbqt}",
            f"center_x = {box['center_x']}", f"center_y = {box['center_y']}", f"center_z = {box['center_z']}",
            f"size_x = {box['size_x']}", f"size_y = {box['size_y']}", f"size_z = {box['size_z']}",
            f"exhaustiveness = {EXHAUST}", f"num_modes = {N_MODES}", f"energy_range = {ENERGY}",
            "cpu = 1", f"seed = {SEED}", f"out = {out}",
        ]) + "\n"
    )
    try:
        proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return {"pair": pair, "target": target, "ligand": lig, "status": "skip",
                "mode1_energy": None, "score_S": None, "reason": f"timeout_{TIMEOUT}s"}
    (HOLD / "logs" / f"{target}_{lig}.log").write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
    if proc.returncode != 0 or not out.exists():
        return {"pair": pair, "target": target, "ligand": lig, "status": "fail",
                "mode1_energy": None, "score_S": None, "reason": (proc.stderr or "")[-300:]}
    split_modes(out, pose)
    e = parse_e(out)
    return {"pair": pair, "target": target, "ligand": lig, "status": "success",
            "mode1_energy": e, "score_S": None if e is None else -e, "reason": ""}


def main() -> int:
    global TIMEOUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--timeout", type=int, default=TIMEOUT, help="per-job Vina timeout seconds; timeout → skip")
    args = ap.parse_args()
    TIMEOUT = args.timeout
    panels = sorted(FROZEN_HOLDOUT_DIR.glob("holdout_panel_HO*_v1.csv"))
    if not panels:
        raise SystemExit(f"no frozen holdout panels in {FROZEN_HOLDOUT_DIR}")
    boxes = {p.stem.replace("_box", ""): json.loads(p.read_text()) for p in (LOCAL / "boxes").glob("*_box.json")}
    jobs = []
    for panel_csv in panels:
        rows = list(csv.DictReader(panel_csv.open()))
        pair = rows[0]["pair"]
        print(f"prep ligands for {pair} n={len(rows)} (frozen IDs; timeout={TIMEOUT}s → skip)", flush=True)
        for r in rows:
            hid = r.get("holdout_id") or r.get("panel_id")
            pdbqt = prep_ligand(hid, r["canonical_smiles"])
            for t in PAIR_TARGETS[pair]:
                jobs.append((pair, t, hid, pdbqt, boxes[t]))
    print(f"holdout dock jobs={len(jobs)}", flush=True)
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(dock_one, *j): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            results.append(r)
            if i % 25 == 0 or r["status"] not in ("success", "exists"):
                print(f"[{i}/{len(jobs)}] {r['status']} {r['pair']} {r['target']} {r['ligand']} {r.get('reason','')}", flush=True)
    tab = HOLD / "tables"
    tab.mkdir(parents=True, exist_ok=True)
    path = tab / "holdout_scores_vina_mode1_v1.csv"
    fields = ["pair", "target", "ligand", "status", "mode1_energy", "score_S", "reason"]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(results, key=lambda x: (x["pair"], x["target"], x["ligand"])))
    ok = sum(1 for r in results if r["status"] in ("success", "exists"))
    skip = sum(1 for r in results if r["status"] == "skip")
    fail = sum(1 for r in results if r["status"] == "fail")
    print(f"done ok={ok} skip={skip} fail={fail} / {len(results)} -> {path}", flush=True)
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
