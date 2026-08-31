#!/usr/bin/env python3
"""T1.4 B3: Build PM110 strict panel, prep, dock E=16, RTM, optional GNINA."""
from __future__ import annotations

import csv
import json
import random
import subprocess
import sys
import urllib.request
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parents[3]
PM48 = REPO / "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv"
FC = REPO / "data/public_pair_selection/pik3ca_mtor_fourclass_chembl_ids.csv"
ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel110_rdkit_v0")
PM48_ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_rdkit_v0")
REPO_OUT = REPO / "data/pik3ca_mtor_panel110_rdkit_v0"

SEED = 20260729
DOCK_SEED = 20260727
EXHAUST = 16
N_MODES = 9
TARGETS = ["4L23", "4JT6"]
VINA = "/home/gwj/miniconda3/bin/vina"
MEKO_PY = "/home/gwj/miniconda3/bin/python"
MAX_WORKERS = 6
QUOTA = {"dual": 30, "A_only": 30, "B_only": 30, "neither": 25}


def strict_class(pA, pB):
    if pA >= 6.5 and pB >= 6.5:
        return "dual"
    if pA >= 6.5 and pB <= 5.5:
        return "A_only"
    if pB >= 6.5 and pA <= 5.5:
        return "B_only"
    if pA <= 5.5 and pB <= 5.5:
        return "neither"
    return None


def murcko(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return ""
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)


def fetch_smiles(cid, cache):
    if cid in cache:
        return cache[cid]
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json"
    req = urllib.request.Request(url, headers={"User-Agent": "dualfourclass-jcim/0.1"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = json.load(resp)
    smi = (data.get("molecule_structures") or {}).get("canonical_smiles")
    cache[cid] = smi
    time.sleep(0.03)
    return smi


def build_panel():
    rng = random.Random(SEED)
    pm48 = pd.read_csv(PM48)
    exclude = set(pm48["molecule_chembl_id"])
    fc = pd.read_csv(FC)
    rows = []
    seen_scaff = set()
    for _, r in fc.iterrows():
        cls = strict_class(r["pchembl_PIK3CA"], r["pchembl_MTOR"])
        if cls is None or r["molecule_chembl_id"] in exclude:
            continue
        rows.append({**r.to_dict(), "class": cls})
    by_cls = {c: [] for c in QUOTA}
    for r in rows:
        by_cls[r["class"]].append(r)
    cache = {}
    panel = []
    # include all PM48 first
    for i, r in pm48.iterrows():
        panel.append({
            "panel_id": r["panel_id"],
            "molecule_chembl_id": r["molecule_chembl_id"],
            "class": r["class"],
            "pchembl_PIK3CA": r["pchembl_PIK3CA"],
            "pchembl_MTOR": r["pchembl_MTOR"],
            "min_pchembl": r["min_pchembl"],
            "smiles": r["smiles"],
            "murcko_scaffold": r.get("murcko_scaffold") or murcko(r["smiles"]),
            "seed": "pm48_carryover",
        })
        seen_scaff.add(panel[-1]["murcko_scaffold"])

    idx = len(pm48)
    for cls, q in QUOTA.items():
        pool = by_cls[cls][:]
        rng.shuffle(pool)
        got = sum(1 for p in panel if p["class"] == cls)
        need = max(0, q - got)
        for r in pool:
            if need <= 0:
                break
            smi = fetch_smiles(r["molecule_chembl_id"], cache)
            if not smi:
                continue
            sc = murcko(smi)
            if sc in seen_scaff:
                continue
            seen_scaff.add(sc)
            idx += 1
            panel.append({
                "panel_id": f"PM110_{idx:03d}",
                "molecule_chembl_id": r["molecule_chembl_id"],
                "class": cls,
                "pchembl_PIK3CA": r["pchembl_PIK3CA"],
                "pchembl_MTOR": r["pchembl_MTOR"],
                "min_pchembl": min(r["pchembl_PIK3CA"], r["pchembl_MTOR"]),
                "smiles": smi,
                "murcko_scaffold": sc,
                "seed": SEED,
            })
            need -= 1
    return panel, cache


def setup_workspace():
    ROOT.mkdir(parents=True, exist_ok=True)
    REPO_OUT.mkdir(parents=True, exist_ok=True)
    for sub in ("receptors", "boxes"):
        dst = ROOT / sub
        if not dst.exists():
            dst.symlink_to(PM48_ROOT / sub)


def prep_ligand(panel_id, smiles):
    from pathlib import Path

    pdbqt = ROOT / "ligands_pdbqt" / f"{panel_id}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = DOCK_SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        AllChem.EmbedMolecule(mol, randomSeed=DOCK_SEED)
    sdf = ROOT / "ligands_sdf" / f"{panel_id}.sdf"
    sdf.parent.mkdir(parents=True, exist_ok=True)
    pdbqt.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(sdf))
    mol.SetProp("_Name", panel_id)
    w.write(mol)
    w.close()
    proc = subprocess.run(
        [MEKO_PY, "/home/gwj/miniconda3/bin/mk_prepare_ligand.py", "-i", str(sdf), "-o", str(pdbqt)],
        capture_output=True, text=True,
    )
    return pdbqt if proc.returncode == 0 else None


def dock_panel(panel):
    import re

    def run_one(target, lig, pdbqt):
        box = json.loads((ROOT / "boxes" / f"{target}_box.json").read_text())
        rec = ROOT / "receptors" / f"{target}_receptor.pdbqt"
        out = ROOT / "logs" / "vina" / f"{target}_{lig}_out.pdbqt"
        log = ROOT / "logs" / "vina" / f"{target}_{lig}.log"
        out.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            VINA, "--receptor", str(rec), "--ligand", str(pdbqt),
            "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]), "--center_z", str(box["center_z"]),
            "--size_x", str(box["size_x"]), "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
            "--exhaustiveness", str(EXHAUST), "--num_modes", str(N_MODES), "--seed", str(DOCK_SEED), "--cpu", "1", "--out", str(out),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        log.write_text(proc.stdout + proc.stderr)
        affs = [float(m.group(1)) for line in (proc.stdout + proc.stderr).splitlines() if (m := re.search(r"^\s*1\s+(-?\d+\.\d+)", line))]
        return target, lig, min(affs) if affs else None, proc.returncode == 0

    jobs = []
    for r in panel:
        pq = prep_ligand(r["panel_id"], r["smiles"])
        if pq:
            for t in TARGETS:
                jobs.append((t, r["panel_id"], pq))

    results = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = [ex.submit(run_one, t, l, p) for t, l, p in jobs]
        for i, fut in enumerate(as_completed(futs), 1):
            results.append(fut.result())
            if i % 20 == 0:
                print(f"docked {i}/{len(jobs)}", flush=True)
    return results


def main():
    setup_workspace()
    panel, _ = build_panel()
    print(f"Panel size: {len(panel)}", flush=True)
    tab = REPO_OUT / "tables"
    tab.mkdir(parents=True, exist_ok=True)
    (ROOT / "tables").mkdir(parents=True, exist_ok=True)
    fields = ["panel_id", "molecule_chembl_id", "class", "pchembl_PIK3CA", "pchembl_MTOR", "min_pchembl", "smiles", "murcko_scaffold", "seed"]
    for path in (tab / "panel_v0_110.csv", ROOT / "tables" / "panel_v0_110.csv"):
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(panel)

    results = dock_panel(panel)
    # write long scores
    rows = []
    for target, lig, aff, ok in results:
        if aff is not None:
            rows.append({"ligand": lig, "target": target, "vina_best": aff, "status": "ok" if ok else "fail"})
    with (tab / "scores_vina_long.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["ligand", "target", "vina_best", "status"])
        w.writeheader()
        w.writerows(rows)
    (ROOT / "tables" / "scores_vina_long.csv").write_text((tab / "scores_vina_long.csv").read_text())
    print(f"Done panel={len(panel)} docked={len(rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
