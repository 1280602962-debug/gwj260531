#!/usr/bin/env python3
"""RDKit ETKDGv3 + meeko PDBQT for Track B panels.

Same prep as frozen DualFourClass ligands (seed 20260727, largest fragment,
MMFF 200 steps). Does not dock. Cloud has no Vina.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

SEED = 20260727
TABLES = Path(__file__).resolve().parents[1] / "tables"
CACHE = Path(__file__).resolve().parents[1] / "cache" / "track_b_ligands"


def resolve_mk() -> str:
    from shutil import which

    for cand in (
        which("mk_prepare_ligand.py"),
        str(Path.home() / ".local/bin/mk_prepare_ligand.py"),
        str(Path(sys.executable).resolve().parent / "mk_prepare_ligand.py"),
    ):
        if cand and Path(cand).exists():
            return cand
    raise SystemExit("mk_prepare_ligand.py not found")


MK = resolve_mk()


def prep_one(panel_id: str, smiles: str, pair: str, cls: str) -> dict:
    sdf = CACHE / "sdf" / f"{panel_id}.sdf"
    pdbqt = CACHE / "pdbqt" / f"{panel_id}.pdbqt"
    sdf.parent.mkdir(parents=True, exist_ok=True)
    pdbqt.parent.mkdir(parents=True, exist_ok=True)
    rec = {
        "panel_id": panel_id,
        "pair": pair,
        "class": cls,
        "status": "fail",
        "reason": "",
        "pdbqt": "",
    }
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        rec["reason"] = "bad_smiles"
        return rec
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        if AllChem.EmbedMolecule(mol, randomSeed=SEED) != 0:
            rec["reason"] = "embed_fail"
            return rec
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    w = Chem.SDWriter(str(sdf))
    mol.SetProp("_Name", panel_id)
    w.write(mol)
    w.close()
    proc = subprocess.run(
        [sys.executable, MK, "-i", str(sdf), "-o", str(pdbqt)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    if proc.returncode != 0 or not pdbqt.exists() or pdbqt.stat().st_size == 0:
        rec["reason"] = (proc.stderr or proc.stdout or "meeko_fail")[-300:]
        return rec
    rec["status"] = "ok"
    rec["pdbqt"] = str(pdbqt.name)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    jobs = []
    for csv_path in sorted((TABLES / "track_b_panels").glob("panel_*_v1.csv")):
        for row in csv.DictReader(csv_path.open()):
            jobs.append((row["panel_id"], row["canonical_smiles"], row["pair"], row["class"]))
    print(f"prepping {len(jobs)} ligands with {args.workers} workers...", flush=True)
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(prep_one, *j) for j in jobs]
        for i, fut in enumerate(as_completed(futs), 1):
            rows.append(fut.result())
            if i % 50 == 0 or i == len(jobs):
                n_ok = sum(1 for r in rows if r["status"] == "ok")
                print(f"  {i}/{len(jobs)} done; ok={n_ok}", flush=True)
    rows.sort(key=lambda r: r["panel_id"])
    out = TABLES / "track_b_ligand_prep_status_v1.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["panel_id", "pair", "class", "status", "reason", "pdbqt"])
        w.writeheader()
        w.writerows(rows)
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_fail = len(rows) - n_ok
    print(f"wrote {out} ok={n_ok} fail={n_fail}")
    if n_fail:
        for r in rows:
            if r["status"] != "ok":
                print(f"FAIL {r['panel_id']} {r['pair']} {r['reason'][:120]}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
