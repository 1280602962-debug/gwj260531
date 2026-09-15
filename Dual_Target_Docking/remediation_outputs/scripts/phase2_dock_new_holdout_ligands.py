#!/usr/bin/env python3
"""Dock 7 new AChE/BChE unused-pool holdout ligands (protocol-identical to HOAB).

Reuses 53 unchanged holdout scores. New ligands use:
  RDKit ETKDG seed 20260727 + MMFF + meeko
  Vina 1.2.7 E=8 n_modes=9 energy_range=3 seed=20260727
  receptors 4EY7/4BDS (hash-identical to ACHE/BCHE)
  boxes 4EY7_box.json / 4BDS_box.json
"""
from __future__ import annotations

import csv
import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
SRC = ROOT / "data/ache_bche_panel_v0"
CMP = ROOT / "remediation_outputs/phase2_ache_bche/holdout_old_vs_corrected.csv"
NEW_HO = ROOT / "remediation_outputs/phase2_ache_bche/holdout_panel_HOAB_corrected.csv"
OUT = ROOT / "remediation_outputs/phase2_ache_bche/vina_holdout_new"
VINA = "/home/gwj/miniconda3/bin/vina"
MEKO_PY = "/home/gwj/miniconda3/bin/python"
MKO = "/home/gwj/miniconda3/bin/mk_prepare_ligand.py"
SEED = 20260727
WORKERS = 2


def parse_mode1(path: Path) -> float | None:
    m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)", path.read_text(errors="replace"))
    return float(m.group(1)) if m else None


def prep(lig: str, smiles: str) -> Path:
    sdf = OUT / "ligands_sdf" / f"{lig}.sdf"
    pdbqt = OUT / "ligands_pdbqt" / f"{lig}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"bad smiles {lig}")
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
    mol.SetProp("_Name", lig)
    w.write(mol)
    w.close()
    proc = subprocess.run([MEKO_PY, MKO, "-i", str(sdf), "-o", str(pdbqt)], capture_output=True, text=True)
    if proc.returncode != 0 or not pdbqt.exists():
        raise RuntimeError(f"meeko {lig}: {proc.stderr[-400:]}")
    return pdbqt


def dock(target: str, rec_name: str, lig: str, ligand: Path, box: dict) -> dict:
    rec = SRC / "receptors" / f"{rec_name}_receptor.pdbqt"
    out = OUT / "out" / f"{rec_name}_{lig}_out.pdbqt"
    log = OUT / "logs" / f"{rec_name}_{lig}.log"
    conf = OUT / "confs" / f"{rec_name}_{lig}.txt"
    for d in (OUT / "out", OUT / "logs", OUT / "confs"):
        d.mkdir(parents=True, exist_ok=True)
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
    if out.exists() and parse_mode1(out) is not None:
        return {"target": rec_name, "ligand": lig, "status": "cached", "affinity": parse_mode1(out)}
    proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True)
    log.write_text((proc.stdout or "") + (proc.stderr or ""))
    aff = parse_mode1(out) if out.exists() else None
    return {
        "target": rec_name,
        "ligand": lig,
        "status": "ok" if aff is not None else f"FAIL rc={proc.returncode}",
        "affinity": aff if aff is not None else "",
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    added = [r["ligand_id"] for r in csv.DictReader(CMP.open()) if r["old_holdout"] == "no" and r["new_holdout"] == "yes"]
    panel = {r["molecule_chembl_id"]: r for r in csv.DictReader(NEW_HO.open())}
    boxes = {
        "4EY7": json.loads((SRC / "boxes/4EY7_box.json").read_text()),
        "4BDS": json.loads((SRC / "boxes/4BDS_box.json").read_text()),
    }
    jobs = []
    for cid in added:
        r = panel[cid]
        lig = r["holdout_id"]
        ligand = prep(lig, r["smiles"])
        jobs.append(("A", "4EY7", lig, ligand, boxes["4EY7"], cid, r["class"]))
        jobs.append(("B", "4BDS", lig, ligand, boxes["4BDS"], cid, r["class"]))
    print(f"jobs={len(jobs)} new_ligands={len(added)}")
    results = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(dock, rec, rec, lig, pdbqt, box): (cid, cls, lig, rec) for _arm, rec, lig, pdbqt, box, cid, cls in jobs}
        # fix: dock() signature is (target, rec_name, lig, ligand, box)
        pass
    # rewrite more clearly
    results = []
    packed = []
    for cid in added:
        r = panel[cid]
        lig = r["holdout_id"]
        ligand = OUT / "ligands_pdbqt" / f"{lig}.pdbqt"
        packed.append((cid, r["class"], lig, ligand))
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = []
        meta = {}
        for cid, cls, lig, ligand in packed:
            for rec_name, box in boxes.items():
                fut = ex.submit(dock, rec_name, rec_name, lig, ligand, box)
                meta[fut] = (cid, cls, lig, rec_name)
                futs.append(fut)
        done = 0
        for fut in as_completed(futs):
            rec = fut.result()
            cid, cls, lig, rec_name = meta[fut]
            rec = dict(rec)
            rec["molecule_chembl_id"] = cid
            rec["class"] = cls
            rec["holdout_id"] = lig
            results.append(rec)
            done += 1
            print(f"[{done}/{len(futs)}] {rec['target']} {lig} {rec['status']} {rec.get('affinity','')}", flush=True)
    dest = OUT / "scores_vina_mode1_new_holdout.csv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)
    n_ok = sum(1 for r in results if r["status"] in {"ok", "cached"})
    print(f"done {n_ok}/{len(results)} -> {dest}")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
