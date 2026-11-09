#!/usr/bin/env python3
"""Prepare all 110 EGFR/HER2 ligands with the Track-B RDKit/Meeko protocol.

ETKDGv3 seed 20260727, MMFF maxIters=200, Meeko 0.7.1.
Largest-fragment rule matches Track-B prep_track_b_ligands_v1.py.
Does not dock. Does not use LigPrep or OpenBabel.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/egfr_her2_uniform_rdkit_v1"
PANEL = ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"
SEED = 20260727
MMFF_ITERS = 200


def resolve_mk() -> str:
    from shutil import which

    for cand in (
        which("mk_prepare_ligand.py"),
        str(Path(sys.executable).resolve().parent / "mk_prepare_ligand.py"),
    ):
        if cand and Path(cand).exists():
            return str(cand)
    raise SystemExit("mk_prepare_ligand.py not found")


def prep_one(ligand_id: str, smiles: str, mk: str) -> dict:
    sdf_dir = OUT / "ligands_sdf"
    pdbqt_dir = OUT / "ligands_pdbqt"
    sdf_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    sdf = sdf_dir / f"{ligand_id}.sdf"
    pdbqt = pdbqt_dir / f"{ligand_id}.pdbqt"
    rec = {
        "ligand_id": ligand_id,
        "smiles": smiles,
        "fragment_status": "",
        "embed_status": "",
        "mmff_status": "",
        "meeko_status": "",
        "pdbqt_status": "",
        "n_frags": "",
        "reason": "",
    }
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        rec.update(
            fragment_status="fail",
            embed_status="not_run",
            mmff_status="not_run",
            meeko_status="not_run",
            pdbqt_status="fail",
            reason="bad_smiles",
        )
        return rec
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    rec["n_frags"] = str(len(frags))
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        rec["fragment_status"] = "largest_organic_or_heaviest_kept"
    else:
        rec["fragment_status"] = "single"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) == 0:
        rec["embed_status"] = "etkdgv3_ok"
    elif AllChem.EmbedMolecule(mol, randomSeed=SEED) == 0:
        rec["embed_status"] = "etkdg_fallback"
    else:
        rec.update(
            embed_status="fail",
            mmff_status="not_run",
            meeko_status="not_run",
            pdbqt_status="fail",
            reason="embed_fail",
        )
        return rec
    try:
        code = AllChem.MMFFOptimizeMolecule(mol, maxIters=MMFF_ITERS)
        rec["mmff_status"] = {0: "ok", 1: "not_converged", -1: "setup_fail"}.get(code, str(code))
    except Exception as exc:
        rec["mmff_status"] = f"exception:{type(exc).__name__}"
    w = Chem.SDWriter(str(sdf))
    mol.SetProp("_Name", ligand_id)
    w.write(mol)
    w.close()
    proc = subprocess.run(
        [mk, "-i", str(sdf), "-o", str(pdbqt)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    if proc.returncode != 0 or not pdbqt.exists() or pdbqt.stat().st_size == 0:
        rec["meeko_status"] = "fail"
        rec["pdbqt_status"] = "fail"
        rec["reason"] = (proc.stderr or proc.stdout or "meeko_fail")[-300:]
        return rec
    rec["meeko_status"] = "ok"
    rec["pdbqt_status"] = "ok"
    return rec


def main() -> int:
    mk = resolve_mk()
    rows_in = list(csv.DictReader(PANEL.open(encoding="utf-8-sig", newline="")))
    if len(rows_in) != 110:
        raise SystemExit(f"expected 110 panel rows, got {len(rows_in)}")
    out_rows = []
    for i, row in enumerate(rows_in, 1):
        rec = prep_one(row["panel_id"], row["smiles"], mk)
        out_rows.append(rec)
        print(f"{i}/110 {rec['ligand_id']} pdbqt={rec['pdbqt_status']} embed={rec['embed_status']}", flush=True)
    dest = OUT / "tables"
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / "ligand_prep_status.csv"
    fields = [
        "ligand_id",
        "smiles",
        "fragment_status",
        "embed_status",
        "mmff_status",
        "meeko_status",
        "pdbqt_status",
        "n_frags",
        "reason",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)
    n_ok = sum(1 for r in out_rows if r["pdbqt_status"] == "ok")
    print(f"wrote {path} ok={n_ok} fail={len(out_rows) - n_ok}")
    return 0 if n_ok == 110 else 1


if __name__ == "__main__":
    raise SystemExit(main())
