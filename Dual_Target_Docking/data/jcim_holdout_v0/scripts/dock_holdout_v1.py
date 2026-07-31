#!/usr/bin/env python3
"""Dock the frozen holdout panels against the already-frozen receptors.

Protocol is NOT re-tuned here: ligand prep (RDKit ETKDGv3 seed=20260727 + MMFF200
+ meeko default PDBQT), receptor PDBQT, box, and exhaustiveness are all reused
byte-for-byte from the frozen main-panel packs (Methods 2.4-2.5, Table 1). Only
the ligand identities are new (drawn by build_holdout_candidate_pool_v1.py from
the ChEMBL pool never touched during panel construction).

Resumable: skips any (pair, receptor, ligand) whose score is already recorded.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
HOLDOUT = ROOT / "data/jcim_holdout_v0"
MK_PREPARE_LIGAND = str(Path.home() / ".local/bin/mk_prepare_ligand.py")
SEED = 20260727  # frozen protocol seed, same as main panels

PAIR_CONFIGS = {
    "HOPM": {
        "pair_label": "PIK3CA/mTOR",
        "receptors": {
            "A": {"name": "4L23", "pdbqt": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4L23_receptor.pdbqt", "box": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/boxes/4L23_box.json"},
            "B": {"name": "4JT6", "pdbqt": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4JT6_receptor.pdbqt", "box": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/boxes/4JT6_box.json"},
        },
        "exhaustiveness": 16,
    },
    "HOAB": {
        "pair_label": "AChE/BChE",
        "receptors": {
            "A": {"name": "4EY7", "pdbqt": ROOT / "data/ache_bche_panel_v0/receptors/4EY7_receptor.pdbqt", "box": ROOT / "data/ache_bche_panel_v0/boxes/4EY7_box.json"},
            "B": {"name": "4BDS", "pdbqt": ROOT / "data/ache_bche_panel_v0/receptors/4BDS_receptor.pdbqt", "box": ROOT / "data/ache_bche_panel_v0/boxes/4BDS_box.json"},
        },
        "exhaustiveness": 8,
    },
    "HOAP": {
        "pair_label": "PIK3CA/PIK3CB",
        "receptors": {
            "A": {"name": "4L23", "pdbqt": ROOT / "data/pik3ca_pik3cb_panel_v0/receptors/4L23_receptor.pdbqt", "box": ROOT / "data/pik3ca_pik3cb_panel_v0/boxes/4L23_box.json"},
            "B": {"name": "2WXF", "pdbqt": ROOT / "data/pik3ca_pik3cb_panel_v0/receptors/2WXF_receptor.pdbqt", "box": ROOT / "data/pik3ca_pik3cb_panel_v0/boxes/2WXF_box.json"},
        },
        "exhaustiveness": 8,
    },
}


def prep_ligand(smiles: str, out_dir: Path, lig_id: str) -> Path | None:
    sdf = out_dir / "ligands_sdf" / f"{lig_id}.sdf"
    pdbqt = out_dir / "ligands_pdbqt" / f"{lig_id}.pdbqt"
    sdf.parent.mkdir(parents=True, exist_ok=True)
    pdbqt.parent.mkdir(parents=True, exist_ok=True)
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        if AllChem.EmbedMolecule(mol, randomSeed=SEED) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    w = Chem.SDWriter(str(sdf))
    mol.SetProp("_Name", lig_id)
    w.write(mol)
    w.close()
    proc = subprocess.run(
        [sys.executable, MK_PREPARE_LIGAND, "-i", str(sdf), "-o", str(pdbqt)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    if proc.returncode != 0 or not pdbqt.exists():
        (out_dir / "logs" / f"prep_fail_{lig_id}.log").write_text(proc.stdout + "\n" + proc.stderr)
        return None
    return pdbqt


def dock_one(args) -> dict:
    prefix, receptor_key, receptor_cfg, lig_id, ligand_pdbqt_str, exhaust, out_dir_str = args
    from vina import Vina  # re-import inside worker process

    out_dir = Path(out_dir_str)
    ligand_pdbqt = Path(ligand_pdbqt_str)
    box = json.loads(Path(receptor_cfg["box"]).read_text())
    pose_dir = out_dir / "poses" / receptor_cfg["name"] / lig_id
    pose_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = pose_dir / "out.pdbqt"

    t0 = time.time()
    try:
        v = Vina(sf_name="vina", cpu=1, seed=SEED, verbosity=0)
        v.set_receptor(str(receptor_cfg["pdbqt"]))
        v.set_ligand_from_file(str(ligand_pdbqt))
        v.compute_vina_maps(
            center=[box["center_x"], box["center_y"], box["center_z"]],
            box_size=[box["size_x"], box["size_y"], box["size_z"]],
        )
        v.dock(exhaustiveness=exhaust, n_poses=9)
        v.write_poses(str(out_pdbqt), n_poses=9, overwrite=True, energy_range=3)
        energies = v.energies(n_poses=9)
        mode1 = float(energies[0][0])
        status = "success"
        reason = ""
    except Exception as exc:  # noqa: BLE001
        mode1 = None
        status = "fail"
        reason = f"{type(exc).__name__}: {exc}"[:300]
    dt = time.time() - t0
    return {
        "prefix": prefix,
        "receptor_key": receptor_key,
        "receptor": receptor_cfg["name"],
        "ligand": lig_id,
        "vina_mode1": mode1,
        "status": status,
        "reason": reason,
        "seconds": round(dt, 1),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, choices=list(PAIR_CONFIGS))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None, help="debug: only first N ligands")
    args = ap.parse_args()

    cfg = PAIR_CONFIGS[args.prefix]
    panel_csv = HOLDOUT / "tables" / f"holdout_panel_{args.prefix}.csv"
    panel = pd.read_csv(panel_csv)
    if args.limit:
        panel = panel.head(args.limit)

    out_dir = HOLDOUT / args.prefix
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)

    scores_path = HOLDOUT / "tables" / f"scores_vina_mode1_{args.prefix}.csv"
    done_keys = set()
    if scores_path.exists():
        prev = pd.read_csv(scores_path)
        done_keys = set(zip(prev["receptor"], prev["ligand"]))
        print(f"resuming: {len(done_keys)} (receptor,ligand) already scored")
    else:
        prev = pd.DataFrame()

    print(f"prepping ligands for {cfg['pair_label']} ({len(panel)} ligands)...", flush=True)
    ligand_pdbqt = {}
    for _, row in panel.iterrows():
        lig_id = row["holdout_id"]
        p = prep_ligand(row["smiles"], out_dir, lig_id)
        if p is not None:
            ligand_pdbqt[lig_id] = str(p)
        else:
            print(f"  PREP FAIL {lig_id}", flush=True)

    jobs = []
    for rk, rcfg in cfg["receptors"].items():
        for lig_id, pdbqt_path in ligand_pdbqt.items():
            if (rcfg["name"], lig_id) in done_keys:
                continue
            jobs.append((args.prefix, rk, rcfg, lig_id, pdbqt_path, cfg["exhaustiveness"], str(out_dir)))

    print(f"docking jobs to run: {len(jobs)} (workers={args.workers})", flush=True)
    results = list(prev.to_dict("records"))
    t_start = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(dock_one, j): j for j in jobs}
        n_done = 0
        for fut in as_completed(futs):
            try:
                res = fut.result()
            except Exception as exc:  # noqa: BLE001
                j = futs[fut]
                res = {
                    "prefix": j[0], "receptor_key": j[1], "receptor": j[2]["name"],
                    "ligand": j[3], "vina_mode1": None, "status": "fail",
                    "reason": f"worker_crash: {exc}", "seconds": None,
                }
            results.append(res)
            n_done += 1
            elapsed = time.time() - t_start
            print(
                f"[{n_done}/{len(jobs)}] {res['status']} {res['receptor']} {res['ligand']} "
                f"score={res['vina_mode1']} ({res['seconds']}s) elapsed={elapsed:.0f}s",
                flush=True,
            )
            if n_done % 10 == 0:
                pd.DataFrame(results).to_csv(scores_path, index=False)

    pd.DataFrame(results).to_csv(scores_path, index=False)
    ok = sum(1 for r in results if r.get("status") == "success")
    print(f"DONE {cfg['pair_label']}: {ok}/{len(results)} success -> {scores_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
