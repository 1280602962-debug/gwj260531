#!/usr/bin/env python3
"""Dock frozen MCL1/Bcl-xL panel96 with Vina on frozen 3WIY/3WIZ receptors.

- Membership: data/jcim_novelty_v0/tables/mcl1_bclxl_chembl_panel96_v1.csv (copied to tables/)
- Skip failed ligands; never stop the panel
- Pose cache: skip if poses/<pdb>/<panel_id>/mode_01.pdbqt exists
- Protocol: seed 20260727, E=8, modes=9, energy_range=3
"""
from __future__ import annotations

import csv
import json
import re
import subprocess
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/mcl1_bclxl_panel_v0"
TAB = OUT / "tables"
LIG_SDF = OUT / "ligands_sdf"
LIG_PQ = OUT / "ligands_pdbqt"
POSE = OUT / "poses"
LOG = OUT / "logs" / "vina"
CONF = OUT / "logs" / "vina_confs"
for d in (LIG_SDF, LIG_PQ, POSE, LOG, CONF):
    d.mkdir(parents=True, exist_ok=True)

PANEL = TAB / "mcl1_bclxl_chembl_panel96_v1.csv"
VINA = "/home/gwj/miniconda3/bin/vina"
PY = "/home/gwj/miniconda3/bin/python"
MK_LIG = "/home/gwj/miniconda3/bin/mk_prepare_ligand.py"
SEED = 20260727
EXHAUST = 8
N_MODES = 9
WORKERS = 6
TIMEOUT = 900

TARGETS = [
    {"name": "MCL1", "pdb": "3WIY"},
    {"name": "BCL2L1", "pdb": "3WIZ"},
]


def utc():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_box(pdb: str) -> dict:
    return json.loads((OUT / "boxes" / f"{pdb}_box.json").read_text())


def prep_ligand(panel_id: str, smiles: str) -> tuple[Path | None, str]:
    pdbqt = LIG_PQ / f"{panel_id}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt, "exists"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "bad_smiles"
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        if AllChem.EmbedMolecule(mol, randomSeed=SEED) != 0:
            return None, "embed_fail"
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    sdf = LIG_SDF / f"{panel_id}.sdf"
    w = Chem.SDWriter(str(sdf))
    mol.SetProp("_Name", panel_id)
    w.write(mol)
    w.close()
    proc = subprocess.run(
        [PY, MK_LIG, "-i", str(sdf), "-o", str(pdbqt)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not pdbqt.exists():
        return None, f"meeko_fail:{(proc.stderr or '')[-200:]}"
    return pdbqt, "ok"


def torsdof(pdbqt: Path) -> int:
    for line in pdbqt.read_text().splitlines():
        if line.startswith("TORSDOF"):
            return int(line.split()[1])
    return 0


def parse_mode1_affinity(log_text: str, out_pdbqt: Path) -> float | None:
    # Prefer REMARK VINA RESULT in out file
    if out_pdbqt.exists():
        for line in out_pdbqt.read_text().splitlines():
            if "VINA RESULT" in line:
                parts = line.split()
                try:
                    return float(parts[3])
                except Exception:
                    pass
            if line.startswith("MODEL"):
                break
    # fallback: table in log
    m = re.search(r"^\s*1\s+(-?\d+\.\d+)", log_text, re.M)
    if m:
        return float(m.group(1))
    return None


def split_modes(out_pdbqt: Path, dest: Path) -> int:
    dest.mkdir(parents=True, exist_ok=True)
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
        (dest / f"mode_{i:02d}.pdbqt").write_text("\n".join(m) + "\n")
    return len(models)


def dock_one(row: dict, tgt: dict) -> dict:
    panel_id = row["panel_id"]
    pdb = tgt["pdb"]
    pose_dir = POSE / pdb / panel_id
    mode1 = pose_dir / "mode_01.pdbqt"
    base = {
        "panel_id": panel_id,
        "molecule_chembl_id": row["molecule_chembl_id"],
        "class": row["class"],
        "target": tgt["name"],
        "pdb": pdb,
        "pchembl_MCL1": row["pchembl_MCL1"],
        "pchembl_BCL2L1": row["pchembl_BCL2L1"],
    }
    if mode1.exists():
        # recover affinity from log if present
        logp = LOG / f"{pdb}_{panel_id}.log"
        aff = None
        if logp.exists():
            aff = parse_mode1_affinity(logp.read_text(), LOG / f"{pdb}_{panel_id}_out.pdbqt")
        return {
            **base,
            "status": "exists",
            "vina_mode1": aff if aff is not None else "",
            "n_modes": len(list(pose_dir.glob("mode_*.pdbqt"))),
            "reason": "",
        }
    lig, st = prep_ligand(panel_id, row["smiles"])
    if lig is None:
        return {**base, "status": "prep_fail", "vina_mode1": "", "n_modes": 0, "reason": st}
    td = torsdof(lig)
    if td >= 28:
        return {
            **base,
            "status": "skip",
            "vina_mode1": "",
            "n_modes": 0,
            "reason": f"torsdof={td}>=28",
        }
    use_e = EXHAUST
    use_to = TIMEOUT
    if td >= 20:
        use_e = min(EXHAUST, 4)
        use_to = max(TIMEOUT, 1200)
    box = load_box(pdb)
    rec = OUT / "receptors" / f"{pdb}_receptor.pdbqt"
    out_pdbqt = LOG / f"{pdb}_{panel_id}_out.pdbqt"
    logp = LOG / f"{pdb}_{panel_id}.log"
    conf = CONF / f"{pdb}_{panel_id}.txt"
    conf.write_text(
        "\n".join(
            [
                f"receptor = {rec}",
                f"ligand = {lig}",
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                f"exhaustiveness = {use_e}",
                f"num_modes = {N_MODES}",
                "energy_range = 3",
                "cpu = 1",
                f"seed = {SEED}",
                f"out = {out_pdbqt}",
            ]
        )
        + "\n"
    )
    try:
        proc = subprocess.run(
            [VINA, "--config", str(conf)],
            capture_output=True,
            text=True,
            timeout=use_to,
        )
        logp.write_text(proc.stdout + "\n" + proc.stderr)
        if proc.returncode != 0 or not out_pdbqt.exists():
            return {
                **base,
                "status": "dock_fail",
                "vina_mode1": "",
                "n_modes": 0,
                "reason": (proc.stderr or proc.stdout)[-300:],
            }
        n = split_modes(out_pdbqt, pose_dir)
        aff = parse_mode1_affinity(proc.stdout, out_pdbqt)
        return {
            **base,
            "status": "success",
            "vina_mode1": aff if aff is not None else "",
            "n_modes": n,
            "reason": "",
            "torsdof": td,
            "exhaustiveness_used": use_e,
        }
    except subprocess.TimeoutExpired:
        return {
            **base,
            "status": "timeout",
            "vina_mode1": "",
            "n_modes": 0,
            "reason": f"timeout_{use_to}s",
        }
    except Exception as e:
        return {
            **base,
            "status": "error",
            "vina_mode1": "",
            "n_modes": 0,
            "reason": f"{type(e).__name__}:{e}",
        }


def write_scores(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    keys = [
        "panel_id",
        "molecule_chembl_id",
        "class",
        "target",
        "pdb",
        "pchembl_MCL1",
        "pchembl_BCL2L1",
        "status",
        "vina_mode1",
        "n_modes",
        "reason",
        "torsdof",
        "exhaustiveness_used",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    rows = list(csv.DictReader(PANEL.open()))
    jobs = [(r, t) for r in rows for t in TARGETS]
    print(f"jobs={len(jobs)} workers={WORKERS} started={utc()}", flush=True)
    results = []
    done = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(dock_one, r, t): (r["panel_id"], t["pdb"]) for r, t in jobs}
        for fut in as_completed(futs):
            pid, pdb = futs[fut]
            try:
                res = fut.result()
            except Exception:
                res = {
                    "panel_id": pid,
                    "pdb": pdb,
                    "status": "error",
                    "vina_mode1": "",
                    "n_modes": 0,
                    "reason": traceback.format_exc()[-300:],
                }
            results.append(res)
            done += 1
            if done % 10 == 0 or done == len(jobs):
                ok = sum(1 for r in results if r["status"] in {"success", "exists"})
                bad = sum(1 for r in results if r["status"] not in {"success", "exists"})
                print(
                    f"[{done}/{len(jobs)}] ok={ok} failed_or_skipped={bad} "
                    f"last={pid}@{pdb} status={res.get('status')}",
                    flush=True,
                )
                write_scores(results, TAB / "vina_scores_MBX_partial_v1.csv")
    write_scores(results, TAB / "vina_scores_MBX_v1.csv")
    skips = [r for r in results if r["status"] not in {"success", "exists"}]
    write_scores(skips, TAB / "vina_skips_MBX_v1.csv")
    meta = {
        "finished_utc": utc(),
        "n_jobs": len(jobs),
        "n_success_or_exists": sum(
            1 for r in results if r["status"] in {"success", "exists"}
        ),
        "n_fail": len(skips),
        "seed": SEED,
        "exhaustiveness": EXHAUST,
        "workers": WORKERS,
    }
    (TAB / "vina_run_meta_v1.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2), flush=True)


if __name__ == "__main__":
    main()
