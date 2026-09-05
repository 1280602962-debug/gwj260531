#!/usr/bin/env python3
"""Independent GNINA docking search for JAK1/TYK2 only.

Formulation-gap rule from the original paper: independent pose generation
was EGFR/HER2 and PIK3CA/mTOR. The only new pair that qualifies is
JAK1/TYK2. Do not run this on F2/F10, JAK1/JAK2, or either PPAR pair.

Not a Vina-pose rescore. Readout = mode-1 minimizedAffinity.
Seed 20260727, E=8, nine modes. Does not replace Table 2.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_track_b_v0"
TABLES = ROOT / "tables"
LIG_DIR = ROOT / "cache" / "track_b_ligands" / "pdbqt"
OUT = LOCAL / "independent_gnina_jak1_tyk2"
GNINA_DEFAULT = Path("/mnt/d/CADD paper exercise/gnina/bin/gnina")
GNINA_LIB_DEFAULT = Path("/mnt/d/CADD paper exercise/gnina/conda_env/lib")
SEED = 20260727
N_MODES = 9
EXHAUSTIVENESS = 8
PAIR = "JAK1/TYK2"
TARGETS = ["6N7A", "3LXP"]
PANEL = TABLES / "track_b_panels" / "panel_JAK1_TYK2_v1.csv"


def parse_mode1_affinity(out_pdbqt: Path) -> float | None:
    text = out_pdbqt.read_text().splitlines()
    in_model = False
    for line in text:
        if line.startswith("MODEL 1") or (line.startswith("MODEL") and "1" in line.split()[1:2]):
            in_model = True
        elif in_model and line.startswith("REMARK minimizedAffinity"):
            return float(line.split()[2])
        elif in_model and line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
        elif line.startswith("ENDMDL") and in_model:
            break
    return None


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


def load_box(target: str) -> dict:
    box = json.loads((LOCAL / "boxes" / f"{target}_box.json").read_text())
    if "center" in box:
        box["center_x"], box["center_y"], box["center_z"] = box["center"]
    if "size" in box:
        box["size_x"], box["size_y"], box["size_z"] = box["size"]
    return box


def receptor_for(target: str) -> Path:
    for name in (f"{target}_protein.pdb", f"{target}_receptor.pdb", f"{target}_receptor.pdbqt"):
        p = LOCAL / "receptors" / name
        if p.exists():
            return p
    raise SystemExit(f"missing receptor for {target}")


def dock_one(target: str, lig: str, gnina: Path, lib: Path) -> dict:
    t0 = time.time()
    pose_dir = OUT / "poses" / target / lig
    out_pdbqt = OUT / "logs" / f"{target}_{lig}_gnina_out.pdbqt"
    log_path = OUT / "logs" / f"{target}_{lig}.log"
    (OUT / "logs").mkdir(parents=True, exist_ok=True)
    if (pose_dir / "mode_01.pdbqt").exists():
        aff = parse_mode1_affinity(pose_dir / "mode_01.pdbqt")
        return {
            "pair": PAIR,
            "target": target,
            "ligand": lig,
            "gnina_mode1": aff,
            "score_S": None if aff is None else -aff,
            "n_modes": len(list(pose_dir.glob("mode_*.pdbqt"))),
            "status": "exists",
            "seconds": 0.0,
        }
    lig_pdbqt = LIG_DIR / f"{lig}.pdbqt"
    if not lig_pdbqt.exists():
        return {
            "pair": PAIR,
            "target": target,
            "ligand": lig,
            "gnina_mode1": None,
            "score_S": None,
            "n_modes": 0,
            "status": "fail",
            "reason": f"missing ligand {lig_pdbqt}",
            "seconds": 0.0,
        }
    box = load_box(target)
    rec = receptor_for(target)
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{lib}:{env.get('LD_LIBRARY_PATH', '')}"
    cmd = [
        str(gnina),
        "--no_gpu",
        "-r",
        str(rec),
        "-l",
        str(lig_pdbqt),
        "--center_x",
        str(box["center_x"]),
        "--center_y",
        str(box["center_y"]),
        "--center_z",
        str(box["center_z"]),
        "--size_x",
        str(box["size_x"]),
        "--size_y",
        str(box["size_y"]),
        "--size_z",
        str(box["size_z"]),
        "--exhaustiveness",
        str(EXHAUSTIVENESS),
        "--num_modes",
        str(N_MODES),
        "--seed",
        str(SEED),
        "--cpu",
        "1",
        "-o",
        str(out_pdbqt),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    log_path.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
    if proc.returncode != 0 or not out_pdbqt.exists():
        return {
            "pair": PAIR,
            "target": target,
            "ligand": lig,
            "gnina_mode1": None,
            "score_S": None,
            "n_modes": 0,
            "status": "fail",
            "reason": (proc.stderr or proc.stdout or "")[-400:],
            "seconds": round(time.time() - t0, 1),
        }
    n = split_modes(out_pdbqt, pose_dir)
    aff = parse_mode1_affinity(out_pdbqt)
    return {
        "pair": PAIR,
        "target": target,
        "ligand": lig,
        "gnina_mode1": aff,
        "score_S": None if aff is None else -aff,
        "n_modes": n,
        "status": "success",
        "seconds": round(time.time() - t0, 1),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--gnina", type=Path, default=GNINA_DEFAULT)
    ap.add_argument("--gnina-lib", type=Path, default=GNINA_LIB_DEFAULT)
    args = ap.parse_args()
    if not args.gnina.exists():
        raise SystemExit(f"GNINA not found: {args.gnina}")
    panel = list(csv.DictReader(PANEL.open()))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "logs").mkdir(parents=True, exist_ok=True)
    jobs = [(t, r["panel_id"]) for r in panel for t in TARGETS]
    scores_path = LOCAL / "tables" / "scores_gnina_independent_jak1_tyk2_v1.csv"
    print(f"{PAIR}: {len(jobs)} jobs workers={args.workers} (do not add other pairs)", flush=True)
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(dock_one, t, lig, args.gnina, args.gnina_lib): (t, lig) for t, lig in jobs}
        n = 0
        for fut in as_completed(futs):
            res = fut.result()
            results.append(res)
            n += 1
            print(
                f"[{n}/{len(jobs)}] {res['status']} {res['target']} {res['ligand']} "
                f"aff={res.get('gnina_mode1')}",
                flush=True,
            )
            if n % 10 == 0:
                _write(results, scores_path)
    _write(results, scores_path)
    ok = sum(1 for r in results if r.get("status") in ("success", "exists"))
    print(f"DONE {PAIR}: {ok}/{len(results)} -> {scores_path}", flush=True)
    return 0 if ok == len(results) else 1


def _write(rows, path: Path):
    if not rows:
        return
    fields = ["pair", "target", "ligand", "gnina_mode1", "score_S", "n_modes", "status", "seconds", "reason"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
