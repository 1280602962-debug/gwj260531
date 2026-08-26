#!/usr/bin/env python3
"""Independent pose generation with GNINA 1.3.x docking search (not Vina-pose rescore).

Frozen protocol: docs/AGENT_COMMAND_INDEPENDENT_POSE_GENERATION_V1.md
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/jcim_independent_dock_v0"
GNINA = Path("/mnt/d/CADD paper exercise/gnina/bin/gnina")
GNINA_LIB = Path("/mnt/d/CADD paper exercise/gnina/conda_env/lib")
SEED = 20260727
N_MODES = 9

PANELS = {
    "EGFR/HER2": {
        "panel_root": ROOT / "results/egfr_her2_panel120_v0",
        "panel_csv": "tables/panel_v0_120.csv",
        "targets": ["3POZ", "3RCD"],
        "exhaustiveness": 8,
    },
    "PIK3CA/mTOR": {
        "panel_root": ROOT / "results/pik3ca_mtor_panel48_rdkit_v0",
        "panel_csv": "tables/panel_v0_48.csv",
        "targets": ["4L23", "4JT6"],
        "exhaustiveness": 16,
    },
}


def load_panel(panel_root: Path, panel_csv: str) -> list[dict]:
    with (panel_root / panel_csv).open() as fh:
        return list(csv.DictReader(fh))


def load_box(panel_root: Path, target: str) -> dict:
    box = json.loads((panel_root / "boxes" / f"{target}_box.json").read_text())
    if "center" in box:
        box["center_x"], box["center_y"], box["center_z"] = box["center"]
    if "size" in box:
        box["size_x"], box["size_y"], box["size_z"] = box["size"]
    return box


def parse_mode1_affinity(out_pdbqt: Path) -> float | None:
    text = out_pdbqt.read_text().splitlines()
    in_model = False
    for line in text:
        if line.startswith("MODEL 1"):
            in_model = True
        elif in_model and line.startswith("REMARK minimizedAffinity"):
            return float(line.split()[2])
        elif in_model and line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
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


def dock_one(
    pair: str,
    target: str,
    lig_id: str,
    panel_root: Path,
    exhaustiveness: int,
) -> dict:
    t0 = time.time()
    lig_pdbqt = panel_root / "ligands_pdbqt" / f"{lig_id}.pdbqt"
    receptor = panel_root / "receptors" / f"{target}_protein.pdb"
    if not receptor.exists():
        receptor = panel_root / "receptors" / f"{target}_receptor.pdbqt"
    box = load_box(panel_root, target)
    pose_dir = OUT / "poses" / pair.replace("/", "_") / target / lig_id
    out_pdbqt = OUT / "logs/gnina_dock" / f"{target}_{lig_id}_gnina_out.pdbqt"
    log_path = OUT / "logs/gnina_dock" / f"{target}_{lig_id}.log"

    if (pose_dir / "mode_01.pdbqt").exists():
        aff = parse_mode1_affinity(pose_dir / "mode_01.pdbqt")
        return {
            "pair": pair,
            "target": target,
            "ligand": lig_id,
            "gnina_mode1": aff,
            "n_modes": len(list(pose_dir.glob("mode_*.pdbqt"))),
            "status": "exists",
            "seconds": 0.0,
        }

    if not lig_pdbqt.exists():
        return {
            "pair": pair,
            "target": target,
            "ligand": lig_id,
            "gnina_mode1": None,
            "n_modes": 0,
            "status": "fail",
            "reason": f"missing ligand {lig_pdbqt}",
            "seconds": round(time.time() - t0, 1),
        }

    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{GNINA_LIB}:{env.get('LD_LIBRARY_PATH', '')}"
    cmd = [
        str(GNINA),
        "--no_gpu",
        "-r",
        str(receptor),
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
        str(exhaustiveness),
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
            "pair": pair,
            "target": target,
            "ligand": lig_id,
            "gnina_mode1": None,
            "n_modes": 0,
            "status": "fail",
            "reason": (proc.stderr or proc.stdout)[-400:],
            "seconds": round(time.time() - t0, 1),
        }

    n_modes = split_modes(out_pdbqt, pose_dir)
    aff = parse_mode1_affinity(out_pdbqt)
    return {
        "pair": pair,
        "target": target,
        "ligand": lig_id,
        "gnina_mode1": aff,
        "n_modes": n_modes,
        "status": "success",
        "seconds": round(time.time() - t0, 1),
    }


def write_scores(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def run_panel(pair: str, cfg: dict, workers: int, limit: int | None) -> None:
    panel_root = cfg["panel_root"]
    panel = load_panel(panel_root, cfg["panel_csv"])
    if limit:
        panel = panel[:limit]
    scores_path = OUT / "tables" / f"gnina_dock_scores_{pair.replace('/', '_')}.csv"
    done = {}
    if scores_path.exists():
        with scores_path.open() as fh:
            for r in csv.DictReader(fh):
                if r.get("status") in ("success", "exists"):
                    done[(r["target"], r["ligand"])] = r

    jobs = []
    for row in panel:
        lig = row.get("panel_id") or row.get("ligand")
        for target in cfg["targets"]:
            if (target, lig) in done:
                continue
            jobs.append((pair, target, lig, panel_root, cfg["exhaustiveness"]))

    print(f"{pair}: {len(jobs)} new jobs, {len(done)} cached, workers={workers}", flush=True)
    results = list(done.values()) if done else []
    # reload all cached rows properly
    if scores_path.exists():
        with scores_path.open() as fh:
            results = list(csv.DictReader(fh))

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {
            ex.submit(dock_one, *j): j for j in jobs
        }
        n = 0
        for fut in as_completed(futs):
            res = fut.result()
            # replace or append
            key = (res["target"], res["ligand"])
            results = [r for r in results if (r["target"], r["ligand"]) != key]
            results.append(res)
            n += 1
            print(
                f"[{n}/{len(jobs)}] {res['status']} {res['target']} {res['ligand']} "
                f"aff={res.get('gnina_mode1')} ({res.get('seconds', 0)}s) "
                f"elapsed={time.time()-t0:.0f}s",
                flush=True,
            )
            if n % 5 == 0:
                write_scores(results, scores_path)

    write_scores(results, scores_path)
    ok = sum(1 for r in results if r.get("status") in ("success", "exists"))
    print(f"DONE {pair}: {ok}/{len(results)} -> {scores_path}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", choices=list(PANELS) + ["all"], default="all")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=None, help="debug: first N ligands only")
    args = ap.parse_args()
    if not GNINA.exists():
        raise SystemExit(f"GNINA not found: {GNINA}")
    pairs = list(PANELS) if args.pair == "all" else [args.pair]
    for pair in pairs:
        run_panel(pair, PANELS[pair], args.workers, args.limit)


if __name__ == "__main__":
    main()
