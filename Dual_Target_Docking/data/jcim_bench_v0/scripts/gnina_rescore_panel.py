#!/usr/bin/env python3
"""GNINA CNN rescore of existing Vina poses (best-of-K over modes)."""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SEED = 20260727
GNINA_ROOT = Path("/mnt/d/CADD paper exercise/gnina")
GNINA_BIN = GNINA_ROOT / "bin" / "gnina"
OBABEL = Path("/home/gwj/miniconda3/envs/cadd_tools/bin/obabel")
if not OBABEL.exists():
    OBABEL = Path("/home/gwj/miniconda3/bin/obabel")


def env_for_gnina():
    env = os.environ.copy()
    lib = str(GNINA_ROOT / "conda_env" / "lib")
    env["LD_LIBRARY_PATH"] = lib + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    env["PATH"] = str(GNINA_ROOT / "bin") + ":" + env.get("PATH", "")
    return env


def pdbqt_to_sdf(pdbqt: Path, sdf: Path) -> bool:
    sdf.parent.mkdir(parents=True, exist_ok=True)
    if sdf.exists() and sdf.stat().st_size > 0:
        return True
    r = subprocess.run(
        [str(OBABEL), str(pdbqt), "-O", str(sdf)],
        capture_output=True,
        text=True,
    )
    return r.returncode == 0 and sdf.exists() and sdf.stat().st_size > 0


def parse_gnina_sdf(sdf: Path):
    """Return best CNNscore / CNNaffinity / minimizedAffinity from multi-pose sdf."""
    if not sdf.exists():
        return None
    text = sdf.read_text(errors="ignore")
    # split molecules by $$$$
    mols = text.split("$$$$")
    best = None
    for block in mols:
        if not block.strip():
            continue
        props = {}
        lines = block.splitlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith("> <"):
                key = lines[i][3:].strip().strip(">")
                val = lines[i + 1].strip() if i + 1 < len(lines) else ""
                props[key] = val
                i += 2
            else:
                i += 1
        cnn = None
        for k in ("CNNscore", "CNN_VS", "CNNaffinity", "minimizedAffinity"):
            if k in props:
                try:
                    cnn = float(props[k])
                    break
                except ValueError:
                    pass
        # prefer CNNscore if present
        score = None
        aff = None
        try:
            if "CNNscore" in props:
                score = float(props["CNNscore"])
            if "CNNaffinity" in props:
                aff = float(props["CNNaffinity"])
            if aff is None and "minimizedAffinity" in props:
                aff = float(props["minimizedAffinity"])
        except ValueError:
            continue
        if score is None and aff is None:
            continue
        # higher CNNscore better; for affinity gnina CNNaffinity higher=better usually
        key = score if score is not None else aff
        if best is None or key > best[0]:
            best = (key, score, aff)
    if best is None:
        return None
    return {"cnn_score": best[1], "cnn_affinity": best[2]}


def rescore_one(args):
    root, target, lig, receptor, mode_path, out_dir, timeout = args
    out_sdf = out_dir / f"{target}_{lig}_{mode_path.stem}.sdf"
    log = out_dir / f"{target}_{lig}_{mode_path.stem}.log"
    if out_sdf.exists() and out_sdf.stat().st_size > 0:
        parsed = parse_gnina_sdf(out_sdf)
        return {
            "target": target,
            "ligand": lig,
            "mode": mode_path.stem,
            "status": "exists",
            **(parsed or {}),
        }
    tmp_sdf = out_dir / "tmp_in" / f"{target}_{lig}_{mode_path.stem}.sdf"
    if not pdbqt_to_sdf(mode_path, tmp_sdf):
        return {
            "target": target,
            "ligand": lig,
            "mode": mode_path.stem,
            "status": "fail_convert",
        }
    cmd = [
        str(GNINA_BIN),
        "--no_gpu",
        "-r",
        str(receptor),
        "-l",
        str(tmp_sdf),
        "--cnn_scoring",
        "rescore",
        "--minimize",
        "--seed",
        str(SEED),
        "--cpu",
        "1",
        "-o",
        str(out_sdf),
    ]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env_for_gnina(),
        )
        log.write_text(proc.stdout + "\n" + proc.stderr)
        if proc.returncode != 0 or not out_sdf.exists():
            return {
                "target": target,
                "ligand": lig,
                "mode": mode_path.stem,
                "status": "fail",
                "reason": (proc.stderr or proc.stdout)[-200:],
            }
    except subprocess.TimeoutExpired:
        log.write_text(f"TIMEOUT {timeout}s\n")
        return {
            "target": target,
            "ligand": lig,
            "mode": mode_path.stem,
            "status": "timeout",
        }
    parsed = parse_gnina_sdf(out_sdf)
    return {
        "target": target,
        "ligand": lig,
        "mode": mode_path.stem,
        "status": "success",
        **(parsed or {}),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--targets", nargs="+", required=True)
    ap.add_argument(
        "--receptor-map",
        nargs="+",
        required=True,
        help="TARGET=/path/to/protein.pdb",
    )
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--modes", default="all", help="all or mode_01")
    args = ap.parse_args()
    root = Path(args.root)
    rec = {}
    for item in args.receptor_map:
        k, v = item.split("=", 1)
        rec[k] = Path(v)

    out_dir = root / "logs" / "gnina"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "tmp_in").mkdir(exist_ok=True)

    jobs = []
    for target in args.targets:
        pose_root = root / "poses" / target
        if not pose_root.exists():
            print("missing poses", pose_root, flush=True)
            continue
        for lig_dir in sorted(pose_root.iterdir()):
            if not lig_dir.is_dir():
                continue
            lig = lig_dir.name
            modes = sorted(lig_dir.glob("mode_*.pdbqt"))
            if args.modes == "mode_01":
                modes = [m for m in modes if m.name == "mode_01.pdbqt"]
            for mp in modes:
                jobs.append((root, target, lig, rec[target], mp, out_dir, args.timeout))

    print(f"GNINA jobs={len(jobs)} workers={args.workers}", flush=True)
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(rescore_one, j): j for j in jobs}
        done = 0
        for fut in as_completed(futs):
            res = fut.result()
            results.append(res)
            done += 1
            if done % 25 == 0 or res["status"] not in ("success", "exists"):
                print(
                    f"[{done}/{len(jobs)}] {res['status']} {res['target']} {res['ligand']} {res['mode']}",
                    flush=True,
                )

    long_csv = root / "tables" / "scores_gnina_long.csv"
    with long_csv.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "target",
                "ligand",
                "mode",
                "status",
                "cnn_score",
                "cnn_affinity",
                "reason",
            ],
        )
        w.writeheader()
        for r in results:
            w.writerow(
                {
                    "target": r.get("target"),
                    "ligand": r.get("ligand"),
                    "mode": r.get("mode"),
                    "status": r.get("status"),
                    "cnn_score": r.get("cnn_score", ""),
                    "cnn_affinity": r.get("cnn_affinity", ""),
                    "reason": r.get("reason", ""),
                }
            )
    print("wrote", long_csv, flush=True)

    # best-of-K per ligand/target by cnn_score (fallback cnn_affinity)
    best = {}
    for r in results:
        if r.get("status") not in ("success", "exists"):
            continue
        key = (r["ligand"], r["target"])
        score = r.get("cnn_score")
        aff = r.get("cnn_affinity")
        rank = score if score is not None else aff
        if rank is None:
            continue
        if key not in best or rank > best[key][0]:
            best[key] = (rank, score, aff)

    # wide table
    ligands = sorted({r["ligand"] for r in results})
    wide_rows = []
    for lig in ligands:
        row = {"ligand": lig}
        for t in args.targets:
            b = best.get((lig, t))
            if b:
                row[f"gnina_cnn_{t}"] = b[1] if b[1] is not None else ""
                row[f"gnina_aff_{t}"] = b[2] if b[2] is not None else ""
            else:
                row[f"gnina_cnn_{t}"] = ""
                row[f"gnina_aff_{t}"] = ""
        wide_rows.append(row)
    wide = root / "tables" / "scores_gnina_best.csv"
    if wide_rows:
        with wide.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(wide_rows[0].keys()))
            w.writeheader()
            w.writerows(wide_rows)
        print("wrote", wide, flush=True)

    ok = sum(1 for r in results if r["status"] in ("success", "exists"))
    print(f"done ok={ok}/{len(results)}", flush=True)
    return 0 if ok > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
