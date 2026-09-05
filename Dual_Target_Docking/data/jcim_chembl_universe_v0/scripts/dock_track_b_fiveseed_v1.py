#!/usr/bin/env python3
"""Five-seed Vina for the five post-census pairs.

Reuses production seed 20260727 scores. Redocks 20260811–20260814 into
local_track_b_v0/multiseed/{seed}/. Ligands stay the frozen ETKDG-20260727
Meeko PDBQTs. Does not overwrite production poses/ or Table 2.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import dock_track_b_production_v1 as prod  # noqa: E402

PRIMARY = 20260727
NEW_SEEDS = [20260811, 20260812, 20260813, 20260814]
ALL_SEEDS = [PRIMARY] + NEW_SEEDS
LOCAL = prod.LOCAL
VINA_DEFAULT = os.environ.get("VINA", prod.VINA)


def write_conf(target, lig, ligand_pdbqt, box, exhaust, seed, pose_root, log_root):
    conf_dir = log_root / "vina_confs"
    conf_dir.mkdir(parents=True, exist_ok=True)
    (log_root / "vina").mkdir(parents=True, exist_ok=True)
    conf = conf_dir / f"{target}_{lig}.txt"
    rec = LOCAL / "receptors" / f"{target}_receptor.pdbqt"
    out = log_root / "vina" / f"{target}_{lig}_out.pdbqt"
    conf.write_text(
        "\n".join(
            [
                f"receptor = {rec}",
                f"ligand = {ligand_pdbqt}",
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                f"exhaustiveness = {exhaust}",
                f"num_modes = {prod.N_MODES}",
                f"energy_range = {prod.ENERGY_RANGE}",
                "cpu = 1",
                f"seed = {seed}",
                f"out = {out}",
            ]
        )
        + "\n"
    )
    return conf, out


def run_one(pair, target, lig, ligand_pdbqt, box, exhaust, timeout_s, seed, pose_root, log_root, vina):
    pose_dir = pose_root / target / lig
    mode1 = pose_dir / "mode_01.pdbqt"
    if mode1.exists() and mode1.stat().st_size > 0:
        e = prod.parse_mode1_energy(mode1)
        return {
            "pair": pair,
            "target": target,
            "ligand": lig,
            "seed": seed,
            "status": "exists",
            "n_modes": len(list(pose_dir.glob("mode_*.pdbqt"))),
            "mode1_energy": e,
            "score_S": None if e is None else -e,
            "reason": "",
            "seconds": 0.0,
        }
    if not ligand_pdbqt.exists() or ligand_pdbqt.stat().st_size == 0:
        return {
            "pair": pair,
            "target": target,
            "ligand": lig,
            "seed": seed,
            "status": "fail",
            "n_modes": 0,
            "mode1_energy": None,
            "score_S": None,
            "reason": "missing_ligand_pdbqt",
            "seconds": 0.0,
        }
    td = prod.torsdof(ligand_pdbqt)
    if td >= prod.SKIP_TORSDOF_GE:
        return {
            "pair": pair,
            "target": target,
            "ligand": lig,
            "seed": seed,
            "status": "skip",
            "n_modes": 0,
            "mode1_energy": None,
            "score_S": None,
            "reason": f"skip_torsdof={td}_ge_{prod.SKIP_TORSDOF_GE}",
            "seconds": 0.0,
        }
    conf, out = write_conf(target, lig, ligand_pdbqt, box, exhaust, seed, pose_root, log_root)
    log = log_root / "vina" / f"{target}_{lig}.log"
    t0 = time.time()
    import subprocess

    try:
        proc = subprocess.run(
            [vina, "--config", str(conf)],
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
        log.write_text(proc.stdout + "\n" + proc.stderr)
        rc = proc.returncode
        err = (proc.stderr or proc.stdout or "")[-300:]
    except subprocess.TimeoutExpired as e:
        elapsed = time.time() - t0
        log.write_text(f"TIMEOUT after {timeout_s}s seed={seed} torsdof={td}\n")
        if out.exists():
            out.unlink()
        return {
            "pair": pair,
            "target": target,
            "ligand": lig,
            "seed": seed,
            "status": "skip",
            "n_modes": 0,
            "mode1_energy": None,
            "score_S": None,
            "reason": f"timeout_{timeout_s}s_torsdof={td}",
            "seconds": round(elapsed, 1),
        }
    elapsed = time.time() - t0
    if rc != 0 or not out.exists() or out.stat().st_size == 0:
        return {
            "pair": pair,
            "target": target,
            "ligand": lig,
            "seed": seed,
            "status": "fail",
            "n_modes": 0,
            "mode1_energy": None,
            "score_S": None,
            "reason": err or f"vina_rc={rc}",
            "seconds": round(elapsed, 1),
        }
    n = prod.split_modes(out, pose_dir)
    e = prod.parse_mode1_energy(out)
    return {
        "pair": pair,
        "target": target,
        "ligand": lig,
        "seed": seed,
        "status": "success",
        "n_modes": n,
        "mode1_energy": e,
        "score_S": None if e is None else -e,
        "reason": "",
        "seconds": round(elapsed, 1),
    }


def copy_primary():
    src = LOCAL / "tables" / "scores_vina_mode1_v1.csv"
    dest_dir = LOCAL / "tables" / "multiseed"
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"scores_vina_mode1_seed{PRIMARY}.csv"
    if not src.exists():
        raise SystemExit(f"missing primary scores: {src}")
    rows = list(csv.DictReader(src.open()))
    with dest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pair", "target", "ligand", "seed", "mode1_energy", "score_S", "status"])
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "pair": r["pair"],
                    "target": r["target"],
                    "ligand": r["ligand"],
                    "seed": PRIMARY,
                    "mode1_energy": r["mode1_energy"],
                    "score_S": r["score_S"],
                    "status": r.get("status", "reused_primary"),
                }
            )
    print(f"reused primary {src} -> {dest} n={len(rows)}", flush=True)
    return dest


def write_seed_tables(results, seed):
    dest_dir = LOCAL / "tables" / "multiseed"
    dest_dir.mkdir(parents=True, exist_ok=True)
    status_path = dest_dir / f"job_status_seed{seed}.csv"
    fields = ["pair", "target", "ligand", "seed", "status", "n_modes", "mode1_energy", "score_S", "reason", "seconds"]
    with status_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in sorted(results, key=lambda x: (x["pair"], x["target"], x["ligand"])):
            w.writerow({k: r.get(k, "") for k in fields})
    score_path = dest_dir / f"scores_vina_mode1_seed{seed}.csv"
    with score_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pair", "target", "ligand", "seed", "mode1_energy", "score_S", "status"])
        w.writeheader()
        for r in sorted(results, key=lambda x: (x["pair"], x["target"], x["ligand"])):
            if r["status"] in ("success", "exists") and r.get("mode1_energy") is not None:
                w.writerow(
                    {
                        "pair": r["pair"],
                        "target": r["target"],
                        "ligand": r["ligand"],
                        "seed": seed,
                        "mode1_energy": r["mode1_energy"],
                        "score_S": r["score_S"],
                        "status": r["status"],
                    }
                )
    return score_path


def run_seed(seed, workers, timeout, exhaust, vina, pair_filter):
    pose_root = LOCAL / "multiseed" / str(seed) / "poses"
    log_root = LOCAL / "multiseed" / str(seed) / "logs"
    pose_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)
    jobs = prod.load_jobs(pair_filter)
    print(f"seed={seed} jobs={len(jobs)} workers={workers} vina={vina}", flush=True)
    results = []
    t_all = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {
            ex.submit(
                run_one,
                pair,
                t,
                lig,
                pdbqt,
                box,
                exhaust,
                timeout,
                seed,
                pose_root,
                log_root,
                vina,
            ): (pair, t, lig)
            for pair, t, lig, pdbqt, box in jobs
        }
        done = 0
        for fut in as_completed(futs):
            res = fut.result()
            results.append(res)
            done += 1
            if done % 25 == 0 or done == len(jobs) or res["status"] not in ("success", "exists"):
                print(
                    f"[seed {seed} {done}/{len(jobs)}] {res['status']} {res['pair']} {res['target']} {res['ligand']}"
                    + (f" reason={res.get('reason', '')[:80]}" if res.get("reason") else ""),
                    flush=True,
                )
            if done % 100 == 0:
                write_seed_tables(results, seed)
    path = write_seed_tables(results, seed)
    ok = sum(1 for r in results if r["status"] in ("success", "exists"))
    print(
        f"seed {seed} done ok={ok}/{len(results)} elapsed_h={(time.time() - t_all) / 3600:.2f} -> {path}",
        flush=True,
    )
    return 0 if all(r["status"] != "fail" for r in results) else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--exhaustiveness", type=int, default=prod.EXHAUSTIVENESS)
    ap.add_argument("--vina", default=VINA_DEFAULT)
    ap.add_argument("--pairs", nargs="*", default=None)
    ap.add_argument("--seeds", nargs="*", type=int, default=None, help="default: all five; 20260727 is reuse-only")
    args = ap.parse_args()
    seeds = args.seeds or ALL_SEEDS
    pair_filter = set(args.pairs) if args.pairs else None
    rc = 0
    if PRIMARY in seeds:
        copy_primary()
    for seed in seeds:
        if seed == PRIMARY:
            continue
        if seed not in NEW_SEEDS:
            raise SystemExit(f"unexpected seed {seed}; locked list is {ALL_SEEDS}")
        rc |= run_seed(seed, args.workers, args.timeout, args.exhaustiveness, args.vina, pair_filter)
    print("five-seed pack finished; do not paste these AUROCs into Table 2", flush=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
