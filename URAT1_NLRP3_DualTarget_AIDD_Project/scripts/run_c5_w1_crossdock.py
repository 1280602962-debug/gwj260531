#!/usr/bin/env python3
"""C5 W1: four-structure URAT1 cross-dock (Job A).

Order:
  1) Gate: benzbromarone @ 9DKA × seeds 42/43/44
     CNNscore-selected Top-1 RMSD vs 9DKA/R75 ≤ 2.0 Å for ALL seeds
  2) If gate passes: remaining new cells (32 total new; 4 reused)

Reuse (do not re-run):
  - lesinurad @ 9DKB seeds 42/43/44
  - benzbromarone @ 9DKB seed42

Settings locked to config/docking_c5_w1.yaml (exh=32, modes=9, rescore).
This machine has no GPU → always --no_gpu (same numeric settings).
Timeouts: skip molecule and continue (recorded as status=timeout).
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from parse_c1_sdf_readouts import parse_sdf_readouts, pose_rmsd, load_poses  # noqa: E402

SEEDS = [42, 43, 44]
TIMEOUT_SEC = 7200
RMSD_GATE = 2.0

LIGANDS = {
    "lesinurad": PROJECT_ROOT
    / "data/campaigns/c1/01_ligand_prep/forced_recovery/pdbqt/lesinurad.pdbqt",
    "benzbromarone": PROJECT_ROOT
    / "data/campaigns/c1/01_ligand_prep/forced_recovery/pdbqt/benzbromarone.pdbqt",
    "TD-3": PROJECT_ROOT
    / "data/campaigns/c5/01_ligand_prep/w1_refs/pdbqt/TD-3.pdbqt",
}

REFS = {
    "lesinurad": PROJECT_ROOT
    / "data/campaigns/c5/01_ligand_prep/w1_crystal_refs/lesinurad_A1AIL_crystal_ref.sdf",
    "benzbromarone": PROJECT_ROOT
    / "data/campaigns/c5/01_ligand_prep/w1_crystal_refs/benzbromarone_R75_crystal_ref.sdf",
    "TD-3": PROJECT_ROOT
    / "data/campaigns/c5/01_ligand_prep/w1_crystal_refs/TD3_A1A45_crystal_ref.sdf",
}

# receptor key in docking_c5_w1.yaml → short dir name
TARGET_DIR = {
    "urat1_9dkb": "9dkb",
    "urat1_9dka": "9dka",
    "urat1_9dkc": "9dkc",
    "urat1_9dk9": "9dk9",
}

REUSE = {
    ("lesinurad", "urat1_9dkb", 42): PROJECT_ROOT
    / "data/campaigns/c1/02_selfdock/urat1_9dkb/seed42/lesinurad_out.sdf",
    ("lesinurad", "urat1_9dkb", 43): PROJECT_ROOT
    / "data/campaigns/c1/02_selfdock/urat1_9dkb/seed43/lesinurad_out.sdf",
    ("lesinurad", "urat1_9dkb", 44): PROJECT_ROOT
    / "data/campaigns/c1/02_selfdock/urat1_9dkb/seed44/lesinurad_out.sdf",
    ("benzbromarone", "urat1_9dkb", 42): PROJECT_ROOT
    / "data/campaigns/c1/03_forced_recovery/urat1_9dkb/seed42/benzbromarone_out.sdf",
}


def out_path(ligand: str, target: str, seed: int) -> Path:
    return (
        PROJECT_ROOT
        / "data/campaigns/c5/01_crossdock"
        / TARGET_DIR[target]
        / f"seed{seed}"
        / f"{ligand}_out.sdf"
    )


def run_gnina(
    gnina: Path,
    receptor: Path,
    ligand: Path,
    center: list[float],
    size: list[float],
    out_sdf: Path,
    seed: int,
    cpu: int,
    timeout: int,
) -> str:
    """Return status: ok|timeout|fail|exists."""
    if out_sdf.exists() and out_sdf.stat().st_size > 0:
        return "exists"
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    log = out_sdf.with_suffix(".log")
    cmd = [
        str(gnina),
        "-r",
        str(receptor),
        "-l",
        str(ligand),
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(size[0]),
        "--size_y",
        str(size[1]),
        "--size_z",
        str(size[2]),
        "--exhaustiveness",
        "32",
        "--num_modes",
        "9",
        "--cpu",
        str(cpu),
        "--cnn_scoring",
        "rescore",
        "--seed",
        str(seed),
        "-o",
        str(out_sdf),
        "--log",
        str(log),
        "--no_gpu",
    ]
    print("RUN:", " ".join(cmd), flush=True)
    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout
        )
    except subprocess.TimeoutExpired:
        (out_sdf.parent / f"{out_sdf.stem}_TIMEOUT.txt").write_text(
            f"timeout_sec={timeout}\n"
        )
        print(f"TIMEOUT after {timeout}s: {out_sdf}", flush=True)
        return "timeout"
    (out_sdf.parent / f"{out_sdf.stem}_stdout.txt").write_text(
        (proc.stdout or "") + "\n" + (proc.stderr or "")
    )
    dt = time.time() - t0
    if proc.returncode != 0 or not out_sdf.exists() or out_sdf.stat().st_size == 0:
        print(f"FAIL rc={proc.returncode} dt={dt:.0f}s: {out_sdf}", flush=True)
        return "fail"
    print(f"OK dt={dt:.0f}s: {out_sdf}", flush=True)
    return "ok"


def ensure_reuse(ligand: str, target: str, seed: int) -> Path:
    src = REUSE[(ligand, target, seed)]
    dst = out_path(ligand, target, seed)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.exists():
        shutil.copy2(src, dst)
        (dst.parent / f"{dst.stem}_REUSED_FROM.txt").write_text(str(src) + "\n")
    return dst


def evaluate_cell(ligand: str, target: str, seed: int, sdf: Path) -> dict:
    ref = REFS[ligand]
    out = {
        "ligand": ligand,
        "target": target,
        "seed": seed,
        "sdf": str(sdf),
        "status": "ok" if sdf.exists() else "missing",
    }
    if not sdf.exists():
        return out
    try:
        readouts = parse_sdf_readouts(sdf)
        out.update(
            {
                "n_poses": readouts.get("n_poses"),
                "CNNscore_star": readouts.get("C1_P0"),
                "CNNaffinity_star": readouts.get("C1_P2star"),
            }
        )
        poses = load_poses(sdf)
        if not poses:
            out["rmsd_error"] = "no_poses"
            return out
        # CNNscore-selected pose index
        scores = []
        for i, m in enumerate(poses):
            try:
                scores.append((float(m.GetProp("CNNscore")), i))
            except Exception:
                continue
        if not scores:
            out["rmsd_error"] = "no_CNNscore"
            return out
        i_star = max(scores)[1]
        ref_poses = load_poses(ref)
        if not ref_poses:
            out["rmsd_error"] = "ref_empty"
            return out
        rmsd_star = pose_rmsd(poses[i_star], ref_poses[0])
        rmsds = []
        for p in poses:
            try:
                rmsds.append(pose_rmsd(p, ref_poses[0]))
            except Exception:
                rmsds.append(float("nan"))
        finite = [x for x in rmsds if x == x]
        out["rmsd_cnnscore_selected"] = float(rmsd_star)
        out["rmsd_best_of_n"] = float(min(finite)) if finite else None
        out["pass_rmsd_gate"] = bool(rmsd_star <= RMSD_GATE)
    except Exception as exc:
        out["rmsd_error"] = str(exc)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["gate", "rest", "all"], default="all")
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=TIMEOUT_SEC)
    args = ap.parse_args()

    cfg = yaml.safe_load((PROJECT_ROOT / "config/docking_c5_w1.yaml").read_text())
    gnina = PROJECT_ROOT / "tools" / "gnina"
    if not gnina.exists():
        raise SystemExit(f"gnina missing: {gnina}")

    # version note (worklist locks 1.3.1; local binary may differ)
    ver = subprocess.run(
        [str(gnina), "--version"], capture_output=True, text=True
    )
    version_line = (ver.stdout or ver.stderr or "").splitlines()[:1]
    meta = {
        "gnina_version_raw": version_line,
        "settings": {
            "exhaustiveness": 32,
            "num_modes": 9,
            "cnn_scoring": "rescore",
            "seeds": SEEDS,
            "no_gpu": True,
            "cpu": args.cpu,
            "timeout_sec": args.timeout,
            "rmsd_gate": RMSD_GATE,
            "pose_selection": "CNNscore",
        },
        "started_unix": time.time(),
    }
    out_root = PROJECT_ROOT / "data/campaigns/c5/01_crossdock"
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "run_meta.json").write_text(json.dumps(meta, indent=2))

    # materialize reused cells
    for key in REUSE:
        ensure_reuse(*key)

    gate_jobs = [("benzbromarone", "urat1_9dka", s) for s in SEEDS]
    all_pairs = [
        (lig, tgt, seed)
        for lig in ("lesinurad", "benzbromarone", "TD-3")
        for tgt in ("urat1_9dkb", "urat1_9dka", "urat1_9dkc", "urat1_9dk9")
        for seed in SEEDS
    ]
    new_jobs = [j for j in all_pairs if j not in REUSE]

    results: list[dict] = []

    def do_job(ligand: str, target: str, seed: int) -> dict:
        dst = out_path(ligand, target, seed)
        if (ligand, target, seed) in REUSE:
            status = "reused"
            ensure_reuse(ligand, target, seed)
        else:
            tcfg = cfg["targets"][target]
            status = run_gnina(
                gnina,
                PROJECT_ROOT / tcfg["prepared_receptor"],
                LIGANDS[ligand],
                tcfg["center"],
                tcfg["size"],
                dst,
                seed,
                args.cpu,
                args.timeout,
            )
        row = evaluate_cell(ligand, target, seed, dst)
        row["run_status"] = status
        print(
            f"EVAL {ligand}@{TARGET_DIR[target]} seed{seed}: "
            f"status={status} rmsd={row.get('rmsd_cnnscore_selected')} "
            f"pass={row.get('pass_rmsd_gate')}",
            flush=True,
        )
        return row

    if args.phase in ("gate", "all"):
        print("==== PHASE GATE: benzbromarone @ 9DKA × 3 ====", flush=True)
        gate_rows = []
        for job in gate_jobs:
            row = do_job(*job)
            gate_rows.append(row)
            results.append(row)
        gate_path = out_root / "gate_benzbromarone_9dka.json"
        gate_pass = all(
            r.get("run_status") in ("ok", "exists", "reused")
            and r.get("pass_rmsd_gate") is True
            for r in gate_rows
        )
        gate_path.write_text(
            json.dumps(
                {
                    "pass": gate_pass,
                    "rule": "CNNscore Top-1 RMSD ≤ 2.0 Å vs 9DKA/R75 for ALL seeds 42/43/44",
                    "rows": gate_rows,
                },
                indent=2,
            )
        )
        pd.DataFrame(gate_rows).to_csv(
            out_root / "gate_benzbromarone_9dka.csv", index=False
        )
        print(f"GATE PASS={gate_pass} → {gate_path}", flush=True)
        if not gate_pass:
            print(
                "GATE FAILED — stop. Do NOT start Task 2/3 or remaining W1.",
                flush=True,
            )
            (out_root / "STOP_GATE_FAILED").write_text(
                json.dumps(gate_rows, indent=2)
            )
            if args.phase == "all":
                raise SystemExit(2)
            return

    if args.phase in ("rest", "all"):
        # remaining new jobs excluding gate jobs already done
        rest = [
            j
            for j in new_jobs
            if j not in gate_jobs or args.phase == "rest"
        ]
        if args.phase == "all":
            rest = [j for j in new_jobs if j not in gate_jobs]
        print(f"==== PHASE REST: {len(rest)} jobs ====", flush=True)
        for job in rest:
            results.append(do_job(*job))

    # also evaluate reused cells for completeness
    for job in REUSE:
        if not any(
            r["ligand"] == job[0] and r["target"] == job[1] and r["seed"] == job[2]
            for r in results
        ):
            results.append(do_job(*job))

    df = pd.DataFrame(results)
    df.to_csv(out_root / "w1_crossdock_metrics.csv", index=False)
    (out_root / "w1_crossdock_metrics.json").write_text(
        json.dumps(results, indent=2)
    )
    print("W1 done →", out_root / "w1_crossdock_metrics.csv", flush=True)


if __name__ == "__main__":
    main()
