#!/usr/bin/env python3
"""Run the URAT1 D1 capability matrix locally.

27 cells = 3 ligands × 3 holo 9DK receptors × 3 seeds.
7 cells are reused; 20 are new GNINA jobs.
Benzbromarone Top-1 > 2 Å does not stop the matrix.
9DK9 apo is not included.

This does not redock the 156-library, reopen rank_track, or change 12/21.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from urat1_capability_lib import (  # noqa: E402
    D1_OUT,
    LIGAND_PDBQT,
    REUSE_SOURCES,
    SEEDS,
    TARGET_PDB,
    cell_sdf,
    d1_new_jobs,
    d1_selfdock,
    evaluate_sdf,
    inventory_rows,
    load_w1_engine_config,
)


def find_gnina(explicit: str | None) -> Path:
    if explicit:
        path = Path(explicit).expanduser()
        if not path.exists():
            raise SystemExit(f"gnina not found: {path}")
        return path
    env = os.environ.get("GNINA")
    if env:
        path = Path(env).expanduser()
        if path.exists():
            return path
    tools = PROJECT_ROOT / "tools" / "gnina"
    if tools.exists():
        return tools
    which = shutil.which("gnina")
    if which:
        return Path(which)
    raise SystemExit(
        "gnina not found. Put it on PATH, set GNINA, use --gnina, or place tools/gnina."
    )


def gnina_version(binary: Path) -> str:
    proc = subprocess.run([str(binary), "--version"], capture_output=True, text=True)
    text = (proc.stdout or proc.stderr or "").strip().splitlines()
    return text[0] if text else "unknown"


def run_gnina(
    binary: Path,
    receptor: Path,
    ligand: Path,
    center: list[float],
    size: list[float],
    out_sdf: Path,
    seed: int,
    cpu: int,
    timeout: int,
    no_gpu: bool,
) -> str:
    if out_sdf.exists() and out_sdf.stat().st_size > 0:
        return "exists"
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    log = out_sdf.with_suffix(".log")
    cmd = [
        str(binary),
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
    ]
    if no_gpu:
        cmd.append("--no_gpu")
    print("RUN:", " ".join(cmd), flush=True)
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        (out_sdf.parent / f"{out_sdf.stem}_TIMEOUT.txt").write_text(f"timeout_sec={timeout}\n")
        print(f"TIMEOUT {timeout}s {out_sdf}", flush=True)
        return "timeout"
    (out_sdf.parent / f"{out_sdf.stem}_stdout.txt").write_text(
        (proc.stdout or "") + "\n" + (proc.stderr or "")
    )
    dt = time.time() - t0
    if proc.returncode != 0 or not out_sdf.exists() or out_sdf.stat().st_size == 0:
        print(f"FAIL rc={proc.returncode} dt={dt:.0f}s {out_sdf}", flush=True)
        return "fail"
    print(f"OK dt={dt:.0f}s {out_sdf}", flush=True)
    return "ok"


def materialize_reuse(ligand: str, target: str, seed: int) -> Path:
    src = REUSE_SOURCES[(ligand, target, seed)]
    dst = cell_sdf(ligand, target, seed)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.exists():
        shutil.copy2(src, dst)
        (dst.parent / f"{dst.stem}_REUSED_FROM.txt").write_text(
            str(src.relative_to(PROJECT_ROOT)) + "\n"
        )
    return dst


def write_tables(rows: list[dict], out_root: Path) -> None:
    slim = []
    for row in rows:
        slim.append({k: v for k, v in row.items() if k not in {"per_mode", "reference"}})
        if "reference" in row:
            slim[-1]["reference_transform"] = row["reference"].get("transform")
            slim[-1]["reference_ca_rmsd_A"] = row["reference"].get("ca_rmsd_A")
    df = pd.DataFrame(slim)
    df.to_csv(out_root / "d1_capability_summary.csv", index=False)
    (out_root / "d1_capability_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    per_mode = []
    for row in rows:
        for mode in row.get("per_mode") or []:
            per_mode.append(
                {
                    "ligand": row["ligand"],
                    "pdb": row.get("pdb"),
                    "seed": row["seed"],
                    "self_dock": row.get("self_dock"),
                    **mode,
                }
            )
    if per_mode:
        pd.DataFrame(per_mode).to_csv(out_root / "d1_capability_per_mode.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true", help="Print the 20 new jobs and exit")
    ap.add_argument("--score-only", action="store_true", help="Evaluate existing SDFs; no GNINA")
    ap.add_argument("--gnina", default=None, help="GNINA binary (else GNINA env, tools/gnina, PATH)")
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=7200)
    ap.add_argument("--no-gpu", action="store_true", help="Pass --no_gpu to GNINA")
    ap.add_argument(
        "--only",
        action="append",
        default=[],
        help="Limit to ligand,target,seed (repeatable). target like urat1_9dkc",
    )
    args = ap.parse_args()

    D1_OUT.mkdir(parents=True, exist_ok=True)
    inv = inventory_rows()
    pd.DataFrame(inv).to_csv(D1_OUT / "d1_job_inventory.csv", index=False)
    (D1_OUT / "d1_job_inventory.json").write_text(json.dumps(inv, indent=2) + "\n")

    wanted = d1_new_jobs()
    if args.only:
        filters = set()
        for item in args.only:
            parts = [p.strip() for p in item.split(",")]
            if len(parts) != 3:
                raise SystemExit(f"--only needs ligand,target,seed; got {item}")
            filters.add((parts[0], parts[1], int(parts[2])))
        wanted = [j for j in wanted if j in filters]
        reuse_wanted = [j for j in REUSE_SOURCES if j in filters]
    else:
        reuse_wanted = list(REUSE_SOURCES)

    print(
        f"D1 inventory: {len(inv)} cells, {len(d1_new_jobs())} new, "
        f"{len(REUSE_SOURCES)} reuse; this invocation new={len(wanted)}",
        flush=True,
    )
    for lig, tgt, seed in wanted:
        print(
            f"  NEW {lig}@{TARGET_PDB[tgt]} seed{seed} "
            f"{'self' if d1_selfdock(lig, tgt) else 'cross'} -> {cell_sdf(lig, tgt, seed)}",
            flush=True,
        )
    if args.dry_run:
        return

    for key in reuse_wanted:
        materialize_reuse(*key)

    rows: list[dict] = []
    if not args.score_only:
        cfg = load_w1_engine_config()
        binary = find_gnina(args.gnina)
        meta = {
            "layer": "d1",
            "gnina": str(binary),
            "gnina_version": gnina_version(binary),
            "exhaustiveness": 32,
            "num_modes": 9,
            "cnn_scoring": "rescore",
            "seeds": list(SEEDS),
            "no_gpu": bool(args.no_gpu),
            "cpu": args.cpu,
            "timeout_sec": args.timeout,
            "native_like_cutoff_A": 2.0,
            "gate": "none_diagnostic_matrix_does_not_stop_on_top1",
            "started_unix": time.time(),
        }
        (D1_OUT / "run_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
        for lig, tgt, seed in wanted:
            tcfg = cfg["targets"][tgt]
            status = run_gnina(
                binary,
                PROJECT_ROOT / tcfg["prepared_receptor"],
                LIGAND_PDBQT[lig],
                tcfg["center"],
                tcfg["size"],
                cell_sdf(lig, tgt, seed),
                seed,
                args.cpu,
                args.timeout,
                args.no_gpu,
            )
            row = evaluate_sdf(lig, tgt, seed, cell_sdf(lig, tgt, seed))
            row["run_status"] = status
            rows.append(row)
            print(
                f"EVAL {lig}@{TARGET_PDB[tgt]} seed{seed}: {status} "
                f"top1={row.get('top1_pose_rmsd_A')} best9={row.get('best_of_9_pose_rmsd_A')} "
                f"rank={row.get('best_rmsd_cnnscore_rank')} class={row.get('capability_class')}",
                flush=True,
            )

    eval_keys = reuse_wanted + (wanted if args.score_only else [])
    if args.score_only:
        eval_keys = reuse_wanted + wanted
    seen = {(r["ligand"], r["target"], r["seed"]) for r in rows}
    for lig, tgt, seed in eval_keys:
        if (lig, tgt, seed) in seen:
            continue
        sdf = cell_sdf(lig, tgt, seed)
        if not sdf.exists() and (lig, tgt, seed) in REUSE_SOURCES:
            materialize_reuse(lig, tgt, seed)
        row = evaluate_sdf(lig, tgt, seed, sdf)
        row["run_status"] = "reused" if (lig, tgt, seed) in REUSE_SOURCES else "scored"
        rows.append(row)

    write_tables(rows, D1_OUT)
    print("Wrote", D1_OUT / "d1_capability_summary.csv", flush=True)


if __name__ == "__main__":
    main()
