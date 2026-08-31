#!/usr/bin/env python3
"""B5: dock the frozen PIK3CA/PIK3CB panel into alternate PIK3CA crystals.

Uses AutoDock Vina 1.2.7 CLI (same as the main-panel dock_panel.py), not the
Python vina binding. Alternate receptors 4JPS / 5DXT are already prepared in
this pack (cognate QC PASS). Pocket B keeps frozen 2WXF scores.

Exhaustiveness = 8 to match the PIK3CA/PIK3CB main panel.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_structure_robust_v0"
PANEL_SCORES = ROOT / "data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv"
SEED = 20260727
EXHAUST = 8
N_MODES = 9
N_BOOT = 2000
BOOT_SEED = 20260729
VINA = shutil.which("vina") or "/home/gwj/miniconda3/bin/vina"
MAIN_SUMMARY_MIN = 0.500  # Table 2 θ = 6.0 pocket-matched Vina, PIK3CA/PIK3CB

LIG_DIR_CANDIDATES = [
    Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_pik3cb_panel_v0/ligands_pdbqt"),
    ROOT / "data/pik3ca_pik3cb_panel_v0/ligands_pdbqt",
]


def find_lig_dir() -> Path:
    for p in LIG_DIR_CANDIDATES:
        if p.is_dir() and any(p.glob("PAB_*.pdbqt")):
            return p
    raise FileNotFoundError("PIK3CA/PIK3CB ligand PDBQTs not found")


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def directional(recs, kda, kdb):
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    da = auroc([r[kda] for r in D], [r[kda] for r in A])
    db = auroc([r[kdb] for r in D], [r[kdb] for r in B])
    return da, db


def boot_min(recs, kda, kdb, n_boot=N_BOOT, seed=BOOT_SEED):
    usable = [r for r in recs if r["cls"] in ("dual", "A_only", "B_only")]
    rng = np.random.default_rng(seed)
    idx = np.arange(len(usable))
    mins = []
    for _ in range(n_boot):
        sub = [usable[i] for i in rng.choice(idx, size=len(idx), replace=True)]
        da, db = directional(sub, kda, kdb)
        if da == da and db == db:
            mins.append(min(da, db))
    lo, hi = np.percentile(mins, [2.5, 97.5])
    return float(lo), float(hi)


def parse_mode1(out_pdbqt: Path) -> float | None:
    if not out_pdbqt.exists():
        return None
    for line in out_pdbqt.read_text(errors="ignore").splitlines():
        if "VINA RESULT" in line:
            parts = line.split()
            for tok in parts:
                try:
                    return float(tok)
                except ValueError:
                    continue
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


def write_conf(conf_path: Path, receptor: Path, ligand: Path, box: dict, out_pdbqt: Path) -> None:
    conf_path.parent.mkdir(parents=True, exist_ok=True)
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    conf_path.write_text(
        "\n".join(
            [
                f"receptor = {receptor}",
                f"ligand = {ligand}",
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                f"exhaustiveness = {EXHAUST}",
                f"num_modes = {N_MODES}",
                "energy_range = 3",
                "cpu = 1",
                f"seed = {SEED}",
                f"out = {out_pdbqt}",
            ]
        )
        + "\n"
    )


def dock_one(args):
    receptor, box_path, lig_id, lig_pdbqt, pose_dir, timeout = args
    pose_dir = Path(pose_dir)
    pose_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = pose_dir / "out.pdbqt"
    log = pose_dir / "vina.log"
    t0 = time.time()
    if (pose_dir / "mode_01.pdbqt").exists() and out_pdbqt.exists():
        aff = parse_mode1(out_pdbqt)
        return {
            "ligand": lig_id,
            "vina_mode1": aff,
            "status": "exists" if aff is not None else "fail",
            "reason": "" if aff is not None else "exists_but_unparseable",
            "seconds": round(time.time() - t0, 1),
        }
    box = json.loads(Path(box_path).read_text())
    conf = pose_dir / "vina.conf"
    write_conf(conf, Path(receptor), Path(lig_pdbqt), box, out_pdbqt)
    try:
        proc = subprocess.run(
            [VINA, "--config", str(conf)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        log.write_text(proc.stdout + "\n" + proc.stderr)
        if proc.returncode != 0 or not out_pdbqt.exists():
            return {
                "ligand": lig_id,
                "vina_mode1": None,
                "status": "fail",
                "reason": (proc.stderr or proc.stdout)[-300:],
                "seconds": round(time.time() - t0, 1),
            }
    except subprocess.TimeoutExpired:
        log.write_text(f"TIMEOUT {timeout}s\n")
        return {
            "ligand": lig_id,
            "vina_mode1": None,
            "status": "timeout",
            "reason": f"timeout_{timeout}s",
            "seconds": round(time.time() - t0, 1),
        }
    split_modes(out_pdbqt, pose_dir)
    aff = parse_mode1(out_pdbqt)
    return {
        "ligand": lig_id,
        "vina_mode1": aff,
        "status": "success" if aff is not None else "fail",
        "reason": "" if aff is not None else "unparseable_affinity",
        "seconds": round(time.time() - t0, 1),
    }


def write_metrics(alt: str, results: list[dict]) -> dict:
    panel = pd.read_csv(PANEL_SCORES)
    alt_map = {
        r["ligand"]: r
        for r in results
        if r.get("status") in ("success", "exists") and r.get("vina_mode1") is not None
    }
    recs = []
    for _, row in panel.iterrows():
        lig = row["ligand"]
        if lig not in alt_map:
            continue
        b_fr = row.get("vina_PIK3CB")
        if pd.isna(b_fr):
            continue
        a = float(alt_map[lig]["vina_mode1"])
        b = float(b_fr)
        recs.append({"ligand": lig, "cls": row["class"], "vina_A": -a, "vina_B": -b, "alt_receptor": alt})
    da, db = directional(recs, "vina_B", "vina_A")
    sm = min(da, db) if da == da and db == db else float("nan")
    lo, hi = boot_min(recs, "vina_B", "vina_A") if recs else (float("nan"), float("nan"))
    metric = {
        "pair": "PIK3CA/PIK3CB",
        "alt_receptor": alt,
        "replaced_pocket": "A",
        "kept_pocket": "B=2WXF frozen",
        "exhaustiveness": EXHAUST,
        "n": len(recs),
        "n_dual": sum(r["cls"] == "dual" for r in recs),
        "n_A_only": sum(r["cls"] == "A_only" for r in recs),
        "n_B_only": sum(r["cls"] == "B_only" for r in recs),
        "auroc_D_vs_A": None if da != da else round(da, 4),
        "auroc_D_vs_B": None if db != db else round(db, 4),
        "summary_min": None if sm != sm else round(sm, 4),
        "summary_min_ci_lo": None if lo != lo else round(lo, 4),
        "summary_min_ci_hi": None if hi != hi else round(hi, 4),
        "main_panel_summary_min": MAIN_SUMMARY_MIN,
        "delta_vs_main": None if sm != sm else round(sm - MAIN_SUMMARY_MIN, 4),
        "note": "pocket A replaced by alt PIK3CA crystal; pocket B = frozen 2WXF scores",
    }
    metric_path = OUT / "tables" / f"pocket_matched_PAB_alt{alt}_v1.csv"
    pd.DataFrame([metric]).to_csv(metric_path, index=False)
    lig_path = OUT / "tables" / f"ligand_scores_PAB_alt{alt}_v1.csv"
    pd.DataFrame(recs).to_csv(lig_path, index=False)
    print(metric, flush=True)
    print(f"wrote {metric_path}", flush=True)
    return metric


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alt", required=True, choices=["4JPS", "5DXT"])
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--analyze-only", action="store_true")
    args = ap.parse_args()
    if not Path(VINA).exists():
        raise SystemExit(f"vina not found at {VINA}")

    lig_dir = find_lig_dir()
    receptor = OUT / "receptors" / f"{args.alt}_receptor.pdbqt"
    box = OUT / "receptors" / f"{args.alt}_box.json"
    if not receptor.exists() or not box.exists():
        raise SystemExit(f"PASS receptor assets missing for {args.alt}")

    scores_path = OUT / "tables" / f"scores_vina_mode1_PAB_alt{args.alt}.csv"
    scores_path.parent.mkdir(parents=True, exist_ok=True)
    pose_root = OUT / "poses" / f"PAB_{args.alt}"
    pose_root.mkdir(parents=True, exist_ok=True)

    done = {}
    if scores_path.exists():
        prev = pd.read_csv(scores_path)
        for _, r in prev.iterrows():
            if r.get("status") in ("success", "exists") and pd.notna(r.get("vina_mode1")):
                done[r["ligand"]] = r.to_dict()
        print(f"resuming: {len(done)} already scored", flush=True)

    panel = pd.read_csv(PANEL_SCORES)
    if args.analyze_only:
        if not done:
            raise SystemExit(f"no scores in {scores_path}")
        write_metrics(args.alt, list(done.values()))
        return 0

    jobs = []
    for _, row in panel.iterrows():
        lig = row["ligand"]
        if lig in done:
            continue
        pdbqt = lig_dir / f"{lig}.pdbqt"
        if not pdbqt.exists():
            print(f"MISSING ligand {pdbqt}", flush=True)
            continue
        jobs.append((str(receptor), str(box), lig, str(pdbqt), str(pose_root / lig), args.timeout))

    print(f"{args.alt}: jobs={len(jobs)} workers={args.workers} vina={VINA}", flush=True)
    results = list(done.values())
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(dock_one, j): j for j in jobs}
        n = 0
        for fut in as_completed(futs):
            res = fut.result()
            results.append(res)
            n += 1
            print(
                f"[{n}/{len(jobs)}] {res['status']} {res['ligand']} "
                f"score={res['vina_mode1']} ({res['seconds']}s) elapsed={time.time()-t0:.0f}s",
                flush=True,
            )
            if n % 10 == 0:
                pd.DataFrame(results).to_csv(scores_path, index=False)

    pd.DataFrame(results).to_csv(scores_path, index=False)
    ok = sum(1 for r in results if r.get("status") in ("success", "exists"))
    print(f"DONE {args.alt}: {ok}/{len(results)} -> {scores_path}", flush=True)
    write_metrics(args.alt, results)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
