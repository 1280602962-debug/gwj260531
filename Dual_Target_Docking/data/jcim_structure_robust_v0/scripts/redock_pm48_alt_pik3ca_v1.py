#!/usr/bin/env python3
"""Redock frozen PM48 ligands to PASS alternative PM receptors.

- 4JPS / 5DXT: replace pocket A (PIK3CA); keep frozen 4JT6 scores for pocket B.
- 4JSX: replace pocket B (mTOR); keep frozen 4L23 scores for pocket A.

Then recomputes pocket-matched directional AUROC / summary_min with the same
bootstrap protocol as holdout/main analyses.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from vina import Vina

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/jcim_structure_robust_v0"
# Ligand PDBQTs live in the local results mirror (not always synced into the git pack).
LIG_DIR_CANDIDATES = [
    Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_rdkit_v0/ligands_pdbqt"),
    ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/ligands_pdbqt",
]
PANEL_SCORES = ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv"
SEED = 20260727
EXHAUST = 16
N_BOOT = 2000
BOOT_SEED = 20260729


def find_lig_dir() -> Path:
    for p in LIG_DIR_CANDIDATES:
        if p.is_dir() and any(p.glob("PM48_*.pdbqt")):
            return p
    raise FileNotFoundError("PM48 ligand PDBQTs not found in expected locations")


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


def dock_one(args):
    receptor_pdbqt, box_path, lig_id, lig_pdbqt, pose_dir, exhaust = args
    pose_dir = Path(pose_dir)
    pose_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = pose_dir / "out.pdbqt"
    box = json.loads(Path(box_path).read_text())
    t0 = time.time()
    try:
        v = Vina(sf_name="vina", cpu=1, seed=SEED, verbosity=0)
        v.set_receptor(str(receptor_pdbqt))
        v.set_ligand_from_file(str(lig_pdbqt))
        v.compute_vina_maps(
            center=[box["center_x"], box["center_y"], box["center_z"]],
            box_size=[box["size_x"], box["size_y"], box["size_z"]],
        )
        v.dock(exhaustiveness=exhaust, n_poses=9)
        v.write_poses(str(out_pdbqt), n_poses=9, overwrite=True, energy_range=3)
        mode1 = float(v.energies(n_poses=9)[0][0])
        return {"ligand": lig_id, "vina_mode1": mode1, "status": "success", "reason": "", "seconds": round(time.time() - t0, 1)}
    except Exception as exc:  # noqa: BLE001
        return {
            "ligand": lig_id,
            "vina_mode1": None,
            "status": "fail",
            "reason": f"{type(exc).__name__}: {exc}"[:300],
            "seconds": round(time.time() - t0, 1),
        }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alt", required=True, choices=["4JPS", "5DXT", "4JSX"])
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    # which pocket is replaced
    replace_pocket = "B" if args.alt == "4JSX" else "A"

    lig_dir = find_lig_dir()
    receptor = OUT / "receptors" / f"{args.alt}_receptor.pdbqt"
    box = OUT / "receptors" / f"{args.alt}_box.json"
    if not receptor.exists() or not box.exists():
        raise SystemExit(f"PASS receptor assets missing for {args.alt}")

    scores_path = OUT / "tables" / f"scores_vina_mode1_PM48_alt{args.alt}.csv"
    scores_path.parent.mkdir(parents=True, exist_ok=True)
    pose_root = OUT / "poses" / args.alt
    pose_root.mkdir(parents=True, exist_ok=True)

    done = {}
    if scores_path.exists():
        prev = pd.read_csv(scores_path)
        for _, r in prev.iterrows():
            if r.get("status") == "success":
                done[r["ligand"]] = r.to_dict()
        print(f"resuming: {len(done)} already scored", flush=True)
    else:
        prev = pd.DataFrame()

    panel = pd.read_csv(PANEL_SCORES)
    jobs = []
    for _, row in panel.iterrows():
        lig = row["ligand"]
        if lig in done:
            continue
        pdbqt = lig_dir / f"{lig}.pdbqt"
        if not pdbqt.exists():
            print(f"MISSING ligand {pdbqt}", flush=True)
            continue
        jobs.append((str(receptor), str(box), lig, str(pdbqt), str(pose_root / lig), EXHAUST))

    print(f"{args.alt}: jobs={len(jobs)} workers={args.workers}", flush=True)
    results = list(done.values())
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
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
    ok = sum(1 for r in results if r.get("status") == "success")
    print(f"DONE {args.alt}: {ok}/{len(results)} -> {scores_path}", flush=True)

    # Assemble pocket-matched metrics with one pocket replaced
    alt_map = {r["ligand"]: r for r in results if r.get("status") == "success" and r.get("vina_mode1") is not None}
    recs = []
    for _, row in panel.iterrows():
        lig = row["ligand"]
        if lig not in alt_map:
            continue
        a_fr = row.get("4L23_affinity")
        b_fr = row.get("4JT6_affinity")
        if pd.isna(a_fr) or pd.isna(b_fr):
            continue
        alt = float(alt_map[lig]["vina_mode1"])
        if replace_pocket == "A":
            a, b = alt, float(b_fr)
            note = "pocket A replaced by alt receptor; pocket B = frozen 4JT6 scores"
        else:
            a, b = float(a_fr), alt
            note = "pocket B replaced by alt receptor; pocket A = frozen 4L23 scores"
        recs.append(
            {
                "ligand": lig,
                "cls": row["class"],
                "vina_A": -a,
                "vina_B": -b,
                "alt_receptor": args.alt,
            }
        )
    da, db = directional(recs, "vina_B", "vina_A")
    sm = min(da, db)
    lo, hi = boot_min(recs, "vina_B", "vina_A")
    metric = {
        "alt_receptor": args.alt,
        "replaced_pocket": replace_pocket,
        "n": len(recs),
        "n_dual": sum(r["cls"] == "dual" for r in recs),
        "n_A_only": sum(r["cls"] == "A_only" for r in recs),
        "n_B_only": sum(r["cls"] == "B_only" for r in recs),
        "auroc_D_vs_A": round(da, 4),
        "auroc_D_vs_B": round(db, 4),
        "summary_min": round(sm, 4),
        "summary_min_ci_lo": round(lo, 4),
        "summary_min_ci_hi": round(hi, 4),
        "main_panel_summary_min": 0.6921,
        "note": note,
    }
    metric_path = OUT / "tables" / f"pocket_matched_PM48_alt{args.alt}_v1.csv"
    pd.DataFrame([metric]).to_csv(metric_path, index=False)
    print(metric, flush=True)
    print(f"wrote {metric_path}", flush=True)


if __name__ == "__main__":
    main()
