#!/usr/bin/env python3
"""RTM best-of-K + Vina scores for a dual-target panel pack."""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_PY = RTM_ROOT / "example" / "rtmscore.py"
MODEL = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"
RTM_PYTHON = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")


def vina_from_pdbqt(path: Path):
    for line in path.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def auroc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    wins = 0.0
    for p in pos:
        wins += np.sum(p > neg) + 0.5 * np.sum(p == neg)
    return float(wins / (len(pos) * len(neg)))


def write_pose_sdfs(root: Path, targets: list[str], panel: pd.DataFrame):
    logs = root / "logs" / "rtmscore"
    logs.mkdir(parents=True, exist_ok=True)
    for target in targets:
        for _, r in panel.iterrows():
            lig = r["panel_id"]
            pose_dir = root / "poses" / target / lig
            modes = sorted(pose_dir.glob("mode_*.pdbqt"))
            if not modes:
                continue
            out_sdf = logs / f"{target}_{lig}_poses.sdf"
            if out_sdf.exists():
                continue
            # convert via obabel if available, else skip (RTM needs sdf/mol2)
            # Prefer writing from pdbqt with openbabel
            all_out = root / "logs" / "vina" / f"{target}_{lig}_out.pdbqt"
            src = all_out if all_out.exists() else modes[0]
            # concatenate modes into one multi-model pdbqt then convert
            if all_out.exists():
                subprocess.run(
                    ["obabel", str(all_out), "-O", str(out_sdf)],
                    capture_output=True,
                )
            else:
                # join modes
                tmp = logs / f"{target}_{lig}_joined.pdbqt"
                parts = []
                for i, m in enumerate(modes, 1):
                    body = [
                        ln
                        for ln in m.read_text().splitlines()
                        if not ln.startswith(("MODEL", "ENDMDL"))
                    ]
                    parts.append(f"MODEL {i}\n" + "\n".join(body) + "\nENDMDL")
                tmp.write_text("\n".join(parts) + "\n")
                subprocess.run(
                    ["obabel", str(tmp), "-O", str(out_sdf)], capture_output=True
                )


def run_rtm(root: Path, target: str, lig: str, receptor_pdb: Path):
    logs = root / "logs" / "rtmscore"
    sdf = logs / f"{target}_{lig}_poses.sdf"
    out_csv = logs / f"{target}_{lig}_rtm.csv"
    if out_csv.exists() and out_csv.stat().st_size > 0:
        return out_csv
    if not sdf.exists():
        return None
    # RTMScore example API: python rtmscore.py -p protein.pdb -l ligand.sdf -m model -o out
    cmd = [
        str(RTM_PYTHON),
        str(RTM_PY),
        "-p",
        str(receptor_pdb),
        "-l",
        str(sdf),
        "-m",
        str(MODEL),
        "-o",
        str(out_csv.with_suffix("")),
        "-gen_pocket",
        "-c",
        "10.0",
    ]
    proc = subprocess.run(cmd, cwd=str(RTM_ROOT), capture_output=True, text=True)
    (logs / f"{target}_{lig}_rtm.log").write_text(proc.stdout + "\n" + proc.stderr)
    # RTM may write .csv automatically
    cands = list(logs.glob(f"{target}_{lig}_rtm*.csv"))
    if not cands and out_csv.with_suffix(".csv").exists():
        return out_csv.with_suffix(".csv")
    return cands[0] if cands else (out_csv if out_csv.exists() else None)


def best_rtm(csv_path: Path):
    if csv_path is None or not csv_path.exists():
        return None
    try:
        df = pd.read_csv(csv_path)
    except Exception:
        return None
    # common column names
    for col in ("score", "rtmscore", "RTMScore", "pred"):
        if col in df.columns:
            return float(df[col].max())
    # fallback: last numeric column
    num = df.select_dtypes(include=[np.number])
    if num.shape[1] == 0:
        return None
    return float(num.iloc[:, -1].max())


def best_vina(root: Path, target: str, lig: str):
    pose_dir = root / "poses" / target / lig
    scores = []
    for m in pose_dir.glob("mode_*.pdbqt"):
        s = vina_from_pdbqt(m)
        if s is not None:
            scores.append(s)
    return (min(scores) if scores else None), len(scores)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--targets", nargs=2, required=True)
    ap.add_argument(
        "--receptor-pdb-map",
        nargs="+",
        help="TARGET=path/to/protein.pdb pairs",
        required=True,
    )
    args = ap.parse_args()
    root = Path(args.root)
    panel = pd.read_csv(args.panel)
    rec_map = {}
    for item in args.receptor_pdb_map:
        k, v = item.split("=", 1)
        rec_map[k] = Path(v)

    write_pose_sdfs(root, args.targets, panel)

    rows = []
    for _, r in panel.iterrows():
        lig = r["panel_id"]
        row = dict(r)
        for t in args.targets:
            v, n = best_vina(root, t, lig)
            row[f"vina_{t}"] = v
            row[f"n_modes_{t}"] = n
            rtm_csv = run_rtm(root, t, lig, rec_map[t])
            row[f"rtm_{t}"] = best_rtm(rtm_csv) if rtm_csv else None
        rows.append(row)

    out = root / "tables" / "scores_vina_rtm.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print("wrote", out)

    # directional on -vina and rtm (higher better for rtm; for vina use -score)
    df = pd.DataFrame(rows)
    a, b = args.targets
    summary = []
    for arm_name, col_a, col_b, flip in [
        ("vina", f"vina_{a}", f"vina_{b}", True),
        ("rtm", f"rtm_{a}", f"rtm_{b}", False),
    ]:
        sa = (-df[col_a].astype(float)) if flip else df[col_a].astype(float)
        sb = (-df[col_b].astype(float)) if flip else df[col_b].astype(float)
        tmp = df.copy()
        tmp["_A"], tmp["_B"] = sa, sb
        for label, arm in [("A", "_A"), ("B", "_B")]:
            d = tmp.loc[tmp["class"] == "dual", arm]
            ao = tmp.loc[tmp["class"] == "A_only", arm]
            bo = tmp.loc[tmp["class"] == "B_only", arm]
            summary.append(
                {
                    "channel": arm_name,
                    "end": label,
                    "auroc_D_vs_A": auroc(d, ao),
                    "auroc_D_vs_B": auroc(d, bo),
                }
            )
    sum_path = root / "tables" / "directional_summary.csv"
    pd.DataFrame(summary).to_csv(sum_path, index=False)
    print("wrote", sum_path)
    print(pd.DataFrame(summary).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
