#!/usr/bin/env python3
"""C5 W4: NLRP3 panel rebuild at 7ALV (Job B).

Requires gate file: data/campaigns/c5/01_crossdock/gate_benzbromarone_9dka.json with pass=true
Requires decoy CSV already locked: data/campaigns/c5/02_nlrp3_panel/w4_decoys_locked.csv

Jobs:
  - 8 non-crystal positives × seeds 42/43/44 = 24
  - ≥40 decoys × 3 seeds
  - REP_07837 @ 7ALV seeds 43/44 = 2
Reuse: NP3-146 × 3; background 20×seed42 + 19×seed43/44 (not re-docked here)

Settings: exh=32, modes=9, rescore, seeds 42/43/44, --no_gpu.
Do NOT call run_nlrp3_structural_panel.py (exh=8).
Timeouts: skip and continue.
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

SEEDS = [42, 43, 44]
TIMEOUT_SEC = 7200
DROP_POSITIVE_IDS = {"CHEMBL3183703"}
GATE = PROJECT_ROOT / "data/campaigns/c5/01_crossdock/gate_benzbromarone_9dka.json"
PANEL_CSV = (
    PROJECT_ROOT
    / "data/campaigns/c1/05_metrics/nlrp3_structural_panel/panel_ligands.csv"
)
PANEL_PDBQT = (
    PROJECT_ROOT
    / "data/campaigns/c1/05_metrics/nlrp3_structural_panel/pdbqt"
)
DECOY_CSV = PROJECT_ROOT / "data/campaigns/c5/02_nlrp3_panel/w4_decoys_locked.csv"
OUT = PROJECT_ROOT / "data/campaigns/c5/02_nlrp3_panel"
REP07837 = (
    PROJECT_ROOT
    / "data/campaigns/c1/01_ligand_prep/acid_clinical_chemistry_pass/pdbqt/REP_07837.pdbqt"
)


def run_gnina(gnina, receptor, ligand, center, size, out_sdf, seed, cpu, timeout) -> str:
    if out_sdf.exists() and out_sdf.stat().st_size > 0:
        return "exists"
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    log = out_sdf.with_suffix(".log")
    cmd = [
        str(gnina), "-r", str(receptor), "-l", str(ligand),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", "32", "--num_modes", "9", "--cpu", str(cpu),
        "--cnn_scoring", "rescore", "--seed", str(seed),
        "-o", str(out_sdf), "--log", str(log), "--no_gpu",
    ]
    print("RUN:", " ".join(cmd), flush=True)
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        (out_sdf.parent / f"{out_sdf.stem}_TIMEOUT.txt").write_text(f"timeout={timeout}\n")
        print(f"TIMEOUT {out_sdf}", flush=True)
        return "timeout"
    (out_sdf.parent / f"{out_sdf.stem}_stdout.txt").write_text(
        (proc.stdout or "") + "\n" + (proc.stderr or "")
    )
    if proc.returncode != 0 or not out_sdf.exists() or out_sdf.stat().st_size == 0:
        print(f"FAIL rc={proc.returncode} {out_sdf}", flush=True)
        return "fail"
    print(f"OK {time.time()-t0:.0f}s {out_sdf}", flush=True)
    return "ok"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=TIMEOUT_SEC)
    ap.add_argument("--skip-gate-check", action="store_true")
    args = ap.parse_args()

    if not args.skip_gate_check:
        if not GATE.exists():
            raise SystemExit(f"missing gate file {GATE}")
        gate = json.loads(GATE.read_text())
        if not gate.get("pass"):
            raise SystemExit("W1 gate failed — W4 forbidden")

    if not DECOY_CSV.exists():
        raise SystemExit(f"lock decoys first: {DECOY_CSV}")

    cfg = yaml.safe_load((PROJECT_ROOT / "config/docking_c5_w1.yaml").read_text())
    tcfg = cfg["targets"]["nlrp3_7alv"]
    receptor = PROJECT_ROOT / tcfg["prepared_receptor"]
    center, size = tcfg["center"], tcfg["size"]
    gnina = PROJECT_ROOT / "tools/gnina"

    panel = pd.read_csv(PANEL_CSV)
    positives = panel[
        panel["role"].isin(
            ["crystal_positive", "tool_positive", "chembl_sulfonylurea_active"]
        )
    ].copy()
    positives = positives[~positives["ligand_id"].isin(DROP_POSITIVE_IDS)]

    # Prepare decoy pdbqts if missing
    decoys = pd.read_csv(DECOY_CSV)
    decoy_prep_csv = OUT / "w4_decoys_for_prep.csv"
    decoys_renamed = decoys.rename(
        columns={"decoy_id": "repurposing_id", "canonical_smiles": "canonical_smiles"}
    )[["repurposing_id", "canonical_smiles"]].copy()
    decoys_renamed.to_csv(decoy_prep_csv, index=False)
    decoy_prep_dir = OUT / "decoy_ligands"
    decoy_pdbqt_dir = decoy_prep_dir / "pdbqt"
    if not decoy_pdbqt_dir.exists() or len(list(decoy_pdbqt_dir.glob("*.pdbqt"))) < len(decoys):
        print("Preparing decoy ligands (Dimorphite-DL → Meeko)...", flush=True)
        subprocess.check_call(
            [
                sys.executable,
                str(PROJECT_ROOT / "scripts/prepare_ligands_c1.py"),
                "--input-csv",
                str(decoy_prep_csv),
                "--output-dir",
                str(decoy_prep_dir),
            ]
        )

    jobs = []
    # positives except NP3-146 (reuse selfdock)
    for _, row in positives.iterrows():
        lid = row["ligand_id"]
        if lid == "NP3-146":
            continue
        pdbqt = PANEL_PDBQT / f"{lid}.pdbqt"
        if not pdbqt.exists():
            raise SystemExit(f"missing pdbqt {pdbqt}")
        for seed in SEEDS:
            jobs.append(("positives", lid, pdbqt, seed))

    # REP_07837 seeds 43/44 only
    for seed in (43, 44):
        jobs.append(("background", "REP_07837", REP07837, seed))

    # decoys
    for _, row in decoys.iterrows():
        did = row["decoy_id"]
        pdbqt = decoy_pdbqt_dir / f"{did}.pdbqt"
        if not pdbqt.exists():
            raise SystemExit(f"missing decoy pdbqt {pdbqt}")
        for seed in SEEDS:
            jobs.append(("decoys", did, pdbqt, seed))

    print(f"W4 new jobs: {len(jobs)}", flush=True)
    status_rows = []
    for subset, lid, pdbqt, seed in jobs:
        out_sdf = OUT / subset / f"seed{seed}" / f"{lid}_out.sdf"
        st = run_gnina(
            gnina, receptor, pdbqt, center, size, out_sdf, seed, args.cpu, args.timeout
        )
        status_rows.append(
            {"subset": subset, "ligand_id": lid, "seed": seed, "status": st, "sdf": str(out_sdf)}
        )
        pd.DataFrame(status_rows).to_csv(OUT / "w4_job_status.csv", index=False)

    # symlink/copy NP3-146 reuse
    for seed in SEEDS:
        src = (
            PROJECT_ROOT
            / f"data/campaigns/c1/02_selfdock/nlrp3_7alv/seed{seed}/NP3-146_out.sdf"
        )
        dst = OUT / "positives" / f"seed{seed}" / "NP3-146_out.sdf"
        dst.parent.mkdir(parents=True, exist_ok=True)
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)
            (dst.parent / "NP3-146_REUSED_FROM.txt").write_text(str(src) + "\n")

    print("W4 docking finished →", OUT / "w4_job_status.csv", flush=True)


if __name__ == "__main__":
    main()
