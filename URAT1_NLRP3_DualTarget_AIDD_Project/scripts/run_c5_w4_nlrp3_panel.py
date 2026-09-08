#!/usr/bin/env python3
"""C5 W4: NLRP3 panel rebuild at 7ALV (Job B).

Authorized after W1 benzbromarone@9DKA gate fail with fork search_ok_selection_fail
(see docs/C5_RANKING_AND_NEXT_DOCKS.md / config/campaign_c5.yaml). W4 has no
scientific dependency on the URAT1 gate.

Requires decoy CSV already locked: data/campaigns/c5/02_nlrp3_panel/w4_decoys_locked.csv

Jobs:
  - 8 non-crystal positives × seeds 42/43/44 = 24
  - 40 locked decoys × 3 seeds = 120
  - REP_07837 @ 7ALV seeds 43/44 = 2
  Total new gnina = 146
Reuse: NP3-146 × 3; background clinical acids except REP_07837 seeds 43/44

Prep lock (2026-09-05): positives AND decoys both go through
scripts/prepare_ligands_c1.py (Dimorphite-DL pH 7.4 → Meeko, embed 0xC0FFEE).
Do NOT reuse the old panel Minimal-Meeko-only PDBQTs (no Dimorphite).

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
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

SEEDS = [42, 43, 44]
TIMEOUT_SEC = 7200
DROP_POSITIVE_IDS = {"CHEMBL3183703", "NP3-146"}  # NP3-146 reused from selfdock
GATE = PROJECT_ROOT / "data/campaigns/c5/01_crossdock/gate_benzbromarone_9dka.json"
CAMPAIGN = PROJECT_ROOT / "config/campaign_c5.yaml"
PANEL_CSV = (
    PROJECT_ROOT
    / "data/campaigns/c1/05_metrics/nlrp3_structural_panel/panel_ligands.csv"
)
DECOY_CSV = PROJECT_ROOT / "data/campaigns/c5/02_nlrp3_panel/w4_decoys_locked.csv"
OUT = PROJECT_ROOT / "data/campaigns/c5/02_nlrp3_panel"
POS_PREP_DIR = OUT / "positive_ligands"
DECOY_PREP_DIR = OUT / "decoy_ligands"
REP07837 = (
    PROJECT_ROOT
    / "data/campaigns/c1/01_ligand_prep/acid_clinical_chemistry_pass/pdbqt/REP_07837.pdbqt"
)


def authorize_w4(strict_gate_pass: bool) -> dict:
    """W4 may run when campaign authorizes it after selection-fail audit."""
    meta = {"gate_path": str(GATE), "authorized": False, "reason": ""}
    if not GATE.exists():
        raise SystemExit(f"missing gate file {GATE}")
    gate = json.loads(GATE.read_text())
    meta["gate_pass"] = bool(gate.get("pass"))
    if gate.get("pass"):
        meta["authorized"] = True
        meta["reason"] = "w1_gate_pass"
        return meta
    if strict_gate_pass:
        raise SystemExit("W1 gate failed — W4 forbidden under --strict-gate-pass")
    camp = yaml.safe_load(CAMPAIGN.read_text()) if CAMPAIGN.exists() else {}
    next_docks = (camp.get("ranking_policy") or {}).get("next_docks") or {}
    w4 = next_docks.get("w4_nlrp3_panel") or {}
    fork = ((camp.get("w1_urat1_crossdock") or {}).get("audit_conclusion_2026-09-05") or {}).get(
        "fork"
    )
    if w4.get("status") == "authorized_must_run" and fork == "search_ok_selection_fail":
        meta["authorized"] = True
        meta["reason"] = "campaign_c5_authorized_after_search_ok_selection_fail"
        meta["fork"] = fork
        meta["n_jobs_new"] = w4.get("n_jobs_new")
        return meta
    raise SystemExit(
        "W1 gate failed and campaign does not authorize W4 "
        f"(fork={fork!r}, w4_status={w4.get('status')!r})"
    )


def prepare_via_c1(rows: list[dict], out_dir: Path, label: str) -> Path:
    """Dimorphite-DL pH 7.4 → Meeko into out_dir/pdbqt. Returns pdbqt dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    prep_csv = out_dir / f"{label}_for_prep.csv"
    pd.DataFrame(rows).to_csv(prep_csv, index=False)
    pdbqt_dir = out_dir / "pdbqt"
    need = len(rows)
    have = len(list(pdbqt_dir.glob("*.pdbqt"))) if pdbqt_dir.exists() else 0
    if have >= need:
        print(f"Reuse {label} pdbqt ({have}/{need}): {pdbqt_dir}", flush=True)
        return pdbqt_dir
    print(f"Preparing {label} (Dimorphite-DL → Meeko)...", flush=True)
    subprocess.check_call(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts/prepare_ligands_c1.py"),
            "--input-csv",
            str(prep_csv),
            "--output-dir",
            str(out_dir),
        ]
    )
    return pdbqt_dir


def _unlink_quiet(path: Path) -> None:
    try:
        if path.exists():
            path.unlink()
    except OSError:
        pass


def run_gnina(gnina, receptor, ligand, center, size, out_sdf, seed, cpu, timeout) -> str:
    # Only treat as done if SDF looks complete (has a mode /$$$$). Partial files from
    # prior timeouts must not block a retry and must not count as success.
    if out_sdf.exists() and out_sdf.stat().st_size > 0:
        txt = out_sdf.read_text(errors="ignore")
        if "$$$$" in txt and ("CNNscore" in txt or "minimizedAffinity" in txt or "mode" in txt.lower()):
            return "exists"
        _unlink_quiet(out_sdf)
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
    # start_new_session so timeout can kill the whole gnina process group
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            start_new_session=True,
        )
    except subprocess.TimeoutExpired as exc:
        # Drop partial SDF so we skip forward and can retry later if wanted
        _unlink_quiet(out_sdf)
        marker = out_sdf.parent / f"{out_sdf.stem}_TIMEOUT.txt"
        marker.write_text(
            f"timeout_sec={timeout}\nelapsed_approx_sec={time.time()-t0:.0f}\n"
            f"policy=skip_and_continue\n"
        )
        # Best-effort: kill orphaned process group if still around
        if exc.pid:
            try:
                import os
                import signal

                os.killpg(exc.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError, OSError):
                pass
        print(f"TIMEOUT skip {out_sdf} (>{timeout}s) — continue next job", flush=True)
        return "timeout"
    (out_sdf.parent / f"{out_sdf.stem}_stdout.txt").write_text(
        (proc.stdout or "") + "\n" + (proc.stderr or "")
    )
    if proc.returncode != 0 or not out_sdf.exists() or out_sdf.stat().st_size == 0:
        _unlink_quiet(out_sdf)
        print(f"FAIL skip rc={proc.returncode} {out_sdf} — continue next job", flush=True)
        return "fail"
    print(f"OK {time.time()-t0:.0f}s {out_sdf}", flush=True)
    return "ok"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", type=int, default=4, help="cpus per gnina process")
    ap.add_argument("--workers", type=int, default=1, help="parallel gnina jobs")
    ap.add_argument("--timeout", type=int, default=TIMEOUT_SEC)
    ap.add_argument(
        "--strict-gate-pass",
        action="store_true",
        help="require gate_benzbromarone_9dka.json pass=true (legacy)",
    )
    ap.add_argument(
        "--prep-only",
        action="store_true",
        help="prepare ligands and write job list; do not dock",
    )
    args = ap.parse_args()

    auth = authorize_w4(strict_gate_pass=args.strict_gate_pass)
    print("W4 authorization:", json.dumps(auth), flush=True)

    if not DECOY_CSV.exists():
        raise SystemExit(f"lock decoys first: {DECOY_CSV}")
    if not REP07837.exists():
        raise SystemExit(f"missing {REP07837}")

    cfg = yaml.safe_load((PROJECT_ROOT / "config/docking_c5_w1.yaml").read_text())
    tcfg = cfg["targets"]["nlrp3_7alv"]
    receptor = PROJECT_ROOT / tcfg["prepared_receptor"]
    if not receptor.exists():
        raise SystemExit(f"missing receptor {receptor}")
    center, size = tcfg["center"], tcfg["size"]
    gnina = PROJECT_ROOT / "tools/gnina"
    if not gnina.exists():
        raise SystemExit(f"missing gnina {gnina}")

    panel = pd.read_csv(PANEL_CSV)
    smile_col = "smiles" if "smiles" in panel.columns else "canonical_smiles"
    # Dock 8 non-crystal positives (NP3-146 reused; CHEMBL3183703 == MCC950 dropped)
    positives = panel[
        panel["role"].isin(["tool_positive", "chembl_sulfonylurea_active"])
    ].copy()
    positives = positives[~positives["ligand_id"].isin(DROP_POSITIVE_IDS)]
    if len(positives) != 8:
        raise SystemExit(f"expected 8 non-crystal positives, got {len(positives)}")

    pos_rows = [
        {"repurposing_id": r["ligand_id"], "canonical_smiles": r[smile_col]}
        for _, r in positives.iterrows()
    ]
    pos_pdbqt_dir = prepare_via_c1(pos_rows, POS_PREP_DIR, "positives")

    decoys = pd.read_csv(DECOY_CSV)
    decoy_rows = [
        {"repurposing_id": r["decoy_id"], "canonical_smiles": r["canonical_smiles"]}
        for _, r in decoys.iterrows()
    ]
    decoy_pdbqt_dir = prepare_via_c1(decoy_rows, DECOY_PREP_DIR, "decoys")

    # verify all pdbqts
    missing = []
    for r in pos_rows:
        p = pos_pdbqt_dir / f"{r['repurposing_id']}.pdbqt"
        if not p.exists() or p.stat().st_size == 0:
            missing.append(str(p))
    for r in decoy_rows:
        p = decoy_pdbqt_dir / f"{r['repurposing_id']}.pdbqt"
        if not p.exists() or p.stat().st_size == 0:
            missing.append(str(p))
    if missing:
        raise SystemExit(f"prep incomplete, missing {len(missing)} pdbqt e.g. {missing[:5]}")

    prep_meta = {
        "protocol": "dimorphite_dl_ph7.4_then_meeko_embed_0xC0FFEE",
        "script": "scripts/prepare_ligands_c1.py",
        "do_not_reuse": "data/campaigns/c1/05_metrics/nlrp3_structural_panel/pdbqt (Meeko-only, no Dimorphite)",
        "positives_n": len(pos_rows),
        "decoys_n": len(decoy_rows),
        "positives_dir": str(pos_pdbqt_dir),
        "decoys_dir": str(decoy_pdbqt_dir),
        "rep_07837_pdbqt": str(REP07837),
        "receptor": str(receptor),
        "box": {"center": center, "size": size},
        "authorization": auth,
        "settings": {
            "exhaustiveness": 32,
            "num_modes": 9,
            "cnn_scoring": "rescore",
            "seeds": SEEDS,
            "no_gpu": True,
            "cpu_per_job": args.cpu,
            "workers": args.workers,
        },
        "started_unix": time.time(),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "w4_prep_meta.json").write_text(json.dumps(prep_meta, indent=2) + "\n")

    jobs = []
    for _, row in positives.iterrows():
        lid = row["ligand_id"]
        pdbqt = pos_pdbqt_dir / f"{lid}.pdbqt"
        for seed in SEEDS:
            jobs.append(("positives", lid, pdbqt, seed))
    for seed in (43, 44):
        jobs.append(("background", "REP_07837", REP07837, seed))
    for _, row in decoys.iterrows():
        did = row["decoy_id"]
        pdbqt = decoy_pdbqt_dir / f"{did}.pdbqt"
        for seed in SEEDS:
            jobs.append(("decoys", did, pdbqt, seed))

    print(f"W4 new jobs: {len(jobs)} (expect 146)", flush=True)
    if len(jobs) != 146:
        print(f"WARNING: expected 146 jobs, got {len(jobs)}", flush=True)

    job_list = [
        {
            "subset": s,
            "ligand_id": lid,
            "seed": seed,
            "pdbqt": str(pdbqt),
            "out_sdf": str(OUT / s / f"seed{seed}" / f"{lid}_out.sdf"),
        }
        for s, lid, pdbqt, seed in jobs
    ]
    pd.DataFrame(job_list).to_csv(OUT / "w4_job_list.csv", index=False)

    if args.prep_only:
        print("prep-only done →", OUT / "w4_prep_meta.json", flush=True)
        return

    status_path = OUT / "w4_job_status.csv"
    status_rows: list[dict] = []
    if status_path.exists():
        status_rows = pd.read_csv(status_path).to_dict(orient="records")

    def _one(item):
        subset, lid, pdbqt, seed = item
        out_sdf = OUT / subset / f"seed{seed}" / f"{lid}_out.sdf"
        st = run_gnina(
            gnina, receptor, pdbqt, center, size, out_sdf, seed, args.cpu, args.timeout
        )
        return {
            "subset": subset,
            "ligand_id": lid,
            "seed": seed,
            "status": st,
            "sdf": str(out_sdf),
        }

    workers = max(1, int(args.workers))
    if workers == 1:
        for item in jobs:
            row = _one(item)
            status_rows.append(row)
            pd.DataFrame(status_rows).to_csv(status_path, index=False)
    else:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_one, item): item for item in jobs}
            for fut in as_completed(futs):
                row = fut.result()
                status_rows.append(row)
                pd.DataFrame(status_rows).to_csv(status_path, index=False)

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

    summary = {
        "n_jobs": len(jobs),
        "status_counts": pd.DataFrame(status_rows)["status"].value_counts().to_dict()
        if status_rows
        else {},
        "status_csv": str(status_path),
        "finished_unix": time.time(),
    }
    (OUT / "w4_run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("W4 docking finished →", status_path, flush=True)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
