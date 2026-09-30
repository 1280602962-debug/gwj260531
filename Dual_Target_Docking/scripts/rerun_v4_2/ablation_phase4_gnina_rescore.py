#!/usr/bin/env python3
"""Phase 4: GNINA score_only on unchanged Vina seed42 poses. No AUROC. No minimize."""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import (  # noqa: E402
    GNINA_BIN,
    GNINA_LIB,
    QA,
    REC,
    RESCORING_WORKERS,
    SEC,
)

WORK = SEC / "gnina_rescore"
MANIFEST = SEC / "rescoring_job_manifest.csv"
LEDGER = QA / "gnina_vina_rescore_ledger.jsonl"
LONG = SEC / "GNINA_VINA_POSE_RESCORE_POSE_LONG.csv"
MASTER = SEC / "GNINA_VINA_POSE_RESCORE_MASTER.csv"

CNN_RE = {
    "CNNscore": re.compile(r"CNNscore[:\s]+([0-9eE.+-]+)"),
    "CNNaffinity": re.compile(r"CNNaffinity[:\s]+([0-9eE.+-]+)"),
    "CNN_VS": re.compile(r"CNN_?VS[:\s]+([0-9eE.+-]+)"),
    "affinity": re.compile(r"(?:Affinity|VINA RESULT:|minimizedAffinity)[:\s]+([0-9eE.+-]+)"),
}


def env() -> dict:
    e = os.environ.copy()
    e["LD_LIBRARY_PATH"] = f"{GNINA_LIB}:{e.get('LD_LIBRARY_PATH', '')}"
    return e


def split_models(pdbqt_text: str) -> list[str]:
    chunks, cur = [], []
    for ln in pdbqt_text.splitlines(keepends=True):
        if ln.startswith("MODEL") and cur:
            chunks.append("".join(cur))
            cur = [ln]
        else:
            cur.append(ln)
    if cur:
        chunks.append("".join(cur))
    cleaned = []
    for chunk in chunks:
        keep = []
        for ln in chunk.splitlines(keepends=True):
            if ln.startswith(("MODEL", "ENDMDL")):
                continue
            keep.append(ln)
        text = "".join(keep)
        if "ATOM" in text or "HETATM" in text:
            cleaned.append(text)
    return cleaned or [pdbqt_text]


def parse_scores(text: str) -> dict:
    out = {k: "" for k in ("CNNscore", "CNNaffinity", "CNN_VS", "empirical_affinity")}
    m = CNN_RE["CNNscore"].search(text)
    if m:
        out["CNNscore"] = m.group(1)
    m = CNN_RE["CNNaffinity"].search(text)
    if m:
        out["CNNaffinity"] = m.group(1)
    m = CNN_RE["CNN_VS"].search(text)
    if m:
        out["CNN_VS"] = m.group(1)
    # last affinity-like number that is not CNN
    affs = CNN_RE["affinity"].findall(text)
    if affs:
        out["empirical_affinity"] = affs[-1]
    return out


def score_pose(receptor: Path, pose: Path, log: Path) -> dict:
    cmd = [
        str(GNINA_BIN), "--receptor", str(receptor), "--ligand", str(pose),
        "--score_only", "--no_gpu", "--scoring", "vina", "--cnn_scoring", "rescore",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env(), timeout=300)
    log.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
    scores = parse_scores((proc.stdout or "") + "\n" + (proc.stderr or ""))
    if pose.is_file():
        extra = parse_scores(pose.read_text(errors="replace"))
        for k, v in extra.items():
            if not scores.get(k) and v:
                scores[k] = v
    scores["rc"] = proc.returncode
    return scores


def run_one(row: dict) -> dict:
    pdb, eid = row["pdb_id"], row["global_ligand_entity_id"]
    jobdir = WORK / f"{pdb}__{eid}"
    jobdir.mkdir(parents=True, exist_ok=True)
    status_p = jobdir / "status.json"
    if status_p.is_file():
        prev = json.loads(status_p.read_text())
        if prev.get("job_status") in {"SUCCESS", "MISSING_NO_VINA_POSE"}:
            return prev
    rec = {
        "pdb_id": pdb,
        "global_ligand_entity_id": eid,
        "pairs": row["pairs"],
        "vina_status": row["vina_status"],
        "job_status": "",
        "n_poses": 0,
        "M1_CNNscore": "",
        "M1b_CNNscore": "",
        "CNNaffinity_mode1": "",
        "CNN_VS_mode1": "",
        "empirical_affinity_mode1": "",
        "reason": "",
        "n_retries": 0,
    }
    if row.get("rescoring_eligible") != "1" or not row.get("pose_path"):
        rec.update(job_status="MISSING_NO_VINA_POSE", reason="no_vina_pose")
        status_p.write_text(json.dumps(rec, indent=2) + "\n")
        rec["_long"] = []
        return rec
    poses = split_models(Path(row["pose_path"]).read_text(errors="replace"))
    rec["n_poses"] = len(poses)
    receptor = REC / f"{pdb}_receptor.pdbqt"
    long_rows = []

    def once() -> list[dict]:
        out = []
        for i, block in enumerate(poses, 1):
            pose_p = jobdir / f"vina_mode{i}.pdbqt"
            safe = "".join(ln for ln in block.splitlines(keepends=True) if not ln.startswith(("MODEL", "ENDMDL")))
            pose_p.write_text(safe, encoding="utf-8")
            sc = score_pose(receptor, pose_p, jobdir / f"vina_mode{i}.gnina.log")
            out.append({"mode": i, **sc})
        return out

    scored = once()
    if any(s.get("rc") != 0 or not s.get("CNNscore") for s in scored):
        rec["n_retries"] = 1
        scored = once()
    ok = all(s.get("CNNscore") not in ("", None) for s in scored)
    cnn = []
    for s in scored:
        try:
            cnn.append(float(s["CNNscore"]))
        except (TypeError, ValueError):
            cnn.append(float("nan"))
        long_rows.append({
            "pdb_id": pdb,
            "global_ligand_entity_id": eid,
            "pairs": row["pairs"],
            "vina_seed": 42,
            "vina_mode": s["mode"],
            "is_vina_mode1": int(s["mode"] == 1),
            "CNNscore": s.get("CNNscore", ""),
            "CNNaffinity": s.get("CNNaffinity", ""),
            "CNN_VS": s.get("CNN_VS", ""),
            "empirical_affinity": s.get("empirical_affinity", ""),
            "gnina_rc": s.get("rc", ""),
        })
    if not ok or not scored:
        rec.update(job_status="TECHNICAL_FAIL", reason="parse_or_rc", _long=long_rows)
    else:
        rec.update(
            job_status="SUCCESS",
            M1_CNNscore=scored[0].get("CNNscore", ""),
            CNNaffinity_mode1=scored[0].get("CNNaffinity", ""),
            CNN_VS_mode1=scored[0].get("CNN_VS", ""),
            empirical_affinity_mode1=scored[0].get("empirical_affinity", ""),
            M1b_CNNscore=f"{max(cnn):.6f}" if cnn else "",
            _long=long_rows,
        )
    slim = {k: v for k, v in rec.items() if k != "_long"}
    status_p.write_text(json.dumps(slim, indent=2) + "\n")
    rec["_long"] = long_rows
    return rec


def main() -> int:
    WORK.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8-sig", newline="")))
    done = []
    pending = []
    for r in rows:
        st = WORK / f"{r['pdb_id']}__{r['global_ligand_entity_id']}" / "status.json"
        if st.is_file():
            prev = json.loads(st.read_text())
            if prev.get("job_status") in {"SUCCESS", "MISSING_NO_VINA_POSE"}:
                # reload long from pose logs later if needed
                done.append(prev)
                continue
        pending.append(r)
    print(f"gnina_rescore total={len(rows)} cached={len(done)} pending={len(pending)} workers={RESCORING_WORKERS}", flush=True)
    long_all = []
    masters = { (d["pdb_id"], d["global_ligand_entity_id"]): d for d in done if "pdb_id" in d }

    # recover long rows from existing success dirs
    for d in done:
        jobdir = WORK / f"{d['pdb_id']}__{d['global_ligand_entity_id']}"
        lp = jobdir / "long.csv"
        if lp.is_file():
            long_all.extend(csv.DictReader(lp.open(encoding="utf-8-sig")))

    with ThreadPoolExecutor(max_workers=RESCORING_WORKERS) as pool:
        futs = {pool.submit(run_one, r): r for r in pending}
        for i, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            long_rows = rec.pop("_long", [])
            jobdir = WORK / f"{rec['pdb_id']}__{rec['global_ligand_entity_id']}"
            if long_rows:
                with (jobdir / "long.csv").open("w", encoding="utf-8", newline="") as h:
                    w = csv.DictWriter(h, fieldnames=list(long_rows[0].keys()), lineterminator="\n")
                    w.writeheader()
                    w.writerows(long_rows)
            long_all.extend(long_rows)
            masters[(rec["pdb_id"], rec["global_ligand_entity_id"])] = rec
            with LEDGER.open("a", encoding="utf-8") as h:
                h.write(json.dumps({k: rec[k] for k in rec}, sort_keys=True) + "\n")
            if i % 20 == 0 or rec.get("job_status") != "SUCCESS":
                print(f"{len(masters)}/{len(rows)} {rec['job_status']} {rec['pdb_id']} {rec['global_ligand_entity_id']}", flush=True)

    master_rows = [masters[k] for k in sorted(masters)]
    fields = [
        "pdb_id", "global_ligand_entity_id", "pairs", "vina_status", "job_status",
        "n_poses", "M1_CNNscore", "M1b_CNNscore", "CNNaffinity_mode1",
        "CNN_VS_mode1", "empirical_affinity_mode1", "reason", "n_retries",
    ]
    with MASTER.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(master_rows)
    if long_all:
        with LONG.open("w", encoding="utf-8", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(long_all[0].keys()), extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            w.writerows(long_all)
    counts = {}
    for r in master_rows:
        counts[r.get("job_status", "?")] = counts.get(r.get("job_status", "?"), 0) + 1
    print(json.dumps({"n": len(master_rows), **counts, "AUROC": "NOT_COMPUTED"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
