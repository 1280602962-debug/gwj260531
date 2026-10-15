#!/usr/bin/env python3
"""Run the stage3b GNINA queue, two jobs at a time."""
import json, os, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/home/gwj/work/stage3b-refdock/URAT1_NLRP3_DualTarget_AIDD_Project")
JOBS = json.loads((ROOT / "results/nlrp3_model_v2/stage3b_local/new32/job_list.json").read_text())
GNINA = "/mnt/d/CADD paper exercise/gnina/bin/gnina"
LIB = "/mnt/d/CADD paper exercise/gnina/conda_env/lib"
TIMEOUT_S = 45 * 60
MAX_JOBS = 2
LEDGER = ROOT / "results/nlrp3_model_v2/stage3b_local/new32/job_ledger.jsonl"

def already_done(job):
    out = Path(job["out"])
    status = out / "status.json"
    sdf = out / "poses_out.sdf"
    if not status.exists() or not sdf.exists() or sdf.stat().st_size == 0:
        return False
    try:
        rec = json.loads(status.read_text())
    except Exception:
        return False
    return rec.get("state") == "completed" and rec.get("n_poses", 0) >= 1 and rec.get("scores_readable") is True

def run_one(job):
    out = Path(job["out"])
    out.mkdir(parents=True, exist_ok=True)
    if already_done(job):
        rec = json.loads((out / "status.json").read_text())
        rec["skipped"] = True
        return rec
    cx, cy, cz = job["center"]
    sx, sy, sz = job["size"]
    cmd = [
        GNINA, "-r", job["receptor"], "-l", job["ligand_pdbqt"],
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", "32", "--num_modes", "9", "--cpu", "4",
        "--cnn_scoring", "rescore", "--seed", str(job["seed"]), "--no_gpu",
        "-o", str(out / "poses_out.sdf"), "--log", str(out / "gnina.log"),
    ]
    (out / "command.sh").write_text(
        "export LD_LIBRARY_PATH=" + LIB + "${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}\n"
        + " ".join(cmd) + "\n"
    )
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = LIB + ((":" + env["LD_LIBRARY_PATH"]) if env.get("LD_LIBRARY_PATH") else "")
    t0 = time.time()
    state = "failed"
    exit_code = None
    n_poses = 0
    scores_readable = False
    err = ""
    try:
        proc = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=TIMEOUT_S, env=env, check=False,
        )
        (out / "stdout.txt").write_bytes(proc.stdout or b"")
        exit_code = proc.returncode
    except subprocess.TimeoutExpired as exc:
        (out / "stdout.txt").write_bytes(exc.stdout or b"")
        state = "timeout"
        err = f"timeout_{TIMEOUT_S}s"
        exit_code = None
    elapsed = time.time() - t0
    sdf = out / "poses_out.sdf"
    if state != "timeout":
        if exit_code != 0:
            state = "failed"
            err = f"exit_{exit_code}"
        elif not sdf.exists() or sdf.stat().st_size == 0:
            state = "failed"
            err = "no_pose_file"
        else:
            # count SDF molecules and whether CNNscore is present; fewer than 9 is still completed
            text = sdf.read_text(errors="replace")
            n_poses = text.count("$$$$")
            scores_readable = "CNNscore" in text and n_poses >= 1
            if scores_readable:
                state = "completed"
            else:
                state = "failed"
                err = "pose_file_unreadable"
    rec = {
        "kind": job["kind"],
        "ligand": job["ligand"],
        "target": job["target"],
        "seed": job["seed"],
        "state": state,
        "n_poses": n_poses,
        "scores_readable": scores_readable,
        "exit_code": exit_code,
        "elapsed_s": round(elapsed, 1),
        "error": err,
        "out": str(out),
        "gnina": "v1.3.2 master:f23dd2b",
        "num_modes_requested": 9,
        "skipped": False,
    }
    (out / "status.json").write_text(json.dumps(rec, indent=2) + "\n")
    with LEDGER.open("a") as f:
        f.write(json.dumps(rec) + "\n")
    print(f"{state} {job['ligand']} {job['target']} seed{job['seed']} poses={n_poses} {elapsed:.0f}s {err}", flush=True)
    return rec

def main():
    pending = [j for j in JOBS if not already_done(j)]
    print(f"queue {len(JOBS)} pending {len(pending)} parallel {MAX_JOBS}", flush=True)
    counts = {}
    with ThreadPoolExecutor(max_workers=MAX_JOBS) as ex:
        futs = [ex.submit(run_one, j) for j in pending]
        for fut in as_completed(futs):
            rec = fut.result()
            counts[rec["state"]] = counts.get(rec["state"], 0) + 1
    print("DONE", json.dumps(counts), flush=True)

if __name__ == "__main__":
    main()
