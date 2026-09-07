#!/usr/bin/env python3
"""Track B production Vina for five pairs (~1100 jobs).

Uses Layer-3-passed receptors/boxes in local_track_b_v0/.
Timeout / TORSDOF>=25 → skip (recorded in job_status). Does not rewrite Table 2.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SEED = 20260727
N_MODES = 9
ENERGY_RANGE = 3
EXHAUSTIVENESS = 8
SKIP_TORSDOF_GE = 25
VINA = "/home/gwj/miniconda3/bin/vina"

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_track_b_v0"
TABLES = ROOT / "tables"
LIG_DIR = ROOT / "cache" / "track_b_ligands" / "pdbqt"

PAIRS = [
    {
        "pair": "F2/F10",
        "panel": TABLES / "track_b_panels" / "panel_F2_F10_v1.csv",
        "targets": ["4UDW", "2JKH"],
    },
    {
        "pair": "JAK1/TYK2",
        "panel": TABLES / "track_b_panels" / "panel_JAK1_TYK2_v1.csv",
        "targets": ["6N7A", "3LXP"],
    },
    {
        "pair": "JAK1/JAK2",
        "panel": TABLES / "track_b_panels" / "panel_JAK1_JAK2_v1.csv",
        "targets": ["6N7A", "8BXH"],
    },
    {
        "pair": "PPARG/PPARA",
        "panel": TABLES / "track_b_panels" / "panel_PPARG_PPARA_v1.csv",
        "targets": ["9V8H", "6LXA"],
    },
    {
        "pair": "PPARA/PPARD",
        "panel": TABLES / "track_b_panels" / "panel_PPARA_PPARD_v1.csv",
        "targets": ["6LXA", "5U3Q"],
    },
]


def torsdof(ligand_pdbqt: Path) -> int:
    for line in ligand_pdbqt.read_text().splitlines():
        if line.startswith("TORSDOF"):
            return int(line.split()[1])
    return 0


def parse_mode1_energy(out_pdbqt: Path) -> float | None:
    """Vina mode-1 affinity from REMARK VINA RESULT on first MODEL."""
    text = out_pdbqt.read_text().splitlines()
    in_model = False
    for line in text:
        if line.startswith("MODEL"):
            in_model = True
            continue
        if in_model and "VINA RESULT" in line:
            # REMARK VINA RESULT:    -8.2      0.000      0.000
            m = re.search(r"VINA RESULT:\s+(-?\d+\.?\d*)", line)
            if m:
                return float(m.group(1))
        if line.startswith("ENDMDL") and in_model:
            break
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


def write_conf(target: str, lig: str, ligand_pdbqt: Path, box: dict, exhaust: int) -> tuple[Path, Path]:
    conf_dir = LOCAL / "logs" / "vina_confs"
    conf_dir.mkdir(parents=True, exist_ok=True)
    (LOCAL / "logs" / "vina").mkdir(parents=True, exist_ok=True)
    conf = conf_dir / f"{target}_{lig}.txt"
    rec = LOCAL / "receptors" / f"{target}_receptor.pdbqt"
    out = LOCAL / "logs" / "vina" / f"{target}_{lig}_out.pdbqt"
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
                f"num_modes = {N_MODES}",
                f"energy_range = {ENERGY_RANGE}",
                "cpu = 1",
                f"seed = {SEED}",
                f"out = {out}",
            ]
        )
        + "\n"
    )
    return conf, out


def run_one(
    pair: str,
    target: str,
    lig: str,
    ligand_pdbqt: Path,
    box: dict,
    exhaust: int,
    timeout_s: int,
) -> dict:
    pose_dir = LOCAL / "poses" / target / lig
    mode1 = pose_dir / "mode_01.pdbqt"
    if mode1.exists() and mode1.stat().st_size > 0:
        e = parse_mode1_energy(mode1)
        return {
            "pair": pair,
            "target": target,
            "ligand": lig,
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
            "status": "fail",
            "n_modes": 0,
            "mode1_energy": None,
            "score_S": None,
            "reason": "missing_ligand_pdbqt",
            "seconds": 0.0,
        }

    td = torsdof(ligand_pdbqt)
    if td >= SKIP_TORSDOF_GE:
        return {
            "pair": pair,
            "target": target,
            "ligand": lig,
            "status": "skip",
            "n_modes": 0,
            "mode1_energy": None,
            "score_S": None,
            "reason": f"skip_torsdof={td}_ge_{SKIP_TORSDOF_GE}",
            "seconds": 0.0,
        }

    conf, out = write_conf(target, lig, ligand_pdbqt, box, exhaust)
    log = LOCAL / "logs" / "vina" / f"{target}_{lig}.log"
    t0 = time.time()
    try:
        proc = subprocess.run(
            [VINA, "--config", str(conf)],
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
        log.write_text(proc.stdout + "\n" + proc.stderr)
        rc = proc.returncode
        err = (proc.stderr or proc.stdout or "")[-300:]
    except subprocess.TimeoutExpired as e:
        elapsed = time.time() - t0
        log.write_text(
            f"TIMEOUT after {timeout_s}s torsdof={td} E={exhaust}\n"
            f"{(e.stdout or b'').decode(errors='ignore') if isinstance(e.stdout, (bytes, bytearray)) else (e.stdout or '')}"
        )
        # Remove partial output so resume does not treat as success.
        if out.exists():
            out.unlink()
        return {
            "pair": pair,
            "target": target,
            "ligand": lig,
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
            "status": "fail",
            "n_modes": 0,
            "mode1_energy": None,
            "score_S": None,
            "reason": err or f"vina_rc={rc}",
            "seconds": round(elapsed, 1),
        }

    n = split_modes(out, pose_dir)
    e = parse_mode1_energy(out)
    return {
        "pair": pair,
        "target": target,
        "ligand": lig,
        "status": "success",
        "n_modes": n,
        "mode1_energy": e,
        "score_S": None if e is None else -e,
        "reason": "",
        "seconds": round(elapsed, 1),
    }


def load_jobs(pair_filter: set[str] | None) -> list[tuple]:
    boxes = {
        p.name.replace("_box.json", ""): json.loads(p.read_text())
        for p in (LOCAL / "boxes").glob("*_box.json")
    }
    jobs = []
    for spec in PAIRS:
        if pair_filter and spec["pair"] not in pair_filter:
            continue
        for t in spec["targets"]:
            if t not in boxes:
                raise SystemExit(f"missing box for {t}")
            rec = LOCAL / "receptors" / f"{t}_receptor.pdbqt"
            if not rec.exists():
                raise SystemExit(f"missing receptor {rec}")
        panel = list(csv.DictReader(spec["panel"].open()))
        for r in panel:
            lig = r["panel_id"]
            pdbqt = LIG_DIR / f"{lig}.pdbqt"
            for t in spec["targets"]:
                jobs.append((spec["pair"], t, lig, pdbqt, boxes[t]))
    # Easy ligands first so flexible ones cannot monopolize the pool.
    jobs.sort(key=lambda j: torsdof(j[3]) if j[3].exists() else 99)
    return jobs


def write_tables(results: list[dict]) -> None:
    tab = LOCAL / "tables"
    tab.mkdir(parents=True, exist_ok=True)
    status_path = tab / "job_status.csv"
    fields = [
        "pair",
        "target",
        "ligand",
        "status",
        "n_modes",
        "mode1_energy",
        "score_S",
        "reason",
        "seconds",
    ]
    with status_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in sorted(results, key=lambda x: (x["pair"], x["target"], x["ligand"])):
            w.writerow({k: r.get(k, "") for k in fields})

    score_path = tab / "scores_vina_mode1_v1.csv"
    with score_path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["pair", "target", "ligand", "mode1_energy", "score_S", "status"],
        )
        w.writeheader()
        for r in sorted(results, key=lambda x: (x["pair"], x["target"], x["ligand"])):
            if r["status"] in ("success", "exists") and r.get("mode1_energy") is not None:
                w.writerow(
                    {
                        "pair": r["pair"],
                        "target": r["target"],
                        "ligand": r["ligand"],
                        "mode1_energy": r["mode1_energy"],
                        "score_S": r["score_S"],
                        "status": r["status"],
                    }
                )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=600, help="per-job Vina timeout seconds; timeout → skip")
    ap.add_argument("--pairs", nargs="*", default=None, help="optional subset, e.g. 'F2/F10'")
    ap.add_argument("--exhaustiveness", type=int, default=EXHAUSTIVENESS)
    args = ap.parse_args()

    if not LIG_DIR.exists():
        raise SystemExit(f"ligand cache missing: {LIG_DIR} — run prep_track_b_ligands_v1.py first")
    n_lig = len(list(LIG_DIR.glob("*.pdbqt")))
    print(f"ligand pdbqt available: {n_lig}", flush=True)

    pair_filter = set(args.pairs) if args.pairs else None
    jobs = load_jobs(pair_filter)
    print(
        f"jobs={len(jobs)} workers={args.workers} timeout={args.timeout}s E={args.exhaustiveness} seed={SEED}",
        flush=True,
    )

    results: list[dict] = []
    t_all = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {
            ex.submit(
                run_one,
                pair,
                t,
                lig,
                pdbqt,
                box,
                args.exhaustiveness,
                args.timeout,
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
                    f"[{done}/{len(jobs)}] {res['status']} {res['pair']} {res['target']} {res['ligand']}"
                    + (f" reason={res.get('reason','')[:80]}" if res.get("reason") else ""),
                    flush=True,
                )
            if done % 100 == 0:
                write_tables(results)

    write_tables(results)
    ok = sum(1 for r in results if r["status"] in ("success", "exists"))
    skip = sum(1 for r in results if r["status"] == "skip")
    fail = sum(1 for r in results if r["status"] == "fail")
    elapsed = time.time() - t_all
    print(
        f"done ok={ok} skip={skip} fail={fail} total={len(results)} elapsed_h={elapsed/3600:.2f}",
        flush=True,
    )
    print(f"wrote {LOCAL / 'tables' / 'job_status.csv'}", flush=True)
    print(f"wrote {LOCAL / 'tables' / 'scores_vina_mode1_v1.csv'}", flush=True)
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
