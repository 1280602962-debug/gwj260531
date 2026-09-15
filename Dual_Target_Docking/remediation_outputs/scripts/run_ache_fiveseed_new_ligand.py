#!/usr/bin/env python3
"""Extra-seed Vina for the replacement AChE A_only ligand AB_056.

Boxes unchanged. Seeds 20260811–20260814; production seed already docked.
"""
from __future__ import annotations

import csv
import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
OUT = ROOT / "remediation_outputs/phase_ache_fiveseed"
REC = ROOT / "data/ache_bche_panel_v0/receptors"
BOX = ROOT / "data/ache_bche_panel_v0/boxes"
LIG = ROOT / "remediation_outputs/phase2_ache_bche/vina_new_ligand/ligands_pdbqt/AB_056.pdbqt"
VINA = "/home/gwj/miniconda3/bin/vina"
SEEDS = [20260811, 20260812, 20260813, 20260814]
TARGETS = (("ACHE", "4EY7"), ("BCHE", "4BDS"))


def parse_mode1(path: Path) -> float | None:
    m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)", path.read_text(errors="replace"))
    return float(m.group(1)) if m else None


def dock(target_key: str, pdb: str, seed: int) -> dict:
    box = json.loads((BOX / f"{target_key}_box.json").read_text())
    rec = REC / f"{target_key}_receptor.pdbqt"
    out_dir = OUT / "out" / f"seed_{seed}"
    log_dir = OUT / "logs" / f"seed_{seed}"
    conf_dir = OUT / "confs" / f"seed_{seed}"
    for d in (out_dir, log_dir, conf_dir):
        d.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{pdb}_AB_056_out.pdbqt"
    log = log_dir / f"{pdb}_AB_056.log"
    conf = conf_dir / f"{pdb}_AB_056.txt"
    if out.exists() and parse_mode1(out) is not None:
        return {
            "pair": "AChE/BChE",
            "seed": seed,
            "ligand": "AB_056",
            "class": "A_only",
            "pdb": pdb,
            "status": "cached",
            "vina_mode1": parse_mode1(out),
        }
    conf.write_text(
        "\n".join(
            [
                f"receptor = {rec}",
                f"ligand = {LIG}",
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                "exhaustiveness = 8",
                "num_modes = 9",
                "energy_range = 3",
                "cpu = 1",
                f"seed = {seed}",
                f"out = {out}",
            ]
        )
        + "\n"
    )
    proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True)
    log.write_text((proc.stdout or "") + (proc.stderr or ""))
    aff = parse_mode1(out) if out.exists() else None
    return {
        "pair": "AChE/BChE",
        "seed": seed,
        "ligand": "AB_056",
        "class": "A_only",
        "pdb": pdb,
        "status": "ok" if aff is not None else f"FAIL rc={proc.returncode}",
        "vina_mode1": aff if aff is not None else "",
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(k, pdb, seed) for seed in SEEDS for k, pdb in TARGETS]
    rows = []
    print(f"ache extra-seed jobs={len(jobs)}", flush=True)
    with ThreadPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(dock, k, pdb, seed) for k, pdb, seed in jobs]
        for fut in as_completed(futs):
            rec = fut.result()
            rows.append(rec)
            print(rec, flush=True)
    dest = OUT / "scores_vina_mode1_AB056_extra_seeds.csv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["seed"], r["pdb"])))
    print("wrote", dest)
    return 0 if all(r["status"] in {"ok", "cached"} for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
