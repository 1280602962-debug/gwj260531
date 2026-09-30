#!/usr/bin/env python3
"""Wait until alternative-receptor production is fully accounted, then start GNINA Phase 10-12.

Also closes the alternative experiment stats (already authorized) so that work is not left hanging.
Does not compute ablation AUROC.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STATUS = ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/alt_production_status.json"
PY = sys.executable


def alt_done() -> bool:
    if not STATUS.is_file():
        return False
    d = json.loads(STATUS.read_text())
    if d.get("PRODUCTION_COMPLETE"):
        return True
    return int(d.get("jobs_accounted", 0)) >= int(d.get("jobs_expected", 734)) and int(d.get("remaining", 1)) == 0


def main() -> int:
    print("watching alternative production...", flush=True)
    while not alt_done():
        d = json.loads(STATUS.read_text()) if STATUS.is_file() else {}
        print(f"alt {d.get('jobs_accounted','?')}/{d.get('jobs_expected','?')} remaining={d.get('remaining','?')}", flush=True)
        time.sleep(60)
    print("ALT_PRODUCTION_ACCOUNTING_COMPLETE launching GNINA independent + alt stats", flush=True)
    # close alt stats in parallel-ish: stats first is minutes; then GNINA
    subprocess.run([PY, str(ROOT / "scripts/rerun_v4_2/alt_score_master_and_stats.py")], check=False)
    gnina_py = Path("/home/gwj/miniconda3/bin/python")
    runner = str(gnina_py if gnina_py.is_file() else PY)
    rc = subprocess.run([runner, str(ROOT / "scripts/rerun_v4_2/ablation_phase10_12_gnina.py")]).returncode
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
