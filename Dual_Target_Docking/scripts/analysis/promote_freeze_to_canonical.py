#!/usr/bin/env python3
"""Copy a verified freeze-rebuild directory onto results/canonical.

Refuses to run unless VERIFICATION_REPORT.md contains PASS.
Does not rebuild the submission pack.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FREEZE = ROOT / "results" / "freeze_rebuild"
CANON = ROOT / "results" / "canonical"
PROCESSED = ROOT / "data" / "processed" / "current_score_master.csv"


def main() -> int:
    report = FREEZE / "VERIFICATION_REPORT.md"
    if not report.is_file() or "result: PASS" not in report.read_text(encoding="utf-8"):
        raise SystemExit("refuse to promote: freeze verification is not PASS")
    CANON.mkdir(parents=True, exist_ok=True)
    copied = []
    for src in sorted(FREEZE.glob("*.csv")):
        dest = CANON / src.name
        shutil.copy2(src, dest)
        copied.append(src.name)
    master = FREEZE / "current_score_master.csv"
    if master.is_file():
        PROCESSED.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(master, PROCESSED)
    md = FREEZE / "score_master_migration_diff.md"
    if md.is_file():
        shutil.copy2(md, CANON / md.name)
    print("promoted", len(copied), "csv files to", CANON)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
