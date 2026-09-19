#!/usr/bin/env python3
"""Replace results/canonical with a verified freeze-rebuild directory.

Refuses to run unless VERIFICATION_REPORT.md contains PASS.
Does not rebuild the submission pack. Does not use checksums.
Canonical CSV filename set must equal the verified rebuild CSV filename set.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.analysis_config import CANONICAL_CSV_NAMES  # noqa: E402

DEFAULT_FREEZE = Path("/tmp/dual_target_freeze_rebuild")
CANON = ROOT / "results" / "canonical"
PROCESSED = ROOT / "data" / "processed" / "current_score_master.csv"


def csv_names(directory: Path) -> set[str]:
    return {p.name for p in directory.glob("*.csv")}


def main() -> int:
    parser = argparse.ArgumentParser(description="Promote a verified freeze directory onto results/canonical.")
    parser.add_argument("--from-dir", default=str(DEFAULT_FREEZE), help="Verified rebuild directory.")
    args = parser.parse_args()
    freeze = Path(args.from_dir)
    if not freeze.is_absolute():
        freeze = (ROOT / freeze).resolve()
    report = freeze / "VERIFICATION_REPORT.md"
    if not report.is_file() or "result: PASS" not in report.read_text(encoding="utf-8"):
        raise SystemExit("refuse to promote: freeze verification is not PASS")
    expected = set(CANONICAL_CSV_NAMES)
    have = csv_names(freeze)
    missing = sorted(expected - have)
    extra = sorted(have - expected)
    if missing:
        raise SystemExit(f"refuse to promote: verified rebuild missing {missing}")
    if extra:
        raise SystemExit(f"refuse to promote: verified rebuild has unexpected CSV {extra}")
    if CANON.exists():
        shutil.rmtree(CANON)
    CANON.mkdir(parents=True)
    copied = []
    for src in sorted(freeze.glob("*.csv")):
        shutil.copy2(src, CANON / src.name)
        copied.append(src.name)
    after = csv_names(CANON)
    if after != have:
        raise SystemExit(f"promotion filename mismatch: canonical={sorted(after)} rebuild={sorted(have)}")
    master = freeze / "current_score_master.csv"
    if master.is_file():
        PROCESSED.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(master, PROCESSED)
    md = freeze / "score_master_migration_diff.md"
    if md.is_file():
        shutil.copy2(md, CANON / md.name)
    print("promoted", len(copied), "csv files to", CANON)
    print("canonical CSV filename set equals verified rebuild")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
