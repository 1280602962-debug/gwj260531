#!/usr/bin/env python3
"""Compare a freeze rebuild directory to committed results/canonical CSVs.

Hash-free. Text/CSV equality only. score_master_migration_diff.md is not a
CSV equality gate. ENV.txt / REBUILD_LOG.txt / VERIFICATION_REPORT.md are
ignored.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
KEY_COLUMNS = (
    "pair",
    "seed",
    "ligand_id",
    "ligand",
    "estimand",
    "contrast",
    "arm",
    "model",
    "engine",
    "replacement",
    "panel",
    "setting",
    "aggregation",
    "protein",
    "pdb",
    "estimator",
    "true_auroc",
    "analysis_set",
    "label_rule",
    "scaled",
)


def fail(msg: str) -> None:
    print("FAIL:", msg)
    raise SystemExit(1)


def csv_names(directory: Path) -> set[str]:
    return {p.name for p in directory.glob("*.csv")}


def read_csv(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        header = list(reader.fieldnames or [])
        rows = list(reader)
    return header, rows


def row_key(header: list[str], row: dict, index: int) -> str:
    parts = [str(row.get(col, "")) for col in KEY_COLUMNS if col in header]
    if parts:
        return "|".join(parts)
    return f"row[{index}]|" + "|".join(str(row.get(col, "")) for col in header[:4])


def compare_file(name: str, rebuilt: Path, canonical: Path, diffs: list[str]) -> None:
    rh, rrows = read_csv(rebuilt / name)
    ch, crows = read_csv(canonical / name)
    if rh != ch:
        diffs.append(f"{name}: header rebuilt={rh} canonical={ch}")
        return
    if len(rrows) != len(crows):
        diffs.append(f"{name}: row count rebuilt={len(rrows)} canonical={len(crows)}")
    n = min(len(rrows), len(crows))
    for i in range(n):
        rr, cr = rrows[i], crows[i]
        key = row_key(ch, cr, i)
        for col in ch:
            rv, cv = rr.get(col, ""), cr.get(col, "")
            if rv != cv:
                diffs.append(
                    f"{name}\tkey={key}\tcolumn={col}\tcanonical={cv!r}\trebuilt={rv!r}"
                )


def main() -> int:
    parser = argparse.ArgumentParser(description="Hash-free CSV equality of rebuild vs canonical.")
    parser.add_argument("--rebuilt", required=True)
    parser.add_argument("--canonical", required=True)
    args = parser.parse_args()
    rebuilt = Path(args.rebuilt)
    canonical = Path(args.canonical)
    if not rebuilt.is_absolute():
        rebuilt = (ROOT / rebuilt).resolve()
    if not canonical.is_absolute():
        canonical = (ROOT / canonical).resolve()
    if not rebuilt.is_dir():
        fail(f"rebuilt directory missing {rebuilt}")
    if not canonical.is_dir():
        fail(f"canonical directory missing {canonical}")

    rnames = csv_names(rebuilt)
    cnames = csv_names(canonical)
    missing = sorted(cnames - rnames)
    extra = sorted(rnames - cnames)
    if missing:
        fail(f"CSV filename set: rebuilt missing {missing}")
    if extra:
        fail(f"CSV filename set: rebuilt extra {extra}")

    diffs: list[str] = []
    for name in sorted(cnames):
        compare_file(name, rebuilt, canonical, diffs)
    if diffs:
        print("FAIL: rebuild CSVs differ from canonical")
        for line in diffs[:200]:
            print(line)
        if len(diffs) > 200:
            print(f"... {len(diffs) - 200} more differences")
        return 1
    print(f"PASS: {len(cnames)}/{len(cnames)} CSV files identical (header, row count, cell text)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
