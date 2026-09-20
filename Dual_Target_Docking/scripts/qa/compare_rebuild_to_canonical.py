#!/usr/bin/env python3
"""Compare a freeze rebuild directory to committed results/canonical CSVs.

Hash-free. Identifiers, memberships, labels, counts, and canonical summary
tables are exact text. The only numeric tolerance is
`ecfp4_oof_predictions.csv:oof_prob` at abs <= 1e-5, for machine-level
LogisticRegression floating-point variation. score_master_migration_diff.md
is not a CSV equality gate. ENV.txt / REBUILD_LOG.txt / VERIFICATION_REPORT.md
are ignored.
"""
from __future__ import annotations

import argparse
import csv
import math
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
    "fold_id",
)
OOF_PROB_FILE = "ecfp4_oof_predictions.csv"
OOF_PROB_COL = "oof_prob"
OOF_PROB_TOL = 1e-5


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


def parse_finite(value: str):
    text = str(value).strip()
    if text == "":
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def compare_file(
    name: str,
    rebuilt: Path,
    canonical: Path,
    diffs: list[str],
    tolerated: list[dict],
) -> tuple[int, int]:
    rh, rrows = read_csv(rebuilt / name)
    ch, crows = read_csv(canonical / name)
    if rh != ch:
        diffs.append(f"{name}: header rebuilt={rh} canonical={ch}")
        return 0, 0
    if len(rrows) != len(crows):
        diffs.append(f"{name}: row count rebuilt={len(rrows)} canonical={len(crows)}")
    n = min(len(rrows), len(crows))
    exact = 0
    for i in range(n):
        rr, cr = rrows[i], crows[i]
        key = row_key(ch, cr, i)
        for col in ch:
            rv, cv = rr.get(col, ""), cr.get(col, "")
            if rv == cv:
                exact += 1
                continue
            if name == OOF_PROB_FILE and col == OOF_PROB_COL:
                rf, cf = parse_finite(rv), parse_finite(cv)
                if rf is not None and cf is not None:
                    abs_diff = abs(rf - cf)
                    if abs_diff <= OOF_PROB_TOL:
                        rec = {
                            "file": name,
                            "key": key,
                            "column": col,
                            "canonical": cv,
                            "rebuilt": rv,
                            "abs_diff": abs_diff,
                        }
                        tolerated.append(rec)
                        print(
                            f"TOLERATED: {name} key={key} {col} "
                            f"canonical={cv} rebuilt={rv} abs_diff={abs_diff:.1e}"
                        )
                        continue
                    diffs.append(
                        f"{name}\tkey={key}\tcolumn={col}\t"
                        f"canonical={cv!r}\trebuilt={rv!r}\tabs_diff={abs_diff}"
                    )
                    continue
            diffs.append(
                f"{name}\tkey={key}\tcolumn={col}\tcanonical={cv!r}\trebuilt={rv!r}"
            )
    return exact, n * len(ch) if ch else 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Hash-free CSV compare of rebuild vs canonical, with OOF-prob tolerance only."
    )
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
    tolerated: list[dict] = []
    exact_cells = 0
    for name in sorted(cnames):
        exact, _ = compare_file(name, rebuilt, canonical, diffs, tolerated)
        exact_cells += exact
    n_tol = len(tolerated)
    max_tol = max((r["abs_diff"] for r in tolerated), default=0.0)
    print(f"exact-match cells: {exact_cells}")
    print(f"tolerated numeric cells: {n_tol}")
    if n_tol:
        print(f"maximum tolerated absolute difference: {max_tol:.1e}")
    else:
        print("maximum tolerated absolute difference: 0")
    if diffs:
        print("FAIL: rebuild CSVs differ from canonical")
        for line in diffs[:200]:
            print(line)
        if len(diffs) > 200:
            print(f"... {len(diffs) - 200} more differences")
        return 1
    print(
        "PASS: all canonical CSVs match; machine-sensitive OOF probabilities "
        "agree within configured tolerance."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
