#!/usr/bin/env python3
"""Remove the withdrawn PIK3CA/PIK3CB pair from paper-facing artifacts.

Does not touch PIK3CA/mTOR. Does not invent replacement BindingDB numbers.
"""
from __future__ import annotations

import csv
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

PAIR_EXACT = {
    "PIK3CA/PIK3CB",
    "PIK3CA_PIK3CB",
    "PIK3CA-PIK3CB",
}
PAIR_RE = re.compile(r"PIK3CA\s*/\s*PIK3CB|PIK3CA_PIK3CB|PIK3CA-PIK3CB")
TARGET_ONLY = re.compile(r"\bPIK3CB\b")
PDB_2WXF = re.compile(r"\b2WXF\b")

SKIP_DIRS = {
    ".git",
    "node_modules",
    "__pycache__",
    ".venv",
}

DELETE_PATHS = [
    "data/pik3ca_pik3cb_panel_v0",
    "data/jcim_novelty_v0/tables/external_panel_PIK3CA_PIK3CB_v1.csv",
    "data/jcim_novelty_v0/tables/external_panel_PIK3CA_PIK3CB_v1.sha256",
    "data/jcim_bench_v0/tables/assembled_PIK3CA_PIK3CB.csv",
    "data/jcim_structure_robust_v0/scripts/redock_pik3cb_alt_pik3ca_v1.py",
    "data/jcim_structure_robust_v0/scripts/run_pik3cb_receptor_swap_local.sh",
    "data/jcim_structure_robust_v0/tables/pocket_matched_PAB_alt4JPS_v1.csv",
    "data/jcim_structure_robust_v0/tables/pocket_matched_PAB_alt5DXT_v1.csv",
    "data/public_pair_selection/mols_PIK3CB.json",
]


def is_pair_value(value: str) -> bool:
    text = (value or "").strip()
    return text in PAIR_EXACT or text.replace("_", "/") in PAIR_EXACT


def row_is_withdrawn(row: dict) -> bool:
    identity_keys = ("pair", "system", "panel", "name", "target_pair")
    for key in identity_keys:
        if key in row and PAIR_RE.search(str(row.get(key) or "")):
            return True
    pdb = str(row.get("pdb") or row.get("pdb_id") or "").upper()
    if pdb == "2WXF":
        return True
    target = str(row.get("target") or row.get("gene") or "")
    protein = str(row.get("protein") or "")
    if target == "PIK3CB" or protein == "PIK3CB":
        pair = str(row.get("pair") or "")
        if "mTOR" not in pair:
            return True
    return False


def filter_csv(path: Path) -> str | None:
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return None
    if "PIK3CA/PIK3CB" not in text and "PIK3CA_PIK3CB" not in text and "2WXF" not in text:
        return None
    try:
        with path.open(encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None:
                return None
            kept = [row for row in reader if not row_is_withdrawn(row)]
            fields = reader.fieldnames
    except Exception:
        return None
    if not fields:
        return None
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(kept)
    return f"csv filtered {path.relative_to(ROOT)} kept={len(kept)}"


def walk_files():
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        yield path


def main() -> int:
    reports = []
    for rel in DELETE_PATHS:
        path = ROOT / rel
        if path.is_dir():
            shutil.rmtree(path)
            reports.append(f"deleted dir {rel}")
        elif path.is_file():
            path.unlink()
            reports.append(f"deleted file {rel}")

    for path in walk_files():
        if path.suffix.lower() == ".csv":
            note = filter_csv(path)
            if note:
                reports.append(note)

    plotted = ROOT / "figures/jcim_article/plotted_values.json"
    if plotted.exists():
        data = json.loads(plotted.read_text(encoding="utf-8"))

        def scrub(obj):
            if isinstance(obj, dict):
                return {
                    k: scrub(v)
                    for k, v in obj.items()
                    if k not in PAIR_EXACT and not PAIR_RE.search(str(k))
                }
            if isinstance(obj, list):
                out = []
                for item in obj:
                    if isinstance(item, str) and (is_pair_value(item) or PAIR_RE.search(item)):
                        continue
                    out.append(scrub(item))
                return out
            return obj

        plotted.write_text(json.dumps(scrub(data), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        reports.append("scrubbed plotted_values.json")

    log = ROOT / "data/jcim_novelty_v0/analysis/PURGE_PIK3CA_PIK3CB_V1.md"
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(
        "# PIK3CA/PIK3CB purge\n\n" + "\n".join(f"- {r}" for r in reports) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {len(reports)} actions to {log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
