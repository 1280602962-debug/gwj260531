#!/usr/bin/env python3
"""
Scan existing Schrödinger Prime MM-GBSA outputs under Docking folder.

Does NOT submit new jobs. Writes inventory TSV for review / downstream ΔΔG.

Usage (plain Python, no Schrödinger required):
  python3 scan_mmgbsa_inventory.py --config jobs_mmgbsa.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from calc_ddg_selectivity import discover_mmgbsa_records, load_config, ligand_matches


def main() -> int:
    parser = argparse.ArgumentParser(description="Inventory existing MM-GBSA result files.")
    parser.add_argument("--config", default="jobs_mmgbsa.json")
    parser.add_argument("--out", default="mmgbsa_results/mmgbsa_inventory.tsv")
    args = parser.parse_args()

    cfg = load_config(Path(args.config).resolve())
    root = Path(cfg.get("root", ".")).resolve()
    include = cfg.get("ligand_include", [])
    exclude = cfg.get("ligand_exclude", ["JNK-IN-8"])
    kinase_map = cfg.get("kinase_by_pdb", {})

    records = discover_mmgbsa_records(root, cfg)
    rows = []
    for rec in records:
        if not ligand_matches(rec.ligand, include, exclude):
            continue
        rows.append(
            {
                "ligand": rec.ligand,
                "pdb_id": rec.pdb_id,
                "kinase": kinase_map.get(rec.pdb_id, ""),
                "dg_bind": "" if rec.energy is None else f"{rec.energy:.3f}",
                "energy_file": str(rec.energy_file.relative_to(root)),
                "source_dir": str(rec.source_dir.relative_to(root)),
                "parse_status": rec.status,
            }
        )

    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = root / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "ligand",
                "pdb_id",
                "kinase",
                "dg_bind",
                "energy_file",
                "source_dir",
                "parse_status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    ok = sum(1 for r in rows if r["parse_status"] == "ok")
    print(f"Inventory: {out_path}")
    print(f"Records  : {len(rows)} total, {ok} with parsed dG Bind")
    if not rows:
        print("WARNING: no MM-GBSA files matched. Edit jobs_mmgbsa.json existing_results paths.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
