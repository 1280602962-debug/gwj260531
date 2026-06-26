#!/usr/bin/env python3
"""
Parse Prime MM-GBSA outputs and compute per-ligand ΔΔG across JNK1/2/3.

Usage:
  %SCHRODINGER%\\run.exe python3 calc_ddg_selectivity.py --config jobs_mmgbsa.json
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
from pathlib import Path


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def normalize_ligand(name: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "_", name.upper()).strip("_")


def parse_stem(stem: str) -> tuple[str, str]:
    base = stem.replace("_prepped", "")
    parts = base.split("_", 1)
    if len(parts) == 2:
        return parts[0].upper(), parts[1]
    return base.upper(), base


def read_text(path: Path) -> str:
    for enc in ("utf-8", "utf-8-sig", "latin-1"):
        try:
            return path.read_text(encoding=enc)
        except UnicodeDecodeError:
            continue
    return path.read_text(errors="ignore")


def parse_energy_from_file(path: Path, column_candidates: list[str]) -> float | None:
    text = read_text(path)
    lower = text.lower()

    # Prime .out / log style: "MMGBSA dG Bind:  -12.34"
    patterns = [
        r"mmgbsa\s+dg\s+bind\s*[:=]\s*([+-]?\d+\.?\d*)",
        r"dg\s*bind\s*[:=]\s*([+-]?\d+\.?\d*)",
        r"binding\s+energy\s*[:=]\s*([+-]?\d+\.?\d*)",
    ]
    for pat in patterns:
        m = re.search(pat, lower, flags=re.IGNORECASE)
        if m:
            return float(m.group(1))

    if path.suffix.lower() in {".csv", ".tsv"}:
        delim = "\t" if path.suffix.lower() == ".tsv" else ","
        with path.open(encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh, delimiter=delim)
            if not reader.fieldnames:
                return None
            fields = {f.lower(): f for f in reader.fieldnames}
            col = None
            for cand in column_candidates:
                key = cand.lower()
                if key in fields:
                    col = fields[key]
                    break
            if col is None:
                for f in reader.fieldnames:
                    if "bind" in f.lower() and "dg" in f.lower():
                        col = f
                        break
            if col is None:
                return None
            for row in reader:
                val = row.get(col, "").strip()
                if val:
                    try:
                        return float(val)
                    except ValueError:
                        continue
    return None


def find_energy_file(job_dir: Path) -> Path | None:
    preferred = [
        "*_mmgbsa.csv",
        "*mmgbsa*.csv",
        "*_out.csv",
        "*.csv",
        "*_mmgbsa.out",
        "*mmgbsa*.out",
        "*.out",
        "*.log",
    ]
    for pattern in preferred:
        hits = sorted(job_dir.glob(pattern))
        if hits:
            return hits[0]
    return None


def collect_rows(root: Path, cfg: dict) -> list[dict]:
    ddg_cfg = cfg.get("ddg", {})
    col_candidates = ddg_cfg.get("energy_column_candidates", [])
    kinase_map = cfg.get("kinase_by_pdb", {})
    out_dir = Path(cfg.get("out_dir", "mmgbsa_results"))
    if not out_dir.is_absolute():
        out_dir = root / out_dir

    per_ligand_pdb: dict[str, dict[str, float]] = {}
    detail_rows: list[dict] = []

    for job_dir in sorted(out_dir.glob("**/*")):
        if not job_dir.is_dir():
            continue
        energy_file = find_energy_file(job_dir)
        if energy_file is None:
            continue

        # job_dir layout: mmgbsa_results/<PDB>/<ligand>/
        try:
            rel = job_dir.relative_to(out_dir)
            pdb_id = rel.parts[0].upper()
            lig_name = rel.parts[1] if len(rel.parts) > 1 else job_dir.name
        except ValueError:
            continue

        energy = parse_energy_from_file(energy_file, col_candidates)
        if energy is None:
            continue

        lig_key = normalize_ligand(lig_name)
        per_ligand_pdb.setdefault(lig_key, {})[pdb_id] = energy
        detail_rows.append(
            {
                "ligand": lig_name,
                "ligand_key": lig_key,
                "pdb_id": pdb_id,
                "kinase": kinase_map.get(pdb_id, ""),
                "dg_bind": f"{energy:.3f}",
                "energy_file": str(energy_file.relative_to(root)),
            }
        )

    return detail_rows, per_ligand_pdb, ddg_cfg


def avg(values: list[float]) -> float | None:
    return statistics.mean(values) if values else None


def fmt(val: float | None) -> str:
    return "" if val is None else f"{val:.3f}"


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute JNK1/2/3 ΔΔG from MM-GBSA outputs.")
    parser.add_argument("--config", default="jobs_mmgbsa.json")
    parser.add_argument("--out-prefix", default="mmgbsa_results/ddg_selectivity")
    args = parser.parse_args()

    cfg = load_config(Path(args.config).resolve())
    root = Path(cfg.get("root", ".")).resolve()
    detail_rows, per_ligand_pdb, ddg_cfg = collect_rows(root, cfg)

    if not detail_rows:
        print("ERROR: no MM-GBSA energies found. Run run_mmgbsa_batch.py first.", file=sys.stderr)
        return 1

    jnk1 = [p.upper() for p in ddg_cfg.get("jnk1_pdbs", ["3ELJ", "4L7F"])]
    jnk2 = [p.upper() for p in ddg_cfg.get("jnk2_pdbs", ["3E7O"])]
    jnk3 = [p.upper() for p in ddg_cfg.get("jnk3_pdbs", ["3TTI", "4WHZ"])]

    summary_rows: list[dict] = []
    for lig_key, pdb_map in sorted(per_ligand_pdb.items()):
        e1 = avg([pdb_map[p] for p in jnk1 if p in pdb_map])
        e2 = avg([pdb_map[p] for p in jnk2 if p in pdb_map])
        e3 = avg([pdb_map[p] for p in jnk3 if p in pdb_map])

        d12 = (e1 - e2) if (e1 is not None and e2 is not None) else None
        d13 = (e1 - e3) if (e1 is not None and e3 is not None) else None
        d23 = (e2 - e3) if (e2 is not None and e3 is not None) else None

        # Display name: first seen ligand string
        lig_display = next(r["ligand"] for r in detail_rows if r["ligand_key"] == lig_key)

        summary_rows.append(
            {
                "ligand": lig_display,
                "dg_bind_jnk1_mean": fmt(e1),
                "dg_bind_jnk2_mean": fmt(e2),
                "dg_bind_jnk3_mean": fmt(e3),
                "ddg_jnk1_minus_jnk2": fmt(d12),
                "ddg_jnk1_minus_jnk3": fmt(d13),
                "ddg_jnk2_minus_jnk3": fmt(d23),
                "n_jnk1_structures": str(sum(1 for p in jnk1 if p in pdb_map)),
                "n_jnk2_structures": str(sum(1 for p in jnk2 if p in pdb_map)),
                "n_jnk3_structures": str(sum(1 for p in jnk3 if p in pdb_map)),
            }
        )

    out_prefix = Path(args.out_prefix)
    if not out_prefix.is_absolute():
        out_prefix = root / out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    detail_path = out_prefix.with_name(out_prefix.name + "_detail.tsv")
    summary_path = out_prefix.with_suffix(".tsv")

    with detail_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["ligand", "ligand_key", "pdb_id", "kinase", "dg_bind", "energy_file"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(detail_rows)

    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "ligand",
                "dg_bind_jnk1_mean",
                "dg_bind_jnk2_mean",
                "dg_bind_jnk3_mean",
                "ddg_jnk1_minus_jnk2",
                "ddg_jnk1_minus_jnk3",
                "ddg_jnk2_minus_jnk3",
                "n_jnk1_structures",
                "n_jnk2_structures",
                "n_jnk3_structures",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Detail : {detail_path}")
    print(f"Summary: {summary_path}")
    print(f"Ligands: {len(summary_rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
