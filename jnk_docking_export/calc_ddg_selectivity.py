#!/usr/bin/env python3
"""
Parse existing Prime MM-GBSA outputs and compute per-ligand ΔΔG across JNK1/2/3.

MM-GBSA jobs are NOT submitted here — only existing Schrödinger results are read.

Usage:
  python3 calc_ddg_selectivity.py --config jobs_mmgbsa.json
  python3 calc_ddg_selectivity.py --config jobs_mmgbsa.json --inventory mmgbsa_results/mmgbsa_inventory.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path

PDB_IDS = ("3ELJ", "4L7F", "3E7O", "3TTI", "4WHZ")


@dataclass
class MmgbsaRecord:
    ligand: str
    pdb_id: str
    energy: float | None
    energy_file: Path
    source_dir: Path
    status: str


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def normalize_ligand(name: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "_", name.upper()).strip("_")


def ligand_matches(name: str, include: list[str], exclude: list[str]) -> bool:
    upper = name.upper()
    for ex in exclude:
        if ex.upper().replace("-", "_") in upper.replace("-", "_"):
            return False
    if not include:
        return True
    norm = upper.replace("-", "_").replace(" ", "_")
    for inc in include:
        token = inc.upper().replace("-", "_").replace(" ", "_")
        if token in norm or norm in token:
            return True
    return False


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

    patterns = [
        r"mmgbsa\s+dg\s+bind\s*[:=]\s*([+-]?\d+\.?\d*)",
        r"dg\s*bind\s*[:=]\s*([+-]?\d+\.?\d*)",
        r"binding\s+energy\s*[:=]\s*([+-]?\d+\.?\d*)",
        r"delta\s+g\s+bind\s*[:=]\s*([+-]?\d+\.?\d*)",
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
                    fl = f.lower()
                    if ("bind" in fl and "dg" in fl) or fl in {"mmgbsa_dg_bind", "dg_bind"}:
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


def find_energy_files(root_dir: Path) -> list[Path]:
    patterns = [
        "*_mmgbsa.csv",
        "*mmgbsa*.csv",
        "*_out.csv",
        "*-out.csv",
        "*_mmgbsa.out",
        "*mmgbsa*.out",
        "*.out",
        "*.log",
        "*.csv",
    ]
    hits: list[Path] = []
    for pattern in patterns:
        hits.extend(root_dir.glob(pattern))
    return sorted(set(hits), key=lambda p: (0 if "mmgbsa" in p.name.lower() else 1, p.name))


def infer_pdb_id(text: str) -> str | None:
    upper = text.upper()
    for pdb in PDB_IDS:
        if pdb in upper:
            return pdb
    return None


def infer_ligand_name(text: str, include: list[str]) -> str | None:
    norm_path = text.upper().replace("-", "_").replace(" ", "_")
    for lig in include:
        token = lig.upper().replace("-", "_").replace(" ", "_")
        if token in norm_path:
            return lig
    return None


def path_should_skip(path: Path, exclude_substrings: list[str]) -> bool:
    low = str(path).lower()
    return any(s.lower() in low for s in exclude_substrings)


def discover_mmgbsa_records(root: Path, cfg: dict) -> list[MmgbsaRecord]:
    existing = cfg.get("existing_results", {})
    ddg_cfg = cfg.get("ddg", {})
    col_candidates = ddg_cfg.get("energy_column_candidates", [])
    include = cfg.get("ligand_include", [])

    search_roots = existing.get("search_roots", [".", "mmgbsa_results"])
    search_globs = existing.get(
        "search_globs",
        ["**/*prime_mmgbsa*", "**/*mmgbsa*", "**/mmgbsa_results/**"],
    )
    exclude_substrings = existing.get("exclude_substrings", ["vsw", "top_5000", "xp_out"])

    candidate_dirs: set[Path] = set()

    out_dir = Path(cfg.get("out_dir", "mmgbsa_results"))
    if not out_dir.is_absolute():
        out_dir = root / out_dir
    if out_dir.exists():
        for job_dir in out_dir.glob("**/*"):
            if job_dir.is_dir() and find_energy_files(job_dir):
                candidate_dirs.add(job_dir)

    for rel_root in search_roots:
        base = Path(rel_root)
        if not base.is_absolute():
            base = root / base
        if not base.exists():
            continue
        for pattern in search_globs:
            for hit in base.glob(pattern):
                if hit.is_dir():
                    candidate_dirs.add(hit)
                elif hit.is_file() and hit.suffix.lower() in {".csv", ".out", ".log", ".tsv"}:
                    candidate_dirs.add(hit.parent)

    records: list[MmgbsaRecord] = []
    seen_energy_files: set[str] = set()

    for job_dir in sorted(candidate_dirs):
        if path_should_skip(job_dir, exclude_substrings):
            continue

        energy_file = None
        for candidate in find_energy_files(job_dir):
            if path_should_skip(candidate, exclude_substrings):
                continue
            energy_file = candidate
            if "mmgbsa" in candidate.name.lower():
                break
        if energy_file is None:
            continue

        key = str(energy_file.resolve())
        if key in seen_energy_files:
            continue
        seen_energy_files.add(key)

        try:
            context = str(job_dir.relative_to(root))
        except ValueError:
            context = str(job_dir)

        pdb_id = infer_pdb_id(context) or infer_pdb_id(energy_file.name) or "UNKNOWN"
        ligand = infer_ligand_name(context, include) or infer_ligand_name(energy_file.name, include)
        if ligand is None:
            ligand = job_dir.name

        energy = parse_energy_from_file(energy_file, col_candidates)
        status = "ok" if energy is not None else "no_energy"
        records.append(
            MmgbsaRecord(
                ligand=ligand,
                pdb_id=pdb_id,
                energy=energy,
                energy_file=energy_file,
                source_dir=job_dir,
                status=status,
            )
        )

    return records


def load_inventory(path: Path, root: Path) -> list[MmgbsaRecord]:
    records: list[MmgbsaRecord] = []
    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            energy_file = root / row["energy_file"]
            source_dir = root / row.get("source_dir", row["energy_file"])
            dg = row.get("dg_bind", "").strip()
            energy = float(dg) if dg else None
            records.append(
                MmgbsaRecord(
                    ligand=row["ligand"],
                    pdb_id=row["pdb_id"].upper(),
                    energy=energy,
                    energy_file=energy_file,
                    source_dir=source_dir,
                    status=row.get("parse_status", "ok"),
                )
            )
    return records


def build_tables(records: list[MmgbsaRecord], cfg: dict) -> tuple[list[dict], dict[str, dict[str, float]], dict]:
    kinase_map = cfg.get("kinase_by_pdb", {})
    include = cfg.get("ligand_include", [])
    exclude = cfg.get("ligand_exclude", ["JNK-IN-8"])

    per_ligand_pdb: dict[str, dict[str, float]] = {}
    detail_rows: list[dict] = []

    for rec in records:
        if rec.energy is None or rec.status != "ok":
            continue
        if not ligand_matches(rec.ligand, include, exclude):
            continue
        if rec.pdb_id == "UNKNOWN":
            continue

        lig_key = normalize_ligand(rec.ligand)
        per_ligand_pdb.setdefault(lig_key, {})[rec.pdb_id] = rec.energy
        detail_rows.append(
            {
                "ligand": rec.ligand,
                "ligand_key": lig_key,
                "pdb_id": rec.pdb_id,
                "kinase": kinase_map.get(rec.pdb_id, ""),
                "dg_bind": f"{rec.energy:.3f}",
                "energy_file": str(rec.energy_file),
                "source_dir": str(rec.source_dir),
            }
        )

    return detail_rows, per_ligand_pdb, cfg.get("ddg", {})


def avg(values: list[float]) -> float | None:
    return statistics.mean(values) if values else None


def fmt(val: float | None) -> str:
    return "" if val is None else f"{val:.3f}"


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute JNK1/2/3 ΔΔG from existing MM-GBSA outputs.")
    parser.add_argument("--config", default="jobs_mmgbsa.json")
    parser.add_argument("--inventory", default="", help="Optional pre-built mmgbsa_inventory.tsv")
    parser.add_argument("--out-prefix", default="mmgbsa_results/ddg_selectivity")
    args = parser.parse_args()

    cfg = load_config(Path(args.config).resolve())
    root = Path(cfg.get("root", ".")).resolve()

    if args.inventory:
        inv_path = Path(args.inventory)
        if not inv_path.is_absolute():
            inv_path = root / inv_path
        records = load_inventory(inv_path, root)
    else:
        records = discover_mmgbsa_records(root, cfg)

    detail_rows, per_ligand_pdb, ddg_cfg = build_tables(records, cfg)

    if not detail_rows:
        print(
            "ERROR: no MM-GBSA energies parsed.\n"
            "  1) Run: python3 scan_mmgbsa_inventory.py --config jobs_mmgbsa.json\n"
            "  2) Edit jobs_mmgbsa.json -> existing_results.search_roots / search_globs\n"
            "  3) Re-run calc_ddg_selectivity.py",
            file=sys.stderr,
        )
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
            fieldnames=[
                "ligand",
                "ligand_key",
                "pdb_id",
                "kinase",
                "dg_bind",
                "energy_file",
                "source_dir",
            ],
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
