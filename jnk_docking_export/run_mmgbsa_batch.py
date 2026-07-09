#!/usr/bin/env python3
"""
Batch Prime MM-GBSA on prepped receptor-ligand complexes.

Requires Schrödinger:
  %SCHRODINGER%\\run.exe python3 run_mmgbsa_batch.py --config jobs_mmgbsa.json

After jobs finish, summarize:
  %SCHRODINGER%\\run.exe python3 calc_ddg_selectivity.py --config jobs_mmgbsa.json
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path


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


def parse_stem(stem: str) -> tuple[str, str]:
    base = stem.replace("_prepped", "")
    parts = base.split("_", 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return base, base


def find_schrodinger() -> Path:
    env = os.environ.get("SCHRODINGER")
    if not env:
        raise RuntimeError("SCHRODINGER environment variable is not set.")
    return Path(env)


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def run_mmgbsa(
    schrod: Path,
    in_file: Path,
    job_dir: Path,
    *,
    ligand_asl: str,
    job_type: str,
    rflexdist: float,
    use_membrane: bool,
    extra_flags: list[str],
    wait: bool,
) -> Path:
    job_dir.mkdir(parents=True, exist_ok=True)
    job_name = in_file.stem
    run_exe = schrod / ("run.exe" if os.name == "nt" else "run")

    cmd = [
        str(run_exe),
        "prime_mmgbsa",
        "-job_type",
        job_type,
        "-ligand",
        ligand_asl,
        "-rflexdist",
        str(rflexdist),
        "-jobname",
        job_name,
        "-odir",
        str(job_dir),
    ]
    if use_membrane:
        cmd.append("-membrane")
    if wait:
        cmd.append("-WAIT")
    else:
        cmd.append("-NOWAIT")
    cmd.extend(extra_flags)
    cmd.append(str(in_file))

    print("  CMD:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return job_dir


def main() -> int:
    parser = argparse.ArgumentParser(description="Batch Prime MM-GBSA for JNK complexes.")
    parser.add_argument("--config", default="jobs_mmgbsa.json")
    parser.add_argument("--complex-dir", default=None, help="Override complexes directory")
    parser.add_argument("--out-dir", default=None, help="Override MM-GBSA output directory")
    parser.add_argument("--no-wait", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    cfg = load_config(Path(args.config).resolve())
    root = Path(cfg.get("root", ".")).resolve()
    complex_dir = Path(args.complex_dir or cfg.get("complex_dir", "complexes_prepped"))
    if not complex_dir.is_absolute():
        complex_dir = root / complex_dir
    out_dir = Path(args.out_dir or cfg.get("out_dir", "mmgbsa_results"))
    if not out_dir.is_absolute():
        out_dir = root / out_dir

    include = cfg.get("ligand_include", [])
    exclude = cfg.get("ligand_exclude", ["JNK-IN-8"])
    mmgbsa_cfg = cfg.get("mmgbsa", {})
    kinase_map = cfg.get("kinase_by_pdb", {})

    files = sorted(complex_dir.glob("**/*_prepped.maegz"))
    if not files:
        files = sorted(complex_dir.glob("**/*.maegz"))
    if not files:
        print(f"ERROR: no complexes under {complex_dir}", file=sys.stderr)
        return 1

    schrod = find_schrodinger()
    rows: list[dict] = []
    errors: list[str] = []
    submitted = 0

    print(f"Complexes: {complex_dir}")
    print(f"Output   : {out_dir}")
    print()

    for path in files:
        pdb_id, lig_name = parse_stem(path.stem)
        if not ligand_matches(lig_name, include, exclude):
            print(f"SKIP {path.name}")
            continue

        kinase = kinase_map.get(pdb_id.upper(), "")
        job_dir = out_dir / pdb_id / sanitize_job(lig_name)
        print(f"MMGBSA {path.name} -> {job_dir.relative_to(out_dir)}")

        if args.dry_run:
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "kinase": kinase,
                    "ligand": lig_name,
                    "input_file": str(path.relative_to(complex_dir)),
                    "job_dir": str(job_dir.relative_to(out_dir)),
                    "status": "dry_run",
                }
            )
            submitted += 1
            continue

        try:
            run_mmgbsa(
                schrod,
                path,
                job_dir,
                ligand_asl=str(cfg.get("ligand_asl", "ligand")),
                job_type=str(mmgbsa_cfg.get("job_type", "REAL_MIN")),
                rflexdist=float(mmgbsa_cfg.get("rflexdist", 5.0)),
                use_membrane=bool(mmgbsa_cfg.get("use_membrane", False)),
                extra_flags=list(mmgbsa_cfg.get("extra_flags", [])),
                wait=not args.no_wait,
            )
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "kinase": kinase,
                    "ligand": lig_name,
                    "input_file": str(path.relative_to(complex_dir)),
                    "job_dir": str(job_dir.relative_to(out_dir)),
                    "status": "submitted",
                }
            )
            submitted += 1
        except subprocess.CalledProcessError as exc:
            msg = f"{path.name}: prime_mmgbsa failed ({exc})"
            print(f"  ERROR: {msg}")
            errors.append(msg)
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "kinase": kinase,
                    "ligand": lig_name,
                    "input_file": str(path.relative_to(complex_dir)),
                    "job_dir": "",
                    "status": "error",
                }
            )

    summary = out_dir / "mmgbsa_jobs.tsv"
    out_dir.mkdir(parents=True, exist_ok=True)
    with summary.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["pdb_id", "kinase", "ligand", "input_file", "job_dir", "status"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print()
    print(f"Submitted: {submitted}")
    print(f"Summary  : {summary}")
    if errors:
        return 2
    return 0 if submitted else 1


def sanitize_job(name: str) -> str:
    return re.sub(r'[<>:"/\\|?*]', "_", name.strip())


if __name__ == "__main__":
    sys.exit(main())
