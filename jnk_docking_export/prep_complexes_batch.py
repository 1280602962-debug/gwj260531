#!/usr/bin/env python3
"""
Run Protein Preparation Wizard on exported receptor-ligand complexes (maegz).

Requires Schrödinger: %SCHRODINGER%\\run.exe python3 prep_complexes_batch.py

Ligand atoms are fixed during prep (default ASL: ligand) so docking poses stay put
while protein receives hydrogens and H-bond network optimization.
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


def sanitize(name: str) -> str:
    return re.sub(r"\s+", "_", name.strip())


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


def find_schrodinger() -> Path:
    env = os.environ.get("SCHRODINGER")
    if not env:
        raise RuntimeError("SCHRODINGER environment variable is not set.")
    root = Path(env)
    prep = root / "prepwizard"
    if os.name == "nt":
        prep = root / "prepwizard.exe"
    if not prep.exists():
        prep = root / "utilities" / "prepwizard"
    if not prep.exists():
        raise RuntimeError(f"prepwizard not found under {root}")
    return root


def run_prepwizard(
    schrod: Path,
    in_file: Path,
    out_file: Path,
    *,
    epik_pH: float,
    fix_asl: str,
    extra_flags: list[str],
    wait: bool,
) -> None:
    run_exe = schrod / "run.exe" if os.name == "nt" else schrod / "run"
    prep_cmd = schrod / "prepwizard"
    if os.name == "nt":
        prep_cmd = schrod / "prepwizard.exe"

    cmd = [
        str(run_exe if run_exe.exists() else prep_cmd),
    ]
    if run_exe.exists():
        cmd.append("prepwizard")
    cmd += [
        str(in_file),
        str(out_file),
        "-epik_pH",
        str(epik_pH),
        "-fix",
        fix_asl,
        "-WAIT" if wait else "-NOWAIT",
    ]
    cmd.extend(extra_flags)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    print("  CMD:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def main() -> int:
    parser = argparse.ArgumentParser(description="Batch Protein Prep on exported complexes.")
    parser.add_argument("--config", default="jobs_step1.json")
    parser.add_argument("--in-dir", default="complexes_mae", help="Raw exported maegz directory")
    parser.add_argument("--out-dir", default="complexes_prepped", help="Prepped output directory")
    parser.add_argument("--summary", default="complexes_prepped/prep_summary.tsv")
    parser.add_argument("--no-wait", action="store_true", help="Submit prep jobs without -WAIT")
    args = parser.parse_args()

    cfg = load_config(Path(args.config).resolve())
    options = cfg.get("options", {})
    prep_cfg = cfg.get("prep", {})
    include = options.get("ligand_include", [])
    exclude = options.get("ligand_exclude", ["JNK-IN-8"])

    root = Path(cfg.get("root", ".")).resolve()
    in_root = Path(args.in_dir)
    if not in_root.is_absolute():
        in_root = root / in_root
    out_root = Path(args.out_dir)
    if not out_root.is_absolute():
        out_root = root / out_root

    schrod = find_schrodinger()
    maegz_files = sorted(in_root.glob("**/*.maegz"))
    if not maegz_files:
        print(f"ERROR: no .maegz under {in_root}", file=sys.stderr)
        return 1

    rows: list[dict] = []
    errors: list[str] = []
    done = 0

    print(f"Input  : {in_root}")
    print(f"Output : {out_root}")
    print(f"Ligands: include={include}, exclude={exclude}")
    print()

    for in_path in maegz_files:
        rel = in_path.relative_to(in_root)
        parts = in_path.stem.split("_", 1)
        pdb_id = parts[0] if parts else in_path.stem
        lig_part = parts[1] if len(parts) > 1 else in_path.stem
        if not ligand_matches(lig_part, include, exclude):
            print(f"SKIP {rel} (ligand filter)")
            continue

        out_path = out_root / rel.parent / f"{in_path.stem}_prepped.maegz"
        print(f"PREP {rel} -> {out_path.relative_to(out_root)}")
        try:
            run_prepwizard(
                schrod,
                in_path,
                out_path,
                epik_pH=float(prep_cfg.get("epik_pH", 7.0)),
                fix_asl=str(prep_cfg.get("fix_ligand_asl", "ligand")),
                extra_flags=list(prep_cfg.get("prepwizard_flags", [])),
                wait=not args.no_wait,
            )
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "ligand": lig_part,
                    "input_file": str(rel),
                    "output_file": str(out_path.relative_to(out_root)),
                    "status": "ok",
                }
            )
            done += 1
        except subprocess.CalledProcessError as exc:
            msg = f"{rel}: prepwizard failed ({exc})"
            print(f"  ERROR: {msg}")
            errors.append(msg)
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "ligand": lig_part,
                    "input_file": str(rel),
                    "output_file": "",
                    "status": "error",
                }
            )

    summary_path = Path(args.summary)
    if not summary_path.is_absolute():
        summary_path = root / summary_path
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["pdb_id", "ligand", "input_file", "output_file", "status"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print()
    print(f"Prepped: {done}/{len(maegz_files)} files (after ligand filter)")
    print(f"Summary: {summary_path}")
    if errors:
        print(f"Errors : {len(errors)}")
        return 2
    return 0 if done else 1


if __name__ == "__main__":
    sys.exit(main())
