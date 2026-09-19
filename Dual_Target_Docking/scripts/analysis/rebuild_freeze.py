#!/usr/bin/env python3
"""Rebuild retained DualFourClass results into a new directory.

Does not patch manuscripts, copy old CSVs, or rebuild the submission pack.
Failure of any step leaves the freeze directory incomplete; do not fill gaps
from results/canonical.

Default output is /tmp/dual_target_freeze_rebuild. Pass --outdir to override.
Adjudication is part of the chain (not a hidden intermediate).
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
PYTHON = sys.executable
DEFAULT_FREEZE = Path("/tmp/dual_target_freeze_rebuild")

STEPS = (
    ("00_adjudicate", ["analysis/adjudicate_activity_records.py"], False),
    ("01_build_master", ["analysis/build_current_score_master.py"], True),
    ("02_canonical_stats", ["analysis/compute_canonical_results.py"], True),
    ("03_descriptors", ["analysis/compute_descriptor_baselines.py"], True),
    ("04_ecfp4", ["analysis/fit_ecfp4_models.py"], True),
    ("05_class_chemistry", ["analysis/compute_class_chemistry.py"], True),
    ("06_leave_one_document", ["analysis/compute_leave_one_document.py"], True),
    ("07_detectable_effect", ["analysis/compute_detectable_effect.py"], True),
)

REQUIRED_AFTER = {
    "00_adjudicate": (),
    "01_build_master": ("current_score_master.csv",),
    "02_canonical_stats": (
        "primary_summary_min.csv",
        "primary_directional_auroc.csv",
        "two_pocket_mean_ranking.csv",
        "top10_operating_points.csv",
        "and_filter_operating_points.csv",
        "fixed_score_negative_class_delta.csv",
        "computational_robustness.csv",
        "receptor_substitution.csv",
        "max_vs_median_sensitivity.csv",
        "holdout_metrics.csv",
        "cognate_rmsd.csv",
        "five_seed_summary_min.csv",
        "five_seed_fixed_membership_sensitivity.csv",
        "protocol_sensitivity.csv",
        "external_eligibility.csv",
    ),
    "03_descriptors": ("descriptor_baselines.csv", "descriptor_nested_scaffold_cv.csv"),
    "04_ecfp4": ("ecfp4_incremental_information.csv", "ecfp4_oof_predictions.csv", "model_fold_assignments.csv"),
    "05_class_chemistry": ("class_chemistry_summary.csv",),
    "06_leave_one_document": ("leave_one_document_delta.csv",),
    "07_detectable_effect": ("detectable_effect_simulation.csv",),
}

ADJUDICATION_OUTPUTS = (
    ROOT / "data/processed/activity_adjudication/ligand_activity_aggregate_v1.csv",
    ROOT / "data/processed/activity_adjudication/excluded_activity_rows_applied_v1.csv",
    ROOT / "data/processed/activity_adjudication/ligand_status_v1.csv",
)


def run(name: str, rel: list[str], uses_io: bool, freeze: Path, log) -> None:
    cmd = [PYTHON, str(SCRIPTS / rel[0])]
    if uses_io:
        cmd.extend(["--outdir", str(freeze), "--master", str(freeze / "current_score_master.csv")])
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SCRIPTS) + os.pathsep + env.get("PYTHONPATH", "")
    print("RUN", name, " ".join(cmd), flush=True)
    log.write(f"\n## {name}\n`{' '.join(cmd)}`\n")
    proc = subprocess.run(cmd, cwd=str(ROOT), env=env, capture_output=True, text=True)
    log.write(proc.stdout)
    if proc.stderr:
        log.write("\nSTDERR\n")
        log.write(proc.stderr)
    log.flush()
    if proc.returncode != 0:
        raise SystemExit(f"{name} failed with exit {proc.returncode}")
    if name == "00_adjudicate":
        missing = [str(p) for p in ADJUDICATION_OUTPUTS if not p.is_file() or p.stat().st_size == 0]
        if missing:
            raise SystemExit(f"adjudicate did not write {missing}")
        return
    for fname in REQUIRED_AFTER[name]:
        path = freeze / fname
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"{name} did not write {fname}; not filling from canonical")


def write_env(path: Path, freeze: Path) -> None:
    import platform

    lines = [
        f"utc={datetime.now(timezone.utc).isoformat()}",
        f"python={sys.version.replace(chr(10), ' ')}",
        f"executable={sys.executable}",
        f"platform={platform.platform()}",
        f"cwd={ROOT}",
        f"freeze_dir={freeze}",
        "note=package versions below are the analysis environment; this file is not a PASS/FAIL gate",
    ]
    for mod in ("numpy", "pandas", "sklearn", "rdkit", "scipy"):
        try:
            m = __import__(mod if mod != "sklearn" else "sklearn")
            if mod == "rdkit":
                from rdkit import rdBase

                lines.append(f"{mod}={rdBase.rdkitVersion}")
            else:
                lines.append(f"{mod}={getattr(m, '__version__', 'import-ok')}")
        except Exception as exc:
            lines.append(f"{mod}=MISSING:{exc}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Rebuild the zero-dock freeze directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_FREEZE), help="Rebuild directory (default: /tmp/dual_target_freeze_rebuild).")
    args = parser.parse_args()
    freeze = Path(args.outdir)
    if not freeze.is_absolute():
        freeze = (ROOT / freeze).resolve()
    if freeze.exists():
        shutil.rmtree(freeze)
    freeze.mkdir(parents=True)
    write_env(freeze / "ENV.txt", freeze)
    log_path = freeze / "REBUILD_LOG.txt"
    with log_path.open("w", encoding="utf-8") as log:
        log.write("# freeze rebuild log\n")
        log.write("No canonical CSV was copied into this directory.\n")
        env = subprocess.run([PYTHON, str(SCRIPTS / "check_analysis_env.py")], cwd=str(ROOT), capture_output=True, text=True)
        log.write("\n## check_analysis_env.py\n")
        log.write(env.stdout)
        log.write(env.stderr)
        if env.returncode != 0:
            raise SystemExit("analysis environment check failed")
        for name, rel, uses_io in STEPS:
            run(name, rel, uses_io, freeze, log)
    print("freeze rebuild complete", freeze)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
