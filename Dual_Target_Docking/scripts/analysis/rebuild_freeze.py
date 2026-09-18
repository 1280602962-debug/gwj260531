#!/usr/bin/env python3
"""Rebuild retained DualFourClass results into a new directory.

Does not patch manuscripts, copy old CSVs, or rebuild the submission pack.
Failure of any step leaves the freeze directory incomplete; do not fill gaps
from results/canonical.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
FREEZE = ROOT / "results" / "freeze_rebuild"
PYTHON = sys.executable

STEPS = (
    ("01_build_master", ["analysis/build_current_score_master.py"]),
    ("02_canonical_stats", ["analysis/compute_canonical_results.py"]),
    ("03_descriptors", ["analysis/compute_descriptor_baselines.py"]),
    ("04_ecfp4", ["analysis/fit_ecfp4_models.py"]),
    ("05_class_chemistry", ["analysis/compute_class_chemistry.py"]),
    ("06_leave_one_document", ["analysis/compute_leave_one_document.py"]),
    ("07_detectable_effect", ["analysis/compute_detectable_effect.py"]),
)

REQUIRED_AFTER = {
    "01_build_master": ("current_score_master.csv",),
    "02_canonical_stats": (
        "primary_summary_min.csv",
        "primary_directional_auroc.csv",
        "two_pocket_mean_ranking.csv",
        "computational_robustness.csv",
        "receptor_substitution.csv",
        "max_vs_median_sensitivity.csv",
        "holdout_metrics.csv",
        "cognate_rmsd.csv",
    ),
    "03_descriptors": ("descriptor_baselines.csv",),
    "04_ecfp4": ("ecfp4_incremental_information.csv", "ecfp4_oof_predictions.csv", "model_fold_assignments.csv"),
    "05_class_chemistry": ("class_chemistry_summary.csv",),
    "06_leave_one_document": ("leave_one_document_delta.csv",),
    "07_detectable_effect": ("detectable_effect_simulation.csv",),
}


def run(name: str, rel: list[str], log) -> None:
    cmd = [PYTHON, str(SCRIPTS / rel[0]), "--outdir", str(FREEZE), "--master", str(FREEZE / "current_score_master.csv")]
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
    for fname in REQUIRED_AFTER[name]:
        path = FREEZE / fname
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"{name} did not write {fname}; not filling from canonical")


def write_env(path: Path) -> None:
    import platform

    lines = [
        f"utc={datetime.now(timezone.utc).isoformat()}",
        f"python={sys.version.replace(chr(10), ' ')}",
        f"executable={sys.executable}",
        f"platform={platform.platform()}",
        f"cwd={ROOT}",
        f"freeze_dir={FREEZE}",
        f"baseline_sha_expected=17435413410290960da3437de640509705f00b51",
    ]
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(ROOT.parent), text=True).strip()
        lines.append(f"git_head={sha}")
    except Exception as exc:
        lines.append(f"git_head_error={exc}")
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
    if FREEZE.exists():
        shutil.rmtree(FREEZE)
    FREEZE.mkdir(parents=True)
    write_env(FREEZE / "ENV.txt")
    log_path = FREEZE / "REBUILD_LOG.txt"
    with log_path.open("w", encoding="utf-8") as log:
        log.write("# freeze rebuild log\n")
        log.write("No canonical CSV was copied into this directory.\n")
        env = subprocess.run([PYTHON, str(SCRIPTS / "check_analysis_env.py")], cwd=str(ROOT), capture_output=True, text=True)
        log.write("\n## check_analysis_env.py\n")
        log.write(env.stdout)
        log.write(env.stderr)
        if env.returncode != 0:
            raise SystemExit("analysis environment check failed")
        for name, rel in STEPS:
            run(name, rel, log)
    print("freeze rebuild complete", FREEZE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
