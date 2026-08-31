#!/usr/bin/env python3
"""Lightweight local environment check for DualFourClass-Bench.

This script distinguishes between:
1. minimal zero-dock analysis dependencies
2. optional heavy docking/rescoring tools

It is intentionally non-invasive: it only reports status and exits non-zero
when required minimal dependencies are missing.
"""
from __future__ import annotations

import importlib
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GNINA_HINT = Path("/mnt/d/CADD paper exercise/gnina/bin/gnina")
VINA_HINT = Path("/home/gwj/miniconda3/bin/vina")
RTM_PY_HINT = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")

REQUIRED_PYTHON = [
    ("numpy", "numpy"),
    ("pandas", "pandas"),
    ("matplotlib", "matplotlib"),
    ("rdkit", "rdkit"),
    ("scipy", "scipy"),
    ("sklearn", "scikit-learn"),
]

OPTIONAL_PYTHON = [
    ("Bio", "biopython"),
]

OPTIONAL_BINARIES = [
    ("vina", str(VINA_HINT)),
    ("gnina", str(GNINA_HINT)),
    ("obabel", "Open Babel in PATH"),
]

OPTIONAL_PATHS = [
    ("RTMScore python", RTM_PY_HINT),
]


def check_module(module_name: str) -> tuple[bool, str]:
    try:
        mod = importlib.import_module(module_name)
    except Exception as exc:  # pragma: no cover - diagnostic path
        return False, str(exc)
    version = getattr(mod, "__version__", None)
    return True, version or "import ok"


def check_binary(name: str) -> tuple[bool, str]:
    path = shutil.which(name)
    if path:
        return True, path
    return False, "not found in PATH"


def check_version(cmd: list[str]) -> str:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
    except Exception as exc:  # pragma: no cover - diagnostic path
        return f"version check failed: {exc}"
    out = (proc.stdout or proc.stderr).strip().splitlines()
    return out[0] if out else "version unavailable"


def main() -> int:
    print("== DualFourClass-Bench local environment check ==")
    print(f"repo: {ROOT}")
    print(f"python: {sys.executable}")
    print()

    missing_required = []

    print("## Required for zero-dock analysis")
    for module_name, pip_name in REQUIRED_PYTHON:
        ok, detail = check_module(module_name)
        mark = "OK" if ok else "MISSING"
        print(f"- {mark:7} {pip_name:14} {detail}")
        if not ok:
            missing_required.append(pip_name)

    print()
    print("## Optional for extended analyses")
    for module_name, pip_name in OPTIONAL_PYTHON:
        ok, detail = check_module(module_name)
        mark = "OK" if ok else "MISSING"
        print(f"- {mark:7} {pip_name:14} {detail}")

    print()
    print("## Optional docking / rescoring tools")
    for name, hint in OPTIONAL_BINARIES:
        ok, detail = check_binary(name)
        mark = "OK" if ok else "MISSING"
        print(f"- {mark:7} {name:14} {detail}")
    for label, path in OPTIONAL_PATHS:
        ok = path.exists()
        mark = "OK" if ok else "MISSING"
        print(f"- {mark:7} {label:14} {path}")

    print()
    print("## Quick version hints")
    print(f"- python      {sys.version.split()[0]}")
    if shutil.which("python3"):
        print(f"- matplotlib  {check_version([sys.executable, '-c', 'import matplotlib; print(matplotlib.__version__)'])}")
        print(f"- pandas      {check_version([sys.executable, '-c', 'import pandas; print(pandas.__version__)'])}")
        print(f"- numpy       {check_version([sys.executable, '-c', 'import numpy; print(numpy.__version__)'])}")

    print()
    if missing_required:
        print("Result: minimal analysis environment is NOT ready.")
        print("Install the missing packages from `requirements-analysis.txt`.")
        return 1

    print("Result: minimal analysis environment is ready.")
    print("Use `bash scripts/run_local_repro.sh` to rebuild the main zero-dock outputs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
