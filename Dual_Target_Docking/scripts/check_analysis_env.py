#!/usr/bin/env python3
"""Analysis-environment gate for DualFourClass-Bench.

Checks Python and the packages pinned in requirements-analysis.txt.
Does not inspect docking binaries (see check_docking_env.py).

Exit 0: PIN_MATCH or COMPATIBLE (required packages import).
Exit 1: required package missing or unusable.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQ = ROOT / "requirements-analysis.txt"

REQUIRED_MODULES = {
    "numpy": "numpy",
    "pandas": "pandas",
}
OPTIONAL_MODULES = {
    "matplotlib": "matplotlib",
    "scipy": "scipy",
    "sklearn": "scikit-learn",
    "rdkit": "rdkit",
    "meeko": "meeko",
    "gemmi": "gemmi",
}


def parse_pins() -> dict[str, str]:
    pins: dict[str, str] = {}
    if not REQ.is_file():
        return pins
    for line in REQ.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "==" in line:
            name, ver = line.split("==", 1)
            pins[name.strip().lower()] = ver.strip()
    return pins


def module_version(mod) -> str:
    return str(getattr(mod, "__version__", "import-ok"))


def check_one(module_name: str, pip_name: str, pins: dict[str, str]) -> dict:
    rec = {
        "package": pip_name,
        "module": module_name,
        "ok": False,
        "version": "",
        "pinned": pins.get(pip_name.lower(), ""),
        "status": "MISSING",
        "error": "",
    }
    try:
        mod = importlib.import_module(module_name)
    except Exception as exc:
        rec["error"] = str(exc)
        rec["status"] = "MISSING"
        return rec
    rec["ok"] = True
    rec["version"] = module_version(mod)
    pinned = rec["pinned"]
    if pinned and rec["version"] == pinned:
        rec["status"] = "PIN_MATCH"
    elif pinned:
        rec["status"] = "COMPATIBLE"
    else:
        rec["status"] = "PRESENT"
    return rec


def main() -> int:
    pins = parse_pins()
    print("== DualFourClass-Bench analysis environment ==")
    print(f"python: {sys.version.split()[0]} ({sys.executable})")
    print(f"pins:   {REQ if REQ.is_file() else 'MISSING requirements-analysis.txt'}")
    print()

    rows = []
    print("## Required")
    missing = []
    for module_name, pip_name in REQUIRED_MODULES.items():
        rec = check_one(module_name, pip_name, pins)
        rows.append(rec)
        print(f"- {rec['status']:12} {pip_name:14} have={rec['version'] or '-'} pin={rec['pinned'] or '-'}")
        if not rec["ok"]:
            missing.append(pip_name)
            if rec["error"]:
                print(f"    {rec['error']}")

    print()
    print("## Optional / extended analysis")
    for module_name, pip_name in OPTIONAL_MODULES.items():
        rec = check_one(module_name, pip_name, pins)
        rows.append(rec)
        print(f"- {rec['status']:12} {pip_name:14} have={rec['version'] or '-'} pin={rec['pinned'] or '-'}")

    print()
    py = sys.version_info
    py_ok = py >= (3, 10)
    print(f"python compatibility: {'OK' if py_ok else 'FAIL'} (need >= 3.10; CI uses 3.12)")
    if not py_ok:
        print("Result: FAIL")
        return 1
    if missing:
        print("Result: FAIL (required packages missing)")
        print("Install from Dual_Target_Docking/requirements-analysis.txt")
        return 1

    statuses = {r["status"] for r in rows if r["ok"] and r["package"] in REQUIRED_MODULES.values()}
    if statuses <= {"PIN_MATCH"}:
        verdict = "PIN_MATCH"
    else:
        verdict = "COMPATIBLE"
    print(f"Result: {verdict}")
    print("Docking binaries are not required for revision-validate; use scripts/check_docking_env.py.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
