#!/usr/bin/env python3
"""Analysis-environment gate for DualFourClass-Bench.

Packages that change zero-dock scientific numbers must match
requirements-analysis.txt exactly. A mismatch is FAIL, not COMPATIBLE.

Exit 0: every scientific pin matches.
Exit 1: a required scientific package is missing or the version differs.
"""
from __future__ import annotations

import importlib
import importlib.metadata
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQ = ROOT / "requirements-analysis.txt"

# These packages change Table 2-adjacent chemistry / ML numbers.
SCIENTIFIC_MODULES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "sklearn": "scikit-learn",
    "rdkit": "rdkit",
}
# Artwork only. Reported, never used to accept a scientific mismatch.
FIGURE_MODULES = {
    "matplotlib": "matplotlib",
    "PIL": "Pillow",
}
OPTIONAL_MODULES = {
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


def installed_version(module_name: str, pip_name: str) -> tuple[bool, str, str]:
    try:
        importlib.import_module(module_name)
    except Exception as exc:
        return False, "", str(exc)
    try:
        ver = importlib.metadata.version(pip_name)
    except importlib.metadata.PackageNotFoundError:
        mod = importlib.import_module(module_name)
        ver = str(getattr(mod, "__version__", "import-ok"))
    return True, ver, ""


def print_row(status: str, pip_name: str, have: str, pin: str) -> None:
    print(f"- {status:12} {pip_name:14} have={have or '-'} pin={pin or '-'}")


def main() -> int:
    pins = parse_pins()
    print("== DualFourClass-Bench analysis environment ==")
    print(f"python: {sys.version.split()[0]} ({sys.executable})")
    print(f"pins:   {REQ if REQ.is_file() else 'MISSING requirements-analysis.txt'}")
    print()
    if not pins:
        print("Result: FAIL (requirements-analysis.txt missing or empty)")
        return 1

    print("## Scientific (must match pin)")
    mismatch = []
    missing = []
    for module_name, pip_name in SCIENTIFIC_MODULES.items():
        pin = pins.get(pip_name.lower(), "")
        ok, ver, err = installed_version(module_name, pip_name)
        if not ok:
            print_row("MISSING", pip_name, "", pin)
            if err:
                print(f"    {err}")
            missing.append(pip_name)
            continue
        if not pin:
            print_row("UNPINNED", pip_name, ver, "")
            mismatch.append(pip_name)
        elif ver == pin:
            print_row("PIN_MATCH", pip_name, ver, pin)
        else:
            print_row("MISMATCH", pip_name, ver, pin)
            mismatch.append(pip_name)

    print()
    print("## Figure render (reported only)")
    for module_name, pip_name in FIGURE_MODULES.items():
        pin = pins.get(pip_name.lower(), "")
        ok, ver, err = installed_version(module_name, pip_name)
        if not ok:
            print_row("MISSING", pip_name, "", pin)
            if err:
                print(f"    {err}")
        elif pin and ver == pin:
            print_row("PIN_MATCH", pip_name, ver, pin)
        elif pin:
            print_row("DIFFERS", pip_name, ver, pin)
        else:
            print_row("PRESENT", pip_name, ver, "-")

    print()
    print("## Optional docking helpers (not used by zero-dock chain)")
    for module_name, pip_name in OPTIONAL_MODULES.items():
        pin = pins.get(pip_name.lower(), "")
        ok, ver, _err = installed_version(module_name, pip_name)
        if not ok:
            print_row("ABSENT", pip_name, "", pin)
        elif pin and ver == pin:
            print_row("PIN_MATCH", pip_name, ver, pin)
        elif pin:
            print_row("DIFFERS", pip_name, ver, pin)
        else:
            print_row("PRESENT", pip_name, ver, "-")

    print()
    py = sys.version_info
    py_ok = py >= (3, 10)
    print(f"python compatibility: {'OK' if py_ok else 'FAIL'} (need >= 3.10; CI uses 3.12)")
    if not py_ok:
        print("Result: FAIL")
        return 1
    if missing or mismatch:
        print("Result: FAIL")
        if missing:
            print("Missing scientific packages:", ", ".join(missing))
        if mismatch:
            print("Scientific pin mismatch:", ", ".join(mismatch))
        print("Install from Dual_Target_Docking/requirements-analysis.txt")
        return 1
    print("Result: PIN_MATCH")
    print("Docking binaries are not required for zero-dock rebuild; use scripts/check_docking_env.py.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
