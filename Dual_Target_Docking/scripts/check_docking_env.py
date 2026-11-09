#!/usr/bin/env python3
"""Docking / rescoring environment check.

Resolves binaries from PATH or environment variables:
  VINA_BIN
  GNINA_BIN
  RTMSCORE_PYTHON
  OBABEL_BIN (optional)

Does not treat /home/gwj or /mnt/d/... as official lookup paths.
This script is optional for GitHub Actions revision-validate.
Missing docking tools => exit 1 only when --require is passed.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def resolve_tool(env_name: str, path_name: str) -> tuple[str, str]:
    env = os.environ.get(env_name, "").strip()
    if env:
        path = Path(env).expanduser()
        if path.is_file() or path.is_dir():
            return str(path), f"{env_name}"
        return "", f"{env_name} set but not found: {env}"
    found = shutil.which(path_name)
    if found:
        return found, "PATH"
    return "", f"not in PATH and {env_name} unset"


def version_line(cmd: list[str]) -> str:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
    except Exception as exc:
        return f"version check failed: {exc}"
    text = (proc.stdout or proc.stderr).strip().splitlines()
    return text[0] if text else "version unavailable"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--require",
        action="store_true",
        help="exit 1 if vina/gnina/rtmscore python are missing",
    )
    args = parser.parse_args()

    print("== DualFourClass-Bench docking environment ==")
    print("Official lookup: PATH or VINA_BIN / GNINA_BIN / RTMSCORE_PYTHON")
    print("Host-specific paths are not official logic.")
    print()

    vina, vina_src = resolve_tool("VINA_BIN", "vina")
    gnina, gnina_src = resolve_tool("GNINA_BIN", "gnina")
    rtm, rtm_src = resolve_tool("RTMSCORE_PYTHON", "python")
    if "RTMSCORE_PYTHON" not in os.environ:
        rtm, rtm_src = "", "RTMSCORE_PYTHON unset (do not assume PATH python is RTMScore)"
    obabel, obabel_src = resolve_tool("OBABEL_BIN", "obabel")

    rows = [
        ("vina", vina, vina_src),
        ("gnina", gnina, gnina_src),
        ("rtmscore_python", rtm, rtm_src),
        ("obabel", obabel, obabel_src),
    ]
    missing_required = []
    for name, path, src in rows:
        ok = bool(path)
        mark = "OK" if ok else "MISSING"
        print(f"- {mark:7} {name:16} {path or src}")
        if ok and name == "vina":
            print(f"         {version_line([path, '--version'])}")
        if name in {"vina", "gnina", "rtmscore_python"} and not ok:
            missing_required.append(name)

    print()
    if args.require and missing_required:
        print("Result: docking environment is NOT ready.")
        return 1
    if missing_required:
        print("Result: OPTIONAL docking tools missing; analysis/CI may still proceed.")
        print("Heavy deposited Vina/GNINA/RTM outputs are treated as frozen inputs; zero-dock validation does not rerun them.")
        return 0
    print("Result: docking environment is ready.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
