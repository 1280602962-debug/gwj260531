#!/usr/bin/env python3
"""LEGACY / PRE-REMEDIATION — DO NOT USE FOR CURRENT SUBMISSION.

This freeze path still calls audit_submission_five_rounds_v1.py and
pack_submission_v1.py. The current freeze entry is
scripts/audit/freeze_submission_postfix_v4.py.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(rel: str, extra: list[str] | None = None) -> None:
    args = [sys.executable, str(ROOT / rel), *(extra or [])]
    print("+", " ".join(args), flush=True)
    proc = subprocess.run(args, cwd=ROOT)
    if proc.returncode != 0:
        raise SystemExit(proc.returncode)


def main() -> int:
    print("LEGACY / PRE-REMEDIATION — DO NOT USE FOR CURRENT SUBMISSION.", file=sys.stderr)
    print("Use scripts/audit/freeze_submission_postfix_v4.py", file=sys.stderr)
    return 1
    run("docs/assemble_manuscript_en.py")
    run("docs/assemble_manuscript_zh.py")
    run("data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py")
    run("figures/jcim_article/scripts/audit_figures_pr32.py")
    run("scripts/audit/audit_submission_five_rounds_v1.py")
    run("scripts/audit/pack_submission_v1.py")
    run("data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py", ["--check"])
    run("data/jcim_novelty_v0/scripts/validate_revision_v1.py")
    print("submission freeze OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
