#!/usr/bin/env python3
"""Current freeze entry for the post-fix V4 submission.

Allowed steps only:
  1. post-fix validation
  2. current figure generation (frozen CSVs; no redocking)
  3. V4 pack build
  4. raw SHA256 pack / canonical checksum validation
  5. final audit

Does not call audit_submission_five_rounds_v1.py or pack_submission_v1.py.
"""
from __future__ import annotations

import argparse
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


def resolve_snapshot(explicit: str) -> str:
    if explicit.strip():
        return explicit.strip()
    readme = ROOT / "submission_pack" / "README.md"
    if readme.is_file():
        for line in readme.read_text(encoding="utf-8").splitlines():
            if "source_snapshot_commit:" in line and "`" in line:
                return line.split("`")[1].strip()
    print(
        "FAIL: pass --source-snapshot <commit used to generate pack>; "
        "do not stamp the packaging HEAD into the pack README.",
        file=sys.stderr,
    )
    raise SystemExit(1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-figures", action="store_true", help="do not regenerate artwork (CI / no Arial)")
    parser.add_argument("--skip-promote", action="store_true")
    parser.add_argument(
        "--source-snapshot",
        default="",
        help="commit used to generate the pack; defaults to existing pack README, never silent HEAD",
    )
    args = parser.parse_args()
    snapshot = resolve_snapshot(args.source_snapshot)

    run("scripts/check_analysis_env.py")
    run("scripts/audit/smoke_postfix_v4.py")
    run("scripts/audit/validate_postfix_v4.py")
    if not args.skip_figures:
        fig_args = ["--source-root", str(ROOT)]
        if args.dry_run:
            fig_args.append("--dry-run")
        run("figures/jcim_article/scripts/update_figures_pr32.py", fig_args)
    pack_args: list[str] = ["--source-snapshot", snapshot]
    if args.dry_run:
        pack_args.append("--dry-run")
    if args.skip_promote:
        pack_args.append("--skip-promote")
    run("scripts/build_submission_pack_postfix_v4.py", pack_args)
    if args.dry_run:
        print("freeze dry-run complete (canonical outputs not replaced)")
        return 0
    run("scripts/audit/check_canonical_checksums_v2.py", ["--write"])
    run("scripts/audit/check_canonical_checksums_v2.py")
    run("scripts/audit/check_submission_pack_v4.py")
    print("postfix V4 freeze OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
