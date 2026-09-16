#!/usr/bin/env python3
"""Clean-room smoke for the post-fix V4 release pipeline.

Uses a fresh git clone (or a worktree copy if the tree is dirty) and does not
re-run Vina/GNINA/RTM. Steps:

  install requirements-analysis (best effort, then env gate)
  post-fix smoke
  canonical validation
  optional figure rebuild from frozen CSVs
  submission-pack rebuild
  checksum + pack hash-parity
  git diff must be clean or limited to expected generated files

Exit 1 on any step failure. Does not modify the source tree.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parent
SKIP_DIRS = {
    ".git",
    "__pycache__",
    ".submission_pack_postfix.staging",
    ".submission_pack.staging",
}
EXPECTED_GENERATED = (
    "Dual_Target_Docking/submission_pack/",
    "Dual_Target_Docking/submission_pack_postfix/",
    "Dual_Target_Docking/docs/SUBMISSION_AUDIT_FIVE_ROUNDS_V2_POSTFIX.md",
    "Dual_Target_Docking/docs/FINAL_POSTFIX_SUBMISSION_AUDIT.md",
    "Dual_Target_Docking/data/manuscript_lock/CANONICAL_CHECKSUMS_V2.csv",
    "Dual_Target_Docking/data/manuscript_lock/HEAVY_OUTPUT_LOCK_V4.csv",
)


def run(cmd: list[str], cwd: Path, env: dict | None = None) -> None:
    print("+", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, cwd=str(cwd), env=env)
    if proc.returncode != 0:
        raise SystemExit(proc.returncode)


def git_out(args: list[str], cwd: Path) -> str:
    return subprocess.check_output(["git", *args], cwd=str(cwd), text=True).strip()


def copy_tree(src: Path, dest: Path) -> None:
    dest.mkdir(parents=True)
    for root, dirs, files in os.walk(src):
        rel = Path(root).relative_to(src)
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS and not d.endswith(".egg-info")]
        (dest / rel).mkdir(parents=True, exist_ok=True)
        for name in files:
            if name.endswith(".pyc"):
                continue
            shutil.copy2(Path(root) / name, dest / rel / name)


def resolve_snapshot() -> str:
    sys.path.insert(0, str(ROOT / "scripts" / "audit"))
    from canonical_paths_v4 import git_rev  # noqa: E402

    snapshot = git_rev()
    if snapshot:
        return snapshot
    readme = ROOT / "submission_pack" / "README.md"
    if readme.is_file():
        for line in readme.read_text(encoding="utf-8").splitlines():
            if "source_snapshot_commit:" in line and "`" in line:
                return line.split("`")[1].strip()
    print("FAIL: cannot resolve source_snapshot_commit for clean-room pack rebuild")
    raise SystemExit(1)


def prepare_clone(tmp: Path) -> Path:
    clone = tmp / "repo"
    print("+ git clone --local", REPO, clone, flush=True)
    proc = subprocess.run(["git", "clone", "--local", str(REPO), str(clone)])
    if proc.returncode != 0:
        raise SystemExit(proc.returncode)
    branch = git_out(["rev-parse", "--abbrev-ref", "HEAD"], REPO)
    if branch and branch != "HEAD":
        subprocess.run(["git", "checkout", branch], cwd=str(clone), check=False)
    dirty = git_out(["status", "--porcelain"], REPO)
    if dirty:
        # Include uncommitted Dual_Target_Docking files so pre-commit runs work.
        work = clone / "Dual_Target_Docking"
        if work.exists():
            shutil.rmtree(work)
        copy_tree(ROOT, work)
        print("WARN: source tree is dirty; clean-room copied the working Dual_Target_Docking tree")
        print("git-diff gate compares pre- vs post-rebuild status, not pre-existing dirty files.")
    return clone


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--skip-install", action="store_true")
    args = parser.parse_args()

    tmp = Path(tempfile.mkdtemp(prefix="postfix-v4-cleanroom-"))
    print(f"clean-room workdir: {tmp}")
    snapshot = resolve_snapshot()
    # Prefer the locked pack snapshot, not the packaging HEAD.
    readme = ROOT / "submission_pack" / "README.md"
    if readme.is_file():
        for line in readme.read_text(encoding="utf-8").splitlines():
            if "source_snapshot_commit:" in line and "`" in line:
                snapshot = line.split("`")[1].strip()
                break
    try:
        clone = prepare_clone(tmp)
        work = clone / "Dual_Target_Docking"
        py = sys.executable
        if not args.skip_install:
            req = work / "requirements-analysis.txt"
            proc = subprocess.run([py, "-m", "pip", "install", "-r", str(req)], cwd=str(work))
            if proc.returncode != 0:
                print("WARN: pinned requirements-analysis install failed; env gate will decide compatibility")
        env = os.environ.copy()
        env["PYTHONPATH"] = str(work / "scripts" / "audit")
        before_status = git_out(["status", "--porcelain"], clone)
        run([py, "scripts/check_analysis_env.py"], work, env)
        run([py, "scripts/audit/smoke_postfix_v4.py"], work, env)
        run([py, "scripts/audit/validate_postfix_v4.py"], work, env)
        if not args.skip_figures:
            run(
                [
                    py,
                    "figures/jcim_article/scripts/update_figures_pr32.py",
                    "--source-root",
                    str(work),
                ],
                work,
                env,
            )
        run(
            [
                py,
                "scripts/build_submission_pack_postfix_v4.py",
                "--skip-promote",
                "--source-snapshot",
                snapshot,
            ],
            work,
            env,
        )
        run([py, "scripts/audit/check_canonical_checksums_v2.py"], work, env)
        run([py, "scripts/audit/check_submission_pack_v4.py"], work, env)

        sys.path.insert(0, str(ROOT / "scripts" / "audit"))
        from canonical_paths_v4 import sha256_raw  # noqa: E402

        mismatches = []
        for rel in [
            "submission_pack/manuscript/MANUSCRIPT_JCIM_EN.md",
            "submission_pack/SI/SUPPORTING_INFORMATION_JCIM_EN_V1.md",
            "submission_pack/figures/plotted_values_postfix.json",
            "submission_pack/tables/source/post_fix_master_metrics.csv",
            "submission_pack/tables/source/3POZ_box_corrected.json",
        ]:
            a = ROOT / rel
            b = work / "submission_pack_postfix" / rel.split("submission_pack/", 1)[1]
            if not a.is_file() or not b.is_file():
                mismatches.append(f"missing {rel}")
                continue
            if sha256_raw(a) != sha256_raw(b):
                mismatches.append(rel)
        if mismatches:
            print("FAIL: clean-room pack hash mismatch:\n  " + "\n".join(mismatches))
            return 1

        after_status = git_out(["status", "--porcelain"], clone)
        before_lines = set(before_status.splitlines()) if before_status else set()
        after_lines = set(after_status.splitlines()) if after_status else set()
        new_lines = sorted(after_lines - before_lines)
        bad = []
        for line in new_lines:
            path = line[3:].strip()
            if " -> " in path:
                path = path.split(" -> ", 1)[1]
            if any(path == rel or path.startswith(rel) for rel in EXPECTED_GENERATED):
                continue
            bad.append(line)
        if bad:
            print("FAIL: rebuild created git diff outside the generated allowlist:")
            print("\n".join(bad))
            return 1
        if new_lines:
            print("clean-room git diff after rebuild is limited to expected generated files:")
            print("\n".join(new_lines))
        elif after_status:
            print("clean-room rebuild left no additional git diff (pre-existing dirty files unchanged)")
        else:
            print("clean-room git working tree is clean")
        print("clean-room smoke: provenance/hash parity PASS")
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
