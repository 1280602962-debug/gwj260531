#!/usr/bin/env python3
"""Verify submission_pack copies against source files using raw SHA256.

Requires submission_pack/RELEASE_MANIFEST.csv with:
  path,role,bytes,sha256_raw,source_path,source_sha256_raw,source_snapshot_commit

Pack copy sha256 must equal source sha256 for manuscript/SI/figures/source tables.
README must record source_snapshot_commit, never a confusing 'commit SHA'.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from canonical_paths_v4 import (  # noqa: E402
    DT_ROOT,
    FIGURE_STEMS,
    PACK_SOURCE_TABLES,
    RELEASE_MANIFEST,
    fail,
    is_legacy_path,
    require_exists,
    sha256_raw,
)

PARITY_ROLES = {"manuscript", "si", "lock", "figure", "caption", "plotted", "table", "box_json", "metrics", "audit"}


def load_manifest() -> list[dict[str, str]]:
    path = require_exists(RELEASE_MANIFEST)
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    needed = {
        "path",
        "role",
        "bytes",
        "sha256_raw",
        "source_path",
        "source_sha256_raw",
        "source_snapshot_commit",
    }
    if not rows or needed - set(rows[0].keys()):
        fail(f"{RELEASE_MANIFEST} missing columns {sorted(needed)}")
    commits = {r["source_snapshot_commit"] for r in rows}
    if len(commits) != 1 or not next(iter(commits)):
        fail(f"source_snapshot_commit is ambiguous: {commits}")
    return rows


def expected_pack_files() -> dict[str, str]:
    mapping = {
        "submission_pack/manuscript/MANUSCRIPT_JCIM_EN.md": "docs/MANUSCRIPT_JCIM_EN.md",
        "submission_pack/manuscript/MANUSCRIPT_JCIM_ZH.md": "docs/MANUSCRIPT_JCIM_ZH.md",
        "submission_pack/manuscript/FIGURE_TABLE_LOCK_POSTFIX_V4.md": "docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md",
        "submission_pack/SI/SUPPORTING_INFORMATION_JCIM_EN_V1.md": "docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md",
        "submission_pack/SI/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md": "docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
        "submission_pack/figures/CAPTIONS.md": "figures/jcim_article/CAPTIONS.md",
        "submission_pack/figures/plotted_values_postfix.json": "figures/jcim_article/plotted_values_postfix.json",
        "submission_pack/audit/POST_FIX_AUDIT_REPORT.md": "remediation_outputs/POST_FIX_AUDIT_REPORT.md",
    }
    for name, stem in FIGURE_STEMS.items():
        exts = ("png", "tif") if name == "TOC" else ("pdf", "png", "tif")
        for ext in exts:
            mapping[f"submission_pack/figures/{name}.{ext}"] = f"figures/jcim_article/{stem}.{ext}"
    for rel in PACK_SOURCE_TABLES:
        mapping[f"submission_pack/tables/source/{Path(rel).name}"] = rel
    return mapping


def main() -> int:
    rows = load_manifest()
    by_path = {r["path"]: r for r in rows}
    mapping = expected_pack_files()
    failed = []
    for pack_rel, src_rel in mapping.items():
        pack = require_exists(pack_rel)
        src = require_exists(src_rel)
        if is_legacy_path(src) or is_legacy_path(pack):
            fail(f"pack/source path is legacy: {pack_rel} / {src_rel}")
        pack_hash = sha256_raw(pack)
        src_hash = sha256_raw(src)
        if pack_hash != src_hash:
            failed.append(f"{pack_rel} != {src_rel}")
            continue
        rec = by_path.get(pack_rel)
        if rec is None:
            failed.append(f"{pack_rel} missing from RELEASE_MANIFEST.csv")
            continue
        if rec["sha256_raw"] != pack_hash or rec["source_sha256_raw"] != src_hash:
            failed.append(f"{pack_rel} manifest hashes stale")
            continue
        if rec["source_path"] != src_rel:
            failed.append(f"{pack_rel} source_path {rec['source_path']} != {src_rel}")
            continue
        if rec["bytes"] != str(pack.stat().st_size):
            failed.append(f"{pack_rel} byte size mismatch")
    if failed:
        fail("pack/source hash parity failed:\n  " + "\n  ".join(failed[:40]))

    readme = require_exists("submission_pack/README.md").read_text(encoding="utf-8")
    if "commit SHA:" in readme:
        fail("submission README still uses 'commit SHA:' which is confused with current HEAD")
    if "source_snapshot_commit:" not in readme:
        fail("submission README missing source_snapshot_commit")
    commit = rows[0]["source_snapshot_commit"]
    if commit not in readme:
        fail("README source_snapshot_commit does not match RELEASE_MANIFEST")
    print(f"check_submission_pack_v4: hash parity PASS; source_snapshot_commit={commit}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
