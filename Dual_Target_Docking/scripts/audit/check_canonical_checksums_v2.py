#!/usr/bin/env python3
"""Raw SHA256 coverage for the V4 canonical publication set.

--write  rebuilds data/manuscript_lock/CANONICAL_CHECKSUMS_V2.csv
--check  (default) requires the committed manifest path set to equal
         canonical_paths_v4.required_paths() and compares bytes + sha256_raw.

normalized_sha256 is stored separately and is never compared to sha256_raw.
Legacy / pre-fix files are not canonical dependencies.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from canonical_paths_v4 import (  # noqa: E402
    CHECKSUM_MANIFEST,
    DT_ROOT,
    HEAVY_OUTPUT_LOCK,
    PACK_SOURCE_TABLES,
    SUPERSEDED_ALLOWED,
    TEXT_SUFFIXES,
    atomic_write_text,
    fail,
    is_legacy_path,
    path_role,
    required_paths,
    sha256_normalized,
    sha256_raw,
)

WARN: list[str] = []

HEAVY_OUTPUTS = list(PACK_SOURCE_TABLES) + [
    "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv",
    "data/jcim_holdout_v0/tables/holdout_ligand_scores_v1.csv",
    "data/jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv",
    "audit_outputs/cv_leakage_audit.csv",
    "remediation_outputs/POST_FIX_AUDIT_REPORT.md",
]


def write_heavy_lock() -> None:
    lines = [
        "path,role,bytes,sha256_raw,note\n",
    ]
    for rel in HEAVY_OUTPUTS:
        path = DT_ROOT / rel
        if not path.is_file():
            fail(f"heavy output missing: {rel}")
        if is_legacy_path(path):
            fail(f"heavy output is legacy: {rel}")
        lines.append(
            ",".join(
                [
                    rel,
                    path_role(rel),
                    str(path.stat().st_size),
                    sha256_raw(path),
                    "frozen_scientific_output; do_not_redock_unless_input_or_protocol_hash_changes",
                ]
            )
            + "\n"
        )
    atomic_write_text(DT_ROOT / HEAVY_OUTPUT_LOCK, "".join(lines))


def current_rows() -> list[dict[str, str]]:
    rows = []
    for rel in required_paths():
        if is_legacy_path(rel):
            fail(f"required path is a legacy dependency: {rel}")
        path = DT_ROOT / rel
        if not path.is_file():
            fail(f"required file missing: {rel}")
        if is_legacy_path(path):
            fail(f"required file resolves under legacy: {rel} -> {path}")
        raw = sha256_raw(path)
        nbytes = path.stat().st_size
        row = {
            "path": rel.replace("\\", "/"),
            "role": path_role(rel),
            "bytes": str(nbytes),
            "sha256_raw": raw,
            "normalized_sha256": sha256_normalized(path) if path.suffix.lower() in TEXT_SUFFIXES else "",
        }
        rows.append(row)
    return rows


def write_manifest(rows: list[dict[str, str]]) -> None:
    lines = ["path,role,bytes,sha256_raw,normalized_sha256\n"]
    for row in rows:
        lines.append(
            ",".join(
                [
                    row["path"],
                    row["role"],
                    row["bytes"],
                    row["sha256_raw"],
                    row["normalized_sha256"],
                ]
            )
            + "\n"
        )
    atomic_write_text(DT_ROOT / CHECKSUM_MANIFEST, "".join(lines))


def load_manifest() -> list[dict[str, str]]:
    path = DT_ROOT / CHECKSUM_MANIFEST
    if not path.is_file():
        fail(f"missing checksum manifest {CHECKSUM_MANIFEST}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    needed = {"path", "role", "bytes", "sha256_raw"}
    if not rows or needed - set(rows[0].keys()):
        fail(f"{CHECKSUM_MANIFEST} missing columns {needed}")
    return rows


def scan_unexpected() -> None:
    required = set(required_paths())
    publication_dirs = [
        DT_ROOT / "docs",
        DT_ROOT / "figures" / "jcim_article",
        DT_ROOT / "submission_pack",
        DT_ROOT / "remediation_outputs" / "phase1_boxes",
        DT_ROOT / "remediation_outputs" / "canonical_tables",
    ]
    unexpected_fail = []
    for root in publication_dirs:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if is_legacy_path(path):
                continue
            rel = path.relative_to(DT_ROOT).as_posix()
            if rel in required or rel in SUPERSEDED_ALLOWED:
                continue
            if rel.startswith("figures/jcim_article/input_snapshot/"):
                continue
            if rel.startswith("figures/jcim_article/scripts/"):
                continue
            if any(part.startswith(".") for part in path.parts):
                continue
            if rel.startswith("submission_pack/"):
                unexpected_fail.append(rel)
                continue
            if path.suffix.lower() in {".pdf", ".png", ".tif", ".tiff"} and "Fig6" in path.name:
                if rel.startswith("submission_pack/"):
                    unexpected_fail.append(rel)
                else:
                    WARN.append(f"non-canonical leftover artwork (not packed): {rel}")
                continue
            if rel.startswith("docs/") and path.suffix.lower() == ".md" and "LOCK" in path.name and "POSTFIX_V4" not in path.name:
                WARN.append(f"superseded publication-facing lock present: {rel}")
    if unexpected_fail:
        fail("unexpected publication-facing file(s):\n  " + "\n  ".join(unexpected_fail[:40]))


def check() -> int:
    expected_rows = load_manifest()
    expected = {r["path"]: r for r in expected_rows}
    current = {r["path"]: r for r in current_rows()}
    expected_set = set(expected)
    current_set = set(current)
    if expected_set != current_set:
        missing = sorted(current_set - expected_set)
        extra = sorted(expected_set - current_set)
        msg = []
        if missing:
            msg.append("manifest missing current required paths: " + ", ".join(missing[:20]))
        if extra:
            msg.append("manifest has paths not in required set: " + ", ".join(extra[:20]))
        fail("; ".join(msg))
    failed = []
    mixed = []
    for rel, exp in expected.items():
        got = current[rel]
        if got["bytes"] != exp["bytes"] or got["sha256_raw"] != exp["sha256_raw"]:
            failed.append(rel)
        norm = exp.get("normalized_sha256") or ""
        if norm and norm == exp["sha256_raw"] and (DT_ROOT / rel).suffix.lower() in TEXT_SUFFIXES:
            # identical is fine for LF-only files; mixing would be using normalized as raw
            pass
        if "normalized" in exp["sha256_raw"]:
            mixed.append(rel)
    if mixed:
        fail("sha256_raw column contains normalized hashes")
    if failed:
        fail("raw SHA256 / byte-size mismatch:\n  " + "\n  ".join(failed[:40]))
    scan_unexpected()
    for w in WARN:
        print("WARN ", w)
    print(f"check_canonical_checksums_v2: {len(current)} paths PASS (raw SHA256 + bytes)")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.write:
        if args.dry_run:
            rows = current_rows()
            print(f"dry-run: would write {HEAVY_OUTPUT_LOCK} and {len(rows)} checksum rows to {CHECKSUM_MANIFEST}")
            return 0
        write_heavy_lock()
        rows = current_rows()
        write_manifest(rows)
        print(f"wrote {CHECKSUM_MANIFEST} ({len(rows)} paths)")
        return 0
    return check()


if __name__ == "__main__":
    raise SystemExit(main())
