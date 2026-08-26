#!/usr/bin/env python3
"""SHA-256 manifest for manuscript-facing DualFourClass tables.

Default: write the manifest. --check: verify committed hashes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
MANIFEST = TAB / "REVISION_CHECKSUM_MANIFEST_v1.csv"

WATCH = [
    "data/jcim_novelty_v0/tables/high_confidence_summary_v1.csv",
    "data/jcim_novelty_v0/tables/high_confidence_labels_v1.csv",
    "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv",
    "data/jcim_novelty_v0/tables/complete_case_usable_pchembl_overlap_v1.csv",
    "data/jcim_novelty_v0/tables/source_document_concentration_v1.csv",
    "data/jcim_novelty_v0/tables/cognate_rank_rmsd_reaudit_v1.csv",
    "data/jcim_novelty_v0/tables/class_chemistry_summary_v1.csv",
    "data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv",
    "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
    "data/jcim_strengthen_t0t1_v0/tables/ligand_ml_scaffold_vs_random_v1.csv",
    "data/jcim_novelty_v0/tables/document_blocked_cv_summary_v1.csv",
    "data/jcim_novelty_v0/tables/document_blocked_cv_methods_v1.csv",
    "data/jcim_novelty_v0/tables/document_cluster_bootstrap_v1.csv",
    "data/jcim_novelty_v0/tables/assay_context_audit.csv",
    "data/jcim_novelty_v0/tables/assay_context_priority_ligands_v1.csv",
    "data/jcim_novelty_v0/tables/time_split_class_counts_v1.csv",
    "data/jcim_novelty_v0/tables/time_split_auroc_v1.csv",
    "data/jcim_novelty_v0/tables/document_year_lookup_v1.csv",
    "data/jcim_novelty_v0/tables/cognate_artifact_inventory_v1.csv",
    "docs/MANUSCRIPT_JCIM_EN.md",
    "docs/MANUSCRIPT_JCIM_ZH.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_rows() -> list[dict]:
    rows = []
    for relative in WATCH:
        path = ROOT / relative
        rows.append(
            {
                "path": relative.replace("\\", "/"),
                "present": int(path.exists()),
                "nbytes": path.stat().st_size if path.exists() else "",
                "sha256": sha256(path) if path.exists() else "",
            }
        )
    return rows


def write_manifest(rows: list[dict]) -> None:
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "present", "nbytes", "sha256"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def check_manifest() -> int:
    if not MANIFEST.exists():
        print("missing checksum manifest", MANIFEST)
        return 1
    expected = {row["path"]: row for row in csv.DictReader(MANIFEST.open(encoding="utf-8"))}
    current = {row["path"]: row for row in build_rows()}
    failed = []
    for path, row in expected.items():
        got = current.get(path)
        if got is None or got["sha256"] != row["sha256"] or str(got["present"]) != str(row["present"]):
            failed.append(path)
    if failed:
        print("checksum mismatch:")
        for path in failed:
            print(" -", path)
        return 1
    print(f"checksum OK ({len(expected)} files)")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        return check_manifest()
    rows = build_rows()
    write_manifest(rows)
    print(f"wrote {MANIFEST} ({len(rows)} files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
