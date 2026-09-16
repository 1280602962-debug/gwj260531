#!/usr/bin/env python3
"""Single-source canonical path allowlist for DualFourClass-Bench V4.

Publication/release scripts must import REQUIRED_PATHS from this module.
Legacy / pre-remediation files are listed only so they can be rejected as
canonical dependencies.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
from pathlib import Path

DT_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = DT_ROOT.parent

PRIMARY_PAIRS = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)
WITHDRAWN_PAIRS = ("PIK3CA/PIK3CB",)
ALLOWED_CLASSES = {"dual", "A_only", "B_only", "neither"}
CLASS_SHORT = {"dual": "D", "A_only": "A", "B_only": "B", "neither": "N"}
THETA_PRIMARY = "theta_6.0"

LEGACY_PATH_MARKERS = (
    "/_legacy_archive/",
    "\\_legacy_archive\\",
    "/submission_pack_pre_v4_archive/",
    "\\submission_pack_pre_v4_archive\\",
)

FIGURE_STEMS = {
    "Figure1": "Fig1_four_state_and_supply",
    "Figure2": "Fig2_negative_class_formulation",
    "Figure3": "Fig3_ligand_chemistry",
    "Figure4": "Fig4_mismatched_pocket",
    "Figure5": "Fig5_computational_realization",
    "FigureS1": "FigS1_ligand_chemistry_detail",
    "FigureS2": "FigS2_protocol_sensitivity",
    "FigureS3": "FigS3_cognate_rmsd",
    "FigureS4": "FigS4_label_source_robustness",
    "FigureS5": "FigS5_external_eligibility",
    "TOC": "TOC_graphic",
}

# Official cognate heavy-atom AABB+5 Å / min 20 Å boxes (machine-readable).
OFFICIAL_BOXES = {
    "3POZ": {
        "center_x": 18.816,
        "center_y": 31.837,
        "center_z": 11.725,
        "size_x": 20.931,
        "size_y": 20.0,
        "size_z": 21.342,
        "construction": "cognate_heavy_atom_AABB_plus_5A_min20",
    },
    "3RCD": {
        "center_x": 12.552,
        "center_y": 2.982,
        "center_z": 28.152,
        "size_x": 21.385,
        "size_y": 22.378,
        "size_z": 20.0,
        "construction": "cognate_heavy_atom_AABB_plus_5A_min20",
    },
}

TEXT_SUFFIXES = {".csv", ".md", ".json", ".txt", ".yaml", ".yml", ".tsv", ".py"}

# (relative-to-DT_ROOT, role)
REQUIRED_SPECS: tuple[tuple[str, str], ...] = (
    ("docs/MANUSCRIPT_JCIM_EN.md", "manuscript"),
    ("docs/MANUSCRIPT_JCIM_ZH.md", "manuscript"),
    ("docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md", "si"),
    ("docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md", "si"),
    ("docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md", "lock"),
    ("figures/jcim_article/CAPTIONS.md", "caption"),
    ("figures/jcim_article/MANUSCRIPT_FIGURE_CAPTIONS.md", "caption"),
    ("figures/jcim_article/plotted_values_postfix.json", "plotted"),
    ("remediation_outputs/canonical_tables/post_fix_master_metrics.csv", "metrics"),
    ("remediation_outputs/POST_FIX_AUDIT_REPORT.md", "audit"),
    ("remediation_outputs/phase1_boxes/3POZ_box_corrected.json", "box_json"),
    ("remediation_outputs/phase1_boxes/3RCD_box_corrected.json", "box_json"),
    ("data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv", "table"),
    ("data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv", "table"),
    ("data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv", "table"),
    ("data/jcim_novelty_v0/tables/eight_pair_ranking_operating_point_v1.csv", "table"),
    ("data/jcim_novelty_v0/tables/review_scored_membership_v1.csv", "table"),
    ("data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv", "table"),
    ("data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv", "table"),
    ("data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv", "table"),
    ("data/jcim_novelty_v0/tables/external_slice_summary_v1.csv", "table"),
    ("data/jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv", "table"),
    ("data/jcim_holdout_v0/tables/holdout_ligand_scores_v1.csv", "table"),
    ("data/jcim_holdout_v0/tables/holdout_panel_HOAB.csv", "table"),
    ("data/jcim_holdout_v0/tables/holdout_panel_HOPM.csv", "table"),
    ("data/ache_bche_panel_v0/tables/panel_v0_strict.csv", "table"),
    ("audit_outputs/cv_leakage_audit.csv", "table"),
    ("requirements-analysis.txt", "env_spec"),
    ("figures/jcim_article/scripts/jcim_figure_style.py", "source"),
)

PACK_SOURCE_TABLES = (
    "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
    "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv",
    "data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv",
    "data/jcim_novelty_v0/tables/eight_pair_ranking_operating_point_v1.csv",
    "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv",
    "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv",
    "data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv",
    "data/jcim_novelty_v0/tables/external_slice_summary_v1.csv",
    "remediation_outputs/canonical_tables/post_fix_master_metrics.csv",
    "remediation_outputs/phase1_boxes/3POZ_box_corrected.json",
    "remediation_outputs/phase1_boxes/3RCD_box_corrected.json",
    "figures/jcim_article/plotted_values_postfix.json",
)

CHECKSUM_MANIFEST = "data/manuscript_lock/CANONICAL_CHECKSUMS_V2.csv"
RELEASE_MANIFEST = "submission_pack/RELEASE_MANIFEST.csv"
HEAVY_OUTPUT_LOCK = "data/manuscript_lock/HEAVY_OUTPUT_LOCK_V4.csv"

# Files that exist as superseded companions; WARN if present, never a canonical dependency.
SUPERSEDED_ALLOWED = {
    "docs/FIGURE_PANEL_LOCK_V3.md",
    "docs/SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md",
    "docs/FIGURE_PLACEMENT_LOCK_V1.md",
    "figures/jcim_article/plotted_values.json",
    "POST_FIX_AUDIT_REPORT.md",
}


def required_paths() -> list[str]:
    paths: list[str] = [rel for rel, _ in REQUIRED_SPECS]
    for stem in FIGURE_STEMS.values():
        exts = ("png", "tif") if stem == "TOC_graphic" else ("pdf", "png", "tif")
        for ext in exts:
            paths.append(f"figures/jcim_article/{stem}.{ext}")
    paths.extend(pack_required_paths())
    paths.append(HEAVY_OUTPUT_LOCK)
    # Keep unique order.
    seen: set[str] = set()
    out: list[str] = []
    for p in paths:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def path_role(rel: str) -> str:
    for spec_rel, role in REQUIRED_SPECS:
        if spec_rel == rel:
            return role
    if rel.startswith("figures/jcim_article/Fig") or rel.startswith("figures/jcim_article/TOC"):
        return "figure"
    if rel.startswith("submission_pack/"):
        return "pack"
    if rel == CHECKSUM_MANIFEST:
        return "checksum_manifest"
    if rel == HEAVY_OUTPUT_LOCK:
        return "heavy_lock"
    return "other"


def pack_required_paths() -> list[str]:
    paths = [
        "submission_pack/README.md",
        "submission_pack/RELEASE_MANIFEST.csv",
        "submission_pack/manuscript/MANUSCRIPT_JCIM_EN.md",
        "submission_pack/manuscript/MANUSCRIPT_JCIM_ZH.md",
        "submission_pack/manuscript/FIGURE_TABLE_LOCK_POSTFIX_V4.md",
        "submission_pack/SI/SUPPORTING_INFORMATION_JCIM_EN_V1.md",
        "submission_pack/SI/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
        "submission_pack/figures/CAPTIONS.md",
        "submission_pack/figures/plotted_values_postfix.json",
        "submission_pack/audit/POST_FIX_AUDIT_REPORT.md",
        "submission_pack/audit/SUBMISSION_AUDIT_FIVE_ROUNDS_V2_POSTFIX.md",
        "submission_pack/audit/FINAL_POSTFIX_SUBMISSION_AUDIT.md",
        "submission_pack/tables/source/checksums.csv",
    ]
    for name, stem in FIGURE_STEMS.items():
        exts = ("png", "tif") if name == "TOC" else ("pdf", "png", "tif")
        for ext in exts:
            paths.append(f"submission_pack/figures/{name}.{ext}")
    for rel in PACK_SOURCE_TABLES:
        paths.append(f"submission_pack/tables/source/{Path(rel).name}")
    return paths


def resolve(rel: str) -> Path:
    return (DT_ROOT / rel).resolve()


def is_legacy_path(path: Path | str) -> bool:
    text = str(path).replace("\\", "/")
    return any(marker.replace("\\", "/") in text for marker in LEGACY_PATH_MARKERS)


def sha256_raw(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_normalized(path: Path) -> str:
    raw = path.read_bytes()
    if path.suffix.lower() in TEXT_SUFFIXES:
        if raw.startswith(b"\xef\xbb\xbf"):
            raw = raw[3:]
        raw = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(raw).hexdigest()


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def require_exists(rel: str) -> Path:
    path = resolve(rel)
    if not path.is_file():
        fail(f"required file missing: {rel}")
    if is_legacy_path(path):
        fail(f"canonical path resolves under legacy archive: {rel} -> {path}")
    return path


def load_csv(rel: str, required_cols: list[str], id_cols: list[str] | None = None) -> list[dict[str, str]]:
    path = require_exists(rel)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            fail(f"{rel}: empty CSV")
        missing = [c for c in required_cols if c not in reader.fieldnames]
        if missing:
            fail(f"{rel}: missing required columns {missing}")
        rows = list(reader)
    if id_cols:
        seen: set[tuple[str, ...]] = set()
        for i, row in enumerate(rows):
            key = tuple(row.get(c, "") for c in id_cols)
            if any(part == "" for part in key):
                fail(f"{rel}: blank identifier at row {i + 2} cols={id_cols}")
            if key in seen:
                fail(f"{rel}: duplicate identifier {key}")
            seen.add(key)
    return rows


def require_finite(rel: str, row: dict[str, str], cols: list[str]) -> None:
    for col in cols:
        raw = row.get(col, "")
        if raw is None or str(raw).strip() == "":
            fail(f"{rel}: blank numeric {col} in {row}")
        try:
            value = float(raw)
        except ValueError:
            fail(f"{rel}: non-numeric {col}={raw!r}")
        if not math.isfinite(value):
            fail(f"{rel}: NaN/inf in {col}={raw!r}")


def one(rows: list[dict[str, str]], rel: str, **keys: str) -> dict[str, str]:
    found = [row for row in rows if all(row.get(k) == v for k, v in keys.items())]
    if len(found) != 1:
        fail(f"{rel}: expected exactly 1 row for {keys}, got {len(found)}")
    return found[0]


def near(actual: float, expected: float, tol: float, label: str) -> None:
    if abs(actual - expected) > tol:
        fail(f"{label}: {actual} not within {tol} of {expected}")


def git_rev(cwd: Path | None = None) -> str:
    import subprocess

    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(cwd or REPO_ROOT),
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return ""


def git_commit_date(commit: str, cwd: Path | None = None) -> str:
    import subprocess

    if not commit:
        return ""
    try:
        return subprocess.check_output(
            ["git", "show", "-s", "--format=%cs", commit],
            cwd=str(cwd or REPO_ROOT),
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return ""


def atomic_replace_dir(staging: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    backup = dest.parent / f".{dest.name}.bak.{os.getpid()}"
    if dest.exists():
        if backup.exists():
            shutil.rmtree(backup)
        dest.rename(backup)
    try:
        staging.rename(dest)
    except Exception:
        if backup.exists() and not dest.exists():
            backup.rename(dest)
        raise
    if backup.exists():
        shutil.rmtree(backup)


def atomic_write_text(path: Path, text: str, encoding: str = "utf-8") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding=encoding, newline="\n") as handle:
            handle.write(text)
        tmp.replace(path)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise


def atomic_write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
        tmp.replace(path)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise


def load_json(rel: str) -> dict:
    path = require_exists(rel)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        fail(f"{rel}: invalid JSON ({exc})")
    if not isinstance(data, dict):
        fail(f"{rel}: expected JSON object")
    return data


def auroc(pos: list[float], neg: list[float]) -> float:
    n1 = len(pos)
    n0 = len(neg)
    if n1 == 0 or n0 == 0:
        fail("AUROC undefined: empty class")
    greater = 0.0
    equal = 0.0
    for p in pos:
        for n in neg:
            if p > n:
                greater += 1
            elif p == n:
                equal += 1
    return (greater + 0.5 * equal) / (n1 * n0)
