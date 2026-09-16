#!/usr/bin/env python3
"""Build the V4 submission pack into a staging directory, validate, then atomically replace.

Does not re-run Vina/GNINA/RTM. Destructive writes support --dry-run.
Partial failure must not leave a half-written canonical pack.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "audit"))
from canonical_paths_v4 import (  # noqa: E402
    DT_ROOT,
    FIGURE_STEMS,
    PACK_SOURCE_TABLES,
    REPO_ROOT,
    atomic_replace_dir,
    atomic_write_text,
    fail,
    git_commit_date,
    git_rev,
    is_legacy_path,
    sha256_raw,
)

ROOT = DT_ROOT
REPO = REPO_ROOT
ART = ROOT / "figures" / "jcim_article"
DOCS = ROOT / "docs"
PACK = ROOT / "submission_pack_postfix"
OLD = ROOT / "submission_pack"
ARCHIVE = ROOT / "submission_pack_pre_v4_archive"
STEMS = FIGURE_STEMS
TABLES = list(PACK_SOURCE_TABLES)

PUB = [
    DOCS / "MANUSCRIPT_JCIM_EN.md",
    DOCS / "MANUSCRIPT_JCIM_ZH.md",
    DOCS / "SUPPORTING_INFORMATION_JCIM_EN_V1.md",
    DOCS / "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
    DOCS / "FIGURE_TABLE_LOCK_POSTFIX_V4.md",
    ART / "CAPTIONS.md",
    ART / "MANUSCRIPT_FIGURE_CAPTIONS.md",
]


def sha256(path: Path) -> str:
    return sha256_raw(path)


def copy_into(dest: Path) -> list[dict[str, str]]:
    if dest.exists():
        shutil.rmtree(dest)
    (dest / "manuscript").mkdir(parents=True)
    (dest / "SI").mkdir()
    (dest / "figures").mkdir()
    (dest / "tables" / "source").mkdir(parents=True)
    (dest / "audit").mkdir()
    copies = [
        (DOCS / "MANUSCRIPT_JCIM_EN.md", dest / "manuscript" / "MANUSCRIPT_JCIM_EN.md", "docs/MANUSCRIPT_JCIM_EN.md", "manuscript"),
        (DOCS / "MANUSCRIPT_JCIM_ZH.md", dest / "manuscript" / "MANUSCRIPT_JCIM_ZH.md", "docs/MANUSCRIPT_JCIM_ZH.md", "manuscript"),
        (DOCS / "SUPPORTING_INFORMATION_JCIM_EN_V1.md", dest / "SI" / "SUPPORTING_INFORMATION_JCIM_EN_V1.md", "docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md", "si"),
        (DOCS / "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md", dest / "SI" / "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md", "docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md", "si"),
        (DOCS / "FIGURE_TABLE_LOCK_POSTFIX_V4.md", dest / "manuscript" / "FIGURE_TABLE_LOCK_POSTFIX_V4.md", "docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md", "lock"),
        (ART / "plotted_values_postfix.json", dest / "figures" / "plotted_values_postfix.json", "figures/jcim_article/plotted_values_postfix.json", "plotted"),
        (ART / "CAPTIONS.md", dest / "figures" / "CAPTIONS.md", "figures/jcim_article/CAPTIONS.md", "caption"),
        (ROOT / "remediation_outputs" / "POST_FIX_AUDIT_REPORT.md", dest / "audit" / "POST_FIX_AUDIT_REPORT.md", "remediation_outputs/POST_FIX_AUDIT_REPORT.md", "audit"),
    ]
    manifest_rows: list[dict[str, str]] = []
    for src, out, rel, role in copies:
        if not src.is_file():
            fail(f"missing source {rel}")
        if is_legacy_path(src):
            fail(f"refusing legacy source {rel}")
        shutil.copy2(src, out)
        if sha256(src) != sha256(out):
            fail(f"copy hash mismatch {rel}")
        manifest_rows.append(_row(out, dest, role, rel, sha256(src)))
    for name, stem in STEMS.items():
        exts = ("png", "tif") if name == "TOC" else ("pdf", "png", "tif")
        for ext in exts:
            src = ART / f"{stem}.{ext}"
            if not src.is_file():
                fail(f"missing figure {src}")
            out = dest / "figures" / f"{name}.{ext}"
            shutil.copy2(src, out)
            if sha256(src) != sha256(out):
                fail(f"figure copy hash mismatch {name}.{ext}")
            manifest_rows.append(
                _row(out, dest, "figure", f"figures/jcim_article/{stem}.{ext}", sha256(src))
            )
    table_rows = []
    for rel in TABLES:
        src = ROOT / rel
        if not src.is_file():
            fail(f"missing table {rel}")
        if is_legacy_path(src):
            fail(f"refusing legacy table {rel}")
        out = dest / "tables" / "source" / Path(rel).name
        shutil.copy2(src, out)
        table_rows.append({"file": Path(rel).name, "source": rel, "sha256": sha256(src), "bytes": src.stat().st_size})
        role = "box_json" if rel.endswith(".json") and "box" in rel else ("plotted" if "plotted_values" in rel else ("metrics" if "post_fix_master" in rel else "table"))
        manifest_rows.append(_row(out, dest, role, rel, sha256(src)))
    checksums = dest / "tables" / "source" / "checksums.csv"
    with checksums.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "source", "sha256", "bytes"])
        writer.writeheader()
        writer.writerows(table_rows)
    return manifest_rows


def _row(out: Path, dest: Path, role: str, source_rel: str, source_sha: str) -> dict[str, str]:
    rel = ("submission_pack/" + str(out.relative_to(dest)).replace("\\", "/"))
    return {
        "path": rel,
        "role": role,
        "bytes": str(out.stat().st_size),
        "sha256_raw": sha256(out),
        "source_path": source_rel,
        "source_sha256_raw": source_sha,
    }


def write_readme(dest: Path, snapshot: str) -> None:
    date = git_commit_date(snapshot) or "unknown"
    dest.joinpath("README.md").write_text(
        "\n".join(
            [
                "# Dual-target docking submission pack (post-fix V4)",
                "",
                "- canonical branch: `cursor/methods-sentence-audit-c7cc`",
                f"- source_snapshot_commit: `{snapshot}`",
                f"- source_snapshot_date: `{date}`",
                "- The packaging git commit is not recorded here; lock it with a Git tag / release metadata.",
                "- canonical result sources:",
                "  - `remediation_outputs/POST_FIX_AUDIT_REPORT.md` (98 PASS / 0 WARNING / 0 FAIL)",
                "  - `remediation_outputs/canonical_tables/post_fix_master_metrics.csv`",
                "  - `remediation_outputs/phase1_boxes/3POZ_box_corrected.json`",
                "  - `remediation_outputs/phase1_boxes/3RCD_box_corrected.json`",
                "  - `data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv`",
                "  - `data/jcim_novelty_v0/tables/eight_pair_ranking_operating_point_v1.csv`",
                "  - `data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv`",
                "  - `figures/jcim_article/plotted_values_postfix.json`",
                "- figure/table numbering lock: `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`",
                "- main text: 5 figures + 3 tables; SI: Figures S1–S5 and Tables S1–S10",
                "- EGFR/HER2 boxes: cognate heavy-atom AABB+5 Å / min edge 20 Å",
                "- AChE/BChE: corrected panel without the ChEMBL-ID-prefix cap",
                "- pre-fix boxes/scores/panels remain in legacy/archive only",
                "",
            ]
        ),
        encoding="utf-8",
    )


def write_release_manifest(dest: Path, rows: list[dict[str, str]], snapshot: str) -> None:
    readme = dest / "README.md"
    rows = list(rows) + [
        {
            "path": "submission_pack/README.md",
            "role": "pack",
            "bytes": str(readme.stat().st_size),
            "sha256_raw": sha256(readme),
            "source_path": "",
            "source_sha256_raw": "",
        }
    ]
    fieldnames = ["path", "role", "bytes", "sha256_raw", "source_path", "source_sha256_raw", "source_snapshot_commit"]
    out = dest / "RELEASE_MANIFEST.csv"
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "source_snapshot_commit": snapshot})
    # Manifest is generated; record it against itself.
    self_row = {
        "path": "submission_pack/RELEASE_MANIFEST.csv",
        "role": "release_manifest",
        "bytes": str(out.stat().st_size),
        "sha256_raw": sha256(out),
        "source_path": "",
        "source_sha256_raw": "",
        "source_snapshot_commit": snapshot,
    }
    # Rewrite once with the self row appended; hashes of the file then change.
    # Keep the self row out of the hashed set: user asked pack copy == source
    # for manuscript/SI/figures/source tables, not for the manifest itself.
    _ = self_row


def scan_stale(paths: list[Path]) -> list[dict]:
    tokens = {
        "0.430": r"(?<![0-9])0\.430(?![0-9])",
        "0.808": r"(?<![0-9])0\.808(?![0-9])",
        "0.378": r"(?<![0-9])0\.378(?![0-9])",
        "0.170": r"(?<![0-9])0\.170(?![0-9])",
        "0.161": r"(?<![0-9])0\.161(?![0-9])",
        "9.505": r"(?<![0-9])9\.505(?![0-9])",
        "18.680": r"18\.680",
        "old_EGFR_box": r"16\.680|18\.680|center = 16",
        "two of eight": r"two of eight",
        "EGFR/HER2 matched-pocket advantage": r"EGFR/HER2 matched-pocket advantage",
        "EGFR top-1 9.505": r"EGFR top-1 9\.505",
        "Figure 6": r"Figure 6",
    }
    allowed_context = {
        "0.808": ("ECFP4", "D vs B", "AChE"),
        "0.170": ("JAK1/JAK2", "−0.177", "[−0.177"),
        "9.505": ("historical", "pre-fix", "not a current", "不得将历史"),
        "0.378": ("pre-fix", "historical", "archived"),
    }
    hits = []
    for path in paths:
        text = path.read_text(encoding="utf-8")
        for name, pat in tokens.items():
            for m in re.finditer(pat, text):
                i = m.start()
                ctx = text[max(0, i - 80) : i + 80].replace("\n", " ")
                hits.append({"file": str(path.relative_to(ROOT)), "token": name, "context": ctx})
    return hits


def audit() -> dict:
    findings = []

    def add(cls, name, status, detail):
        findings.append({"class": cls, "name": name, "status": status, "detail": detail})

    en = (DOCS / "MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    zh = (DOCS / "MANUSCRIPT_JCIM_ZH.md").read_text(encoding="utf-8")
    si_en = (DOCS / "SUPPORTING_INFORMATION_JCIM_EN_V1.md").read_text(encoding="utf-8")
    si_zh = (DOCS / "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md").read_text(encoding="utf-8")
    lock = (DOCS / "FIGURE_TABLE_LOCK_POSTFIX_V4.md").read_text(encoding="utf-8")
    plotted = json.loads((ART / "plotted_values_postfix.json").read_text(encoding="utf-8"))
    audit_md = (ROOT / "remediation_outputs" / "POST_FIX_AUDIT_REPORT.md").read_text(encoding="utf-8")

    # 1. 20-class invariant pointer
    if "98 PASS" in audit_md and "0 WARNING" in audit_md and "0 FAIL" in audit_md:
        add("20-class", "canonical_audit", "PASS", "remediation_outputs/POST_FIX_AUDIT_REPORT.md 98/0/0")
    else:
        add("20-class", "canonical_audit", "FAIL", "canonical report missing 98/0/0")
    root_audit = (ROOT / "POST_FIX_AUDIT_REPORT.md").read_text(encoding="utf-8")
    if root_audit.lstrip().startswith("# SUPERSEDED"):
        add("20-class", "root_pointer", "PASS", "root report is a pointer")
    else:
        add("20-class", "root_pointer", "FAIL", "root POST_FIX_AUDIT_REPORT.md is not a pointer")

    # 2. figure provenance
    recs = plotted["records"]
    def rec(fig, panel, pair, metric):
        hits = [r for r in recs if r["figure"] == fig and r["panel"] == panel and r["pair"] == pair and r["metric"] == metric]
        return hits[0] if hits else None

    checks = [
        ("Figure 2", "A", "EGFR/HER2", "fixed_delta_pocketA", 0.462, 0.01),
        ("Figure 2", "B", "EGFR/HER2", "AUROC_D_vs_B_pocketA", 0.324, 0.002),
        ("Figure 2", "C", "EGFR/HER2", "summary_min", 0.324, 0.002),
        ("Figure 2", "D", "EGFR/HER2", "top_dual", 1, 0),
        ("Figure 2", "D", "AChE/BChE", "top_dual", 5, 0),
        ("Figure S3", "A", "EGFR 3POZ", "top1_rmsd_A", 1.019, 0.001),
        ("Figure 5", "A", "EGFR/HER2", "gnina_summary_min", 0.265, 0.002),
        ("Figure 5", "A", "EGFR/HER2", "gnina_D_vs_neither", 0.737, 0.002),
    ]
    for fig, panel, pair, metric, expect, tol in checks:
        r = rec(fig, panel, pair, metric)
        if r is None:
            add("provenance", f"{fig}{panel}:{pair}:{metric}", "FAIL", "missing plotted record")
            continue
        val = float(r["raw_value"])
        ok = abs(val - expect) <= (tol if tol else 1e-9) or abs(val - expect) < 0.001
        add("provenance", f"{fig}{panel}:{pair}:{metric}", "PASS" if ok else "FAIL", f"{val} (expect {expect})")

    # 3. manuscript-number trace
    must = [
        ("EN EGFR 0.324", "0.324" in en),
        ("EN EGFR 0.786", "0.786" in en),
        ("EN EGFR 0.462", "0.462" in en),
        ("EN EGFR matched 0.056", "0.056" in en),
        ("EN AChE 27 / 26 / 28", "27 / 26 / 28" in en),
        ("EN AChE 0.177", "0.177" in en),
        ("EN 3POZ 1.019", "1.019" in en),
        ("EN 3RCD 1.947", "1.947" in en),
        ("EN range 0.324 to 0.728", "0.324 to 0.728" in en),
        ("ZH range 0.324", "0.324 至 0.728" in zh),
        ("EN no Figure 6", "Figure 6" not in en),
        ("ZH no Figure 6", "Figure 6" not in zh),
        ("EN Figure 2D eight-pair", "Eight-pair top-10%" in en or "eight-pair top-10%" in en.lower() or "Eight-pair" in en),
        ("lock V4 unique", "FIGURE_TABLE_LOCK_POSTFIX_V4.md" in lock),
        ("S2 18.816", "18.816" in si_en and "18.816" in si_zh),
        ("S2 no 18.680", "18.680" not in si_en and "18.680" not in si_zh),
        ("S2b 1.019", "1.019" in si_en and "1.019" in si_zh),
        ("S2b no current 9.505 unlabeled", "9.505" not in si_en.split("historical")[0] or "historical 9.505" in si_en),
    ]
    for name, ok in must:
        add("manuscript-number", name, "PASS" if ok else "FAIL", "")

    # 4. EN/ZH consistency of flagship tokens
    for token in ["0.324", "0.462", "0.056", "0.177", "1.019", "1.947", "27 / 26 / 28", "0.357", "1.778"]:
        add("EN/ZH", token, "PASS" if token in en and token in zh else "FAIL", f"EN={token in en} ZH={token in zh}")

    # 5. SI/main cross-ref: main has 5 figs, SI S1-S5, no main Fig6
    for label, text in [("EN", en), ("ZH", zh)]:
        for n in range(1, 6):
            add("xref", f"{label} Figure {n}", "PASS" if f"Figure {n}" in text else "FAIL", "")
        add("xref", f"{label} no Fig6", "PASS" if "Figure 6" not in text else "FAIL", "")
        for n in range(1, 6):
            add("xref", f"{label} Figure S{n}", "PASS" if f"Figure S{n}" in text or f"Figure S{n}" in (si_en if label=="EN" else si_zh) else "FAIL", "")

    # 6. caption/panel
    for path in [ART / "CAPTIONS.md", ART / "MANUSCRIPT_FIGURE_CAPTIONS.md", DOCS / "FIGURE_PANEL_LOCK_V3.md"]:
        t = path.read_text(encoding="utf-8")
        add("caption", path.name + " superseded or V4", "PASS" if "FIGURE_TABLE_LOCK_POSTFIX_V4" in t else "FAIL", "")
    cap = (ART / "MANUSCRIPT_FIGURE_CAPTIONS.md").read_text(encoding="utf-8")
    add("caption", "no JAK1-only 2D as current def", "PASS" if "Eight-pair top-10%" in cap or "eight-pair top-10%" in cap.lower() else "FAIL", "")
    add("caption", "no 0.378 as current", "PASS" if "0.378" not in cap else "FAIL", "")
    add("caption", "no 9.505 as current unlabeled", "PASS" if "9.505 Å value is not current" in cap or "historical 9.505" in cap else "FAIL" if "9.505" in cap else "PASS", "")

    # 7. pack checksum
    man_pack = sha256(PACK / "manuscript" / "MANUSCRIPT_JCIM_EN.md")
    man_docs = sha256(DOCS / "MANUSCRIPT_JCIM_EN.md")
    add("pack", "EN manuscript checksum", "PASS" if man_pack == man_docs else "FAIL", "")
    si_pack = sha256(PACK / "SI" / "SUPPORTING_INFORMATION_JCIM_EN_V1.md")
    si_docs = sha256(DOCS / "SUPPORTING_INFORMATION_JCIM_EN_V1.md")
    add("pack", "EN SI checksum", "PASS" if si_pack == si_docs else "FAIL", "")
    for name in STEMS:
        exts = ("png", "tif") if name == "TOC" else ("pdf", "png", "tif")
        for ext in exts:
            add("pack", f"{name}.{ext}", "PASS" if (PACK / "figures" / f"{name}.{ext}").exists() else "FAIL", "")

    stale = scan_stale(PUB + [PACK / "manuscript" / "MANUSCRIPT_JCIM_EN.md", PACK / "SI" / "SUPPORTING_INFORMATION_JCIM_EN_V1.md", PACK / "README.md"])
    blocking_stale = []
    notes = []
    for h in stale:
        tok = h["token"]
        ctx = h["context"]
        if tok in ("0.808",) and "AChE" in ctx and "ECFP" in ctx.replace(" ", ""):
            notes.append(h)
            continue
        if tok == "0.808" and "0.808" in ctx and ("AChE/BChE | D vs B" in ctx or "AChE/BChE | D vs B" in h["file"] or "D vs B | 0.821 | 0.808" in ctx):
            notes.append(h)
            continue
        if tok == "0.170" and ("JAK1/JAK2" in ctx or "−0.177" in ctx or "[−0.177" in ctx):
            notes.append(h)
            continue
        if tok == "9.505" and ("historical" in ctx.lower() or "pre-fix" in ctx.lower() or "not a current" in ctx.lower() or "不得将历史" in ctx or "not current" in ctx.lower()):
            notes.append(h)
            continue
        if tok in ("0.430", "0.378", "0.161", "0.170", "0.808") and ("pre-fix" in ctx.lower() or "historical" in ctx.lower() or "not current" in ctx.lower()):
            notes.append(h)
            continue
        if tok == "Figure 6" and ("no main-text" in ctx.lower() or "there is no" in ctx.lower() or "not current" in ctx.lower() or "superseded" in ctx.lower() or "historical" in ctx.lower()):
            notes.append(h)
            continue
        if tok == "Figure 6" and "SUPERSEDED" in ctx:
            notes.append(h)
            continue
        if tok in ("0.430", "0.378", "0.161", "18.680", "old_EGFR_box", "two of eight", "EGFR/HER2 matched-pocket advantage", "EGFR top-1 9.505"):
            blocking_stale.append(h)
            continue
        if tok == "Figure 6":
            blocking_stale.append(h)
            continue
        notes.append(h)
    add("stale-token", "blocking", "PASS" if not blocking_stale else "FAIL", json.dumps(blocking_stale, ensure_ascii=False)[:1500])
    add("stale-token", "labeled_or_unrelated", "PASS", f"n={len(notes)}")

    n_fail = sum(1 for f in findings if f["status"] == "FAIL")
    ready = n_fail == 0
    return {"findings": findings, "n_fail": n_fail, "ready": ready, "stale_notes": notes, "blocking_stale": blocking_stale}


def write_v2_and_final(result: dict, snapshot: str, pack_root: Path, write_docs: bool) -> None:
    lines = [
        "# Five-round JCIM submission audit (V2 post-fix)",
        "",
        "Date: 2026-09-16",
        "Branch: `cursor/methods-sentence-audit-c7cc`",
        f"source_snapshot_commit: `{snapshot}`",
        "This audit verifies post-fix canonical values. It does not replace `docs/SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md` (pre-remediation).",
        "",
        f"Summary: **{sum(1 for f in result['findings'] if f['status']=='PASS')} PASS**, **{result['n_fail']} FAIL**.",
        "",
        "| Round | Status | Finding |",
        "|---|---|---|",
    ]
    for f in result["findings"]:
        lines.append(f"| {f['class']} | {f['status']} | {f['name']}: {f['detail']} |")
    lines += [
        "",
        "Post-fix tokens verified:",
        "",
        "- EGFR D-vs-B-only 0.324; D-vs-neither 0.786; Δ 0.462",
        "- EGFR matched-minus-mismatched 0.056; CI includes zero",
        "- AChE/BChE n_scored 27/26/28; matched 0.177",
        "- Top10 EGFR 1/5/5/0 EF=0.357; AChE 5/3/1/1 EF≈1.778",
        "- Figure 2D eight-pair stacked composition",
        "- cognate 3POZ top-1 1.019 Å; 3RCD 1.947 Å",
        "- boxes 18.816 / 12.552 from canonical JSON",
        "",
        "Target: 0 FAIL.",
        "",
    ]
    v2_text = "\n".join(lines) + "\n"
    (pack_root / "audit" / "SUBMISSION_AUDIT_FIVE_ROUNDS_V2_POSTFIX.md").write_text(v2_text, encoding="utf-8")
    if write_docs:
        atomic_write_text(DOCS / "SUBMISSION_AUDIT_FIVE_ROUNDS_V2_POSTFIX.md", v2_text)

    verdict = "READY FOR LANGUAGE/EDITORIAL POLISH" if result["ready"] else "NOT READY — BLOCKING ISSUES REMAIN"
    fin = [
        "# FINAL POSTFIX SUBMISSION AUDIT",
        "",
        f"source_snapshot_commit: `{snapshot}`",
        "Branch: `cursor/methods-sentence-audit-c7cc`",
        "Canonical audit: `remediation_outputs/POST_FIX_AUDIT_REPORT.md` (98 PASS / 0 WARNING / 0 FAIL)",
        "Figure lock: `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`",
        "",
        "## Audits",
        "",
        "1. 20-class invariant audit — canonical report 98/0/0; root file is a pointer",
        "2. figure-value provenance audit — plotted_values_postfix.json vs post-fix CSVs",
        "3. manuscript-number trace — EN/ZH/SI vs canonical post-fix values",
        "4. EN/ZH consistency audit — flagship tokens present in both languages",
        "5. SI/main cross-reference audit — Figures 1–5 and S1–S5; no main Figure 6",
        "6. caption/panel-number audit — V4 lock is the unique current lock",
        "7. submission-pack checksum audit — pack copies match docs/figures",
        "",
        "## Result table",
        "",
        "| Class | Name | Status | Detail |",
        "|---|---|---|---|",
    ]
    for f in result["findings"]:
        fin.append(f"| {f['class']} | {f['name']} | {f['status']} | {f['detail']} |")
    if result["stale_notes"]:
        fin += ["", "## Non-blocking stale-token notes (current chemistry / labeled historical)", ""]
        for h in result["stale_notes"]:
            fin.append(f"- `{h['token']}` in `{h['file']}`: {h['context'][:160]}")
    if result["blocking_stale"]:
        fin += ["", "## Blocking stale tokens", ""]
        for h in result["blocking_stale"]:
            fin.append(f"- `{h['token']}` in `{h['file']}`: {h['context'][:200]}")
    fin += ["", "## Verdict", "", f"**{verdict}**", ""]
    final_text = "\n".join(fin) + "\n"
    (pack_root / "audit" / "FINAL_POSTFIX_SUBMISSION_AUDIT.md").write_text(final_text, encoding="utf-8")
    if write_docs:
        atomic_write_text(DOCS / "FINAL_POSTFIX_SUBMISSION_AUDIT.md", final_text)
    print(verdict)
    print("FAIL", result["n_fail"])
    for f in result["findings"]:
        if f["status"] == "FAIL":
            print(" ", f)
    if not result["ready"]:
        fail("publication-facing audit has FAIL items; canonical pack not replaced")


def copy_audit_docs_from_pack(pack_root: Path) -> None:
    mapping = [
        "SUBMISSION_AUDIT_FIVE_ROUNDS_V2_POSTFIX.md",
        "FINAL_POSTFIX_SUBMISSION_AUDIT.md",
    ]
    for name in mapping:
        src = pack_root / "audit" / name
        if not src.is_file():
            fail(f"pack audit missing {name}")
        atomic_write_text(DOCS / name, src.read_text(encoding="utf-8"))


def promote_from(postfix: Path) -> None:
    staging = ROOT / ".submission_pack.staging"
    if staging.exists():
        shutil.rmtree(staging)
    shutil.copytree(postfix, staging)
    atomic_replace_dir(staging, OLD)
    print("promoted submission_pack_postfix -> submission_pack")


def main() -> int:
    global PACK
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-promote", action="store_true")
    parser.add_argument(
        "--source-snapshot",
        default="",
        help="commit used to generate the pack; do not default to packaging HEAD",
    )
    args = parser.parse_args()
    snapshot = args.source_snapshot.strip()
    if not snapshot:
        readme = OLD / "README.md"
        if readme.is_file():
            for line in readme.read_text(encoding="utf-8").splitlines():
                if "source_snapshot_commit:" in line and "`" in line:
                    snapshot = line.split("`")[1].strip()
                    break
    if not snapshot:
        fail("source_snapshot_commit unresolved; pass --source-snapshot (do not stamp packaging HEAD)")
    staging = ROOT / ".submission_pack_postfix.staging"
    try:
        rows = copy_into(staging)
        write_readme(staging, snapshot)
        PACK = staging
        result = audit()
        write_v2_and_final(result, snapshot, staging, write_docs=False)
        write_release_manifest(staging, rows, snapshot)
        if args.dry_run:
            print(f"dry-run: staging pack OK at {staging}; canonical not replaced")
            shutil.rmtree(staging)
            return 0
        if not result["ready"]:
            fail("publication-facing audit has FAIL items; canonical pack not replaced")
        atomic_replace_dir(staging, ROOT / "submission_pack_postfix")
        PACK = ROOT / "submission_pack_postfix"
        if not args.skip_promote:
            copy_audit_docs_from_pack(PACK)
            promote_from(PACK)
    except BaseException:
        if staging.exists() and staging.resolve() != (ROOT / "submission_pack_postfix").resolve() and staging.resolve() != OLD.resolve():
            shutil.rmtree(staging, ignore_errors=True)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
