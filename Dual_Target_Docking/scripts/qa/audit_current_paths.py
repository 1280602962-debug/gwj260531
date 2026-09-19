#!/usr/bin/env python3
"""Scan current-facing docs for source/protocol/script paths. Existence only."""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT_MD = ROOT / "docs/CURRENT_PATH_AUDIT.md"
OUT_CSV = ROOT / "results/qa/current_path_audit.csv"

SCAN_DIRS = [
    ROOT / "docs",
    ROOT / "submission_pack",
    ROOT / "figures" / "jcim_article",
]
SCAN_FILES = [
    ROOT / "README.md",
    ROOT / "data/jcim_chembl_universe_v0/README.md",
    ROOT / "data/jcim_chembl_universe_v0/tables/track_b_local_run_v1.yaml",
    ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/README.md",
    ROOT / "requirements-analysis.txt",
]
SKIP_NAMES = {
    "FINAL_REVIEWER_RED_TEAM_AUDIT.md",
    "BACKWARD_PROVENANCE_AUDIT.md",
    "STALE_CURRENT_TRUTH_AUDIT.md",
    "PROTOCOL_CONTRADICTION_AUDIT.md",
    "EGFR_LIGAND_PREP_PROVENANCE_AUDIT.md",
    "EGFR_UNIFORM_PREP_IMPACT.md",
    "CURRENT_PATH_AUDIT.md",
    "FINAL_WRITING_FREEZE_AUDIT.md",
}
HISTORICAL_MARKERS = (
    "git history",
    "git-history",
    "historical path",
    "historical / unavailable",
    "historical/unavailable",
    "archived",
    "unavailable",
    "not a current instruction",
    "not current instruction",
    "superseded",
    "supersedes",
    "deleted",
    "must not exist",
    "must be absent",
    "not written",
    "leftover",
    "not a repository",
    "not current sources",
    "not a current source",
    "duplicate authoritative",
    "cleanup notes",
)
KNOWN_FALSE_POS = {
    "results/canonical",
    "docs/",
    "submission_pack/",
    "Table 1",
    "Table 2",
    "Table S1",
    "pair/channel",
    "B = 2000",
    "theta=6.0",
    "θ=6.0",
}
PATH_RE = re.compile(
    r"`((?:data|docs|scripts|results|figures|submission_pack|requirements)[^`\n]{3,220})`"
)
MD_LINK_RE = re.compile(r"\[[^\]]+\]\(([^)]+)\)")


def is_pathish(text: str) -> bool:
    t = text.strip().rstrip(").,;")
    if t in KNOWN_FALSE_POS:
        return False
    if t.startswith("http") or t.startswith("mailto:"):
        return False
    if any(t.endswith(ext) for ext in (".csv", ".md", ".py", ".yaml", ".yml", ".json", ".txt", ".pdbqt", ".pdb", ".png")):
        return True
    if t.startswith(("data/", "docs/", "scripts/", "results/", "figures/", "submission_pack/")):
        return True
    return False


def exists(rel: str) -> bool:
    rel = rel.split()[0].rstrip(").,;")
    rel = rel.replace("\\_", "_")
    if any(ch in rel for ch in "*{}"):
        # glob pattern: pass if any match
        return bool(list(ROOT.glob(rel)))
    return (ROOT / rel).exists()


def nearby_historical(text: str, start: int) -> bool:
    window = text[max(0, start - 180) : start + 180].lower()
    return any(m in window for m in HISTORICAL_MARKERS)


def scan_file(path: Path) -> list[dict]:
    text = path.read_text(encoding="utf-8", errors="replace")
    rel_src = str(path.relative_to(ROOT))
    rows = []
    seen = set()
    for match in PATH_RE.finditer(text):
        cand = match.group(1).strip()
        if not is_pathish(cand):
            continue
        key = (rel_src, cand)
        if key in seen:
            continue
        seen.add(key)
        hist = nearby_historical(text, match.start())
        ok = exists(cand.split("`")[0])
        status = "ok" if ok else ("historical_or_unavailable" if hist else "broken_current")
        rows.append(
            {
                "source_file": rel_src,
                "cited_path": cand,
                "exists": int(ok),
                "historical_marker_nearby": int(hist),
                "status": status,
            }
        )
    return rows


def main() -> int:
    files = []
    for d in SCAN_DIRS:
        if not d.is_dir():
            continue
        files.extend(p for p in d.rglob("*") if p.is_file() and p.suffix.lower() in {".md", ".yaml", ".yml", ".txt", ".csv"})
    files.extend(p for p in SCAN_FILES if p.is_file())
    rows = []
    for path in sorted(set(files)):
        if path.name in SKIP_NAMES:
            continue
        if "archive" in path.parts or "legacy" in path.parts:
            continue
        rows.extend(scan_file(path))
    broken = [r for r in rows if r["status"] == "broken_current"]
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else ["source_file"], lineterminator="\n")
        if rows:
            w.writeheader()
            w.writerows(rows)
    lines = [
        "# Current-path audit",
        "",
        "Date: 2026-09-18",
        "Scanner: `scripts/qa/audit_current_paths.py`",
        "Rule: every current-facing source/protocol/script path must exist, or be marked historical / archived / unavailable.",
        "Checksum manifests are integrity records, not scientific PASS/FAIL authority.",
        "",
        f"- files scanned: {len(set(r['source_file'] for r in rows))}",
        f"- path citations: {len(rows)}",
        f"- exists: {sum(r['exists'] for r in rows)}",
        f"- historical_or_unavailable: {sum(1 for r in rows if r['status']=='historical_or_unavailable')}",
        f"- unresolved broken current paths: **{len(broken)}**",
        "",
    ]
    if broken:
        lines.append("## Broken current paths")
        lines.append("")
        for r in broken:
            lines.append(f"- `{r['source_file']}` → `{r['cited_path']}`")
        lines.append("")
        verdict = "NOT_READY_FOR_WRITING_FREEZE (broken current paths > 0)"
    else:
        lines.append("No unresolved broken current-facing paths.")
        lines.append("")
        verdict = "path audit unresolved = 0"
    lines.append(f"## Verdict")
    lines.append("")
    lines.append(verdict)
    lines.append("")
    lines.append("Machine table: `results/qa/current_path_audit.csv`.")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote", OUT_MD, "broken", len(broken))
    return 0 if not broken else 1


if __name__ == "__main__":
    raise SystemExit(main())
