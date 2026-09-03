#!/usr/bin/env python3
"""Copy KEEP_FOR_RELEASE=1 assets into release/JCIM_submission_v1/.

Does not delete history in the working tree. The release directory is gitignored.
"""
from __future__ import annotations

import csv
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "data" / "manuscript_lock" / "ARTICLE_ASSET_INDEX_v1.csv"
OUT = ROOT / "release" / "JCIM_submission_v1"

PATH_RE = re.compile(r"((?:docs|data|figures|scripts)/[^\s]+)")
BRACE = re.compile(r"\{([^{}]+)\}")


def brace_expand(token: str) -> list[str]:
    m = BRACE.search(token)
    if not m:
        return [token]
    return [BRACE.sub(ext.strip(), token, count=1) for ext in m.group(1).split(",")]


def paths_from_cell(cell: str) -> list[str]:
    out: list[str] = []
    for match in PATH_RE.finditer(cell.replace("`", "")):
        token = match.group(1).rstrip(".,;:)`")
        out.extend(brace_expand(token))
    return out


def copy_file(src: Path, seen: set[str]) -> bool:
    if not src.is_file():
        return False
    key = str(src.relative_to(ROOT))
    if key in seen:
        return False
    dest = OUT / key
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    seen.add(key)
    return True


def main() -> int:
    if not INDEX.exists():
        print(f"missing {INDEX}", file=sys.stderr)
        return 1
    rows = list(csv.DictReader(INDEX.open(encoding="utf-8")))
    if not rows or "keep_for_release" not in rows[0]:
        print("ARTICLE_ASSET_INDEX_v1.csv needs keep_for_release column", file=sys.stderr)
        return 1
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    copied = 0
    missing: list[str] = []
    seen: set[str] = set()
    for r in rows:
        if str(r.get("keep_for_release", "")).strip() not in {"1", "true", "TRUE"}:
            continue
        for rel in paths_from_cell(r.get("canonical_source", "")):
            if any(ch in rel for ch in "*?["):
                hits = list(ROOT.glob(rel))
                if not hits:
                    missing.append(rel)
                    continue
                for src in hits:
                    if copy_file(src, seen):
                        copied += 1
                continue
            src = ROOT / rel
            if copy_file(src, seen):
                copied += 1
            else:
                missing.append(rel)
    extras = [
        ROOT / "docs" / "STATISTICAL_LOCK_V1.md",
        ROOT / "docs" / "MANUSCRIPT_JCIM_EN.md",
        ROOT / "docs" / "SUPPORTING_INFORMATION_JCIM_EN_V1.md",
        ROOT / "docs" / "REFERENCES_JCIM.md",
        ROOT / "README.md",
    ]
    extras.extend(sorted((ROOT / "figures" / "jcim_article").glob("*")))
    for src in extras:
        if copy_file(src, seen):
            copied += 1
    manifest = OUT / "RELEASE_MANIFEST.txt"
    unique_missing = sorted(set(missing))
    manifest.write_text(
        f"copied={copied}\nmissing={len(unique_missing)}\n"
        + "".join(f"MISSING {m}\n" for m in unique_missing),
        encoding="utf-8",
    )
    print(f"release bundle: copied {copied} files -> {OUT}")
    if unique_missing:
        print(f"warning: {len(unique_missing)} listed paths not found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
