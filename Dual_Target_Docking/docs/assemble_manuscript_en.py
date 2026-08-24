#!/usr/bin/env python3
"""Assemble the English JCIM main text from canonical section files.

The English submission text remains `MANUSCRIPT_JCIM_EN.md`.
Chinese working assembly: `python3 docs/assemble_manuscript_zh.py` → `MANUSCRIPT_JCIM_ZH.md`.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "MANUSCRIPT_JCIM_EN.md"

SECTIONS = [
    ("Title and Abstract", "TITLE_AND_ABSTRACT_JCIM_EN_V1.md"),
    ("Introduction", "INTRODUCTION_SECTION_JCIM_EN_V1.md"),
    ("Methods", "METHODS_SECTION_JCIM_EN_V1.md"),
    ("Results", "RESULTS_SECTION_JCIM_EN_V1.md"),
    ("Discussion", "DISCUSSION_SECTION_JCIM_EN_V1.md"),
    ("Conclusions", "CONCLUSIONS_SECTION_JCIM_EN_V1.md"),
]


def strip_bom(text: str) -> str:
    return text.lstrip("\ufeff")


def main() -> None:
    parts = [
        "# DualFourClass-Bench — assembled English manuscript\n",
        "\n",
        "> Working assembly for a JCIM Articles submission. **Canonical sources are the section files listed below.** Edit those files and re-run `python3 docs/assemble_manuscript_en.py`. Do not add `_V2` / `_FINAL` copies of individual sections.\n",
        ">\n",
        "> Claim ceiling: [`../data/jcim_bench_v0/CLAIM_CEILING.md`](../data/jcim_bench_v0/CLAIM_CEILING.md). Figure files: [`../figures/jcim_article/`](../figures/jcim_article/).\n",
        "\n",
        "---\n",
        "\n",
    ]
    for title, name in SECTIONS:
        path = ROOT / name
        body = strip_bom(path.read_text())
        parts.append(f"<!-- BEGIN {name} -->\n")
        parts.append(body.rstrip())
        parts.append(f"\n<!-- END {name} -->\n\n---\n\n")
    OUT.write_text("".join(parts))
    print("wrote", OUT, "bytes", OUT.stat().st_size)


if __name__ == "__main__":
    main()
