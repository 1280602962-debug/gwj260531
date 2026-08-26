#!/usr/bin/env python3
"""Assemble a clean English JCIM manuscript from canonical section files."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "MANUSCRIPT_JCIM_EN.md"

SECTIONS = [
    ("INTRODUCTION_SECTION_JCIM_EN_V1.md", "## 1."),
    ("METHODS_SECTION_JCIM_EN_V1.md", "## 2."),
    ("RESULTS_SECTION_JCIM_EN_V1.md", "## 3."),
    ("DISCUSSION_SECTION_JCIM_EN_V1.md", "## 4."),
    ("CONCLUSIONS_SECTION_JCIM_EN_V1.md", "## 5."),
    ("DATA_AND_SOFTWARE_AVAILABILITY_JCIM_EN.md", "## Data and Software Availability"),
]


def read(name: str) -> str:
    return (ROOT / name).read_text(encoding="utf-8").lstrip("\ufeff")


def title_abstract_keywords() -> tuple[str, str, str]:
    text = read("TITLE_AND_ABSTRACT_JCIM_EN_V1.md")
    title = text.split("## Title", 1)[1].split("## Abstract", 1)[0].strip().strip("*")
    abstract_block = text.split("## Abstract", 1)[1]
    abstract = abstract_block.split("## Keywords", 1)[0].strip()
    keywords = abstract_block.split("## Keywords", 1)[1].strip()
    return title, abstract, keywords


def section_body(name: str, marker: str) -> str:
    text = read(name)
    if marker not in text:
        raise ValueError(f"{name}: missing section marker {marker!r}")
    body = marker + text.split(marker, 1)[1]
    if name.startswith("INTRODUCTION_"):
        body = body.split("\n---\n\n## References", 1)[0]
    body = re.sub(r"\^\(([\d,–-]+)\)", r"[\1]", body)
    return body.strip()


def main() -> None:
    title, abstract, keywords = title_abstract_keywords()
    parts = [f"# {title}\n\n## Abstract\n\n{abstract}\n\n**Keywords:** {keywords}\n\n"]
    parts.extend(section_body(name, marker) + "\n\n" for name, marker in SECTIONS)
    parts.append(read("REFERENCES_JCIM.md").strip() + "\n")
    manuscript = "".join(parts)
    forbidden = ("Companion to", "Do not use:", "Working draft", "not typeset", "写法说明")
    leaked = [term for term in forbidden if term in manuscript]
    if leaked:
        raise ValueError(f"internal author notes leaked into manuscript: {leaked}")
    OUT.write_text(manuscript, encoding="utf-8")
    print("wrote", OUT, "bytes", OUT.stat().st_size)


if __name__ == "__main__":
    main()
