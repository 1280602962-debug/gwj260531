#!/usr/bin/env python3
"""Assemble a clean Chinese JCIM working manuscript from canonical sections."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "MANUSCRIPT_JCIM_ZH.md"

SECTIONS = [
    ("INTRODUCTION_DRAFT_ZH_JCIM_V1.md", "## 1."),
    ("METHODS_DRAFT_ZH_JCIM_V1.md", "## 2."),
    ("RESULTS_DRAFT_ZH_JCIM_V1.md", "## 3."),
    ("DISCUSSION_DRAFT_ZH_JCIM_V1.md", "## 4."),
    ("CONCLUSIONS_DRAFT_ZH_JCIM_V1.md", "## 5."),
    ("DATA_AND_SOFTWARE_AVAILABILITY_JCIM_ZH.md", "## 数据与软件可用性"),
]


def read(name: str) -> str:
    return (ROOT / name).read_text().lstrip("\ufeff")


def title_abstract_keywords() -> tuple[str, str, str]:
    text = read("TITLE_AND_ABSTRACT_JCIM_ZH_V1.md")
    title = text.split("## 题名", 1)[1].split("## 摘要", 1)[0].strip().strip("*")
    abstract_block = text.split("## 摘要", 1)[1]
    abstract = abstract_block.split("## 关键词", 1)[0].strip()
    keywords = abstract_block.split("## 关键词", 1)[1].strip()
    return title, abstract, keywords


def section_body(name: str, marker: str) -> str:
    text = read(name)
    if marker not in text:
        raise ValueError(f"{name}: missing section marker {marker!r}")
    body = marker + text.split(marker, 1)[1]
    if name.startswith("INTRODUCTION_"):
        body = body.split("\n---\n\n[^1]:", 1)[0]
    body = re.sub(r"\[\^(\d+)\]", r"[\1]", body)
    body = re.sub(r"\^\(([\d,–-]+)\)", r"[\1]", body)
    return body.strip()


def main() -> None:
    title, abstract, keywords = title_abstract_keywords()
    parts = [f"# {title}\n\n## 摘要\n\n{abstract}\n\n**关键词：** {keywords}\n\n"]
    parts.extend(section_body(name, marker) + "\n\n" for name, marker in SECTIONS)
    parts.append(read("REFERENCES_JCIM.md").replace("## References", "## 参考文献", 1).strip() + "\n")
    manuscript = "".join(parts)
    forbidden = ("供阅读与内部核对", "写法说明（不进正文）", "投稿以英文为准", "不要使用：")
    leaked = [term for term in forbidden if term in manuscript]
    if leaked:
        raise ValueError(f"internal author notes leaked into manuscript: {leaked}")
    OUT.write_text(manuscript)
    print("wrote", OUT, "bytes", OUT.stat().st_size)


if __name__ == "__main__":
    main()
