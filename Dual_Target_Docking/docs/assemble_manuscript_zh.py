#!/usr/bin/env python3
"""Assemble the Chinese JCIM working manuscript from canonical ZH section files.

Edit the section files, then re-run this script. Do not create _V2/_FINAL copies.
The English submission text remains `MANUSCRIPT_JCIM_EN.md`.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "MANUSCRIPT_JCIM_ZH.md"

SECTIONS = [
    ("题名与摘要", "TITLE_AND_ABSTRACT_JCIM_ZH_V1.md"),
    ("引言", "INTRODUCTION_DRAFT_ZH_JCIM_V1.md"),
    ("方法", "METHODS_DRAFT_ZH_JCIM_V1.md"),
    ("结果", "RESULTS_DRAFT_ZH_JCIM_V1.md"),
    ("讨论", "DISCUSSION_DRAFT_ZH_JCIM_V1.md"),
    ("结论", "CONCLUSIONS_DRAFT_ZH_JCIM_V1.md"),
]

CUT_MARKERS = {
    "METHODS_DRAFT_ZH_JCIM_V1.md": "## 写法说明（不进正文）",
}


def strip_bom(text: str) -> str:
    return text.lstrip("\ufeff")


def apply_cut(name: str, body: str) -> str:
    marker = CUT_MARKERS.get(name)
    if not marker:
        return body
    idx = body.find(marker)
    if idx < 0:
        raise SystemExit(f"{name}: expected cut marker {marker!r} not found")
    return body[:idx].rstrip()


def main() -> None:
    parts = [
        "# DualFourClass-Bench — 组装中文工作稿\n",
        "\n",
        "> 供阅读与内部核对的中文主稿。**投稿以英文 [`MANUSCRIPT_JCIM_EN.md`](MANUSCRIPT_JCIM_EN.md) 为准。** 规范来源是下列中文章节文件；修改那些文件后运行 `python3 docs/assemble_manuscript_zh.py`。不要另开 `_V2` / `_FINAL` 分文件。\n",
        ">\n",
        "> 主张边界：[`../data/jcim_bench_v0/CLAIM_CEILING.md`](../data/jcim_bench_v0/CLAIM_CEILING.md)。主图：[`../figures/jcim_article/`](../figures/jcim_article/)。Methods 文末“写法说明”已在组装时去掉。\n",
        "\n",
        "---\n",
        "\n",
    ]
    for _title, name in SECTIONS:
        path = ROOT / name
        body = apply_cut(name, strip_bom(path.read_text()))
        parts.append(f"<!-- BEGIN {name} -->\n")
        parts.append(body.rstrip())
        parts.append(f"\n<!-- END {name} -->\n\n---\n\n")
    OUT.write_text("".join(parts))
    print("wrote", OUT, "bytes", OUT.stat().st_size)


if __name__ == "__main__":
    main()
