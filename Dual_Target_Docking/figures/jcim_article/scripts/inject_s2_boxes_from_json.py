#!/usr/bin/env python3
"""Rewrite SI Table S2a EGFR/HER2 box rows from canonical JSON.

Do not hard-code coordinates. Read:
  remediation_outputs/phase1_boxes/3POZ_box_corrected.json
  remediation_outputs/phase1_boxes/3RCD_box_corrected.json
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BOX_DIR = ROOT / "remediation_outputs" / "phase1_boxes"
SI_FILES = [
    ROOT / "docs" / "SUPPORTING_INFORMATION_JCIM_EN_V1.md",
    ROOT / "docs" / "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
]


def load_box(pdb: str) -> tuple[str, str]:
    d = json.loads((BOX_DIR / f"{pdb}_box_corrected.json").read_text(encoding="utf-8"))
    center = f"{d['center_x']:.3f}, {d['center_y']:.3f}, {d['center_z']:.3f}"
    size = f"{float(d['size_x']):.3f}, {float(d['size_y']):.3f}, {float(d['size_z']):.3f}"
    return center, size


def row(protein: str, pdb: str, cognate: str, resolution: str) -> str:
    center, size = load_box(pdb)
    return f"| {protein} | {pdb} | {cognate} | {resolution} | {center} | {size} |"


def patch(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    replacements = {
        r"\| EGFR \| 3POZ \| 03P \| 1\.50 \| [^|]+ \| [^|]+ \|": row("EGFR", "3POZ", "03P", "1.50"),
        r"\| HER2 \| 3RCD \| 03P \| 3\.21 \| [^|]+ \| [^|]+ \|": row("HER2", "3RCD", "03P", "3.21"),
    }
    for pat, new in replacements.items():
        text2, n = re.subn(pat, new, text, count=1)
        if n != 1:
            raise SystemExit(f"{path.name}: failed to replace {pat!r} (n={n})")
        text = text2
    path.write_text(text, encoding="utf-8")
    print("patched", path.relative_to(ROOT))
    print("  3POZ", load_box("3POZ"))
    print("  3RCD", load_box("3RCD"))


def main() -> None:
    for path in SI_FILES:
        patch(path)


if __name__ == "__main__":
    main()
