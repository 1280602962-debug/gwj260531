#!/usr/bin/env python3
"""
Build a demo screening library from curated JNK datasets.

Creates data/libraries/screening_demo.smi with unique SMILES from all
isoform curated tables (suitable for end-to-end pipeline testing).

Usage:
    python scripts/build_demo_library.py
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
OUTPUT = ROOT / "data" / "libraries" / "screening_demo.smi"

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def main() -> None:
    smiles: set[str] = set()
    for iso in ["jnk1", "jnk2", "jnk3"]:
        path = PROCESSED / f"{iso}_curated.csv"
        if path.exists():
            df = pd.read_csv(path)
            smiles.update(df["canonical_smiles"].dropna().unique())

    paired = PROCESSED / "paired_set.csv"
    if paired.exists():
        df = pd.read_csv(paired)
        smiles.update(df["canonical_smiles"].dropna().unique())

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT, "w") as f:
        f.write("# Demo screening library built from curated JNK datasets\n")
        for smi in sorted(smiles):
            f.write(f"{smi}\n")

    logger.info("Wrote %d unique SMILES → %s", len(smiles), OUTPUT)


if __name__ == "__main__":
    main()
