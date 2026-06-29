#!/usr/bin/env python3
"""
STAD-AIDD Stage 5: Retrospective benchmark validation + ablation.

Must recover: lesinurad, benzbromarone, MCC950, GDC-2394 in top-K.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_PATH = PROJECT_ROOT / "data" / "benchmarks" / "literature_benchmarks.csv"


def load_benchmarks() -> list[dict]:
    if not BENCHMARK_PATH.exists():
        return []
    with open(BENCHMARK_PATH) as f:
        return list(csv.DictReader(f))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, help="Ranked candidates CSV from screening")
    parser.add_argument("--top-k", type=int, default=500)
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "results" / "validation")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    benchmarks = load_benchmarks()
    must_recover = [b["compound_name"] for b in benchmarks if b.get("validation_role") == "retrospective_must_recover"]

    report = {
        "status": "skeleton",
        "benchmark_compounds": must_recover,
        "metrics": ["Recall@100", "Recall@500", "MRR"],
        "ablations": [
            "Abl-1: single PDB docking",
            "Abl-2: no foundation model",
            "Abl-3: no SLC22 transfer",
            "Abl-4: single-target intersection",
            "Abl-5: no RL generation",
        ],
        "pass_criteria": {
            "urat1_drugs_in_top500": ">= 2/4",
            "nlrp3_tools_in_top500": ">= 1/2",
        },
    }
    out = args.output / "validation_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Validation template: {out}")
    print(f"Benchmark compounds to recover: {must_recover}")


if __name__ == "__main__":
    main()
