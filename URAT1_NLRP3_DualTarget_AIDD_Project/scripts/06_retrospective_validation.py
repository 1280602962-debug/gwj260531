#!/usr/bin/env python3
"""
TAPE-GATE Stage 6: Retrospective validation + ablations (incl. PLK1-style baseline).

Reports Recall@K separately for Path A, Path B, and merged union.
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


def load_must_recover_names(benchmarks: list[dict]) -> list[str]:
    """Positive benchmark compounds from literature_benchmarks.csv."""
    names: list[str] = []
    seen: set[str] = set()
    for b in benchmarks:
        role = b.get("validation_role", "")
        if "must_recover" in role and b.get("compound_name") not in seen:
            names.append(b["compound_name"])
            seen.add(b["compound_name"])
    return names


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, help="Fused ranked candidates CSV")
    parser.add_argument("--top-k", type=int, default=500)
    parser.add_argument("--run-plk1-baseline", action="store_true")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "results" / "validation")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    benchmarks = load_benchmarks()
    must_recover = load_must_recover_names(benchmarks)

    report = {
        "status": "skeleton",
        "framework": "TAPE-GATE",
        "benchmark_compounds": must_recover,
        "metrics_by_path": {
            "path_a_library": ["Recall@100", "Recall@500"],
            "path_b_generative": ["Recall@100", "Recall@500", "mean_novelty_Tc"],
            "union": ["Recall@100", "Recall@500", "MRR"],
        },
        "ablations": [
            "Abl-1: no S_trap (single PDB URAT1)",
            "Abl-2: NLRP3 anchor similarity (PLK1-style) vs assay-conditioned",
            "Abl-3: fixed 0.5/0.5 fusion vs reliability Pareto",
            "Abl-4: Path A only (no generative)",
            "Abl-5: Path B only (no library)",
            "Abl-6: MTL vs independent dual models",
            "Abl-7: no SLC22 transfer",
        ],
        "plk1_style_baseline": {
            "enabled": args.run_plk1_baseline,
            "method": "SVR(URAT1) + 5-anchor ECFP similarity(NLRP3) + 0.5 fusion",
            "purpose": "negative control — not TAPE-GATE primary method",
        },
        "pass_criteria": {
            "note": "Criteria apply after full pipeline implementation; not auto-passed in skeleton mode",
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
