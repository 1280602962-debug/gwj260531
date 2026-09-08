#!/usr/bin/env python3
"""OAT-transfer ablation for URAT1 ML (not docking protocol P2).

Compares URAT1 training with vs without OAT1/OAT3 sequential pretrain.
OAT transfer is not a manuscript main result.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT_ROOT / "scripts"
RESULTS = PROJECT_ROOT / "results" / "training"


def _run(cmd: list[str]) -> None:
    print(f"\n>>> {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=SCRIPTS)


def _load_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def main() -> None:
    parser = argparse.ArgumentParser(description="OAT transfer ablation (baseline vs transfer)")
    parser.add_argument("--output-root", type=Path, default=RESULTS)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    baseline_dir = args.output_root / "ablation_no_oat"
    transfer_dir = args.output_root / "ablation_oat_transfer"
    benchmark_root = PROJECT_ROOT / "results" / "benchmark_backtest"

    train = SCRIPTS / "02_train_asymmetric_models.py"
    bench = SCRIPTS / "07_benchmark_backtest.py"

    _run([
        sys.executable, str(train),
        "--no-oat-transfer",
        "--output", str(baseline_dir),
        "--n-splits", str(args.n_splits),
        "--seed", str(args.seed),
    ])
    _run([
        sys.executable, str(bench),
        "--model-dir", str(baseline_dir),
        "--output", str(benchmark_root / "ablation_no_oat"),
    ])

    _run([
        sys.executable, str(train),
        "--oat-transfer",
        "--output", str(transfer_dir),
        "--n-splits", str(args.n_splits),
        "--seed", str(args.seed),
    ])
    _run([
        sys.executable, str(bench),
        "--model-dir", str(transfer_dir),
        "--output", str(benchmark_root / "ablation_oat_transfer"),
    ])

    base_train = _load_json(baseline_dir / "training_report.json")
    xfer_train = _load_json(transfer_dir / "training_report.json")
    base_bench = _load_json(benchmark_root / "ablation_no_oat" / "benchmark_backtest_report.json")
    xfer_bench = _load_json(benchmark_root / "ablation_oat_transfer" / "benchmark_backtest_report.json")

    def urat1_summary(train_report: dict, bench_report: dict) -> dict:
        cv = train_report["urat1"]["cv_metrics"]
        bt = bench_report["urat1_backtest"]
        return {
            "spearman_oof": cv["spearman"],
            "r2_oof": cv["r2"],
            "rmse_oof": cv["rmse"],
            "roc_auc_p7": cv.get("roc_auc_p7"),
            "ef_5pct_p7": cv.get("ef_5pct_p7"),
            "screening_suitable": train_report["urat1"]["screening_assessment"]["suitable_for_screening"],
            "benchmark_must_recover_pass": f"{bt['must_recover_binary_pass']}/{bt['must_recover_count']}",
            "overall_verdict": bench_report["overall_verdict"]["verdict"],
        }

    base_u = urat1_summary(base_train, base_bench)
    xfer_u = urat1_summary(xfer_train, xfer_bench)

    def delta(key: str) -> float | None:
        b, t = base_u.get(key), xfer_u.get(key)
        if isinstance(b, (int, float)) and isinstance(t, (int, float)):
            return float(t - b)
        return None

    ablation = {
        "experiment": "oat_transfer_ablation",
        "baseline": {
            "label": "no_oat_transfer",
            "model_dir": str(baseline_dir),
            "benchmark_dir": str(benchmark_root / "ablation_no_oat"),
            "urat1": base_u,
        },
        "oat_transfer": {
            "label": "oat1_oat3_sequential_finetune",
            "model_dir": str(transfer_dir),
            "benchmark_dir": str(benchmark_root / "ablation_oat_transfer"),
            "urat1": xfer_u,
        },
        "delta_transfer_minus_baseline": {
            "spearman_oof": delta("spearman_oof"),
            "r2_oof": delta("r2_oof"),
            "rmse_oof": delta("rmse_oof"),
            "roc_auc_p7": delta("roc_auc_p7"),
            "ef_5pct_p7": delta("ef_5pct_p7"),
        },
        "interpretation": (
            "Positive delta on spearman/r2/roc_auc_p7 suggests OAT auxiliary pretrain helps URAT1. "
            "Benchmark must-recover pass (x/4) is the primary go/no-go for library screening."
        ),
    }

    out_path = args.output_root / "oat_transfer_ablation.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(ablation, f, indent=2)

    print("\n=== OAT Transfer Ablation (URAT1) ===")
    print(f"  Baseline  Spearman={base_u['spearman_oof']:.3f}  R2={base_u['r2_oof']:.3f}  "
          f"bench={base_u['benchmark_must_recover_pass']}  verdict={base_u['overall_verdict']}")
    print(f"  Transfer  Spearman={xfer_u['spearman_oof']:.3f}  R2={xfer_u['r2_oof']:.3f}  "
          f"bench={xfer_u['benchmark_must_recover_pass']}  verdict={xfer_u['overall_verdict']}")
    print(f"\nReport: {out_path}")


if __name__ == "__main__":
    main()
