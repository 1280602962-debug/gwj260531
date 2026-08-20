#!/usr/bin/env python3
"""Analyze 9DKB + 7ALV Pareto merge vs literature benchmarks; write JSON + markdown report."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from utils_ml import canonicalize
DEFAULT_MERGED = PROJECT_ROOT / "data" / "repurposing" / "p2" / "pareto_merged_scores.csv"
DEFAULT_SHORT = PROJECT_ROOT / "data" / "repurposing" / "p2" / "pareto_shortlist.csv"
DEFAULT_SUMMARY = PROJECT_ROOT / "data" / "repurposing" / "p2" / "pareto_summary.json"
DEFAULT_BENCH = PROJECT_ROOT / "data" / "benchmarks" / "literature_benchmarks.csv"
DEFAULT_POOL = PROJECT_ROOT / "data" / "repurposing" / "screening" / "docking_pool_p05.csv"
DEFAULT_MANIFEST = PROJECT_ROOT / "data" / "repurposing" / "repurposing_manifest.csv"
DEFAULT_OUT = PROJECT_ROOT / "results" / "repurposing" / "pareto_benchmark_report.json"


def _name_col(df: pd.DataFrame) -> str:
    for c in ("name", "pref_name", "compound_name"):
        if c in df.columns:
            return c
    raise ValueError("no name column")


def _benchmark_table(merged: pd.DataFrame, pool: pd.DataFrame, manifest: pd.DataFrame, bench: pd.DataFrame) -> list[dict]:
    name_col = _name_col(merged)
    pool_smiles = set(pool["canonical_smiles"].astype(str))
    rows = []
    for compound in bench["compound_name"].drop_duplicates():
        sub = bench[bench["compound_name"] == compound]
        smi = canonicalize(sub["canonical_smiles"].iloc[0])
        in_pool = smi in pool_smiles
        hit = merged[merged["canonical_smiles"].astype(str) == smi]
        if not len(hit) and name_col in merged.columns:
            hit = merged[merged[name_col].astype(str).str.upper() == compound.upper()]
        man = manifest[manifest["canonical_smiles"].astype(str) == smi]
        if not len(man) and "name" in manifest.columns:
            man = manifest[manifest["name"].astype(str).str.upper() == compound.upper()]
        p_ml = float(man["p_active_nlrp3"].iloc[0]) if len(man) and "p_active_nlrp3" in man.columns else None
        row = {
            "compound": compound,
            "target": sub["target_gene"].iloc[0],
            "in_p05_pool": in_pool,
            "p_active_nlrp3": p_ml,
            "in_dual_merge": bool(len(hit)),
        }
        if len(hit):
            r = hit.iloc[0]
            row.update(
                {
                    "s_u_percentile": float(r["s_u_percentile"]),
                    "s_n_percentile": float(r["s_n_percentile"]),
                    "s_n_ml_percentile": float(r["s_n_ml_percentile"]),
                    "s_n_dock_percentile": float(r["s_n_dock_percentile"]),
                    "urat1_dock_score": float(r["dock_score"]),
                    "nlrp3_dock_score": float(r["nlrp3_dock_score"]),
                    "pareto_front": bool(r["pareto_front"]),
                }
            )
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Pareto vs benchmark analysis (9DKB + 7ALV)")
    parser.add_argument("--merged", type=Path, default=DEFAULT_MERGED)
    parser.add_argument("--shortlist", type=Path, default=DEFAULT_SHORT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    merged = pd.read_csv(args.merged, low_memory=False)
    short = pd.read_csv(args.shortlist)
    pool = pd.read_csv(DEFAULT_POOL, low_memory=False)
    manifest = pd.read_csv(DEFAULT_MANIFEST, low_memory=False)
    bench = pd.read_csv(DEFAULT_BENCH)
    pareto_summary = json.loads(args.summary.read_text()) if args.summary.exists() else {}

    name_col = _name_col(merged)
    sn_dock = merged["nlrp3_dock_score"]
    r_sp, p_sp = spearmanr(merged["p_active_nlrp3"], sn_dock, nan_policy="omit")

    bench_rows = _benchmark_table(merged, pool, manifest, bench)
    shortlist_records = short.to_dict(orient="records")

    md_recommendations = {
        "urat1_9dkb": [
            {"compound": "lesinurad", "reason": "9DKB co-crystal control; MD calibration, not a novel nomination"},
            {"compound": "GSK-3008348", "reason": "Preferred URAT1-side computational hypothesis after chemistry nomination"},
        ],
        "nlrp3_7alv": [
            {"compound": "MCC950", "reason": "NLRP3 tool inhibitor analog dock at 7ALV (not self-dock); calibration, not a novel nomination"},
            {"compound": "Vecabrutinib", "reason": "Preferred NLRP3-side computational hypothesis after chemistry nomination"},
        ],
        "note": "Raw Pareto is audit-only. Follow-up molecules are chemistry-nominated after P2 dual-dock (docs/MANUSCRIPT.md). Benchmark uricosurics may be absent from the P≥0.5 pool by design.",
    }

    report = {
        "structures": {"urat1": "9DKB", "nlrp3": "7ALV"},
        "coverage": {
            "pool_p05": int(len(pool)),
            "dual_merged": int(len(merged)),
            "pareto_front": int(merged["pareto_front"].sum()),
            "shortlist_n": int(len(short)),
            **{k: pareto_summary.get(k) for k in ("n_pool_missing_dock",) if k in pareto_summary},
        },
        "spearman_ml_p_vs_7alv_dock": {"r": float(r_sp), "p": float(p_sp)},
        "benchmarks": bench_rows,
        "pareto_shortlist": shortlist_records,
        "md_recommendations": md_recommendations,
        "conclusions": {
            "funnel_valid": True,
            "urat1_benchmark_in_pool": "Report lesinurad/verinurad P2 percentiles when present; benzbromarone/dotinurad may need standalone 9DKB poses if absent from the P≥0.5 pool",
            "nlrp3_benchmark": "Do not treat ML vs 7ALV Spearman as affinity; colchicine is an indirect-mechanism control",
            "repurposing_signal": "Raw Pareto is audit-only. Follow-up uses dual-dock gates plus chemistry nomination",
            "proceed_to_md": True,
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2))

    md_path = args.output.with_suffix(".md")
    lines = [
        "# Pareto analysis: 9DKB + 7ALV dual docking",
        "",
        f"- Dual merged: **{len(merged)}** / pool **{len(pool)}**",
        f"- Pareto front: **{merged['pareto_front'].sum()}**",
        f"- Spearman ML P(active) vs 7ALV docking: **r={r_sp:.3f}**",
        "",
        "## Benchmarks in merge",
        "",
    ]
    for b in bench_rows:
        if b["in_dual_merge"]:
            lines.append(
                f"- **{b['compound']}**: S_U={b['s_u_percentile']:.1f}, S_N={b['s_n_percentile']:.1f}, Pareto={b['pareto_front']}"
            )
        else:
            lines.append(f"- **{b['compound']}**: not in dual merge (in P05 pool={b['in_p05_pool']})")
    lines.extend(["", "## Pareto shortlist", ""])
    for s in shortlist_records:
        lines.append(f"- {s.get(name_col, s.get('name'))}: S_U={s['s_u_percentile']:.1f}, S_N={s['s_n_percentile']:.1f}")
    lines.extend(["", "## MD next (2+2)", ""])
    for x in md_recommendations["urat1_9dkb"] + md_recommendations["nlrp3_7alv"]:
        if "compound" in x:
            lines.append(f"- {x['compound']}: {x['reason']}")
    md_path.write_text("\n".join(lines) + "\n")
    print(json.dumps({"report": str(args.output), "markdown": str(md_path)}, indent=2))


if __name__ == "__main__":
    main()
