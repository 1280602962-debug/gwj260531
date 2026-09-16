#!/usr/bin/env python3
"""Trace figure annotations to source CSVs after hard-coded result literals were removed."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
OUT = ROOT / "remediation_outputs/figure_number_trace_post_fix.csv"

ROWS = [
    {
        "figure": "Figure 2B",
        "panel": "neither n= annotation",
        "annotation": "n={n_neg}",
        "source_file": "data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv (D_vs_neither_mean) or five_pair table2 n_neither",
        "source_row": "primary_row(pair)['n_neg']",
        "raw_value": "CSV n_neither / n_neg",
        "display_value": "f-string from CSV; PIK3CA/mTOR historically n=4 is no longer a literal",
    },
    {
        "figure": "Figure 3",
        "panel": "fixed-score ΔAUROC EGFR/HER2 and JAK1/TYK2",
        "annotation": "Δ from equal-score tables",
        "source_file": "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv; five_pair_stack_v1/equal_score_negative_s34_v1.csv",
        "source_row": "contrast D_vs_B_or_neither_pocketA",
        "raw_value": "CSV delta (historically 0.378 / 0.444)",
        "display_value": "read from CSV; audit_figures_pr32.py asserts vs CSV",
    },
    {
        "figure": "Figure 2A",
        "panel": "directional AUROC / summary_min",
        "annotation": "points and CIs",
        "source_file": "unified_threshold_sensitivity_v2.csv; five_pair table2_comparable_theta6_v1.csv",
        "source_row": "primary_row(pair)",
        "raw_value": "auroc_D_vs_A / auroc_D_vs_B / pocket_matched_summary_min",
        "display_value": "plotted from CSV",
    },
    {
        "figure": "Figure 4",
        "panel": "ECFP4 vs incremental",
        "annotation": "AUROC bars",
        "source_file": "incremental_information_v1.csv; ecfp4_incremental_s20s24_v1.csv",
        "source_row": "pair, contrast",
        "raw_value": "CSV AUROC",
        "display_value": "plotted from CSV",
    },
    {
        "figure": "Figure 5",
        "panel": "matched-minus-mismatched / holdout",
        "annotation": "delta and CI",
        "source_file": "wrong_pocket_paired_delta_bootstrap_v1.csv; holdout_pocket_matched_v1.csv; five_pair wrong_pocket_by_channel",
        "source_row": "pair, set/variant",
        "raw_value": "CSV delta / CI",
        "display_value": "plotted from CSV",
    },
    {
        "figure": "scripts",
        "panel": "hard-coded experimental results",
        "annotation": "REMOVED 0.378, 0.444, n=4",
        "source_file": "figures/jcim_article/scripts/plot_jcim_article_figures_v3.py; audit_figures_pr32.py",
        "source_row": "n/a",
        "raw_value": "none remaining for official results",
        "display_value": "axis limits / 0 / 0.5 reference lines only",
    },
]


def main() -> int:
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ROWS[0].keys()))
        w.writeheader()
        w.writerows(ROWS)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
