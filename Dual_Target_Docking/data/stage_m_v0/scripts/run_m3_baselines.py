#!/usr/bin/env python3
"""PART C — M3 trivial baselines vs docking arms."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ANALYSIS, BASELINE_ARMS, DOCK_ARMS, TABLES  # noqa: E402

# summary rule frozen: min(D/A, D/B)
SUMMARY = "summary_min"


def main():
    m1 = pd.read_csv(TABLES / "m1_directional_auroc.csv")
    m1 = m1[m1.subset == "all"].copy()
    rows = []
    pair_gates = {}

    for pair in ("EGFR_HER2", "PIK3CA_mTOR"):
        sub = m1[m1.pair == pair]
        vol = sub[sub.arm.isin(["heavy_atoms", "MW", "cLogP", "TPSA"])]
        best_vol_arm = vol.loc[vol[SUMMARY].idxmax(), "arm"]
        best_vol = float(vol.loc[vol[SUMMARY].idxmax(), SUMMARY])
        beat_any = False
        for _, r in sub.iterrows():
            is_dock = r.arm in DOCK_ARMS
            is_base = r.arm in BASELINE_ARMS
            fail = False
            if is_dock:
                fail = float(r[SUMMARY]) <= best_vol
                if not fail:
                    beat_any = True
            rows.append(
                {
                    "pair": pair,
                    "arm": r.arm,
                    "family": "docking" if is_dock else ("baseline" if is_base else "other"),
                    "auroc_D_vs_A": r.auroc_D_vs_A,
                    "auroc_D_vs_B": r.auroc_D_vs_B,
                    "summary_min": r.summary_min,
                    "summary_mean": r.summary_mean,
                    "best_volume_arm": best_vol_arm,
                    "best_volume_summary_min": round(best_vol, 4),
                    "fail_baseline": bool(fail) if is_dock else "",
                    "beats_best_volume": (not fail) if is_dock else "",
                }
            )
        # pair gate: docking arm beats volume on summary_min
        pair_gates[pair] = "Go" if beat_any else "No-Go"

    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "m3_baselines_vs_arms.csv", index=False)
    pd.DataFrame(
        [{"pair": k, "m3_gate": v} for k, v in pair_gates.items()]
    ).to_csv(TABLES / "m3_gate_summary.csv", index=False)

    md = []
    md.append("# M3 — Trivial baselines vs docking arms\n\n")
    md.append("## Rule (frozen)\n\n")
    md.append(
        "- Mandatory baselines: `heavy_atoms`, `MW`, `cLogP`, `TPSA` "
        "(plus optional `morgan_dual_medsim` = LOO median Tanimoto to other duals).\n"
        "- Arm comparison uses **summary = min(D/A, D/B)**.\n"
        "- Docking arm **fail_baseline** if `summary_min ≤ best_volume_summary_min` "
        "on that pair.\n"
        "- Report pairs **separately**; do not average EGFR with PIK3CA/mTOR.\n\n"
    )
    md.append("## Results\n\n")
    show = out[
        [
            "pair",
            "arm",
            "family",
            "auroc_D_vs_A",
            "auroc_D_vs_B",
            "summary_min",
            "best_volume_arm",
            "best_volume_summary_min",
            "fail_baseline",
        ]
    ]
    md.append(show.to_markdown(index=False) + "\n\n")
    md.append("## Gate (per pair)\n\n")
    for pair, g in pair_gates.items():
        md.append(f"- **{pair}: M3 = {g}**\n")
    md.append(
        "\nOverall note: EGFR expected No-Go (volume ≥ docking on weak end); "
        "PIK3CA/mTOR expected Go if docking exceeds volume on both directions.\n"
    )
    (ANALYSIS / "M3_BASELINES.md").write_text("".join(md))
    print(pair_gates)


if __name__ == "__main__":
    main()
