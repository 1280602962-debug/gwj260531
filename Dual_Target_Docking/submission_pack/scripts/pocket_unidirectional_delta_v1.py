#!/usr/bin/env python3
"""Unidirectional matched-vs-wrong pocket deltas next to Delta summary_min.

Uses already-locked CSVs. Does not redock.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "local_track_b_v0" / "tables" / "pocket_unidirectional_delta_v1.csv"

ORIG = ROOT.parent / "jcim_strengthen_t0t1_v0" / "tables" / "wrong_pocket_paired_delta_bootstrap_v1.csv"
TRACKB = ROOT / "local_track_b_v0" / "tables" / "five_pair_local_channels_v1" / "wrong_pocket_by_channel_v1.csv"

KEEP_CHANNELS = {
    "vina_20260727": "main_panel",
    "holdout_vina_20260727": "unused_pool_holdout",
}


def weaker(a: float, b: float) -> str:
    if a < b:
        return "D_vs_A"
    if b < a:
        return "D_vs_B"
    return "tie"


def rows_from(path: Path, set_col: str | None, channel_filter: dict[str, str] | None):
    out = []
    for r in csv.DictReader(path.open()):
        if channel_filter is not None:
            ch = r.get("channel", "")
            if ch not in channel_filter:
                continue
            dataset = channel_filter[ch]
        else:
            dataset = r[set_col]
        da_m = float(r["matched_D_vs_A"])
        db_m = float(r["matched_D_vs_B"])
        da_w = float(r["wrong_D_vs_A"])
        db_w = float(r["wrong_D_vs_B"])
        out.append(
            {
                "pair": r["pair"],
                "dataset": dataset,
                "n_dual": r["n_dual"],
                "n_A_only": r["n_A_only"],
                "n_B_only": r["n_B_only"],
                "matched_D_vs_A": f"{da_m:.4f}",
                "wrong_D_vs_A": f"{da_w:.4f}",
                "delta_D_vs_A": f"{da_m - da_w:.4f}",
                "matched_D_vs_B": f"{db_m:.4f}",
                "wrong_D_vs_B": f"{db_w:.4f}",
                "delta_D_vs_B": f"{db_m - db_w:.4f}",
                "matched_summary_min": f"{float(r['matched_summary_min']):.4f}",
                "wrong_summary_min": f"{float(r['wrong_summary_min']):.4f}",
                "delta_summary_min": f"{float(r['delta_matched_minus_wrong']):.4f}",
                "delta_summary_min_ci": f"[{float(r['delta_ci_lo']):.3f}, {float(r['delta_ci_hi']):.3f}]",
                "weaker_matched": weaker(da_m, db_m),
                "weaker_wrong": weaker(da_w, db_w),
                "weaker_switched": "yes" if weaker(da_m, db_m) != weaker(da_w, db_w) else "no",
                "weaker_switched_definition": "matched_vs_mismatched_weaker_arm",
                "unidirectional_delta_is": "point_estimate_only",
                "ci_applies_to": "delta_summary_min_only",
                "derivation": "rearranged_from_locked_wrong_pocket_tables",
            }
        )
    return out


def main() -> None:
    rows = rows_from(ORIG, "set", None) + rows_from(TRACKB, None, KEEP_CHANNELS)
    order = [
        "EGFR/HER2",
        "AChE/BChE",
        "PIK3CA/mTOR",
        "F2/F10",
        "JAK1/TYK2",
        "JAK1/JAK2",
        "PPARG/PPARA",
        "PPARA/PPARD",
    ]
    set_order = {"main_panel": 0, "unused_pool_holdout": 1}
    rows.sort(key=lambda r: (order.index(r["pair"]), set_order[r["dataset"]]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {OUT}")


if __name__ == "__main__":
    main()
