#!/usr/bin/env python3
"""C11 — In-silico opportunity-cost table: unbought 2231 vs purchased 690/2157."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "c11_2231_comparison"
OUT.mkdir(parents=True, exist_ok=True)

SHORT = ROOT / "data/shortlist/md_shortlist_final.csv"
TABLE27 = ROOT / "docs/popular_science/data_tables/27_MD16_选择性排序与报价.csv"
PURCHASE = ROOT / "docs/popular_science/data_tables/11_采购清单purchase_after_md.csv"


def main():
    short = pd.read_csv(SHORT)
    t27 = pd.read_csv(TABLE27)
    # normalize id column
    t27["compound_id"] = t27["库内ID"].astype(str)

    ids = ["690", "2157", "2231"]
    rows = []
    for cid in ids:
        s = short[short["compound_id"].astype(str) == cid]
        if s.empty:
            s = short[short["title"].astype(str) == cid]
        r = s.iloc[0]
        m = t27[t27["compound_id"] == cid]
        mrow = m.iloc[0] if len(m) else None
        delta = float(r["score_JNK2"] - r["score_JNK1"])  # not exact archived Δsel; report archived if present
        rows.append(
            {
                "compound_id": cid,
                "HIT_ID": None if mrow is None else mrow.get("HIT_ID"),
                "purchased": cid in {"690", "2157"},
                "group": r.get("group"),
                "smiles": r.get("smiles"),
                "score_JNK1": r.get("score_JNK1"),
                "score_JNK2": r.get("score_JNK2"),
                "score_JNK3": r.get("score_JNK3"),
                "delta_sel_dock_table27": None if mrow is None else mrow.get("delta_sel_dock"),
                "mmgbsa_JNK1": r.get("mmgbsa_JNK1"),
                "md_jnk1_bias_score": None if mrow is None else mrow.get("md_jnk1_bias_score"),
                "hinge_JNK1": None if mrow is None else mrow.get("hinge_JNK1"),
                "hinge_JNK2": None if mrow is None else mrow.get("hinge_JNK2"),
                "hinge_JNK3": None if mrow is None else mrow.get("hinge_JNK3"),
                "RMSD_JNK1": None if mrow is None else mrow.get("RMSD_JNK1"),
                "RMSD_JNK2": None if mrow is None else mrow.get("RMSD_JNK2"),
                "RMSD_JNK3": None if mrow is None else mrow.get("RMSD_JNK3"),
                "pass_md_overall": None if mrow is None else mrow.get("pass_md_overall"),
                "pose_grade": None if mrow is None else mrow.get("pose_grade"),
                "chemotype_sim": r.get("chemotype_sim"),
                "AD_maxTc": r.get("AD_maxTc"),
                "price_CNY_table27": None if mrow is None else mrow.get("小计_元"),
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "c11_690_2157_vs_2231.csv", index=False)

    narrative = {
        "decision": "Purchased 690+2157; did not purchase 2231",
        "why_2231_looked_attractive": (
            "Highest MD JNK1-bias score among shortlist; strongest hinge asymmetry "
            "(hinge JNK1≈0.91 vs JNK2≈0); best Glide score_JNK1 (−11.22)."
        ),
        "why_not_purchased_for_RQ_A": (
            "pass_md_overall = NO / pose_grade C — fails project MD overall gate "
            "(requires JNK1 pass AND (JNK2 OR JNK3) pass). Purchase set prioritized "
            "pose-credible family binders (grade A), not the strongest MD-bias hypothesis."
        ),
        "tradeoff_for_RQ_B": (
            "Omitting 2231 weakens prospective test of MD-predicted JNK1 preference; "
            "2157 is only a secondary bias hypothesis and has anti-JNK1 Δsel_dock (−1.05)."
        ),
        "how_to_use_in_paper": (
            "SI/Discussion: document opportunity cost; do not claim 2231 inactivity; "
            "optional future buy if budget allows."
        ),
    }
    (OUT / "c11_narrative.json").write_text(json.dumps(narrative, indent=2), encoding="utf-8")

    md = [
        "# C11 In-silico Comparison: 2231 (unbought) vs 690/2157 (purchased)",
        "",
        "## Table",
        "",
        df.to_markdown(index=False),
        "",
        "## Narrative (locked with Option A)",
        "",
        f"- **Attractive for RQ-B:** {narrative['why_2231_looked_attractive']}",
        f"- **Excluded for RQ-A gate:** {narrative['why_not_purchased_for_RQ_A']}",
        f"- **Tradeoff:** {narrative['tradeoff_for_RQ_B']}",
        f"- **Paper use:** {narrative['how_to_use_in_paper']}",
        "",
    ]
    (OUT / "C11_2231_COMPARISON.md").write_text("\n".join(md), encoding="utf-8")
    print(df[["compound_id", "purchased", "pose_grade", "pass_md_overall", "md_jnk1_bias_score", "delta_sel_dock_table27", "score_JNK1"]].to_string(index=False))


if __name__ == "__main__":
    main()
