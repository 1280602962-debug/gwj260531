#!/usr/bin/env python3
"""C11 — Purchase rationale table for 690 + 2231 (updated buy set)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "c11_2231_comparison"
OUT.mkdir(parents=True, exist_ok=True)

SHORT = ROOT / "data/shortlist/md_shortlist_final.csv"
TABLE27 = ROOT / "docs/popular_science/data_tables/27_MD16_选择性排序与报价.csv"


def main():
    short = pd.read_csv(SHORT)
    t27 = pd.read_csv(TABLE27)
    t27["compound_id"] = t27["库内ID"].astype(str)

    ids = ["690", "2231", "2157"]  # 2157 retained as not-purchased comparator
    rows = []
    for cid in ids:
        s = short[short["compound_id"].astype(str) == cid]
        if s.empty:
            s = short[short["title"].astype(str) == cid]
        r = s.iloc[0]
        m = t27[t27["compound_id"] == cid]
        mrow = m.iloc[0] if len(m) else None
        rows.append(
            {
                "compound_id": cid,
                "HIT_ID": None if mrow is None else mrow.get("HIT_ID"),
                "purchased": cid in {"690", "2231"},
                "role": (
                    "activity/pan-leaning anchor"
                    if cid == "690"
                    else (
                        "MD JNK1-bias hypothesis"
                        if cid == "2231"
                        else "not purchased (former alt)"
                    )
                ),
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
                "price_CNY_table27": None if mrow is None else mrow.get("小计_元"),
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "c11_690_2231_purchase_panel.csv", index=False)

    narrative = {
        "decision": "Purchased 690 + 2231 (replaces prior 690+2157 plan)",
        "why_690": (
            "G1 grade-A, pass_md_overall, pan-leaning hinge (high on all isoforms) — "
            "RQ-A family-activity anchor with credible pose QC."
        ),
        "why_2231": (
            "Strongest MD JNK1-bias score; hinge J1≫J2; favorable Δsel_dock (+3.37) and "
            "best score_JNK1 (−11.22). Purchased to give RQ-B a real prospective test."
        ),
        "acknowledged_risk": (
            "2231 pose_grade C and pass_md_overall = NO under archived MD gate. "
            "Extended 200 ns used ligand restraints — must not be cited as unrestrained proof. "
            "Week 2–3 priority: unrestrained MD replicas (C3) before over-interpreting bias."
        ),
        "why_not_2157": (
            "Budget/slot tradeoff; 2157 had MD bias #2 but anti-JNK1 Δsel_dock (−1.05). "
            "Kept as in-silico comparator only."
        ),
        "option_A_framing": (
            "Core contribution remains selectivity-predictor failure + family pipeline. "
            "2231 upgrades secondary RQ-B; null/mixture outcomes remain publishable."
        ),
    }
    (OUT / "c11_narrative.json").write_text(json.dumps(narrative, indent=2), encoding="utf-8")

    md = [
        "# C11 Purchase Panel: 690 + 2231 (updated)",
        "",
        "## Table",
        "",
        df.to_markdown(index=False),
        "",
        "## Narrative",
        "",
        f"- **690:** {narrative['why_690']}",
        f"- **2231:** {narrative['why_2231']}",
        f"- **Risk:** {narrative['acknowledged_risk']}",
        f"- **2157:** {narrative['why_not_2157']}",
        f"- **Option A:** {narrative['option_A_framing']}",
        "",
    ]
    (OUT / "C11_2231_COMPARISON.md").write_text("\n".join(md), encoding="utf-8")
    print(df[["compound_id", "purchased", "role", "pose_grade", "pass_md_overall", "md_jnk1_bias_score", "delta_sel_dock_table27"]].to_string(index=False))


if __name__ == "__main__":
    main()
