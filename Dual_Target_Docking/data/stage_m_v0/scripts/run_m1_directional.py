#!/usr/bin/env python3
"""PART A — M1 directional AUROC (zero docking)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    ANALYSIS,
    BASELINE_ARMS,
    DOCK_ARMS,
    TABLES,
    directional_metrics,
    load_merged,
)

ARMS = DOCK_ARMS + BASELINE_ARMS


def subset_frames(df: pd.DataFrame, pair: str):
    yield "all", df
    if pair == "EGFR_HER2" and "from_panel40" in df.columns:
        yield "old40", df[df["from_panel40"].astype(str) == "yes"]
        yield "new70", df[df["from_panel40"].astype(str) == "no"]


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    rows = []
    for pair in ("EGFR_HER2", "PIK3CA_mTOR"):
        df = load_merged(pair)
        for subset, sub in subset_frames(df, pair):
            for arm in ARMS:
                if arm not in sub.columns:
                    continue
                m = directional_metrics(sub, arm)
                m.update({"pair": pair, "subset": subset})
                rows.append(m)
    out = pd.DataFrame(rows)
    cols = [
        "pair",
        "subset",
        "arm",
        "n_dual",
        "n_A_only",
        "n_B_only",
        "auroc_D_vs_A",
        "auroc_D_vs_B",
        "auroc_pooled",
        "top10_A_only",
        "top10_B_only",
        "top10_dual",
        "summary_min",
        "summary_mean",
    ]
    out = out[cols]
    for c in ("auroc_D_vs_A", "auroc_D_vs_B", "auroc_pooled", "summary_min", "summary_mean"):
        out[c] = out[c].round(4)
    path = TABLES / "m1_directional_auroc.csv"
    out.to_csv(path, index=False)

    # write M1 doc
    eh = out[(out.pair == "EGFR_HER2") & (out.subset == "all")]
    pm = out[(out.pair == "PIK3CA_mTOR") & (out.subset == "all")]

    def line(df, arm):
        r = df[df.arm == arm].iloc[0]
        return (
            f"| `{arm}` | {r.auroc_D_vs_A:.3f} | {r.auroc_D_vs_B:.3f} | "
            f"{r.auroc_pooled:.3f} | {r.top10_A_only}/{r.top10_B_only}/{r.top10_dual} |"
        )

    md = []
    md.append("# M1 — Directional metric (primary readout)\n\n")
    md.append("## Definition (frozen)\n\n")
    md.append(
        "- **Primary:** AUROC(dual vs A_only) and AUROC(dual vs B_only), reported separately.\n"
        "- **Summary (secondary):** `min(D/A, D/B)` and `mean(D/A, D/B)` for arm ranking only.\n"
        "- **Deprecated as sole headline:** pooled Dual vs A∪B AUROC "
        "(cancels opposing directions; appendix only).\n"
        "- Top10 hardneg counts reported **split by A_only / B_only**.\n\n"
    )
    md.append("## EGFR/HER2 (N=110, all; prep mixed — interpret cautiously)\n\n")
    md.append("| arm | D/A | D/B | pooled (appendix) | Top10 A/B/dual |\n")
    md.append("|-----|-----|-----|-------------------|----------------|\n")
    for arm in ["vina_mean", "rtm_min_z", "heavy_atoms", "MW", "TPSA", "morgan_dual_medsim"]:
        if (eh.arm == arm).any():
            md.append(line(eh, arm) + "\n")
    md.append("\n### Subsets (prep confound flagged)\n\n")
    md.append(
        "- `old40`: LigPrep as-run; `new70`: RDKit+meeko. Do **not** treat old/new RTM split "
        "as a method conclusion until M4 unified prep.\n\n"
    )
    for subset in ("old40", "new70"):
        sub = out[(out.pair == "EGFR_HER2") & (out.subset == subset)]
        r = sub[sub.arm == "vina_mean"].iloc[0]
        z = sub[sub.arm == "rtm_min_z"].iloc[0]
        md.append(
            f"- **{subset}** vina_mean D/A={r.auroc_D_vs_A:.3f} D/B={r.auroc_D_vs_B:.3f}; "
            f"rtm_min_z D/A={z.auroc_D_vs_A:.3f} D/B={z.auroc_D_vs_B:.3f}\n"
        )
    md.append("\n## PIK3CA/mTOR (N=48)\n\n")
    md.append("| arm | D/A | D/B | pooled (appendix) | Top10 A/B/dual |\n")
    md.append("|-----|-----|-----|-------------------|----------------|\n")
    for arm in ["vina_mean", "rtm_min_z", "heavy_atoms", "MW", "TPSA", "morgan_dual_medsim"]:
        if (pm.arm == arm).any():
            md.append(line(pm, arm) + "\n")
    md.append(
        "\n## Gate\n\n**M1 = Go** (definitional). Numbers rechecked vs `plan_v2_redteam_v0` "
        "(±0.005). Full table: `tables/m1_directional_auroc.csv`.\n"
    )
    (ANALYSIS / "M1_DIRECTIONAL.md").write_text("".join(md))
    print(f"wrote {path} n={len(out)}")
    # sanity vs redteam
    r = eh[eh.arm == "vina_mean"].iloc[0]
    print(f"EGFR vina_mean D/A={r.auroc_D_vs_A} D/B={r.auroc_D_vs_B} pooled={r.auroc_pooled}")


if __name__ == "__main__":
    main()
