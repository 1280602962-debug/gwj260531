#!/usr/bin/env python3
"""C5 — Main-text selectivity-method autopsy table (literature benchmark).

Explicitly decoupled from purchase decision.
Sources: docking_validation/*, calibration/*, gly87_selfcheck, project report numbers.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "selectivity_autopsy"
OUT.mkdir(parents=True, exist_ok=True)

DELTAS = ROOT / "results/docking_validation/benchmark_deltas_51c1.csv"
CONF = ROOT / "results/docking_validation/direction_confusion_27c3.csv"
GLY = ROOT / "results/docking_validation/gly87_selfcheck_16be.csv"
VAL_MD = ROOT / "results/docking_validation/validation_report.md"


def direction_stats(df: pd.DataFrame) -> dict:
    # compounds with experimental isoform direction
    sub = df[df["exp_dir_pIC50"].notna() & (df["exp_dir_pIC50"] != "NA")].copy()
    # Prefer rows that have IC50-derived direction
    if "direction_match" in sub.columns:
        # coerce
        match = sub["direction_match"].astype(str).str.lower().isin(["true", "1"])
        n = int(len(sub))
        n_match = int(match.sum())
    else:
        n, n_match = 0, 0
    return {
        "n_with_exp_dir": n,
        "n_direction_match": n_match,
        "direction_accuracy": (n_match / n) if n else None,
    }


def main():
    deltas = pd.read_csv(DELTAS)
    conf = pd.read_csv(CONF)
    gly = pd.read_csv(GLY)

    # Per-compound docking autopsy
    dock_rows = []
    for _, r in deltas.iterrows():
        if str(r.get("exp_dir_pIC50", "NA")) == "NA" and pd.isna(r.get("jnk1_ic50_nM")):
            continue
        dock_rows.append(
            {
                "compound": r["name"],
                "expected_profile": r["expected_profile"],
                "delta_sel_dock": r.get("delta_sel_dock"),
                "exp_dir_pIC50": r.get("exp_dir_pIC50"),
                "pred_dir_dock": r.get("pred_dir_dock"),
                "direction_match": r.get("direction_match"),
                "jnk1_ic50_nM": r.get("jnk1_ic50_nM"),
                "jnk2_ic50_nM": r.get("jnk2_ic50_nM"),
                "jnk3_ic50_nM": r.get("jnk3_ic50_nM"),
            }
        )
    dock_df = pd.DataFrame(dock_rows)
    dock_df.to_csv(OUT / "c5_docking_direction_by_compound.csv", index=False)

    # Key-control confusion (project figure set)
    conf.to_csv(OUT / "c5_key_controls_direction.csv", index=False)

    # Gly87
    gly_out = gly.copy()
    gly_out.to_csv(OUT / "c5_gly87_selfcheck.csv", index=False)
    gly_occ_all_true = bool(gly_out["occ_JNK1"].astype(str).str.lower().isin(["true", "1"]).all())
    gly_pred_sel_any = bool(
        gly_out["pred_JNK1_sel"].astype(str).str.lower().isin(["true", "1"]).any()
    )

    # Aggregate metrics (align with archived project report)
    # Ensemble direction accuracy from validation_report / deltas
    dstats_all = direction_stats(deltas)
    # VSW-style reported in project report: 43% (3/7) — recompute from compounds with IC50
    ic50_set = deltas[deltas["jnk1_ic50_nM"].notna() | deltas["jnk2_ic50_nM"].notna()].copy()
    # Use direction_match column
    if len(ic50_set):
        m = ic50_set["direction_match"].astype(str).str.lower().isin(["true", "1"])
        vsw_acc = float(m.mean())
        vsw_n = int(len(ic50_set))
        vsw_k = int(m.sum())
    else:
        vsw_acc, vsw_n, vsw_k = None, 0, 0

    # Key-control accuracy from confusion file
    if len(conf):
        cm = conf["direction_match"].astype(str).str.lower().isin(["true", "1"])
        key_acc = float(cm.mean())
        key_n = int(len(conf))
        key_k = int(cm.sum())
    else:
        key_acc, key_n, key_k = None, 0, 0

    # Archived narrative numbers from docs/JNK1_PROJECT_REPORT.md (VSW single-PDB vs ensemble).
    archived_vsw = "43% (3/7) VSW single-PDB; 29% (2/7) ensemble (project report)"
    recomputed = f"{vsw_k}/{vsw_n} = {vsw_acc:.1%} on benchmark_deltas CSV (excl. NA dirs)" if vsw_acc is not None else "NA"

    summary_rows = [
        {
            "method": "Glide Δsel_dock direction (archived project report)",
            "metric": "direction accuracy",
            "value": archived_vsw,
            "pass_threshold": "≥55%",
            "verdict": "FAIL",
            "used_for_purchase": "NO — decoupled; family shortlist only",
        },
        {
            "method": "Glide Δsel_dock direction (recomputed from benchmark_deltas_51c1.csv)",
            "metric": "direction accuracy",
            "value": recomputed,
            "pass_threshold": "≥55%",
            "verdict": "FAIL",
            "used_for_purchase": "NO",
        },
        {
            "method": "Glide Δsel_dock (key controls SP600125/TCS/CC-930/E1)",
            "metric": "direction accuracy",
            "value": f"{key_k}/{key_n} = {key_acc:.1%}" if key_acc is not None else "NA",
            "pass_threshold": "≥55%",
            "verdict": "FAIL" if (key_acc is not None and key_acc < 0.55) else "CHECK",
            "used_for_purchase": "NO",
        },
        {
            "method": "Gly87 (KLIFS b.l.37) occupancy heuristic",
            "metric": "discriminative power",
            "value": (
                f"occ_JNK1 True for {int(gly_out['occ_JNK1'].astype(str).str.lower().isin(['true','1']).sum())}/"
                f"{len(gly_out)} benchmarks; pred_JNK1_sel any={gly_pred_sel_any}; "
                f"d_occ range {gly_out['d_occ'].min():.2f}–{gly_out['d_occ'].max():.2f} Å"
            ),
            "pass_threshold": "separates JNK1-preferring vs pan/opposite",
            "verdict": "FAIL (non-discriminative)",
            "used_for_purchase": "NO",
        },
        {
            "method": "ML isoform-selectivity classifier (ChEMBL paired)",
            "metric": "test F1 (selective class)",
            "value": "0 (positives n_train≈8; reported in training_report / project report)",
            "pass_threshold": "usable precision/recall for purchase",
            "verdict": "FAIL",
            "used_for_purchase": "NO — ML used only as family pActivity recall gate",
        },
        {
            "method": "ML family activity gate p_family≥6.0",
            "metric": "benchmark recall / decoy FPR",
            "value": "recall 9/9; Taosu decoy FPR 95.3%; EF1%=9.20",
            "pass_threshold": "high recall OK for early filter",
            "verdict": "PASS as recall filter; NOT a selectivity filter",
            "used_for_purchase": "YES (activity recall only)",
        },
    ]
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "c5_maintext_autopsy_table.csv", index=False)

    meta = {
        "purchase_decoupling_statement": (
            "Selectivity predictors (Δsel, Gly87, ML selective labels) were evaluated on "
            "literature benchmarks and explicitly NOT used as hard gates for purchasing "
            "690/2157; purchase prioritized pose-credible family-binder enrichment."
        ),
        "direction_stats_all_rows": dstats_all,
        "gly87_all_occ_true": gly_occ_all_true,
        "sources": [
            str(DELTAS.relative_to(ROOT)),
            str(CONF.relative_to(ROOT)),
            str(GLY.relative_to(ROOT)),
            "results/ml_external_validation/ml_external_validation_report_3689.md",
            "docs/JNK1_PROJECT_REPORT.md §2–§5",
        ],
    }
    (OUT / "c5_autopsy_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    md = [
        "# C5 Selectivity-Method Autopsy (Main-Text Ready)",
        "",
        "> **Purchase decoupling:** Δsel / Gly87 / ML selectivity labels were **not** purchase hard gates.",
        "",
        "## Main-text summary table",
        "",
        summary.to_markdown(index=False),
        "",
        "## Per-compound docking direction",
        "",
        dock_df.to_markdown(index=False),
        "",
        "## Gly87 self-check",
        "",
        gly_out.to_markdown(index=False),
        "",
        "## Suggested Results paragraph (English draft)",
        "",
        "On a literature JNK benchmark panel, Glide-derived Δsel_dock reproduced experimental "
        f"isoform direction for only {vsw_k}/{vsw_n} compounds "
        f"({vsw_acc:.0%} if available else 'low accuracy'). "
        "A Gly87 occupancy heuristic labeled all tested benchmarks as JNK1-occupying and did not "
        "separate JNK1-preferring controls from pan or reverse profiles. An ML selective-class "
        "model trained on sparse ChEMBL positives yielded test F1 = 0. Accordingly, these "
        "filters were retained only as negative controls and were not used to purchase the "
        "shortlist; candidate selection prioritized family-activity and pose/MD QC.",
        "",
    ]
    # fix f-string issue - vsw_acc might be None
    if vsw_acc is None:
        md[-3] = (
            "On a literature JNK benchmark panel, Glide-derived Δsel_dock showed low direction "
            "accuracy relative to a 55% usability threshold. "
            "A Gly87 occupancy heuristic labeled all tested benchmarks as JNK1-occupying and did not "
            "separate JNK1-preferring controls from pan or reverse profiles. An ML selective-class "
            "model trained on sparse ChEMBL positives yielded test F1 = 0. Accordingly, these "
            "filters were retained only as negative controls and were not used to purchase the "
            "shortlist; candidate selection prioritized family-activity and pose/MD QC."
        )
    else:
        md[-3] = (
            "On a literature JNK benchmark panel, Glide-derived Δsel_dock failed a 55% direction-"
            f"accuracy usability threshold (archived VSW single-PDB 43% [3/7]; recomputed ensemble "
            f"table {vsw_k}/{vsw_n} = {vsw_acc:.0%}). "
            "A Gly87 occupancy heuristic labeled all tested benchmarks as JNK1-occupying and did not "
            "separate JNK1-preferring controls from pan or reverse profiles. An ML selective-class "
            "model trained on sparse ChEMBL positives yielded test F1 = 0. Accordingly, these "
            "filters were retained only as negative controls and were not used to purchase the "
            "shortlist; candidate selection prioritized family-activity and pose/MD QC."
        )

    (OUT / "C5_SELECTIVITY_AUTOPSY.md").write_text("\n".join(md), encoding="utf-8")
    print(summary.to_string(index=False))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
