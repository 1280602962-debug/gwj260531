#!/usr/bin/env python3
"""Finalize BindingDB external-slice feasibility as SI flow tables (no new docking)."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data/jcim_novelty_v0/tables"
AN = ROOT / "data/jcim_novelty_v0/analysis"
SI = ROOT / "docs"
AN.mkdir(parents=True, exist_ok=True)

FLOW_IN = TAB / "external_candidate_flow.csv"
SUM_IN = TAB / "external_slice_summary_v1.csv"


def main():
    flow = list(csv.DictReader(FLOW_IN.open(encoding="utf-8")))
    summary = list(csv.DictReader(SUM_IN.open(encoding="utf-8")))

    # Wide SI flow table: one row per pair with funnel stages
    stages = [
        ("native_paired_theta6", "raw / θ=6.0 paired BindingDB"),
        ("drop_shared_literature", "after literature exclusion"),
        ("drop_shared_structure", "after structure exclusion"),
        ("drop_neighbors_lt_0.70", "after ECFP4 < 0.70 exclusion"),
    ]
    by_pair = {}
    for r in flow:
        by_pair.setdefault(r["pair"], {})[r["layer"]] = r

    wide_rows = []
    for pair, layers in by_pair.items():
        row = {"pair": pair}
        for key, _label in stages:
            L = layers.get(key, {})
            row[f"{key}_n"] = L.get("n_ligands", "")
            row[f"{key}_dual"] = L.get("n_dual", "")
            row[f"{key}_A"] = L.get("n_A_only", "")
            row[f"{key}_B"] = L.get("n_B_only", "")
            row[f"{key}_neither"] = L.get("n_neither", "")
        s = next((x for x in summary if x["pair"] == pair), {})
        row["gate"] = s.get("gate", "")
        row["packaged_as_external_evaluation"] = s.get("packaged_as_external_evaluation", "0")
        row["n_sources_dual"] = s.get("n_sources_dual", "")
        row["n_sources_A_only"] = s.get("n_sources_A_only", "")
        row["n_sources_B_only"] = s.get("n_sources_B_only", "")
        row["top_doc_frac_dual"] = s.get("top_doc_frac_dual", "")
        row["note"] = s.get("note", "")
        wide_rows.append(row)

    out_wide = TAB / "bindingdb_external_feasibility_flow_v1.csv"
    with out_wide.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(wide_rows[0].keys()))
        w.writeheader()
        w.writerows(wide_rows)

    # Long SI-friendly funnel
    long_rows = []
    for pair, layers in by_pair.items():
        for i, (key, label) in enumerate(stages, start=1):
            L = layers.get(key, {})
            long_rows.append(
                {
                    "pair": pair,
                    "step": i,
                    "stage_key": key,
                    "stage_label": label,
                    "n_ligands": L.get("n_ligands", ""),
                    "n_dual": L.get("n_dual", ""),
                    "n_A_only": L.get("n_A_only", ""),
                    "n_B_only": L.get("n_B_only", ""),
                    "n_neither": L.get("n_neither", ""),
                }
            )
        long_rows.append(
            {
                "pair": pair,
                "step": 5,
                "stage_key": "eligible_pairs",
                "stage_label": "pairs meeting primary external gate",
                "n_ligands": "0",
                "n_dual": "0",
                "n_A_only": "0",
                "n_B_only": "0",
                "n_neither": "0",
            }
        )
    out_long = TAB / "bindingdb_external_feasibility_funnel_long_v1.csv"
    with out_long.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(long_rows[0].keys()))
        w.writeheader()
        w.writerows(long_rows)

    md = [
        "# BindingDB external-slice feasibility (final SI packaging)",
        "",
        "**Status:** completed negative feasibility result. **No docking.**",
        "",
        "Pre-frozen filters (literature → structure → ECFP4 < 0.70 → class/source gates)",
        "yielded **0** target pairs eligible for primary external evaluation.",
        "`packaged_as_external_evaluation = 0` for all audited pairs.",
        "",
        "## Funnel (machine-readable)",
        "",
        f"- Wide: `{out_wide.relative_to(ROOT)}`",
        f"- Long: `{out_long.relative_to(ROOT)}`",
        f"- Per-layer source: `{FLOW_IN.relative_to(ROOT)}`",
        f"- Gate summary: `{SUM_IN.relative_to(ROOT)}`",
        "",
        "## Condensed flow",
        "",
        "```",
        "raw BindingDB (θ=6.0 paired)",
        "        ↓",
        "literature exclusion",
        "        ↓",
        "structure exclusion",
        "        ↓",
        "ECFP4 max-sim < 0.70",
        "        ↓",
        "class / source concentration gates",
        "        ↓",
        "0 eligible pairs  →  no docking",
        "```",
        "",
        "## Claim ceiling",
        "",
        "- Allowed: strict database-external bidirectional hard negatives are scarce under this freeze.",
        "- Forbidden: relaxing thresholds after seeing counts; calling remaining ligands an external set;",
        "  presenting BindingDB as external validation.",
        "",
        "| pair | after ECFP n | dual/A/B/neither | gate | docked |",
        "|---|---:|---|---|---|",
    ]
    for s in summary:
        md.append(
            f"| {s['pair']} | {s.get('after_ecfp_lt_0.70','')} | "
            f"{s.get('n_dual','')}/{s.get('n_A_only','')}/{s.get('n_B_only','')}/{s.get('n_neither','')} | "
            f"{s.get('gate','')} | no |"
        )
    (AN / "BINDINGDB_EXTERNAL_FEASIBILITY_SI_V1.md").write_text("\n".join(md) + "\n")

    # Patch EN SI outline section if present
    si_path = SI / "SUPPORTING_INFORMATION_JCIM_EN_V1.md"
    if si_path.exists():
        text = si_path.read_text(encoding="utf-8")
        block = (
            "\n## S9. BindingDB supply freeze (negative feasibility)\n\n"
            "| SI table | Content | Role |\n"
            "|---------|---------|------|\n"
            "| S48–S49 / flow | `bindingdb_external_feasibility_flow_v1.csv`, funnel long table | Prespecified supply freeze |\n"
            "| Verdict | `BINDINGDB_EXTERNAL_FEASIBILITY_SI_V1.md` | 0 eligible pairs; no docking |\n\n"
            "Do not interpret remaining ligands as an external evaluation set.\n"
        )
        if "S9. BindingDB supply freeze" not in text:
            # replace or append
            if "## S9." in text:
                pass
            else:
                text = text.rstrip() + "\n" + block
                si_path.write_text(text + "\n", encoding="utf-8")
    print("wrote", out_wide)
    print("wrote", out_long)
    print("wrote", AN / "BINDINGDB_EXTERNAL_FEASIBILITY_SI_V1.md")


if __name__ == "__main__":
    main()
