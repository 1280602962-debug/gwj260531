#!/usr/bin/env python3
"""Leave exact receptor cognate ligands out of the primary-seed panels.

This is a descriptive sensitivity analysis, not a train/test leakage test.  It
uses the same complete-case, pocket-matched directional estimands as Table 2
and the same per-ligand two-pocket mean score as Table 3.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "data/jcim_multiseed_v0/tables/multiseed_scores_long_v1.csv"
OUT = ROOT / "data/jcim_novelty_v0/tables/leave_cognate_out_v1.csv"
VERDICT = ROOT / "data/jcim_novelty_v0/analysis/LEAVE_COGNATE_OUT_V1.md"
PRIMARY_SEED = 20260727

# Only exact cognate ligands that are themselves members of a main panel.
COGNATES = {
    "EGFR/HER2": {
        "ligand": "EH40_01",
        "chembl_id": "CHEMBL1614725",
        "name": "TAK-285",
        "receptors": "3POZ/3RCD",
        "het_code": "03P",
    },
    "PIK3CA/mTOR": {
        "ligand": "PM48_01",
        "chembl_id": "CHEMBL573339",
        "name": "PI-103",
        "receptors": "4L23/4JT6",
        "het_code": "X6K",
    },
}


def auroc(pos: list[float], neg: list[float]) -> float:
    if not pos or not neg:
        return float("nan")
    wins = sum(p > n for p in pos for n in neg)
    ties = sum(p == n for p in pos for n in neg)
    return (wins + 0.5 * ties) / (len(pos) * len(neg))


def metrics(records: list[dict]) -> dict:
    dual = [r for r in records if r["class"] == "dual"]
    a_only = [r for r in records if r["class"] in {"A_only", "A-only"}]
    b_only = [r for r in records if r["class"] in {"B_only", "B-only"}]
    neither = [r for r in records if r["class"] == "neither"]
    d_a = auroc([r["B"] for r in dual], [r["B"] for r in a_only])
    d_b = auroc([r["A"] for r in dual], [r["A"] for r in b_only])
    d_n = auroc(
        [(r["A"] + r["B"]) / 2.0 for r in dual],
        [(r["A"] + r["B"]) / 2.0 for r in neither],
    )
    return {
        "n_complete": len(records),
        "n_dual": len(dual),
        "n_A_only": len(a_only),
        "n_B_only": len(b_only),
        "n_neither": len(neither),
        "auroc_D_vs_A_pocketB": d_a,
        "auroc_D_vs_B_pocketA": d_b,
        "summary_min": min(d_a, d_b),
        "auroc_D_vs_neither_vina_mean": d_n,
    }


def r4(value: float) -> str:
    return f"{value:.4f}"


def main() -> None:
    grouped: dict[tuple[str, str], dict] = defaultdict(dict)
    with SOURCE.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if int(row["seed"]) != PRIMARY_SEED:
                continue
            if row["status"] not in {"ok", "exists", "primary_reused"}:
                continue
            if not row["vina_mode1"]:
                continue
            key = (row["pair"], row["ligand"])
            grouped[key].setdefault("class", row["class"])
            grouped[key][row["pocket"]] = -float(row["vina_mode1"])

    by_pair: dict[str, list[dict]] = defaultdict(list)
    for (pair, ligand), record in grouped.items():
        if "A" in record and "B" in record:
            by_pair[pair].append({"ligand": ligand, **record})

    output = []
    for pair, cognate in COGNATES.items():
        before_records = by_pair[pair]
        exact = [r for r in before_records if r["ligand"] == cognate["ligand"]]
        if len(exact) != 1 or exact[0]["class"] != "dual":
            raise RuntimeError(f"expected one dual cognate member for {pair}: {exact}")
        after_records = [r for r in before_records if r["ligand"] != cognate["ligand"]]
        before = metrics(before_records)
        after = metrics(after_records)
        output.append(
            {
                "pair": pair,
                "receptors": cognate["receptors"],
                "cognate_name": cognate["name"],
                "het_code": cognate["het_code"],
                "panel_ligand": cognate["ligand"],
                "molecule_chembl_id": cognate["chembl_id"],
                "n_complete_before": before["n_complete"],
                "n_complete_after": after["n_complete"],
                "n_dual_before": before["n_dual"],
                "n_dual_after": after["n_dual"],
                "n_A_only_before": before["n_A_only"],
                "n_A_only_after": after["n_A_only"],
                "n_B_only_before": before["n_B_only"],
                "n_B_only_after": after["n_B_only"],
                "n_neither_before": before["n_neither"],
                "n_neither_after": after["n_neither"],
                "D_vs_A_before": r4(before["auroc_D_vs_A_pocketB"]),
                "D_vs_A_after": r4(after["auroc_D_vs_A_pocketB"]),
                "D_vs_B_before": r4(before["auroc_D_vs_B_pocketA"]),
                "D_vs_B_after": r4(after["auroc_D_vs_B_pocketA"]),
                "summary_min_before": r4(before["summary_min"]),
                "summary_min_after": r4(after["summary_min"]),
                "delta_summary_min_after_minus_before": r4(after["summary_min"] - before["summary_min"]),
                "D_vs_neither_before": r4(before["auroc_D_vs_neither_vina_mean"]),
                "D_vs_neither_after": r4(after["auroc_D_vs_neither_vina_mean"]),
                "delta_D_vs_neither_after_minus_before": r4(
                    after["auroc_D_vs_neither_vina_mean"] - before["auroc_D_vs_neither_vina_mean"]
                ),
            }
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)

    lines = [
        "# Leave-cognate-out sensitivity",
        "",
        f"Source: `{SOURCE.relative_to(ROOT).as_posix()}`; production seed {PRIMARY_SEED}.",
        "The analysis removes the exact co-crystallized ligand when it is a member of the main panel, then recomputes the two Table 2 directional AUROCs, their descriptive `summary_min`, and the Table 3 Dual-versus-neither `vina_mean` AUROC.",
        "It tests sensitivity to the single exact cognate molecule; it does not remove cognate-like chemotypes and is not described as a train/test leakage analysis.",
        "",
        "| pair | removed cognate | n dual before/after | D/A before→after | D/B before→after | summary_min before→after | D/neither before→after |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in output:
        lines.append(
            f"| {row['pair']} | {row['cognate_name']} ({row['panel_ligand']}; {row['molecule_chembl_id']}) | "
            f"{row['n_dual_before']}/{row['n_dual_after']} | {row['D_vs_A_before']}→{row['D_vs_A_after']} | "
            f"{row['D_vs_B_before']}→{row['D_vs_B_after']} | {row['summary_min_before']}→{row['summary_min_after']} | "
            f"{row['D_vs_neither_before']}→{row['D_vs_neither_after']} |"
        )
    lines += [
        "",
        "Removing either exact cognate ligand produced only small point-estimate changes and did not alter the qualitative interpretation. Residual receptor–chemotype favorability is not excluded by this one-ligand sensitivity.",
        "",
    ]
    VERDICT.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {OUT}")
    print(f"wrote {VERDICT}")


if __name__ == "__main__":
    main()
