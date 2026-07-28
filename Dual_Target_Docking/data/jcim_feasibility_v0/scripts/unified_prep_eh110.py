#!/usr/bin/env python3
"""Assemble a fully uniform-prep EGFR/HER2 panel120 and recompute M1 metrics.

panel120 mixed two ligand preparations: the 40 legacy ligands reused LigPrep
poses while the 70 new ligands used RDKit ETKDG + meeko. M4-min re-docked the
legacy 40 under RDKit + meeko, so an all-RDKit EH110 can be assembled from
existing score tables with **no new docking**.

RTM z-scores are recomputed on the assembled panel, since the stored
`rtm_min_z` was standardised within its original (mixed) panel.
"""

from __future__ import annotations

import csv
import statistics as st
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "tables"

PANEL = ROOT / "data" / "egfr_her2_panel120_v0" / "tables" / "panel_v0_120.csv"
MIXED = ROOT / "data" / "egfr_her2_panel120_v0" / "tables" / "ablation_ligand_scores.csv"
REPREP = (
    ROOT / "data" / "egfr_her2_panel40_reprep_rdkit_v0" / "tables" / "ablation_ligand_scores.csv"
)

ARMS = ["vina_mean", "vina_min", "rtm_min", "rtm_mean", "rtm_min_z"]
BASELINES = ["heavy_atoms", "MW", "cLogP", "TPSA"]


def auroc(pos: list[float], neg: list[float]) -> float:
    if not pos or not neg:
        return float("nan")
    total = 0
    wins = 0.0
    for p in pos:
        for q in neg:
            total += 1
            if p > q:
                wins += 1.0
            elif p == q:
                wins += 0.5
    return wins / total


def build() -> list[dict]:
    panel = {r["panel_id"]: r for r in csv.DictReader(PANEL.open())}
    mixed = list(csv.DictReader(MIXED.open()))
    reprep = {r["ligand"]: r for r in csv.DictReader(REPREP.open())}

    rows = []
    for row in mixed:
        src = reprep[row["ligand"]] if row["from_panel40"] == "yes" else row
        mol = Chem.MolFromSmiles(panel[row["ligand"]]["smiles"])
        rows.append(
            {
                "ligand": row["ligand"],
                "class": row["class"],
                "prep": "rdkit_meeko",
                "vina_mean": float(src["vina_mean"]),
                "vina_min": float(src["vina_min"]),
                "rtm_mean": float(src["rtm_mean"]),
                "rtm_min": float(src["rtm_min"]),
                "rtm_3POZ": float(src["rtm_3POZ"]),
                "rtm_3RCD": float(src["rtm_3RCD"]),
                "heavy_atoms": float(mol.GetNumHeavyAtoms()),
                "MW": Descriptors.MolWt(mol),
                "cLogP": Descriptors.MolLogP(mol),
                "TPSA": Descriptors.TPSA(mol),
            }
        )

    for end in ("rtm_3POZ", "rtm_3RCD"):
        vals = [r[end] for r in rows]
        mu, sd = st.mean(vals), st.pstdev(vals)
        for r in rows:
            r[f"{end}_z"] = (r[end] - mu) / sd
    for r in rows:
        r["rtm_min_z"] = min(r["rtm_3POZ_z"], r["rtm_3RCD_z"])
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = build()

    scores = OUT / "eh110_unified_prep_scores.csv"
    fields = [k for k in rows[0] if not k.startswith("_")]
    with scores.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    by_class = {c: [r for r in rows if r["class"] == c] for c in ("dual", "A_only", "B_only")}
    results = []
    for arm in ARMS + BASELINES:
        get = lambda rs: [r[arm] for r in rs]  # noqa: E731
        d, a, b = get(by_class["dual"]), get(by_class["A_only"]), get(by_class["B_only"])
        da, db = auroc(d, a), auroc(d, b)
        results.append(
            {
                "pair": "EGFR_HER2",
                "subset": "all_110_unified_rdkit",
                "arm": arm,
                "family": "docking" if arm in ARMS else "baseline",
                "auroc_D_vs_A": round(da, 4),
                "auroc_D_vs_B": round(db, 4),
                "summary_min": round(min(da, db), 4),
            }
        )

    best_baseline = max(
        (r for r in results if r["family"] == "baseline"), key=lambda r: r["summary_min"]
    )
    for r in results:
        if r["family"] == "docking":
            r["fail_baseline"] = r["summary_min"] <= best_baseline["summary_min"]
        else:
            r["fail_baseline"] = ""
        r["best_baseline_arm"] = best_baseline["arm"]
        r["best_baseline_summary_min"] = best_baseline["summary_min"]

    metrics = OUT / "eh110_unified_prep_directional.csv"
    with metrics.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)

    print(f"wrote {scores}\nwrote {metrics}\n")
    print(f"n = {len(rows)} (all RDKit ETKDG + meeko; no new docking)")
    print(f"{'arm':14s} {'family':>9s} {'D/A':>7s} {'D/B':>7s} {'min':>7s}  fail_baseline")
    for r in results:
        print(
            f"{r['arm']:14s} {r['family']:>9s} {r['auroc_D_vs_A']:7.3f} "
            f"{r['auroc_D_vs_B']:7.3f} {r['summary_min']:7.3f}  {r['fail_baseline']}"
        )
    print(
        f"\nbest non-docking baseline: {best_baseline['arm']} "
        f"(summary_min = {best_baseline['summary_min']})"
    )


if __name__ == "__main__":
    main()
