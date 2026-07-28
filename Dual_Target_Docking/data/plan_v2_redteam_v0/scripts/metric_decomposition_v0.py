#!/usr/bin/env python3
"""Measurement audit for the DualFourClass task (no new docking).

Recomputes, from already-frozen score tables:
  1. directional AUROC decomposition: dual-vs-A_only and dual-vs-B_only
     reported separately instead of the pooled dual-vs-(A u B) readout;
  2. trivial non-docking baselines (MW, cLogP, heavy-atom count) on the
     same splits, as required negative controls;
  3. label-margin census: how many ligands of each class sit within a
     given distance of the pChEMBL activity cutoff.

Exploration-pool analysis only. Does not change any frozen protocol.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

DUAL_ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "tables"
CUTOFF = 6.0
MARGIN = 0.5


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


def load_pair(pair_dir: Path, panel_csv: str, score_csv: str, id_col: str):
    panel = {r[id_col]: r for r in csv.DictReader((pair_dir / panel_csv).open())}
    scores = list(csv.DictReader((pair_dir / score_csv).open()))
    for row in scores:
        meta = panel[row["ligand"]]
        mol = Chem.MolFromSmiles(meta["smiles"])
        row["MW"] = Descriptors.MolWt(mol)
        row["cLogP"] = Descriptors.MolLogP(mol)
        row["heavy_atoms"] = float(mol.GetNumHeavyAtoms())
        row["_meta"] = meta
    return panel, scores


def directional_table(rows, arms, pair, subset):
    out = []
    by_class = {c: [r for r in rows if r["class"] == c] for c in ("dual", "A_only", "B_only")}
    for arm in arms:
        get = lambda rs: [float(r[arm]) for r in rs]  # noqa: E731
        d = get(by_class["dual"])
        a = get(by_class["A_only"])
        b = get(by_class["B_only"])
        out.append(
            {
                "pair": pair,
                "subset": subset,
                "arm": arm,
                "n_dual": len(d),
                "n_A_only": len(a),
                "n_B_only": len(b),
                "auroc_dual_vs_A_only": round(auroc(d, a), 4),
                "auroc_dual_vs_B_only": round(auroc(d, b), 4),
                "auroc_dual_vs_pooled": round(auroc(d, a + b), 4),
            }
        )
    return out


def margin_table(panel_rows, pair, end_cols):
    counts = {}
    for row in panel_rows:
        cls = row["class"]
        vals = []
        for col in end_cols:
            try:
                vals.append(float(row[col]))
            except (TypeError, ValueError):
                continue
        near = any(abs(v - CUTOFF) <= MARGIN for v in vals)
        rec = counts.setdefault(cls, {"pair": pair, "class": cls, "n": 0, "n_within_margin": 0})
        rec["n"] += 1
        rec["n_within_margin"] += int(near)
    for rec in counts.values():
        rec["frac_within_margin"] = round(rec["n_within_margin"] / rec["n"], 3)
        rec["cutoff"] = CUTOFF
        rec["margin"] = MARGIN
    return list(counts.values())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    directional: list[dict] = []
    margins: list[dict] = []

    eh_dir = DUAL_ROOT / "data" / "egfr_her2_panel120_v0"
    eh_panel, eh_scores = load_pair(
        eh_dir, "tables/panel_v0_120.csv", "tables/ablation_ligand_scores.csv", "panel_id"
    )
    eh_arms = ["vina_mean", "vina_min", "rtm_mean", "rtm_min_z", "MW", "cLogP", "heavy_atoms"]
    directional += directional_table(eh_scores, eh_arms, "EGFR_HER2", "all_110")
    directional += directional_table(
        [r for r in eh_scores if r["from_panel40"] == "yes"], eh_arms, "EGFR_HER2", "old_40_ligprep"
    )
    directional += directional_table(
        [r for r in eh_scores if r["from_panel40"] == "no"], eh_arms, "EGFR_HER2", "new_70_rdkit"
    )
    margins += margin_table(
        eh_panel.values(), "EGFR_HER2", ["pchembl_EGFR", "pchembl_HER2"]
    )

    pm_dir = DUAL_ROOT / "data" / "pik3ca_mtor_panel48_v0"
    pm_panel, pm_scores = load_pair(
        pm_dir, "tables/panel_v0_48.csv", "tables/ablation_ligand_scores.csv", "panel_id"
    )
    pm_arms = ["vina_mean", "vina_min", "rtm_mean", "rtm_min_z", "MW", "cLogP", "heavy_atoms"]
    directional += directional_table(pm_scores, pm_arms, "PIK3CA_mTOR", "all_48")
    margins += margin_table(
        pm_panel.values(), "PIK3CA_mTOR", ["pchembl_PIK3CA", "pchembl_MTOR"]
    )

    dpath = OUT / "directional_auroc_v0.csv"
    with dpath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(directional[0].keys()))
        w.writeheader()
        w.writerows(directional)

    mpath = OUT / "label_margin_v0.csv"
    with mpath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(margins[0].keys()))
        w.writeheader()
        w.writerows(margins)

    print(json.dumps({"directional": str(dpath), "margins": str(mpath)}, indent=2))
    for row in directional:
        print(
            f"{row['pair']:12s} {row['subset']:16s} {row['arm']:12s} "
            f"D/A={row['auroc_dual_vs_A_only']:.3f} "
            f"D/B={row['auroc_dual_vs_B_only']:.3f} "
            f"pooled={row['auroc_dual_vs_pooled']:.3f}"
        )
    for row in margins:
        print(
            f"{row['pair']:12s} {row['class']:8s} n={row['n']:3d} "
            f"within±{MARGIN} of {CUTOFF}: {row['n_within_margin']:3d} "
            f"({row['frac_within_margin']:.2f})"
        )


if __name__ == "__main__":
    main()
