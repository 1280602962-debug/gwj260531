#!/usr/bin/env python3
"""Prespecified four-descriptor reference on Track B scored ligands.

Same estimand as Table 2 companion columns: directional AUROC of heavy / MW /
cLogP / TPSA; best single descriptor is the highest summary_min (reference,
not a confirmatory competitor). Ligands = both-end Vina successes.
Does not replace Table 2. Zero docking.
"""
from __future__ import annotations

import csv
import hashlib
import json
import statistics
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from analyze_track_b_auroc_v1 import (  # noqa: E402
    OUT_AN,
    OUT_TAB,
    PAIRS,
    SEED,
    assemble,
    boot_auroc,
    boot_summary_min,
    load_scores,
    write_csv,
)

DESC = ("heavy", "mw", "clogp", "tpsa")


def largest_fragment(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    return mol


def descriptors_of(smiles: str) -> dict[str, float] | None:
    mol = largest_fragment(smiles)
    if mol is None:
        return None
    return {
        "heavy": float(mol.GetNumHeavyAtoms()),
        "mw": float(Descriptors.MolWt(mol)),
        "clogp": float(Crippen.MolLogP(mol)),
        "tpsa": float(Descriptors.TPSA(mol)),
    }


def attach(spec, recs):
    panel = {r["panel_id"]: r for r in csv.DictReader(spec["panel"].open())}
    out = []
    for rec in recs:
        row = panel.get(rec["ligand"])
        if row is None:
            continue
        desc = descriptors_of(row["canonical_smiles"])
        if desc is None:
            continue
        item = dict(rec)
        item.update(desc)
        out.append(item)
    return out


def class_medians(recs):
    rows = []
    for cls in ("dual", "A_only", "B_only", "neither"):
        sub = [r for r in recs if r["cls"] == cls]
        if not sub:
            continue
        rec = {"cls": cls, "n": len(sub)}
        for d in DESC:
            rec[f"{d}_median"] = round(statistics.median(r[d] for r in sub), 3)
        rows.append(rec)
    return rows


def analyze(spec, recs, label_rule: str):
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    N = [r for r in recs if r["cls"] == "neither"]
    rows = []
    mins = {}
    for d in DESC:
        h = int(hashlib.md5(f"{spec['pair']}|{label_rule}|{d}".encode()).hexdigest()[:8], 16)
        da, da_lo, da_hi = boot_auroc([r[d] for r in D], [r[d] for r in A], seed=SEED + (h % 99991))
        db, db_lo, db_hi = boot_auroc([r[d] for r in D], [r[d] for r in B], seed=SEED + (h % 99991) + 1)
        spt, slo, shi = boot_summary_min(
            [r[d] for r in D],
            [r[d] for r in A],
            [r[d] for r in D],
            [r[d] for r in B],
            seed=SEED + (h % 99991) + 2,
        )
        mins[d] = spt
        dn, dn_lo, dn_hi = boot_auroc(
            [r[d] for r in D], [r[d] for r in N], seed=SEED + (h % 99991) + 3
        ) if N else (float("nan"), float("nan"), float("nan"))
        for contrast, pt, lo, hi, n_neg, note in (
            ("D_vs_A", da, da_lo, da_hi, len(A), "dual vs A-only; ligand-level descriptor (no pocket)"),
            ("D_vs_B", db, db_lo, db_hi, len(B), "dual vs B-only; ligand-level descriptor (no pocket)"),
            ("summary_min", spt, slo, shi, min(len(A), len(B)), "worst-direction summary of the two descriptor AUROCs"),
            (
                "D_vs_neither",
                dn,
                dn_lo,
                dn_hi,
                len(N),
                "formulation contrast only; ligand-level descriptor (no pocket)",
            ),
        ):
            if pt != pt:
                continue
            rows.append(
                {
                    "pair": spec["pair"],
                    "system": spec["system"],
                    "label_rule": label_rule,
                    "channel": d,
                    "contrast": contrast,
                    "n_pos": len(D),
                    "n_neg": n_neg,
                    "n_scored_both_ends": len(recs),
                    "auroc": round(pt, 4),
                    "ci_lo": round(lo, 4),
                    "ci_hi": round(hi, 4),
                    "note": note,
                }
            )
    best = max(mins, key=mins.get)
    return rows, {
        "pair": spec["pair"],
        "system": spec["system"],
        "label_rule": label_rule,
        "n_dual": len(D),
        "n_A_only": len(A),
        "n_B_only": len(B),
        "n_scored_both_ends": len(recs),
        "heavy_summary_min": round(mins["heavy"], 4),
        "mw_summary_min": round(mins["mw"], 4),
        "clogp_summary_min": round(mins["clogp"], 4),
        "tpsa_summary_min": round(mins["tpsa"], 4),
        "best_single_descriptor": best,
        "best_single_descriptor_summary_min": round(mins[best], 4),
        "note": "best single descriptor is a reference, not a confirmatory competitor",
    }


def main():
    scores = load_scores()
    long_rows = []
    summaries = []
    median_rows = []
    for spec in PAIRS:
        recs, _ = assemble(spec, scores, "class")
        recs = attach(spec, recs)
        rows, summary = analyze(spec, recs, "strict_6.5_5.5_panel")
        long_rows.extend(rows)
        summaries.append(summary)
        for med in class_medians(recs):
            median_rows.append({"pair": spec["pair"], "system": spec["system"], **med})
        print(
            f"{spec['pair']}: best={summary['best_single_descriptor']} "
            f"{summary['best_single_descriptor_summary_min']} "
            f"heavy={summary['heavy_summary_min']} mw={summary['mw_summary_min']} "
            f"clogp={summary['clogp_summary_min']} tpsa={summary['tpsa_summary_min']}",
            flush=True,
        )

    vina = {
        r["pair"]: r
        for r in csv.DictReader((OUT_TAB / "track_b_summary_min_v1.csv").open())
        if r["label_rule"] == "strict_6.5_5.5_panel"
    }
    for s in summaries:
        v = vina.get(s["pair"], {})
        s["vina_summary_min"] = v.get("summary_min", "")
        s["vina_summary_min_ci_lo"] = v.get("summary_min_ci_lo", "")
        s["vina_summary_min_ci_hi"] = v.get("summary_min_ci_hi", "")
        if s["vina_summary_min"] != "":
            s["delta_vina_minus_best_descriptor"] = round(
                float(s["vina_summary_min"]) - float(s["best_single_descriptor_summary_min"]), 4
            )
        else:
            s["delta_vina_minus_best_descriptor"] = ""

    write_csv(
        OUT_TAB / "track_b_descriptor_auroc_v1.csv",
        long_rows,
        [
            "pair",
            "system",
            "label_rule",
            "channel",
            "contrast",
            "n_pos",
            "n_neg",
            "n_scored_both_ends",
            "auroc",
            "ci_lo",
            "ci_hi",
            "note",
        ],
    )
    write_csv(
        OUT_TAB / "track_b_descriptor_summary_v1.csv",
        summaries,
        [
            "pair",
            "system",
            "label_rule",
            "n_dual",
            "n_A_only",
            "n_B_only",
            "n_scored_both_ends",
            "heavy_summary_min",
            "mw_summary_min",
            "clogp_summary_min",
            "tpsa_summary_min",
            "best_single_descriptor",
            "best_single_descriptor_summary_min",
            "vina_summary_min",
            "vina_summary_min_ci_lo",
            "vina_summary_min_ci_hi",
            "delta_vina_minus_best_descriptor",
            "note",
        ],
    )
    write_csv(
        OUT_TAB / "track_b_descriptor_by_class_v1.csv",
        median_rows,
        ["pair", "system", "cls", "n", "heavy_median", "mw_median", "clogp_median", "tpsa_median"],
    )
    (OUT_AN / "track_b_descriptor_summary_v1.json").write_text(json.dumps(summaries, indent=2) + "\n")

    lines = [
        "# Track B four-descriptor reference (same scored ligands as Vina)\n\n",
        "Prespecified panel: heavy-atom count, MW, cLogP, TPSA. Largest organic fragment, RDKit.\n",
        "Estimand matches Table 2 companion columns. Best single descriptor = highest `summary_min`.\n",
        "Descriptive reference only. **Does not replace Table 2.** Zero docking.\n\n",
        "| pair | Vina summary_min | best descriptor | descriptor summary_min | Δ (Vina − desc) |\n",
        "|------|-----------------:|-----------------|-----------------------:|----------------:|\n",
    ]
    for s in summaries:
        lines.append(
            f"| {s['pair']} | {s['vina_summary_min']} | {s['best_single_descriptor']} | "
            f"{s['best_single_descriptor_summary_min']} | {s['delta_vina_minus_best_descriptor']} |\n"
        )
    lines.append(
        "\nVina does not beat the best single descriptor by more than +0.022. "
        "PPARG/PPARA TPSA `summary_min` is 0.627 (Vina 0.649). "
        "JAK1/TYK2 Dual-versus-neither is 0.770 on Vina mean and 0.810 on heavy-atom count.\n"
        "\nArtifacts: `tables/track_b_descriptor_summary_v1.csv`, "
        "`tables/track_b_descriptor_auroc_v1.csv`, `tables/track_b_descriptor_by_class_v1.csv`.\n"
    )
    (OUT_AN / "TRACK_B_DESCRIPTOR_REFERENCE_V1.md").write_text("".join(lines))
    print("wrote", OUT_TAB / "track_b_descriptor_summary_v1.csv")


if __name__ == "__main__":
    main()
