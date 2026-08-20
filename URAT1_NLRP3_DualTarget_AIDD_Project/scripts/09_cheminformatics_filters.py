#!/usr/bin/env python3
"""
Non-docking computational module B — medicinal-chemistry structural alerts.

Annotates the EXISTING docking/Pareto data with PAINS (A/B/C), Brenk and NIH
structural-alert flags plus a transparent colloidal-aggregation heuristic.
This is a pure downstream annotation: it never re-docks and never changes Pareto membership.
Run only on the production P2 merge from `results/repurposing/` after `run_funnel_p2.sh`.

Inputs (read-only):
  results/repurposing/pareto_shortlist.csv
  results/repurposing/pareto_merged_scores.csv

Outputs:
  results/cheminformatics/filters_shortlist.csv
  results/cheminformatics/filters_pool.csv
  results/cheminformatics/filters_summary.json

Usage:
  python3 scripts/09_cheminformatics_filters.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PARETO_DIR = PROJECT_ROOT / "results" / "repurposing"
OUT_DIR = PROJECT_ROOT / "results" / "cheminformatics"

SMILES_COL = "canonical_smiles"
NAME_COL = "name"


def build_catalogs():
    """One FilterCatalog per alert family so we can report which family fires."""
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

    families = {
        "PAINS_A": FilterCatalogParams.FilterCatalogs.PAINS_A,
        "PAINS_B": FilterCatalogParams.FilterCatalogs.PAINS_B,
        "PAINS_C": FilterCatalogParams.FilterCatalogs.PAINS_C,
        "BRENK": FilterCatalogParams.FilterCatalogs.BRENK,
        "NIH": FilterCatalogParams.FilterCatalogs.NIH,
    }
    catalogs = {}
    for label, cat in families.items():
        params = FilterCatalogParams()
        params.AddCatalog(cat)
        catalogs[label] = FilterCatalog(params)
    return catalogs


def aggregation_heuristic(mol) -> bool:
    """Transparent, literature-inspired colloidal-aggregation risk flag.

    Not a substitute for Aggregator Advisor; flags the classic high-logP /
    poly-aromatic / low-polar-fraction profile that correlates with promiscuous
    aggregation in biochemical assays. Reported as a heuristic, not ground truth.
    """
    from rdkit.Chem import Descriptors

    logp = Descriptors.MolLogP(mol)
    arom = Descriptors.NumAromaticRings(mol)
    tpsa = Descriptors.TPSA(mol)
    return bool(logp >= 3.5 and arom >= 3 and tpsa < 75)


def annotate(df: pd.DataFrame, catalogs) -> pd.DataFrame:
    from rdkit import Chem

    records = []
    for _, row in df.iterrows():
        smi = row.get(SMILES_COL)
        rec = {c: row.get(c) for c in df.columns if c in (
            "repurposing_id", "chembl_id", NAME_COL, SMILES_COL,
            "max_phase", "s_u_percentile", "s_n_percentile", "pareto_front",
        )}
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if mol is None:
            rec.update({k: None for k in (
                "pains_any", "pains_families", "brenk", "nih",
                "aggregation_risk_heuristic", "alert_descriptions",
            )})
            rec["parse_ok"] = False
            records.append(rec)
            continue

        pains_fam = []
        descriptions = []
        for label in ("PAINS_A", "PAINS_B", "PAINS_C"):
            m = catalogs[label].GetFirstMatch(mol)
            if m is not None:
                pains_fam.append(label)
                descriptions.append(f"{label}:{m.GetDescription()}")
        brenk_m = catalogs["BRENK"].GetFirstMatch(mol)
        nih_m = catalogs["NIH"].GetFirstMatch(mol)
        if brenk_m is not None:
            descriptions.append(f"BRENK:{brenk_m.GetDescription()}")
        if nih_m is not None:
            descriptions.append(f"NIH:{nih_m.GetDescription()}")

        rec.update({
            "parse_ok": True,
            "pains_any": len(pains_fam) > 0,
            "pains_families": ";".join(pains_fam) if pains_fam else "",
            "brenk": brenk_m is not None,
            "nih": nih_m is not None,
            "aggregation_risk_heuristic": aggregation_heuristic(mol),
            "alert_descriptions": " | ".join(descriptions) if descriptions else "",
        })
        records.append(rec)
    return pd.DataFrame(records)


def summarize(name: str, df: pd.DataFrame) -> dict:
    n = int(len(df))
    parsed = df[df["parse_ok"] == True]  # noqa: E712
    return {
        "set": name,
        "n": n,
        "n_parsed": int(len(parsed)),
        "n_pains": int(parsed["pains_any"].sum()),
        "n_brenk": int(parsed["brenk"].sum()),
        "n_nih": int(parsed["nih"].sum()),
        "n_aggregation_heuristic": int(parsed["aggregation_risk_heuristic"].sum()),
        "n_any_alert": int(
            (parsed["pains_any"] | parsed["brenk"] | parsed["nih"]).sum()
        ),
        "pct_any_alert": round(
            100.0 * (parsed["pains_any"] | parsed["brenk"] | parsed["nih"]).mean(), 2
        ) if len(parsed) else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Medicinal-chemistry structural alerts (non-docking module B)")
    parser.add_argument("--shortlist", type=Path, default=PARETO_DIR / "pareto_shortlist.csv")
    parser.add_argument("--pool", type=Path, default=PARETO_DIR / "pareto_merged_scores.csv")
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    catalogs = build_catalogs()
    summaries = []

    shortlist = pd.read_csv(args.shortlist)
    ann_short = annotate(shortlist, catalogs)
    ann_short.to_csv(args.output_dir / "filters_shortlist.csv", index=False)
    summaries.append(summarize("pareto_shortlist", ann_short))
    print(f"Shortlist annotated: {len(ann_short)} rows -> filters_shortlist.csv")

    pool = pd.read_csv(args.pool)
    ann_pool = annotate(pool, catalogs)
    ann_pool.to_csv(args.output_dir / "filters_pool.csv", index=False)
    summaries.append(summarize("dual_docked_pool", ann_pool))
    print(f"Pool annotated: {len(ann_pool)} rows -> filters_pool.csv")

    summary = {"module": "B_cheminformatics_filters", "sets": summaries}
    with open(args.output_dir / "filters_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== Structural-alert summary ===")
    for s in summaries:
        print(
            f"  {s['set']}: n={s['n']} parsed={s['n_parsed']} "
            f"PAINS={s['n_pains']} Brenk={s['n_brenk']} NIH={s['n_nih']} "
            f"aggregation~={s['n_aggregation_heuristic']} any={s['n_any_alert']} ({s['pct_any_alert']}%)"
        )

    print("\nShortlist detail (alerts):")
    cols = [NAME_COL, "max_phase", "pains_any", "pains_families", "brenk", "aggregation_risk_heuristic"]
    print(ann_short[[c for c in cols if c in ann_short.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
