#!/usr/bin/env python3
"""C7 — PAINS / Brenk / physicochemical risk re-check for purchase set."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
try:
    from rdkit import RDConfig
    RDDataDir = RDConfig.RDDataDir
except Exception:
    RDDataDir = ""

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "purchase_risk"
OUT.mkdir(parents=True, exist_ok=True)

MOLECULES = {
    "690": "Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1",
    "2157": "Cc1cnc(NCc2cccc3c2OCCCO3)nc1C",
    "2231": "COc1nc(NCc2ccccc2CN2CCCC2=O)ncc1F",  # not purchased; opportunity-cost comparator
    "E1": "COc1ccc(Nc2ncc(Br)c(Nc3ccc(OCCN(C)C)cc3C(N)=O)n2)cc1",
    "CC-90001": "C[C@@H]1CC[C@@H](Nc2nc(NC(C)(C)C)ncc2C(N)=O)C[C@H]1O",
}


def build_catalog(flag) -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(flag)
    return FilterCatalog(params)


def hits(catalog: FilterCatalog, mol) -> list[str]:
    return [e.GetDescription() for e in catalog.GetMatches(mol)]


def main():
    pains = build_catalog(FilterCatalogParams.FilterCatalogs.PAINS)
    # Brenk / NIH / ZINC often available under FilterCatalogs
    extra_flags = []
    for name in ["BRENK", "NIH", "ZINC"]:
        if hasattr(FilterCatalogParams.FilterCatalogs, name):
            extra_flags.append((name, build_catalog(getattr(FilterCatalogParams.FilterCatalogs, name))))

    rows = []
    for cid, smi in MOLECULES.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(cid)
        pains_h = hits(pains, mol)
        row = {
            "compound_id": cid,
            "role": "purchased" if cid in {"690", "2157"} else ("comparator_unbought" if cid == "2231" else "positive_control"),
            "smiles": smi,
            "MW": round(Descriptors.MolWt(mol), 2),
            "logP": round(Descriptors.MolLogP(mol), 2),
            "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol),
            "TPSA": round(Descriptors.TPSA(mol), 2),
            "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
            "qed": round(QED.qed(mol), 3),
            "lipinski_violations": int(
                (Descriptors.MolWt(mol) > 500)
                + (Descriptors.MolLogP(mol) > 5)
                + (Lipinski.NumHDonors(mol) > 5)
                + (Lipinski.NumHAcceptors(mol) > 10)
            ),
            "PAINS_hit": bool(pains_h),
            "PAINS_families": ";".join(pains_h) if pains_h else "",
        }
        for name, cat in extra_flags:
            h = hits(cat, mol)
            row[f"{name}_hit"] = bool(h)
            row[f"{name}_families"] = ";".join(h[:8]) if h else ""
        # simple aggregator-ish flags (heuristic, not definitive)
        row["risk_notes"] = []
        if row["logP"] > 5:
            row["risk_notes"].append("high_logP")
        if row["TPSA"] < 20:
            row["risk_notes"].append("very_low_TPSA")
        if row["PAINS_hit"]:
            row["risk_notes"].append("PAINS")
        if row.get("BRENK_hit"):
            row["risk_notes"].append("BRENK")
        row["risk_notes"] = ";".join(row["risk_notes"]) if row["risk_notes"] else "none_flagged"
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "c7_purchase_risk_filters.csv", index=False)

    md = [
        "# C7 Purchase Assay-Risk Filters",
        "",
        "RDKit FilterCatalog PAINS (+ Brenk/NIH/ZINC when available).",
        "These are **alerts**, not proof of assay artifact.",
        "",
        df.drop(columns=[c for c in df.columns if c.endswith("_families") and c != "PAINS_families"]).to_markdown(index=False),
        "",
        "## Detail: PAINS / Brenk families",
        "",
    ]
    for _, r in df.iterrows():
        md.append(f"- **{r['compound_id']}**: PAINS=`{r['PAINS_families'] or 'none'}`; notes=`{r['risk_notes']}`")
    md.append("")
    (OUT / "C7_PURCHASE_RISK_REPORT.md").write_text("\n".join(md), encoding="utf-8")
    (OUT / "c7_meta.json").write_text(
        json.dumps({"n": len(df), "rdkit_data": RDDataDir, "out": str(OUT)}, indent=2),
        encoding="utf-8",
    )
    print(df[["compound_id", "role", "PAINS_hit", "qed", "lipinski_violations", "risk_notes"]].to_string(index=False))


if __name__ == "__main__":
    main()
