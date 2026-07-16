#!/usr/bin/env python3
"""C1 — Chemotype novelty audit for purchased candidates 690 & 2157.

Compares ECFP4 (Morgan r=2, 2048-bit) Tanimoto and Bemis–Murcko scaffolds
against literature benchmark JNK ligands and curated ChEMBL JNK sets.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "chemotype_novelty"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PURCHASED = {
    "690": "Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1",
    "2231": "COc1nc(NCc2ccccc2CN2CCCC2=O)ncc1F",
}

# Canonical literature references (same SMILES as project benchmark table)
REFS = {
    "SP600125": "O=C1c2ccccc2-c2n[nH]c3cccc1c23",
    "CC-90001": "C[C@@H]1CC[C@@H](Nc2nc(NC(C)(C)C)ncc2C(N)=O)C[C@H]1O",
    "CC-930": "O[C@H]1CC[C@H](Nc2ncc3nc(Nc4c(F)cc(F)cc4F)n([C@H]4CCOC4)c3n2)CC1",
    "E1": "COc1ccc(Nc2ncc(Br)c(Nc3ccc(OCCN(C)C)cc3C(N)=O)n2)cc1",
    "Q63": "CC(C)(C)CNc1ncc(Br)c(Nc2ccccc2C(N)=O)n1",
    "TCS_JNK_6O": "CCOc1nc(NC(=O)Cc2cc(OC)ccc2OC)cc(N)c1C#N",
    "AS602801": "N#CC(c1ccnc(OCc2ccc(CN3CCOCC3)cc2)n1)c1nc2ccccc2s1",
    "JNK-IN-8": "Cc1cc(NC(=O)c2cccc(NC(=O)/C=C/CN(C)C)c2)ccc1Nc1nccc(-c2cccnc2)n1",
}


def mol_from_smiles(smi: str):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"Bad SMILES: {smi}")
    return m


def fp(m):
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)


def murcko_smiles(m) -> str:
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(m)
        return Chem.MolToSmiles(scaf) if scaf is not None else ""
    except Exception:
        return ""


def max_tc_against(query_fp, library_fps, library_ids):
    best_tc, best_id = -1.0, None
    for lid, lfp in zip(library_ids, library_fps):
        tc = DataStructs.TanimotoSimilarity(query_fp, lfp)
        if tc > best_tc:
            best_tc, best_id = tc, lid
    return best_tc, best_id


def load_chembl_pool() -> pd.DataFrame:
    frames = []
    for iso, path in [
        ("JNK1", ROOT / "data/processed/jnk1_curated.csv"),
        ("JNK2", ROOT / "data/processed/jnk2_curated.csv"),
        ("JNK3", ROOT / "data/processed/jnk3_curated.csv"),
    ]:
        df = pd.read_csv(path)
        df = df.rename(columns={"canonical_smiles": "smiles", "molecule_chembl_id": "chembl_id"})
        df["source_isoform"] = iso
        frames.append(df[["chembl_id", "smiles", "source_isoform"]])
    paired = pd.read_csv(ROOT / "data/processed/paired_set.csv")
    paired = paired.rename(columns={"canonical_smiles": "smiles", "molecule_chembl_id": "chembl_id"})
    paired["source_isoform"] = "paired"
    frames.append(paired[["chembl_id", "smiles", "source_isoform"]])
    pool = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["chembl_id", "smiles"])
    mols, keep = [], []
    for i, row in pool.iterrows():
        m = Chem.MolFromSmiles(str(row["smiles"]))
        if m is None:
            continue
        mols.append(m)
        keep.append(i)
    pool = pool.loc[keep].reset_index(drop=True)
    pool["mol"] = mols
    pool["fp"] = [fp(m) for m in mols]
    pool["murcko"] = [murcko_smiles(m) for m in mols]
    return pool


def main():
    ref_mols = {k: mol_from_smiles(v) for k, v in REFS.items()}
    ref_fps = {k: fp(m) for k, m in ref_mols.items()}
    ref_murcko = {k: murcko_smiles(m) for k, m in ref_mols.items()}

    pool = load_chembl_pool()
    pool_fps = list(pool["fp"])
    pool_ids = list(pool["chembl_id"].astype(str))

    vs_ref_rows = []
    vs_chembl_rows = []
    summary = []

    for cid, smi in PURCHASED.items():
        m = mol_from_smiles(smi)
        qfp = fp(m)
        qmurcko = murcko_smiles(m)

        for rname, rfp in ref_fps.items():
            tc = DataStructs.TanimotoSimilarity(qfp, rfp)
            vs_ref_rows.append(
                {
                    "compound_id": cid,
                    "hit_smiles": smi,
                    "query_murcko": qmurcko,
                    "ref_name": rname,
                    "ref_murcko": ref_murcko[rname],
                    "ecfp4_tanimoto": round(tc, 4),
                    "same_murcko": qmurcko == ref_murcko[rname] and qmurcko != "",
                }
            )

        max_tc, nearest = max_tc_against(qfp, pool_fps, pool_ids)
        nearest_row = pool.loc[pool["chembl_id"].astype(str) == nearest].iloc[0]
        vs_chembl_rows.append(
            {
                "compound_id": cid,
                "hit_smiles": smi,
                "query_murcko": qmurcko,
                "nearest_chembl_id": nearest,
                "nearest_smiles": nearest_row["smiles"],
                "nearest_murcko": nearest_row["murcko"],
                "nearest_source_isoform": nearest_row["source_isoform"],
                "max_ecfp4_tanimoto_chembl_jnk": round(float(max_tc), 4),
                "same_murcko_as_nearest": qmurcko == nearest_row["murcko"] and qmurcko != "",
            }
        )

        max_ref_name = max(ref_fps, key=lambda k: DataStructs.TanimotoSimilarity(qfp, ref_fps[k]))
        max_ref_tc = DataStructs.TanimotoSimilarity(qfp, ref_fps[max_ref_name])
        summary.append(
            {
                "compound_id": cid,
                "smiles": smi,
                "murcko": qmurcko,
                "maxTc_vs_literature_refs": round(float(max_ref_tc), 4),
                "nearest_literature_ref": max_ref_name,
                "maxTc_vs_chembl_jnk_pool": round(float(max_tc), 4),
                "nearest_chembl_id": nearest,
                "interpretation": (
                    "near-neighbor / known-like"
                    if max(max_ref_tc, max_tc) >= 0.55
                    else (
                        "moderate similarity"
                        if max(max_ref_tc, max_tc) >= 0.35
                        else "chemically distant from curated JNK set (ECFP4)"
                    )
                ),
            }
        )

    df_ref = pd.DataFrame(vs_ref_rows)
    df_chembl = pd.DataFrame(vs_chembl_rows)
    df_sum = pd.DataFrame(summary)

    df_ref.to_csv(OUT_DIR / "c1_vs_literature_refs.csv", index=False)
    df_chembl.to_csv(OUT_DIR / "c1_vs_chembl_jnk_nearest.csv", index=False)
    df_sum.to_csv(OUT_DIR / "c1_novelty_summary.csv", index=False)

    report = {
        "fingerprint": "Morgan/ECFP4 radius=2 nBits=2048",
        "scaffold": "Bemis-Murcko (RDKit)",
        "chembl_pool_n": int(len(pool)),
        "purchased": PURCHASED,
        "summary": summary,
        "note": (
            "AD_maxTc in shortlist CSV is vs training actives (project filter); "
            "this audit is independent and focused on literature refs + ChEMBL JNK pool."
        ),
    }
    (OUT_DIR / "c1_novelty_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    md = [
        "# C1 Chemotype Novelty Audit",
        "",
        f"- Fingerprint: Morgan/ECFP4 r=2, 2048-bit",
        f"- ChEMBL/JNK pool size: **{len(pool)}** unique SMILES",
        "",
        "## Summary",
        "",
        df_sum.to_markdown(index=False),
        "",
        "## vs literature references (full matrix)",
        "",
        df_ref.to_markdown(index=False),
        "",
        "## Nearest ChEMBL JNK neighbor",
        "",
        df_chembl.to_markdown(index=False),
        "",
        "## Interpretation guide",
        "",
        "- Tc ≥ 0.55: treat as known-like / near-neighbor risk for novelty claim.",
        "- 0.35–0.55: moderate; discuss scaffold relationship carefully.",
        "- < 0.35: ECFP4-distant from curated set; still may share hinge-binder pharmacophore.",
        "- Same Murcko as E1/CC-90001/SP600125 would be a red flag even at modest Tc.",
        "",
    ]
    (OUT_DIR / "C1_NOVELTY_REPORT.md").write_text("\n".join(md), encoding="utf-8")
    print(df_sum.to_string(index=False))
    print(f"Wrote outputs under {OUT_DIR}")


if __name__ == "__main__":
    main()
