#!/usr/bin/env python3
"""Sample ≥40 NLRP3 W4 decoys BEFORE docking (C5 Job B).

Locks (docs/C5_DOCKING_WORKLIST.md):
  - Source: data/benchmarks/urat1_true_decoy/true_decoys.csv
  - Do NOT use experimental_inactives.csv
  - Do NOT use clinical 156
  - Property windows vs 9 NLRP3 positives (same keys as TrueDecoy summary)
  - Morgan ECFP4 r=2 2048-bit max TC to ANY NLRP3 positive ≤ 0.5
  - Sample seed 0xC5DEC0
  - Write CSV first; docking is a separate step
  - Drop CHEMBL3183703 (SMILES duplicate of MCC950) from positives
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import DataStructs

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PANEL_CSV = (
    PROJECT_ROOT
    / "data/campaigns/c1/05_metrics/nlrp3_structural_panel/panel_ligands.csv"
)
TRUE_DECOYS = (
    PROJECT_ROOT / "data/benchmarks/urat1_true_decoy/true_decoys.csv"
)
SUMMARY = (
    PROJECT_ROOT / "data/benchmarks/urat1_true_decoy/summary.json"
)
OUT_DIR = PROJECT_ROOT / "data/campaigns/c5/02_nlrp3_panel"
SEED = 0xC5DEC0
N_DECOYS = 40
MAX_TC = 0.5
DROP_POSITIVE_IDS = {"CHEMBL3183703"}

PROP_KEYS = [
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
]


def fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def props(mol):
    return {
        "MolWt": Descriptors.MolWt(mol),
        "MolLogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=N_DECOYS)
    ap.add_argument("--seed", type=lambda x: int(x, 0), default=SEED)
    ap.add_argument("--max-tc", type=float, default=MAX_TC)
    args = ap.parse_args()

    windows = json.loads(SUMMARY.read_text())["property_windows"]
    panel = pd.read_csv(PANEL_CSV)
    pos = panel[
        panel["role"].isin(
            ["crystal_positive", "tool_positive", "chembl_sulfonylurea_active"]
        )
    ].copy()
    pos = pos[~pos["ligand_id"].isin(DROP_POSITIVE_IDS)].reset_index(drop=True)
    if len(pos) != 9:
        raise SystemExit(f"expected 9 positives after dedup, got {len(pos)}")

    # Freeze positive list BEFORE any clinical name is needed for decoy sampling
    pos_out = OUT_DIR / "positives_locked.csv"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pos.to_csv(pos_out, index=False)

    pos_mols = []
    pos_fps = []
    pos_props = []
    for _, row in pos.iterrows():
        m = Chem.MolFromSmiles(row["smiles"])
        if m is None:
            raise SystemExit(f"bad positive SMILES {row['ligand_id']}")
        pos_mols.append(m)
        pos_fps.append(fp(m))
        pos_props.append(props(m))

    # envelopes: union of windows around each positive
    # candidate must fall within at least one positive's window for ALL properties
    # (TrueDecoy-style matching envelope)
    lo = {k: min(p[k] for p in pos_props) - windows[k] for k in PROP_KEYS}
    hi = {k: max(p[k] for p in pos_props) + windows[k] for k in PROP_KEYS}

    td = pd.read_csv(TRUE_DECOYS)
    # Prefer property-matched commercial decoys; exclude experimental_inactive rows
    # that were folded into true_decoys.csv (worklist forbids experimental_inactives.csv).
    if "decoy_source" in td.columns:
        td = td[td["decoy_source"].astype(str) != "experimental_inactive"].copy()
    print(f"true_decoy pool after source filter: {len(td)}", flush=True)

    rng = np.random.default_rng(args.seed)
    order = rng.permutation(len(td))

    kept = []
    rejected = {"prop": 0, "tc": 0, "parse": 0}
    for idx in order:
        row = td.iloc[int(idx)]
        smi = row["canonical_smiles"]
        m = Chem.MolFromSmiles(smi)
        if m is None:
            rejected["parse"] += 1
            continue
        pr = props(m)
        if any(pr[k] < lo[k] or pr[k] > hi[k] for k in PROP_KEYS):
            rejected["prop"] += 1
            continue
        f = fp(m)
        max_tc = max(DataStructs.TanimotoSimilarity(f, pf) for pf in pos_fps)
        if max_tc > args.max_tc:
            rejected["tc"] += 1
            continue
        kept.append(
            {
                "decoy_id": f"C5W4D_{len(kept)+1:03d}",
                "canonical_smiles": smi,
                "max_tc_nlrp3_positive": float(max_tc),
                **pr,
                "source_row": int(idx),
                "decoy_source": row.get("decoy_source", ""),
            }
        )
        if len(kept) >= args.n:
            break

    if len(kept) < args.n:
        raise SystemExit(
            f"only sampled {len(kept)}/{args.n} decoys; rejected={rejected}"
        )

    decoy_csv = OUT_DIR / "w4_decoys_locked.csv"
    pd.DataFrame(kept).to_csv(decoy_csv, index=False)
    meta = {
        "n_decoys": len(kept),
        "seed": args.seed,
        "max_tc": args.max_tc,
        "property_windows": windows,
        "envelope_lo": lo,
        "envelope_hi": hi,
        "positives_locked": str(pos_out),
        "positives": pos["ligand_id"].tolist(),
        "rejected": rejected,
        "pool": str(TRUE_DECOYS),
        "note": "CSV locked before docking; do not rename decoys by hand",
    }
    (OUT_DIR / "w4_decoy_sample_meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))
    print("wrote", decoy_csv)


if __name__ == "__main__":
    main()
