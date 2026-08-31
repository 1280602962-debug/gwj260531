#!/usr/bin/env python3
"""Build strict-quota panels for AChE/BChE and PIK3CA/PIK3CB (no docking)."""
from __future__ import annotations

import json
import random
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/home/gwj/repos/gwj260531/Dual_Target_Docking")
SRC = ROOT / "data/public_pair_selection"
SEED = 20260729
HI, LO = 6.5, 5.5
MURCKO_CAP = 5
QUOTA = {"dual": 28, "A_only": 28, "B_only": 28, "neither": 16}  # ~100; shrink if supply short


def load(t):
    return {k: float(v) for k, v in json.loads((SRC / f"mols_{t}.json").read_text()).items()}


def murcko(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except Exception:
        return None


def build_pair(name_a, name_b, out_csv: Path, prefix: str):
    a, b = load(name_a), load(name_b)
    both = sorted(set(a) & set(b))
    rows = []
    for cid in both:
        pa, pb = a[cid], b[cid]
        if pa >= HI and pb >= HI:
            cls = "dual"
        elif pa >= HI and pb <= LO:
            cls = "A_only"
        elif pb >= HI and pa <= LO:
            cls = "B_only"
        elif pa <= LO and pb <= LO:
            cls = "neither"
        else:
            continue  # gray excluded from main panel
        rows.append(
            {
                "molecule_chembl_id": cid,
                "class": cls,
                f"pchembl_{name_a}": pa,
                f"pchembl_{name_b}": pb,
                "min_pchembl": min(pa, pb),
            }
        )
    df = pd.DataFrame(rows)
    rng = random.Random(SEED)
    picked = []
    caps = {}
    for cls, need in QUOTA.items():
        pool = df[df["class"] == cls].to_dict("records")
        rng.shuffle(pool)
        got = 0
        for r in pool:
            # SMILES deferred — use chembl id only for now; fetch later if needed
            # Murcko cap requires smiles; without network use chembl-id diversity only
            key = (cls, r["molecule_chembl_id"][:8])  # weak diversity proxy
            if caps.get(key, 0) >= MURCKO_CAP:
                continue
            caps[key] = caps.get(key, 0) + 1
            picked.append(r)
            got += 1
            if got >= need:
                break
        if got < need:
            print(f"WARN {name_a}/{name_b} class {cls}: got {got}/{need} (supply limit)")

    out = pd.DataFrame(picked)
    out.insert(0, "panel_id", [f"{prefix}_{i:03d}" for i in range(1, len(out) + 1)])
    out["label_rule"] = "strict_6.5_5.5"
    out["gray_excluded"] = True
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)
    print(out_csv, len(out), out["class"].value_counts().to_dict())
    return out


def main():
    build_pair(
        "ACHE",
        "BCHE",
        ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict.csv",
        "AB",
    )
    build_pair(
        "PIK3CA",
        "PIK3CB",
        ROOT / "data/pik3ca_pik3cb_panel_v0/tables/panel_v0_strict.csv",
        "PAB",
    )


if __name__ == "__main__":
    main()
