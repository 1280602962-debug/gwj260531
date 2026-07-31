#!/usr/bin/env python3
"""Build a frozen, zero-fabrication external holdout pool for DualFourClass-Bench.

Design (must match CLAIM_CEILING.md — no new label rule, no post-hoc tuning):
  1. Recompute the SAME strict label rule (pChEMBL >= 6.5 / <= 5.5) used to build the
     three frozen main panels (PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB), from the same
     cached ChEMBL max-pChEMBL dictionaries (`data/public_pair_selection/mols_*.json`).
  2. Remove every ChEMBL ID already present in the corresponding frozen panel CSV
     (the panel that produced the reported Table 2 numbers). What remains is a
     genuine "never touched during panel construction or protocol tuning" pool.
  3. From that remaining pool, draw a fixed-size holdout quota with a NEW seed
     (frozen in this file, chosen before any holdout score exists) and a Murcko
     scaffold cap to avoid a single chemotype dominating any one class.
  4. Fetch SMILES for the sampled holdout IDs only, from the live ChEMBL API.

This script performs NO docking and NO scoring. It only freezes the ligand list.
Docking must reuse the already-frozen receptors/boxes/exhaustiveness/prep protocol
exactly as in Table 1 (Methods 2.4-2.5) with no changes made after seeing results.
"""
from __future__ import annotations

import json
import random
import time
import urllib.request
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "data/public_pair_selection"
OUT = ROOT / "data/jcim_holdout_v0/tables"
OUT.mkdir(parents=True, exist_ok=True)

HI, LO = 6.5, 5.5
HOLDOUT_SEED = 20260731  # frozen at plan time; distinct from every panel-building seed used so far
MURCKO_CAP = 3
QUOTA = {"dual": 20, "A_only": 20, "B_only": 20}

PAIRS = [
    {
        "name_a": "PIK3CA",
        "name_b": "MTOR",
        "pair_label": "PIK3CA/mTOR",
        "used_panel_csv": ROOT / "data/pik3ca_mtor_panel110_rdkit_v0/tables/panel_v0_110.csv",
        "prefix": "HOPM",
    },
    {
        "name_a": "ACHE",
        "name_b": "BCHE",
        "pair_label": "AChE/BChE",
        "used_panel_csv": ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict_with_smiles.csv",
        "prefix": "HOAB",
    },
    {
        "name_a": "PIK3CA",
        "name_b": "PIK3CB",
        "pair_label": "PIK3CA/PIK3CB",
        "used_panel_csv": ROOT / "data/pik3ca_pik3cb_panel_v0/tables/panel_v0_strict_with_smiles.csv",
        "prefix": "HOAP",
    },
]


def load_mols(target: str) -> dict[str, float]:
    return {k: float(v) for k, v in json.loads((SRC / f"mols_{target}.json").read_text()).items()}


def label(pa: float, pb: float) -> str | None:
    if pa >= HI and pb >= HI:
        return "dual"
    if pa >= HI and pb <= LO:
        return "A_only"
    if pb >= HI and pa <= LO:
        return "B_only"
    if pa <= LO and pb <= LO:
        return "neither"
    return None


def fetch_smiles(chembl_id: str, retries: int = 3) -> str | None:
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=15) as resp:
                data = json.loads(resp.read())
            struct = data.get("molecule_structures") or {}
            return struct.get("canonical_smiles")
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return None


def murcko(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except Exception:
        return None


def build_pair(cfg: dict) -> None:
    a, b = cfg["name_a"], cfg["name_b"]
    ma, mb = load_mols(a), load_mols(b)
    both = sorted(set(ma) & set(mb))

    rows = []
    for cid in both:
        pa, pb = ma[cid], mb[cid]
        cls = label(pa, pb)
        if cls is None or cls == "neither":
            continue
        rows.append({"molecule_chembl_id": cid, "class": cls, f"pchembl_{a}": pa, f"pchembl_{b}": pb})
    pool = pd.DataFrame(rows)

    used_ids: set[str] = set()
    if cfg["used_panel_csv"].exists():
        used_ids = set(pd.read_csv(cfg["used_panel_csv"])["molecule_chembl_id"].astype(str))
    else:
        print(f"WARN: used-panel CSV missing for {cfg['pair_label']}: {cfg['used_panel_csv']}")

    pool["used_in_frozen_panel"] = pool["molecule_chembl_id"].isin(used_ids)
    pool_path = OUT / f"strict_pool_full_{cfg['prefix']}.csv"
    pool.to_csv(pool_path, index=False)

    remaining = pool[~pool["used_in_frozen_panel"]].copy()
    rng = random.Random(HOLDOUT_SEED)

    picked_rows = []
    scaffold_caps: dict[tuple[str, str], int] = {}
    smiles_cache: dict[str, str | None] = {}

    for cls, need in QUOTA.items():
        cls_pool = remaining[remaining["class"] == cls].to_dict("records")
        rng.shuffle(cls_pool)
        got = 0
        for r in cls_pool:
            if got >= need:
                break
            cid = r["molecule_chembl_id"]
            if cid not in smiles_cache:
                smiles_cache[cid] = fetch_smiles(cid)
            smi = smiles_cache[cid]
            if not smi:
                continue
            scaf = murcko(smi)
            key = (cls, scaf or cid)
            if scaffold_caps.get(key, 0) >= MURCKO_CAP:
                continue
            scaffold_caps[key] = scaffold_caps.get(key, 0) + 1
            r = dict(r)
            r["smiles"] = smi
            r["murcko_scaffold"] = scaf
            picked_rows.append(r)
            got += 1
        if got < need:
            print(f"WARN {cfg['pair_label']} class {cls}: got {got}/{need} from unused pool")

    holdout = pd.DataFrame(picked_rows)
    if len(holdout):
        holdout.insert(0, "holdout_id", [f"{cfg['prefix']}_{i:03d}" for i in range(1, len(holdout) + 1)])
    holdout["label_rule"] = "strict_6.5_5.5"
    holdout["holdout_seed"] = HOLDOUT_SEED
    holdout["source"] = "chembl_unused_pool_post_panel_freeze"

    holdout_path = OUT / f"holdout_panel_{cfg['prefix']}.csv"
    holdout.to_csv(holdout_path, index=False)

    print(f"== {cfg['pair_label']} ==")
    print(f"  strict pool total: {len(pool)}  used_in_panel: {pool['used_in_frozen_panel'].sum()}")
    print(f"  remaining pool by class: {remaining['class'].value_counts().to_dict()}")
    print(f"  holdout drawn: {len(holdout)}  by class: {holdout['class'].value_counts().to_dict() if len(holdout) else {}}")
    print(f"  -> {holdout_path}")


def main() -> None:
    for cfg in PAIRS:
        build_pair(cfg)


if __name__ == "__main__":
    main()
