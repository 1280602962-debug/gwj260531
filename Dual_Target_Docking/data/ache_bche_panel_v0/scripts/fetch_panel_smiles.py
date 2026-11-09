#!/usr/bin/env python3
"""Fetch ChEMBL SMILES for strict panels (API recovered)."""
from __future__ import annotations

import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path("/home/gwj/repos/gwj260531/Dual_Target_Docking")
CACHE = ROOT / "data/ache_bche_panel_v0/tables/smiles_cache.json"
# shared cache for both packs
CACHE.parent.mkdir(parents=True, exist_ok=True)


def fetch_smiles(cid: str) -> str | None:
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json"
    req = urllib.request.Request(url, headers={"User-Agent": "dualfourclass-jcim/0.1"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = json.load(resp)
    structs = data.get("molecule_structures") or {}
    return structs.get("canonical_smiles")


def enrich(panel_csv: Path, cache: dict):
    df = pd.read_csv(panel_csv)
    smiles = []
    for i, cid in enumerate(df["molecule_chembl_id"]):
        if cid in cache and cache[cid]:
            smiles.append(cache[cid])
            continue
        try:
            s = fetch_smiles(cid)
            cache[cid] = s
            smiles.append(s)
            print(f"[{i+1}/{len(df)}] {cid} ok={bool(s)}", flush=True)
            time.sleep(0.05)
        except Exception as e:
            print(f"FAIL {cid}: {e}", flush=True)
            smiles.append(None)
            time.sleep(0.5)
    df["smiles"] = smiles
    out = panel_csv.with_name(panel_csv.stem + "_with_smiles.csv")
    df.to_csv(out, index=False)
    # also overwrite main panel for docking
    df.to_csv(panel_csv, index=False)
    print("wrote", out, "n_smiles", df["smiles"].notna().sum())
    return df


def main():
    cache = {}
    if CACHE.exists():
        cache = json.loads(CACHE.read_text())
    enrich(ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict.csv", cache)
    enrich(ROOT / "data/pik3ca_pik3cb_panel_v0/tables/panel_v0_strict.csv", cache)
    CACHE.write_text(json.dumps(cache, indent=2))
    # mirror cache
    (ROOT / "data/pik3ca_pik3cb_panel_v0/tables/smiles_cache.json").write_text(
        json.dumps(cache, indent=2)
    )


if __name__ == "__main__":
    main()
