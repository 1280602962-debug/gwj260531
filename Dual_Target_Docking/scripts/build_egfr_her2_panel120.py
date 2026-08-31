#!/usr/bin/env python3
"""Build EGFR/HER2 expanded panel (~100–120) for Stage-1 S1 gate.

Rules (from AGENT_COMMAND_STAGE1_EGFR_EXPAND.md + panel40):
  - paired only (measured both ends in mols_EGFR ∩ mols_HER2)
  - theta pChEMBL = 6.0; untested ≠ inactive
  - keep all panel40 members
  - prioritize A_only / B_only hard negatives
  - Murcko scaffold quota ≤5 per class
  - MW 180–750, heavy atoms ≤55
  - architecture NOT used as selection filter
"""
from __future__ import annotations

import csv
import json
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "public_pair_selection"
PANEL40 = ROOT / "data" / "egfr_her2_panel40_v0" / "tables" / "panel_v0_40.csv"
OUT_DIR = ROOT / "data" / "egfr_her2_panel120_v0"
THETA = 6.0
TARGET_N = 110
MAX_PER_SCAFFOLD = 5
MW_MIN, MW_MAX, HA_MAX = 180.0, 750.0, 55
# Prefer hardneg-heavy mix; dual/neither as balance
QUOTA = {"dual": 28, "A_only": 38, "B_only": 32, "neither": 12}  # sum=110
SEED = 20260727


def classify(pa: float, pb: float) -> str:
    if pa >= THETA and pb >= THETA:
        return "dual"
    if pa >= THETA and pb < THETA:
        return "A_only"
    if pb >= THETA and pa < THETA:
        return "B_only"
    return "neither"


def chembl_molecule(cid: str, retries: int = 4) -> dict:
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json"
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.load(r)
        except Exception as e:
            last = e
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"ChEMBL fetch failed {cid}: {last}")


def murcko(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)


def ok_physchem(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mw = Descriptors.MolWt(mol)
    ha = mol.GetNumHeavyAtoms()
    return (MW_MIN <= mw <= MW_MAX) and (ha <= HA_MAX)


def fp(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "tables").mkdir(exist_ok=True)

    egfr = json.loads((CACHE / "mols_EGFR.json").read_text())
    her2 = json.loads((CACHE / "mols_HER2.json").read_text())
    paired = sorted(set(egfr) & set(her2))

    panel40 = list(csv.DictReader(PANEL40.open()))
    keep_ids = {r["molecule_chembl_id"] for r in panel40}
    print("panel40", len(panel40), Counter(r["class"] for r in panel40))

    # Build candidate pool by class
    by_class: dict[str, list[str]] = defaultdict(list)
    for cid in paired:
        by_class[classify(float(egfr[cid]), float(her2[cid]))].append(cid)
    for k in by_class:
        # hardneg: prefer larger activity gap
        if k == "A_only":
            by_class[k].sort(key=lambda c: (float(egfr[c]) - float(her2[c])), reverse=True)
        elif k == "B_only":
            by_class[k].sort(key=lambda c: (float(her2[c]) - float(egfr[c])), reverse=True)
        elif k == "dual":
            by_class[k].sort(key=lambda c: min(float(egfr[c]), float(her2[c])), reverse=True)
        else:
            by_class[k].sort(key=lambda c: max(float(egfr[c]), float(her2[c])))

    # Start with panel40 rows (preserve IDs EH40_xx)
    selected: list[dict] = []
    scaffold_count = defaultdict(int)  # (class, scaffold) -> n
    selected_fps = defaultdict(list)  # class -> fps
    chembl_to_row = {}

    smiles_cache_path = OUT_DIR / "tables" / "smiles_cache.json"
    smiles_cache = {}
    if smiles_cache_path.exists():
        smiles_cache = json.loads(smiles_cache_path.read_text())

    def ensure_smiles(cid: str) -> tuple[str, str, str]:
        if cid in smiles_cache:
            d = smiles_cache[cid]
            return d["smiles"], d.get("pref_name") or "", d.get("inchi_key") or ""
        molj = chembl_molecule(cid)
        smi = (molj.get("molecule_structures") or {}).get("canonical_smiles") or ""
        pref = molj.get("pref_name") or ""
        inch = (molj.get("molecule_structures") or {}).get("standard_inchi_key") or ""
        smiles_cache[cid] = {"smiles": smi, "pref_name": pref, "inchi_key": inch}
        if len(smiles_cache) % 20 == 0:
            smiles_cache_path.write_text(json.dumps(smiles_cache))
        return smi, pref, inch

    # Add panel40 first
    for r in panel40:
        cid = r["molecule_chembl_id"]
        smi = r["smiles"]
        scaf = murcko(smi)
        cls = r["class"]
        selected.append(
            {
                "panel_id": r["panel_id"],
                "molecule_chembl_id": cid,
                "class": cls,
                "pchembl_EGFR": float(r["pchembl_EGFR"]),
                "pchembl_HER2": float(r["pchembl_HER2"]),
                "min_pchembl": float(r["min_pchembl"]),
                "pref_name": r.get("pref_name") or "",
                "smiles": smi,
                "inchi_key": r.get("inchi_key") or "",
                "max_phase": r.get("max_phase") or "",
                "murcko_scaffold": scaf,
                "from_panel40": "yes",
                "role_note": r.get("role_note") or "panel40 retained",
            }
        )
        scaffold_count[(cls, scaf)] += 1
        f = fp(smi)
        if f is not None:
            selected_fps[cls].append(f)
        chembl_to_row[cid] = selected[-1]
        smiles_cache[cid] = {
            "smiles": smi,
            "pref_name": r.get("pref_name") or "",
            "inchi_key": r.get("inchi_key") or "",
        }

    def can_add(cid: str, cls: str, smi: str) -> bool:
        if cid in chembl_to_row:
            return False
        if not ok_physchem(smi):
            return False
        scaf = murcko(smi)
        if not scaf:
            return False
        if scaffold_count[(cls, scaf)] >= MAX_PER_SCAFFOLD:
            return False
        f = fp(smi)
        if f is None:
            return False
        # soft diversity: avoid near-identical to already selected in class
        for g in selected_fps[cls]:
            if DataStructs.TanimotoSimilarity(f, g) >= 0.90:
                return False
        return True

    # Fill quotas
    rng_order = {c: list(by_class[c]) for c in QUOTA}
    for cls, quota in QUOTA.items():
        have = sum(1 for r in selected if r["class"] == cls)
        need = max(0, quota - have)
        print(f"fill {cls}: have {have}, need {need}, pool {len(rng_order[cls])}")
        added = 0
        for cid in rng_order[cls]:
            if added >= need:
                break
            if cid in keep_ids or cid in chembl_to_row:
                continue
            try:
                smi, pref, inch = ensure_smiles(cid)
            except Exception as e:
                print("skip fetch", cid, e)
                continue
            if not smi or not can_add(cid, cls, smi):
                continue
            scaf = murcko(smi)
            selected.append(
                {
                    "panel_id": "",  # assign later
                    "molecule_chembl_id": cid,
                    "class": cls,
                    "pchembl_EGFR": float(egfr[cid]),
                    "pchembl_HER2": float(her2[cid]),
                    "min_pchembl": min(float(egfr[cid]), float(her2[cid])),
                    "pref_name": pref,
                    "smiles": smi,
                    "inchi_key": inch,
                    "max_phase": "",
                    "murcko_scaffold": scaf,
                    "from_panel40": "no",
                    "role_note": "stage1 expand hardneg-priority",
                }
            )
            scaffold_count[(cls, scaf)] += 1
            selected_fps[cls].append(fp(smi))
            chembl_to_row[cid] = selected[-1]
            added += 1
        print(f"  added {added}")

    # If under TARGET_N, top up with dual then neither then A_only
    for cls in ("dual", "neither", "A_only", "B_only"):
        if len(selected) >= TARGET_N:
            break
        for cid in rng_order[cls]:
            if len(selected) >= TARGET_N:
                break
            if cid in chembl_to_row:
                continue
            try:
                smi, pref, inch = ensure_smiles(cid)
            except Exception:
                continue
            if not smi or not can_add(cid, cls, smi):
                continue
            scaf = murcko(smi)
            selected.append(
                {
                    "panel_id": "",
                    "molecule_chembl_id": cid,
                    "class": cls,
                    "pchembl_EGFR": float(egfr[cid]),
                    "pchembl_HER2": float(her2[cid]),
                    "min_pchembl": min(float(egfr[cid]), float(her2[cid])),
                    "pref_name": pref,
                    "smiles": smi,
                    "inchi_key": inch,
                    "max_phase": "",
                    "murcko_scaffold": scaf,
                    "from_panel40": "no",
                    "role_note": "stage1 expand top-up",
                }
            )
            scaffold_count[(cls, scaf)] += 1
            selected_fps[cls].append(fp(smi))
            chembl_to_row[cid] = selected[-1]

    # Assign panel IDs: keep EH40_*; new get EH120_041+
    next_i = 41
    for r in selected:
        if r["from_panel40"] == "yes":
            continue
        r["panel_id"] = f"EH120_{next_i:03d}"
        next_i += 1

    # Sort: panel40 first by id, then new
    def sort_key(r):
        if r["from_panel40"] == "yes":
            return (0, r["panel_id"])
        return (1, r["panel_id"])

    selected.sort(key=sort_key)

    out_csv = OUT_DIR / "tables" / "panel_v0_120.csv"
    fields = [
        "panel_id",
        "molecule_chembl_id",
        "class",
        "pchembl_EGFR",
        "pchembl_HER2",
        "min_pchembl",
        "pref_name",
        "smiles",
        "inchi_key",
        "max_phase",
        "murcko_scaffold",
        "from_panel40",
        "role_note",
    ]
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in selected:
            w.writerow({k: r.get(k, "") for k in fields})

    smiles_cache_path.write_text(json.dumps(smiles_cache, indent=0))

    counts = Counter(r["class"] for r in selected)
    scaf_max = max(scaffold_count.values()) if scaffold_count else 0
    manifest = OUT_DIR / "MANIFEST.md"
    manifest.write_text(
        f"""# MANIFEST — egfr_her2_panel120_v0 (Stage-1 expand)

- theta_pchembl: **{THETA}**
- N: **{len(selected)}**
- class counts: {dict(counts)}
- Murcko max per (class,scaffold): **{MAX_PER_SCAFFOLD}** (observed max={scaf_max})
- retained panel40: {sum(1 for r in selected if r['from_panel40']=='yes')}
- new ligands: {sum(1 for r in selected if r['from_panel40']=='no')}
- docking protocol: 3POZ/3RCD, E=8, seed=20260727, n_modes=9, RTM best-of-9
- architecture: not used as selection filter
- warning flags: diagnostic only (not gated into score)
"""
    )
    print("wrote", out_csv, "N=", len(selected), dict(counts))


if __name__ == "__main__":
    main()
