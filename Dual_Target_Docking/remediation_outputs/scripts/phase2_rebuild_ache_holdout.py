#!/usr/bin/env python3
"""Rebuild AChE/BChE unused-pool holdout after the no-ID-prefix-cap main panel.

Same rule as data/jcim_holdout_v0/scripts/build_holdout_candidate_pool_v1.py:
  unused pool = strict 6.5/5.5 D/A/B minus NEW main IDs
  seed 20260731, quota dual/A_only/B_only = 20/20/20, Murcko cap 3
  neither not sampled
  skip ligands without SMILES

SMILES are taken from local caches first; ChEMBL API only if missing.
Does not consult docking scores or AUROC.
"""
from __future__ import annotations

import csv
import json
import random
import time
import urllib.request
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
SRC = ROOT / "data/public_pair_selection"
OUT = ROOT / "remediation_outputs/phase2_ache_bche"
NEW_PANEL = OUT / "panel_v0_strict_no_id_prefix_cap.csv"
OLD_HOLDOUT = ROOT / "data/jcim_holdout_v0/tables/holdout_panel_HOAB.csv"
OLD_MAIN = ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict_with_smiles.csv"
HI, LO = 6.5, 5.5
HOLDOUT_SEED = 20260731
MURCKO_CAP = 3
QUOTA = {"dual": 20, "A_only": 20, "B_only": 20}


def load_mols(t: str) -> dict[str, float]:
    return {k: float(v) for k, v in json.loads((SRC / f"mols_{t}.json").read_text()).items()}


def classify(pa: float, pb: float) -> str | None:
    if pa >= HI and pb >= HI:
        return "dual"
    if pa >= HI and pb <= LO:
        return "A_only"
    if pb >= HI and pa <= LO:
        return "B_only"
    if pa <= LO and pb <= LO:
        return "neither"
    return None


def murcko(smiles: str | None) -> str:
    if not smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except Exception:
        return ""


def load_smiles() -> dict[str, str]:
    out: dict[str, str] = {}
    caches = [
        ROOT / "data/ache_bche_panel_v0/tables/smiles_cache.json",
        ROOT / "data/jcim_strengthen_t0t1_v0/tables/enrichment_smiles_cache.json",
        ROOT / "data/pik3ca_pik3cb_panel_v0/tables/smiles_cache.json",
    ]
    for p in caches:
        if p.exists():
            for k, v in json.loads(p.read_text()).items():
                if v:
                    out[str(k)] = v
    csvs = [
        ROOT / "data/jcim_novelty_v0/tables/chembl_smiles_cache_theta6_v1.csv",
        OLD_HOLDOUT,
        OLD_MAIN,
        NEW_PANEL,
        ROOT / "data/jcim_holdout_v0/tables/holdout_ligand_scores_v1.csv",
    ]
    for p in csvs:
        if not p.exists():
            continue
        with p.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                cid = row.get("molecule_chembl_id") or row.get("chembl_id")
                smi = row.get("smiles") or row.get("canonical_smiles")
                if cid and smi and cid not in out:
                    out[cid] = smi
    return out


def fetch_smiles(chembl_id: str) -> str | None:
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as resp:
                data = json.loads(resp.read())
            struct = data.get("molecule_structures") or {}
            return struct.get("canonical_smiles")
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return None


def main() -> int:
    a, b = load_mols("ACHE"), load_mols("BCHE")
    both = sorted(set(a) & set(b))
    new_main = {r["molecule_chembl_id"] for r in csv.DictReader(NEW_PANEL.open())}
    old_main = {r["molecule_chembl_id"] for r in csv.DictReader(OLD_MAIN.open())}
    old_hold = list(csv.DictReader(OLD_HOLDOUT.open()))
    old_hold_ids = {r["molecule_chembl_id"] for r in old_hold}
    smiles = load_smiles()

    pool = []
    for cid in both:
        cls = classify(a[cid], b[cid])
        if cls is None or cls == "neither":
            continue
        pool.append(
            {
                "molecule_chembl_id": cid,
                "class": cls,
                "pchembl_ACHE": a[cid],
                "pchembl_BCHE": b[cid],
                "used_in_new_main": cid in new_main,
            }
        )
    remaining = [r for r in pool if not r["used_in_new_main"]]
    rng = random.Random(HOLDOUT_SEED)
    picked = []
    caps: dict[tuple[str, str], int] = {}
    n_fetch = 0
    n_skip_smi = 0
    for cls, need in QUOTA.items():
        cls_pool = [r for r in remaining if r["class"] == cls]
        rng.shuffle(cls_pool)
        got = 0
        for r in cls_pool:
            if got >= need:
                break
            cid = r["molecule_chembl_id"]
            smi = smiles.get(cid)
            if not smi:
                smi = fetch_smiles(cid)
                n_fetch += 1
                if smi:
                    smiles[cid] = smi
            if not smi:
                n_skip_smi += 1
                continue
            scaf = murcko(smi)
            key = (cls, scaf or cid)
            if caps.get(key, 0) >= MURCKO_CAP:
                continue
            caps[key] = caps.get(key, 0) + 1
            rec = dict(r)
            rec["smiles"] = smi
            rec["murcko_scaffold"] = scaf
            rec["holdout_seed"] = HOLDOUT_SEED
            rec["label_rule"] = "strict_6.5_5.5"
            rec["source"] = "unused_pool_after_no_id_prefix_cap_main"
            picked.append(rec)
            got += 1
        if got < need:
            print(f"WARN class {cls}: {got}/{need}")

    for i, r in enumerate(picked, 1):
        r["holdout_id"] = f"HOAB_{i:03d}"

    dest = OUT / "holdout_panel_HOAB_corrected.csv"
    fields = [
        "holdout_id",
        "molecule_chembl_id",
        "class",
        "pchembl_ACHE",
        "pchembl_BCHE",
        "used_in_new_main",
        "smiles",
        "murcko_scaffold",
        "label_rule",
        "holdout_seed",
        "source",
    ]
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(picked)

    new_ids = {r["molecule_chembl_id"] for r in picked}
    overlap_id = new_ids & new_main
    from rdkit.Chem import MolFromSmiles
    from rdkit.Chem import MolToSmiles

    def canon(s: str) -> str:
        m = MolFromSmiles(s)
        return MolToSmiles(m) if m is not None else s

    main_smi = {canon(r["smiles"]) for r in csv.DictReader(NEW_PANEL.open()) if r.get("smiles")}
    hold_smi = {canon(r["smiles"]) for r in picked if r.get("smiles")}
    overlap_smi = main_smi & hold_smi

    cmp_path = OUT / "holdout_old_vs_corrected.csv"
    all_ids = sorted(old_hold_ids | new_ids)
    old_map = {r["molecule_chembl_id"]: r for r in old_hold}
    new_map = {r["molecule_chembl_id"]: r for r in picked}
    with cmp_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["ligand_id", "old_holdout", "new_holdout", "class", "reason"],
        )
        w.writeheader()
        for cid in all_ids:
            in_o, in_n = cid in old_hold_ids, cid in new_ids
            cls = (new_map.get(cid) or old_map.get(cid) or {}).get("class", "")
            if in_o and in_n:
                reason = "kept"
            elif in_n:
                reason = "added after main-panel unused-pool rebuild"
            else:
                reason = "removed after main-panel unused-pool rebuild"
            w.writerow(
                {
                    "ligand_id": cid,
                    "old_holdout": "yes" if in_o else "no",
                    "new_holdout": "yes" if in_n else "no",
                    "class": cls,
                    "reason": reason,
                }
            )

    summary = OUT / "holdout_rebuild_summary.txt"
    n_add = sum(1 for cid in new_ids if cid not in old_hold_ids)
    n_rm = sum(1 for cid in old_hold_ids if cid not in new_ids)
    n_keep = len(new_ids & old_hold_ids)
    # new main must not appear in holdout; old-main ligand CHEMBL173607 may now be eligible
    summary.write_text(
        "\n".join(
            [
                "AChE/BChE unused-pool holdout rebuild",
                f"seed {HOLDOUT_SEED} Murcko_cap {MURCKO_CAP} quota {QUOTA}",
                f"new main n={len(new_main)} old main n={len(old_main)}",
                f"main ID delta added={len(new_main-old_main)} removed={len(old_main-new_main)}",
                f"unused remaining n={len(remaining)}",
                f"holdout n={len(picked)} kept={n_keep} added={n_add} removed={n_rm}",
                f"ligand ID overlap main∩holdout = {len(overlap_id)}",
                f"canonical SMILES overlap main∩holdout = {len(overlap_smi)}",
                f"API fetches {n_fetch} skipped_no_smiles {n_skip_smi}",
                f"CHEMBL3960861 in holdout: {'yes' if 'CHEMBL3960861' in new_ids else 'no'}",
                f"CHEMBL173607 in holdout: {'yes' if 'CHEMBL173607' in new_ids else 'no'}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(summary.read_text())
    print("wrote", dest, cmp_path)
    if overlap_id or overlap_smi:
        print("ERROR overlap not zero")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
