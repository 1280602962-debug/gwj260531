#!/usr/bin/env python3
"""Rebuild AChE/BChE primary panel without the ChEMBL-ID prefix cap.

Rule (deterministic, no docking/AUROC peeking):
  identical to data/ache_bche_panel_v0/scripts/build_strict_panels.py
  EXCEPT the (class, molecule_chembl_id[:8]) cap is removed.

Kept from the original builder:
  - candidate pool = mols_ACHE.json ∩ mols_BCHE.json
  - strict 6.5 / 5.5 class assignment (gray excluded)
  - quota dual/A_only/B_only/neither = 28/28/28/16
  - seed 20260729
  - sorted ChEMBL ID intersection, then Random(seed).shuffle per class
  - take first N that remain after (now empty) extra caps
  - panel_id AB_001... in draw order

This matches the 'no additional scaffold cap' principle used for the other
strict-quota pairs. It does NOT add Track B small-molecule filters that were
never in the AChE builder.
"""
from __future__ import annotations

import csv
import json
import random
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
OUT = ROOT / "remediation_outputs"
SRC = ROOT / "data/public_pair_selection"
OLD_PANEL = ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict_with_smiles.csv"
SMILES_CACHE = ROOT / "data/ache_bche_panel_v0/tables/smiles_cache.json"
THETA6_CACHE = ROOT / "data/jcim_novelty_v0/tables/chembl_smiles_cache_theta6_v1.csv"
SEED = 20260729
HI, LO = 6.5, 5.5
QUOTA = {"dual": 28, "A_only": 28, "B_only": 28, "neither": 16}


def load_mols(t: str) -> dict[str, float]:
    return {k: float(v) for k, v in json.loads((SRC / f"mols_{t}.json").read_text()).items()}


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
    if SMILES_CACHE.exists():
        for k, v in json.loads(SMILES_CACHE.read_text()).items():
            if v:
                out[str(k)] = v
    if THETA6_CACHE.exists():
        with THETA6_CACHE.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                cid = row.get("chembl_id") or row.get("molecule_chembl_id")
                smi = row.get("smiles") or row.get("canonical_smiles")
                if cid and smi and cid not in out:
                    out[cid] = smi
    if OLD_PANEL.exists():
        with OLD_PANEL.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                cid = row["molecule_chembl_id"]
                smi = row.get("smiles") or ""
                if cid and smi and cid not in out:
                    out[cid] = smi
    return out


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


def main() -> int:
    a, b = load_mols("ACHE"), load_mols("BCHE")
    both = sorted(set(a) & set(b))
    rows = []
    supply = {"dual": 0, "A_only": 0, "B_only": 0, "neither": 0, "gray": 0}
    for cid in both:
        pa, pb = a[cid], b[cid]
        cls = classify(pa, pb)
        if cls is None:
            supply["gray"] += 1
            continue
        supply[cls] += 1
        rows.append(
            {
                "molecule_chembl_id": cid,
                "class": cls,
                "pchembl_ACHE": pa,
                "pchembl_BCHE": pb,
                "min_pchembl": min(pa, pb),
            }
        )
    pool = rows
    rng = random.Random(SEED)
    picked = []
    for cls, need in QUOTA.items():
        cls_pool = [r for r in pool if r["class"] == cls]
        rng.shuffle(cls_pool)
        take = cls_pool[:need]
        if len(take) < need:
            print(f"WARN class {cls}: {len(take)}/{need} (supply={supply[cls]})")
        picked.extend(take)

    smiles_map = load_smiles()
    new_ids = {r["molecule_chembl_id"] for r in picked}
    old_rows = list(csv.DictReader(OLD_PANEL.open(newline="", encoding="utf-8")))
    old_ids = {r["molecule_chembl_id"] for r in old_rows}
    old_class = {r["molecule_chembl_id"]: r["class"] for r in old_rows}

    new_by_id = {}
    for i, r in enumerate(picked, 1):
        cid = r["molecule_chembl_id"]
        smi = smiles_map.get(cid, "")
        rec = {
            "panel_id": f"AB_{i:03d}",
            "molecule_chembl_id": cid,
            "class": r["class"],
            "pchembl_ACHE": r["pchembl_ACHE"],
            "pchembl_BCHE": r["pchembl_BCHE"],
            "min_pchembl": r["min_pchembl"],
            "label_rule": "strict_6.5_5.5",
            "gray_excluded": True,
            "smiles": smi,
            "scaffold_id": murcko(smi) if smi else "",
            "reason_selected": (
                f"strict_quota_sample seed={SEED} class={r['class']} "
                f"quota={QUOTA[r['class']]}; no ID-prefix cap; no Murcko cap"
            ),
            "source_document": "",
            "smiles_available": bool(smi),
        }
        new_by_id[cid] = rec

    cmp_rows = []
    all_ids = sorted(old_ids | new_ids)
    for cid in all_ids:
        cls = new_by_id.get(cid, {}).get("class") or old_class.get(cid, "")
        smi = smiles_map.get(cid, new_by_id.get(cid, {}).get("smiles", ""))
        in_old = cid in old_ids
        in_new = cid in new_ids
        if in_old and in_new:
            reason = "kept: selected by both original ID-prefix-capped draw and no-cap draw"
        elif in_new:
            reason = "added: selected after removing ChEMBL-ID prefix cap"
        else:
            reason = "removed: previously selected only because of ID-prefix-cap skipping earlier shuffle items"
        cmp_rows.append(
            {
                "ligand_id": cid,
                "old_panel": "yes" if in_old else "no",
                "new_panel": "yes" if in_new else "no",
                "class": cls,
                "reason_selected": reason,
                "scaffold_id": murcko(smi) if smi else "",
                "source_document": "",
                "old_panel_id": next((r["panel_id"] for r in old_rows if r["molecule_chembl_id"] == cid), ""),
                "new_panel_id": new_by_id.get(cid, {}).get("panel_id", ""),
            }
        )

    dest = OUT / "phase2_ache_bche"
    dest.mkdir(parents=True, exist_ok=True)
    new_path = dest / "panel_v0_strict_no_id_prefix_cap.csv"
    with new_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(next(iter(new_by_id.values())).keys()))
        w.writeheader()
        w.writerows(new_by_id[r["molecule_chembl_id"]] for r in picked)

    cmp_path = OUT / "ache_bche_old_vs_corrected_panel.csv"
    with cmp_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cmp_rows[0].keys()))
        w.writeheader()
        w.writerows(cmp_rows)

    n_add = sum(1 for r in cmp_rows if r["old_panel"] == "no" and r["new_panel"] == "yes")
    n_rm = sum(1 for r in cmp_rows if r["old_panel"] == "yes" and r["new_panel"] == "no")
    n_keep = sum(1 for r in cmp_rows if r["old_panel"] == "yes" and r["new_panel"] == "yes")
    n_smi = sum(1 for r in picked if smiles_map.get(r["molecule_chembl_id"]))
    summary = dest / "panel_reconstruction_summary.txt"
    summary.write_text(
        "\n".join(
            [
                "AChE/BChE panel reconstruction WITHOUT ChEMBL ID prefix cap",
                f"pool intersection n={len(both)}",
                f"strict supply D/A/B/N/gray={supply['dual']}/{supply['A_only']}/{supply['B_only']}/{supply['neither']}/{supply['gray']}",
                f"quota {QUOTA}",
                f"seed {SEED}",
                f"old n={len(old_ids)} new n={len(new_ids)} kept={n_keep} added={n_add} removed={n_rm}",
                f"new ligands with cached SMILES {n_smi}/{len(picked)}",
                "sampling: sorted intersection -> per-class Random(20260729).shuffle -> first N",
                "no docking scores consulted",
                "no AUROC consulted",
                "Track B small-mol filter NOT applied (was not in original AChE builder)",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(summary.read_text())
    print("wrote", cmp_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
