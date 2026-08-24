#!/usr/bin/env python3
"""A4: max vs median pChEMBL on frozen DualFourClass panels. No new docking.

Re-fetches assay-level pchembl_value from the ChEMBL Web API for every
ligand–target pair in the four frozen panels, then relabels under θ = 6.0
and recomputes pocket-matched directional AUROC on the frozen Vina scores.
"""
from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_novelty_v0"
TAB = OUT / "tables"
CACHE = OUT / "cache" / "chembl_activity"
TAB.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

THETA = 6.0
BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

TARGETS = {
    "EGFR/HER2": ("CHEMBL203", "CHEMBL1824"),
    "AChE/BChE": ("CHEMBL220", "CHEMBL1914"),
    "PIK3CA/PIK3CB": ("CHEMBL4005", "CHEMBL3145"),
    "PIK3CA/mTOR": ("CHEMBL4005", "CHEMBL2842"),
}

PANELS = {
    "EGFR/HER2": dict(
        scores="data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        panel="data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        vina_a="3POZ_affinity",
        vina_b="3RCD_affinity",
        p_a="pchembl_EGFR",
        p_b="pchembl_HER2",
        id_from_panel=True,
    ),
    "AChE/BChE": dict(
        scores="data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        panel=None,
        vina_a="vina_ACHE",
        vina_b="vina_BCHE",
        p_a="pchembl_ACHE",
        p_b="pchembl_BCHE",
        id_from_panel=False,
    ),
    "PIK3CA/PIK3CB": dict(
        scores="data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        panel=None,
        vina_a="vina_PIK3CA",
        vina_b="vina_PIK3CB",
        p_a="pchembl_PIK3CA",
        p_b="pchembl_PIK3CB",
        id_from_panel=False,
    ),
    "PIK3CA/mTOR": dict(
        scores="data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        panel="data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        vina_a="4L23_affinity",
        vina_b="4JT6_affinity",
        p_a="pchembl_PIK3CA",
        p_b="pchembl_MTOR",
        id_from_panel=True,
    ),
}


def fnum(v):
    try:
        return float(v) if v not in ("", None) else None
    except (TypeError, ValueError):
        return None


def load_csv(p: Path):
    with p.open() as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def boot_auroc(pos, neg, n_boot=2000, seed=20260729):
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    vals = [auroc(rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True)) for _ in range(n_boot)]
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(auroc(pos, neg)), float(lo), float(hi)


def classify(pa, pb, theta=THETA):
    if pa is None or pb is None:
        return None
    a_on, b_on = pa >= theta, pb >= theta
    if a_on and b_on:
        return "dual"
    if a_on and not b_on:
        return "A_only"
    if b_on and not a_on:
        return "B_only"
    return "neither"


def fetch_pchembl(mol_id: str, target_id: str) -> dict:
    key = f"{mol_id}__{target_id}"
    cache = CACHE / f"{key}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    values = []
    offset = 0
    n_rows = 0
    while True:
        q = urllib.parse.urlencode(
            {
                "molecule_chembl_id": mol_id,
                "target_chembl_id": target_id,
                "limit": 200,
                "offset": offset,
            }
        )
        url = f"{BASE}?{q}"
        last_err = None
        for attempt in range(5):
            try:
                req = urllib.request.Request(url, headers={"Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=60) as resp:
                    payload = json.loads(resp.read().decode())
                last_err = None
                break
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as e:
                last_err = e
                time.sleep(2.0 * (attempt + 1))
        if last_err is not None:
            rec = {
                "molecule_chembl_id": mol_id,
                "target_chembl_id": target_id,
                "n_activity_rows": None,
                "n_pchembl": 0,
                "api_max": None,
                "api_median": None,
                "values": [],
                "fetch_error": str(last_err),
            }
            cache.write_text(json.dumps(rec))
            return rec
        acts = payload.get("activities") or []
        n_rows += len(acts)
        for a in acts:
            v = a.get("pchembl_value")
            if v not in (None, ""):
                try:
                    values.append(float(v))
                except (TypeError, ValueError):
                    pass
        page_meta = payload.get("page_meta") or {}
        nxt = page_meta.get("next")
        if not nxt or not acts:
            break
        offset += 200
        time.sleep(0.05)
    rec = {
        "molecule_chembl_id": mol_id,
        "target_chembl_id": target_id,
        "n_activity_rows": n_rows,
        "n_pchembl": len(values),
        "api_max": max(values) if values else None,
        "api_median": float(np.median(values)) if values else None,
        "values": values,
    }
    cache.write_text(json.dumps(rec))
    return rec


def assemble_jobs():
    jobs = []
    ligands = []
    for pair, cfg in PANELS.items():
        scores = load_csv(ROOT / cfg["scores"])
        smimap = {}
        classmap = {}
        pmap = {}
        if cfg["panel"]:
            for r in load_csv(ROOT / cfg["panel"]):
                pid = r.get("panel_id")
                smimap[pid] = r.get("molecule_chembl_id")
                classmap[pid] = r.get("class")
                pmap[pid] = (fnum(r.get(cfg["p_a"])), fnum(r.get(cfg["p_b"])))
        ta, tb = TARGETS[pair]
        for r in scores:
            a, b = fnum(r.get(cfg["vina_a"])), fnum(r.get(cfg["vina_b"]))
            if a is None or b is None:
                continue
            lig = r.get("ligand") or r.get("panel_id")
            mol = r.get("molecule_chembl_id") or smimap.get(lig)
            if not mol:
                continue
            cached_a, cached_b = pmap.get(lig, (fnum(r.get(cfg["p_a"])), fnum(r.get(cfg["p_b"]))))
            cls = r.get("class") or classmap.get(lig)
            rec = {
                "pair": pair,
                "ligand": lig,
                "molecule_chembl_id": mol,
                "frozen_class": cls,
                "cached_pA": cached_a,
                "cached_pB": cached_b,
                "vina_A": -a,
                "vina_B": -b,
                "target_A": ta,
                "target_B": tb,
            }
            ligands.append(rec)
            jobs.append((mol, ta))
            jobs.append((mol, tb))
    uniq = sorted(set(jobs))
    return ligands, uniq


def main():
    ligands, uniq = assemble_jobs()
    print(f"ligands with both scores: {len(ligands)}; unique ligand-target queries: {len(uniq)}")
    fetched = {}
    for i, (mol, tgt) in enumerate(uniq, 1):
        rec = fetch_pchembl(mol, tgt)
        fetched[(mol, tgt)] = rec
        if i % 25 == 0 or i == len(uniq):
            print(f"  fetched {i}/{len(uniq)}")

    per_lt = []
    flips = []
    for lig in ligands:
        ra = fetched[(lig["molecule_chembl_id"], lig["target_A"])]
        rb = fetched[(lig["molecule_chembl_id"], lig["target_B"])]
        max_a, med_a = ra.get("api_max"), ra.get("api_median")
        max_b, med_b = rb.get("api_max"), rb.get("api_median")
        cls_max = classify(max_a, max_b)
        cls_med = classify(med_a, med_b)
        row = {
            **{k: lig[k] for k in (
                "pair", "ligand", "molecule_chembl_id", "frozen_class",
                "cached_pA", "cached_pB",
            )},
            "api_max_A": max_a,
            "api_median_A": med_a,
            "n_pchembl_A": ra.get("n_pchembl"),
            "n_activity_A": ra.get("n_activity_rows"),
            "cache_matches_api_max_A": int(
                lig["cached_pA"] is not None and max_a is not None and abs(lig["cached_pA"] - max_a) < 1e-6
            ),
            "api_max_B": max_b,
            "api_median_B": med_b,
            "n_pchembl_B": rb.get("n_pchembl"),
            "n_activity_B": rb.get("n_activity_rows"),
            "cache_matches_api_max_B": int(
                lig["cached_pB"] is not None and max_b is not None and abs(lig["cached_pB"] - max_b) < 1e-6
            ),
            "class_max_theta6": cls_max,
            "class_median_theta6": cls_med,
            "flip_max_to_median": int(cls_max != cls_med),
            "max_ne_median_A": int(max_a is not None and med_a is not None and abs(max_a - med_a) > 1e-6),
            "max_ne_median_B": int(max_b is not None and med_b is not None and abs(max_b - med_b) > 1e-6),
        }
        per_lt.append(row)
        if cls_max != cls_med:
            flips.append(row)

    write_csv(TAB / "assay_max_vs_median_ligand_v1.csv", per_lt)

    summary = []
    auroc_rows = []
    lig_by_pair = {}
    for lig in ligands:
        lig_by_pair.setdefault(lig["pair"], []).append(lig)

    rec_by_key = {(r["pair"], r["ligand"]): r for r in per_lt}

    for pair, recs in lig_by_pair.items():
        n = len(recs)
        n_flip = sum(rec_by_key[(pair, r["ligand"])]["flip_max_to_median"] for r in recs)
        n_max_ne_med = sum(
            rec_by_key[(pair, r["ligand"])]["max_ne_median_A"] or rec_by_key[(pair, r["ligand"])]["max_ne_median_B"]
            for r in recs
        )
        n_cache_ok = sum(
            rec_by_key[(pair, r["ligand"])]["cache_matches_api_max_A"]
            and rec_by_key[(pair, r["ligand"])]["cache_matches_api_max_B"]
            for r in recs
        )
        summary.append(
            {
                "pair": pair,
                "n_ligands_scored": n,
                "n_cache_matches_both_max": n_cache_ok,
                "n_any_end_max_ne_median": n_max_ne_med,
                "n_class_flip_theta6": n_flip,
                "frac_class_flip_theta6": round(n_flip / n, 4) if n else "",
            }
        )
        for agg, cls_key in (("max", "class_max_theta6"), ("median", "class_median_theta6")):
            dual = [r for r in recs if rec_by_key[(pair, r["ligand"])][cls_key] == "dual"]
            aonly = [r for r in recs if rec_by_key[(pair, r["ligand"])][cls_key] == "A_only"]
            bonly = [r for r in recs if rec_by_key[(pair, r["ligand"])][cls_key] == "B_only"]
            da_pt, da_lo, da_hi = boot_auroc([r["vina_B"] for r in dual], [r["vina_B"] for r in aonly])
            db_pt, db_lo, db_hi = boot_auroc([r["vina_A"] for r in dual], [r["vina_A"] for r in bonly])
            sm = min(da_pt, db_pt) if da_pt == da_pt and db_pt == db_pt else float("nan")
            auroc_rows.append(
                {
                    "pair": pair,
                    "aggregation": agg,
                    "n_dual": len(dual),
                    "n_A_only": len(aonly),
                    "n_B_only": len(bonly),
                    "auroc_D_vs_A": None if da_pt != da_pt else round(da_pt, 4),
                    "ci_lo_D_vs_A": None if da_lo != da_lo else round(da_lo, 4),
                    "ci_hi_D_vs_A": None if da_hi != da_hi else round(da_hi, 4),
                    "auroc_D_vs_B": None if db_pt != db_pt else round(db_pt, 4),
                    "ci_lo_D_vs_B": None if db_lo != db_lo else round(db_lo, 4),
                    "ci_hi_D_vs_B": None if db_hi != db_hi else round(db_hi, 4),
                    "summary_min": None if sm != sm else round(sm, 4),
                    "underpowered": int(min(len(dual), len(aonly), len(bonly)) < 8),
                }
            )

    write_csv(TAB / "assay_max_vs_median_summary_v1.csv", summary)
    write_csv(TAB / "assay_max_vs_median_auroc_v1.csv", auroc_rows)
    write_csv(TAB / "assay_max_vs_median_flips_v1.csv", flips)
    print("summary:")
    for r in summary:
        print(r)
    print("auroc:")
    for r in auroc_rows:
        print(r)


if __name__ == "__main__":
    main()
