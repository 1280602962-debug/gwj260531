#!/usr/bin/env python3
"""Pocket-matched directional diagnostics + confound controls (no new docking).

Motivation: dual vs A_only differ ONLY in measured activity at target B, so the
decision-relevant score for that contrast is the pocket-B score, not a pooled
mean over both pockets. This script therefore reports, per pair:

  1. pooled score directional AUROC (the metric used in v1 tables)
  2. pocket-matched directional AUROC (D/A via pocket B, D/B via pocket A)
  3. wrong-pocket control (D/A via pocket A, D/B via pocket B)
     -> if this is far from 0.5, the separation is molecule-level confound
        (size / potency / chemotype), not pocket-specific recognition
  4. ligand-efficiency normalised scores (score / heavy atoms)
  5. worst-pocket aggregation
  6. size-stratified (heavy-atom tertile) pocket-matched AUROC with cell counts
  7. bootstrap CIs for the pocket-matched summary and the confound gap

Outputs go to data/jcim_bench_v0/tables/ with prefix pocket_matched_.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_bench_v0" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
SEED = 20260729

# pair -> (score table, pocket-A vina col, pocket-B vina col, rtm A, rtm B,
#          smiles source table or None if smiles inline, key column)
SPEC = {
    "EGFR/HER2": dict(
        scores="data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        vina_a="3POZ_affinity",
        vina_b="3RCD_affinity",
        rtm_a="rtm_3POZ",
        rtm_b="rtm_3RCD",
        panel="data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        panel_key="panel_id",
    ),
    "AChE/BChE": dict(
        scores="data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        vina_a="vina_ACHE",
        vina_b="vina_BCHE",
        rtm_a="rtm_ACHE",
        rtm_b="rtm_BCHE",
        panel=None,
        panel_key=None,
    ),
    "PIK3CA/PIK3CB": dict(
        scores="data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        vina_a="vina_PIK3CA",
        vina_b="vina_PIK3CB",
        rtm_a="rtm_PIK3CA",
        rtm_b="rtm_PIK3CB",
        panel=None,
        panel_key=None,
    ),
    "PIK3CA/mTOR": dict(
        scores="data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        vina_a="4L23_affinity",
        vina_b="4JT6_affinity",
        rtm_a="rtm_4L23",
        rtm_b="rtm_4JT6",
        panel="data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        panel_key="panel_id",
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


def assemble(pair: str, cfg: dict) -> list[dict]:
    rows = load_csv(ROOT / cfg["scores"])
    smimap = {}
    if cfg["panel"]:
        for r in load_csv(ROOT / cfg["panel"]):
            smimap[r[cfg["panel_key"]]] = r.get("smiles")
    out = []
    for r in rows:
        a, b = fnum(r.get(cfg["vina_a"])), fnum(r.get(cfg["vina_b"]))
        if a is None or b is None:
            continue
        lig = r.get("ligand") or r.get("panel_id")
        smi = r.get("smiles") or smimap.get(lig)
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        ha = mol.GetNumHeavyAtoms()
        ra, rb = fnum(r.get(cfg["rtm_a"])), fnum(r.get(cfg["rtm_b"]))
        # vina affinities are negative kcal/mol; flip so higher = better
        rec = {
            "pair": pair,
            "ligand": lig,
            "cls": r.get("class"),
            "heavy": float(ha),
            "vina_A": -a,
            "vina_B": -b,
            "vina_mean": -(a + b) / 2.0,
            "vina_worst": min(-a, -b),
            "le_A": -a / ha,
            "le_B": -b / ha,
            "le_mean": (-(a + b) / 2.0) / ha,
            "le_worst": min(-a, -b) / ha,
        }
        if ra is not None and rb is not None:
            rec["rtm_A"] = ra
            rec["rtm_B"] = rb
            rec["rtm_worst"] = min(ra, rb)
        out.append(rec)
    return out


def directional(recs, key_for_DA, key_for_DB):
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    da = auroc(
        [r[key_for_DA] for r in D if key_for_DA in r],
        [r[key_for_DA] for r in A if key_for_DA in r],
    )
    db = auroc(
        [r[key_for_DB] for r in D if key_for_DB in r],
        [r[key_for_DB] for r in B if key_for_DB in r],
    )
    return da, db


def boot_ci(recs, key_da, key_db, n_boot=N_BOOT, seed=SEED):
    usable = [r for r in recs if r["cls"] in ("dual", "A_only", "B_only") and key_da in r and key_db in r]
    if len(usable) < 8:
        return None
    rng = np.random.default_rng(seed)
    idx = np.arange(len(usable))
    mins = []
    for _ in range(n_boot):
        ii = rng.choice(idx, size=len(idx), replace=True)
        sub = [usable[i] for i in ii]
        da, db = directional(sub, key_da, key_db)
        if da != da or db != db:
            continue
        mins.append(min(da, db))
    if len(mins) < n_boot // 2:
        return None
    lo, hi = np.percentile(mins, [2.5, 97.5])
    return float(lo), float(hi)


def main():
    packs = {pair: assemble(pair, cfg) for pair, cfg in SPEC.items()}

    main_rows = []
    strat_rows = []
    for pair, recs in packs.items():
        n_d = sum(r["cls"] == "dual" for r in recs)
        n_a = sum(r["cls"] == "A_only" for r in recs)
        n_b = sum(r["cls"] == "B_only" for r in recs)

        variants = [
            ("pooled_vina_mean", "vina_mean", "vina_mean", "pooled score for both contrasts (v1 main metric)"),
            ("pocket_matched_vina", "vina_B", "vina_A", "D/A via pocket B; D/B via pocket A"),
            ("wrong_pocket_control_vina", "vina_A", "vina_B", "control: uninformative pocket; ~0.5 expected if no confound"),
            ("worst_pocket_vina", "vina_worst", "vina_worst", "min over pockets"),
            ("le_pocket_matched", "le_B", "le_A", "ligand-efficiency normalised (score / heavy atoms)"),
            ("le_wrong_pocket_control", "le_A", "le_B", "LE control"),
            ("le_worst_pocket", "le_worst", "le_worst", "LE worst pocket"),
            ("heavy_atoms_baseline", "heavy", "heavy", "trivial size baseline"),
        ]
        if any("rtm_A" in r for r in recs):
            variants += [
                ("pocket_matched_rtm", "rtm_B", "rtm_A", "RTM pocket-matched"),
                ("wrong_pocket_control_rtm", "rtm_A", "rtm_B", "RTM control"),
            ]

        for name, kda, kdb, note in variants:
            da, db = directional(recs, kda, kdb)
            if da != da or db != db:
                continue
            ci = boot_ci(recs, kda, kdb, seed=SEED + abs(hash((pair, name))) % 99991)
            main_rows.append(
                {
                    "pair": pair,
                    "variant": name,
                    "n_dual": n_d,
                    "n_A_only": n_a,
                    "n_B_only": n_b,
                    "auroc_D_vs_A": round(da, 4),
                    "auroc_D_vs_B": round(db, 4),
                    "summary_min": round(min(da, db), 4),
                    "summary_min_ci_lo": round(ci[0], 4) if ci else "",
                    "summary_min_ci_hi": round(ci[1], 4) if ci else "",
                    "note": note,
                }
            )

        # size-stratified pocket-matched
        has = sorted(r["heavy"] for r in recs)
        t1, t2 = has[len(has) // 3], has[2 * len(has) // 3]
        for lab, sel in (
            ("small", lambda h: h <= t1),
            ("mid", lambda h: t1 < h <= t2),
            ("large", lambda h: h > t2),
        ):
            sub = [r for r in recs if sel(r["heavy"])]
            nd = sum(r["cls"] == "dual" for r in sub)
            na = sum(r["cls"] == "A_only" for r in sub)
            nb = sum(r["cls"] == "B_only" for r in sub)
            da, db = directional(sub, "vina_B", "vina_A")
            strat_rows.append(
                {
                    "pair": pair,
                    "stratum": lab,
                    "heavy_cut_lo": t1,
                    "heavy_cut_hi": t2,
                    "n": len(sub),
                    "n_dual": nd,
                    "n_A_only": na,
                    "n_B_only": nb,
                    "underpowered": int(min(nd, na, nb) < 5),
                    "auroc_D_vs_A": round(da, 4) if da == da else "",
                    "auroc_D_vs_B": round(db, 4) if db == db else "",
                    "summary_min": round(min(da, db), 4) if da == da and db == db else "",
                }
            )

    write_csv(TAB / "pocket_matched_directional_v1.csv", main_rows)
    write_csv(TAB / "pocket_matched_size_strata_v1.csv", strat_rows)

    # confound gap summary: pocket-matched vs wrong-pocket
    gaps = []
    for pair in packs:
        pm = next((r for r in main_rows if r["pair"] == pair and r["variant"] == "pocket_matched_vina"), None)
        wp = next((r for r in main_rows if r["pair"] == pair and r["variant"] == "wrong_pocket_control_vina"), None)
        le = next((r for r in main_rows if r["pair"] == pair and r["variant"] == "le_pocket_matched"), None)
        hv = next((r for r in main_rows if r["pair"] == pair and r["variant"] == "heavy_atoms_baseline"), None)
        if not (pm and wp and le and hv):
            continue
        gaps.append(
            {
                "pair": pair,
                "pocket_matched_min": pm["summary_min"],
                "wrong_pocket_min": wp["summary_min"],
                "specificity_gap": round(pm["summary_min"] - wp["summary_min"], 4),
                "wrong_pocket_dev_from_0p5": round(abs(wp["summary_min"] - 0.5), 4),
                "le_normalised_min": le["summary_min"],
                "size_baseline_min": hv["summary_min"],
                "survives_le_normalisation": bool(le["summary_min"] > hv["summary_min"]),
            }
        )
    write_csv(TAB / "pocket_specificity_gap_v1.csv", gaps)

    meta = {
        "rationale": "dual vs A_only is defined only by activity at B; the decision-relevant score is the pocket-B score",
        "n_boot": N_BOOT,
        "seed": SEED,
        "gaps": gaps,
    }
    (TAB / "pocket_matched_meta_v1.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
