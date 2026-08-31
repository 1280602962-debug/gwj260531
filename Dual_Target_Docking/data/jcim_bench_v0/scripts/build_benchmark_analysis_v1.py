#!/usr/bin/env python3
"""Build the DualFourClass-Bench analysis pack (no new docking).

Produces, under data/jcim_bench_v0/:
  - assembled per-pair ligand tables (scores + physchem + GNINA + pChEMBL)
  - directional AUROC point estimates vs trivial baselines
  - ligand-bootstrap 95% CIs (B=2000) for summary_min and Δ(dock−baseline)
  - Top10 hard-negative composition (+ bootstrap CI on hardneg count)
  - continuous Spearman (score vs min_pChEMBL)
  - threshold sensitivity (θ ∈ {5.5, 6.0, 6.5}) where both-end pChEMBL exist
  - AChE/BChE TPSA confound diagnostics + all-pair descriptor-by-class
  - PIK3CA/mTOR LigPrep vs RDKit prep sensitivity copy
  - pooled-vs-directional contrast + asymmetry table

Exploration / evaluation pool only. Does not invent decision arms.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors


def stable_offset(*parts, modulus=100000):
    payload = "|".join(map(str, parts)).encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:16], 16) % modulus

try:
    from scipy.stats import spearmanr
except ImportError:  # pragma: no cover
    spearmanr = None

RDLogger.DisableLog("rdApp.*")

# scripts/ -> jcim_bench_v0/ -> data/ -> Dual_Target_Docking/
ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_bench_v0"
TAB = OUT / "tables"
AN = OUT / "analysis"
TAB.mkdir(parents=True, exist_ok=True)
AN.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
SEED = 20260729
ARMS = ["vina_mean", "rtm_min_z", "gnina_cnn_min", "heavy", "mw", "clogp", "tpsa"]
DOCK = {"vina_mean", "rtm_min_z", "gnina_cnn_min"}
BASE = {"heavy", "mw", "clogp", "tpsa"}
CUTOFFS = (5.5, 6.0, 6.5)


def fnum(x):
    try:
        if x is None or x == "":
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def load_csv(p: Path):
    with p.open() as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None):
    if not rows:
        path.write_text("")
        return
    fields = fields or list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
        )
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def auroc(pos: list[float], neg: list[float]) -> float:
    if not pos or not neg:
        return float("nan")
    pos_a = np.asarray(pos, dtype=float)
    neg_a = np.asarray(neg, dtype=float)
    diff = pos_a[:, None] - neg_a[None, :]
    wins = (diff > 0).sum() + 0.5 * (diff == 0).sum()
    return float(wins / (len(pos_a) * len(neg_a)))


def assign_fourclass(pA, pB, cut: float) -> str | None:
    if pA is None or pB is None or (isinstance(pA, float) and math.isnan(pA)) or (
        isinstance(pB, float) and math.isnan(pB)
    ):
        return None
    a_act = pA >= cut
    b_act = pB >= cut
    if a_act and b_act:
        return "dual"
    if a_act and not b_act:
        return "A_only"
    if b_act and not a_act:
        return "B_only"
    return "neither"


def phys(smi: str) -> dict | None:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return {
        "heavy": float(mol.GetNumHeavyAtoms()),
        "mw": float(Descriptors.MolWt(mol)),
        "clogp": float(Descriptors.MolLogP(mol)),
        "tpsa": float(Descriptors.TPSA(mol)),
    }


def gnina_min(row: dict) -> float | None:
    vals = []
    for k, v in row.items():
        if k.startswith("gnina_cnn_") and not k.endswith("_min"):
            x = fnum(v)
            if x is not None:
                vals.append(x)
    return min(vals) if vals else None


def directional_point(scores: np.ndarray, labels: np.ndarray) -> tuple[float, float, float] | None:
    pos = scores[labels == "dual"]
    a = scores[labels == "A_only"]
    b = scores[labels == "B_only"]
    if len(pos) == 0 or len(a) == 0 or len(b) == 0:
        return None
    da, db = auroc(pos.tolist(), a.tolist()), auroc(pos.tolist(), b.tolist())
    return da, db, min(da, db)


def pooled_auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    pos = scores[labels == "dual"]
    neg = scores[np.isin(labels, ["A_only", "B_only"])]
    return auroc(pos.tolist(), neg.tolist())


def bootstrap_directional(rows: list[dict], arm: str, n_boot: int = N_BOOT, seed: int = SEED):
    usable = [
        r
        for r in rows
        if arm in r and r.get(arm) is not None and r.get("cls") in ("dual", "A_only", "B_only")
    ]
    if len(usable) < 8:
        return None
    y = np.array([r["cls"] for r in usable])
    s = np.array([r[arm] for r in usable], dtype=float)
    idx = np.arange(len(usable))
    rng = np.random.default_rng(seed)

    point = directional_point(s, y)
    if point is None:
        return None
    das, dbs, mns, pools = [], [], [], []
    for _ in range(n_boot):
        ii = rng.choice(idx, size=len(idx), replace=True)
        m = directional_point(s[ii], y[ii])
        if m is None:
            continue
        das.append(m[0])
        dbs.append(m[1])
        mns.append(m[2])
        pools.append(pooled_auroc(s[ii], y[ii]))
    if len(mns) < n_boot // 2:
        return None

    def ci(xs):
        lo, hi = np.percentile(xs, [2.5, 97.5])
        return float(np.mean(xs)), float(lo), float(hi)

    da_m, da_lo, da_hi = ci(das)
    db_m, db_lo, db_hi = ci(dbs)
    mn_m, mn_lo, mn_hi = ci(mns)
    po_m, po_lo, po_hi = ci(pools)
    return {
        "n": len(usable),
        "n_dual": int((y == "dual").sum()),
        "n_A_only": int((y == "A_only").sum()),
        "n_B_only": int((y == "B_only").sum()),
        "n_boot_ok": len(mns),
        "auroc_D_vs_A": point[0],
        "auroc_D_vs_A_boot_mean": da_m,
        "auroc_D_vs_A_ci_lo": da_lo,
        "auroc_D_vs_A_ci_hi": da_hi,
        "auroc_D_vs_B": point[1],
        "auroc_D_vs_B_boot_mean": db_m,
        "auroc_D_vs_B_ci_lo": db_lo,
        "auroc_D_vs_B_ci_hi": db_hi,
        "summary_min": point[2],
        "summary_min_boot_mean": mn_m,
        "summary_min_ci_lo": mn_lo,
        "summary_min_ci_hi": mn_hi,
        "auroc_pooled": pooled_auroc(s, y),
        "auroc_pooled_boot_mean": po_m,
        "auroc_pooled_ci_lo": po_lo,
        "auroc_pooled_ci_hi": po_hi,
    }


def top10_hardneg(rows: list[dict], arm: str, n_boot: int = N_BOOT, seed: int = SEED):
    usable = [
        r
        for r in rows
        if arm in r and r.get(arm) is not None and r.get("cls") in ("dual", "A_only", "B_only")
    ]
    if len(usable) < 10:
        return None
    y = np.array([r["cls"] for r in usable])
    s = np.array([r[arm] for r in usable], dtype=float)
    # higher score = better rank
    order = np.argsort(-s)[:10]
    labs = y[order]
    point = {
        "n_top10_dual": int((labs == "dual").sum()),
        "n_top10_A_only": int((labs == "A_only").sum()),
        "n_top10_B_only": int((labs == "B_only").sum()),
        "n_top10_hardneg": int(np.isin(labs, ["A_only", "B_only"]).sum()),
    }
    rng = np.random.default_rng(seed)
    idx = np.arange(len(usable))
    hn_counts = []
    for _ in range(n_boot):
        ii = rng.choice(idx, size=len(idx), replace=True)
        order_b = np.argsort(-s[ii])[:10]
        labs_b = y[ii][order_b]
        hn_counts.append(int(np.isin(labs_b, ["A_only", "B_only"]).sum()))
    lo, hi = np.percentile(hn_counts, [2.5, 97.5])
    return {
        **point,
        "n_top10_hardneg_boot_mean": float(np.mean(hn_counts)),
        "n_top10_hardneg_ci_lo": float(lo),
        "n_top10_hardneg_ci_hi": float(hi),
        "n": len(usable),
    }


def assemble_egfr() -> list[dict]:
    panel = {r["panel_id"]: r for r in load_csv(ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv")}
    mixed = load_csv(ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv")
    reprep = {
        r["ligand"]: r
        for r in load_csv(ROOT / "data/egfr_her2_panel40_reprep_rdkit_v0/tables/ablation_ligand_scores.csv")
    }
    gn = {
        r["ligand"]: r
        for r in load_csv(ROOT / "data/egfr_her2_panel120_v0/tables/scores_gnina_best.csv")
    }
    rows = []
    for r in mixed:
        src = reprep[r["ligand"]] if r.get("from_panel40") == "yes" and r["ligand"] in reprep else r
        meta = panel.get(r["ligand"], {})
        smi = meta.get("smiles")
        if not smi:
            continue
        ph = phys(smi)
        if not ph:
            continue
        rec = {
            "pair": "EGFR/HER2",
            "ligand": r["ligand"],
            "cls": r["class"],
            "smiles": smi,
            "label_rule": "theta_6_panel120",
            "prep": "rdkit_meeko_unified",
            "pA": fnum(meta.get("pchembl_EGFR")),
            "pB": fnum(meta.get("pchembl_HER2")),
            "min_pchembl": fnum(meta.get("min_pchembl")),
            **ph,
        }
        for k in ("vina_mean", "vina_min", "rtm_min", "rtm_mean"):
            v = fnum(src.get(k))
            if v is not None:
                rec[k] = v
        ra, rb = fnum(src.get("rtm_3POZ")), fnum(src.get("rtm_3RCD"))
        if ra is not None and rb is not None:
            rec["_rtm_A"], rec["_rtm_B"] = ra, rb
        g = gn.get(r["ligand"])
        if g:
            gm = gnina_min(g)
            if gm is not None:
                rec["gnina_cnn_min"] = gm
        if "vina_mean" not in rec:
            continue
        rows.append(rec)

    a_vals = [r["_rtm_A"] for r in rows if "_rtm_A" in r]
    b_vals = [r["_rtm_B"] for r in rows if "_rtm_B" in r]
    if a_vals and b_vals:
        mu_a, sd_a = st.mean(a_vals), st.pstdev(a_vals)
        mu_b, sd_b = st.mean(b_vals), st.pstdev(b_vals)
        for r in rows:
            if "_rtm_A" in r and "_rtm_B" in r and sd_a > 0 and sd_b > 0:
                r["rtm_min_z"] = min((r["_rtm_A"] - mu_a) / sd_a, (r["_rtm_B"] - mu_b) / sd_b)
    for r in rows:
        r.pop("_rtm_A", None)
        r.pop("_rtm_B", None)
    return rows


def assemble_generic(
    pair: str,
    score_csv: Path,
    panel_csv: Path,
    gnina_csv: Path,
    label_rule: str,
    pA_key: str,
    pB_key: str,
) -> list[dict]:
    panel = {}
    for r in load_csv(panel_csv):
        for k in ("ligand", "panel_id"):
            if r.get(k):
                panel[r[k]] = r
    gn = {r["ligand"]: r for r in load_csv(gnina_csv)} if gnina_csv.exists() else {}
    out = []
    for r in load_csv(score_csv):
        lig = r.get("ligand") or r.get("panel_id")
        meta = panel.get(lig, {})
        smi = r.get("smiles") or meta.get("smiles")
        if not smi:
            continue
        ph = phys(smi)
        if not ph:
            continue
        cls = r.get("class") or meta.get("class")
        pA = fnum(r.get(pA_key) or meta.get(pA_key))
        pB = fnum(r.get(pB_key) or meta.get(pB_key))
        mn = fnum(r.get("min_pchembl") or meta.get("min_pchembl"))
        if mn is None and pA is not None and pB is not None:
            mn = min(pA, pB)
        rec = {
            "pair": pair,
            "ligand": lig,
            "cls": cls,
            "smiles": smi,
            "label_rule": label_rule,
            "prep": r.get("prep") or "rdkit_meeko",
            "pA": pA,
            "pB": pB,
            "min_pchembl": mn,
            **ph,
        }
        for k in ("vina_mean", "vina_min", "rtm_min", "rtm_mean", "rtm_min_z"):
            v = fnum(r.get(k))
            if v is not None:
                rec[k] = v
        if "rtm_min_z" not in rec:
            zs = [
                fnum(v)
                for k, v in r.items()
                if k.startswith("rtm_") and k.endswith("_z") and k != "rtm_min_z"
            ]
            zs = [z for z in zs if z is not None]
            if len(zs) >= 2:
                rec["rtm_min_z"] = min(zs)
        g = gn.get(lig)
        if g:
            gm = gnina_min(g)
            if gm is not None:
                rec["gnina_cnn_min"] = gm
        if "vina_mean" not in rec:
            continue
        out.append(rec)
    return out


def descriptor_by_class(pair: str, rows: list[dict]) -> list[dict]:
    by = defaultdict(list)
    for r in rows:
        by[r["cls"]].append(r)
    summary = []
    for cls, rs in sorted(by.items()):
        for arm in ("tpsa", "clogp", "heavy", "mw", "vina_mean", "rtm_min_z", "gnina_cnn_min", "min_pchembl"):
            vals = [r[arm] for r in rs if arm in r and r[arm] is not None]
            if not vals:
                continue
            summary.append(
                {
                    "pair": pair,
                    "class": cls,
                    "feature": arm,
                    "n": len(vals),
                    "mean": round(st.mean(vals), 4),
                    "median": round(st.median(vals), 4),
                    "stdev": round(st.pstdev(vals), 4) if len(vals) > 1 else 0.0,
                }
            )
    return summary


def ache_feature_vs_hardneg(rows: list[dict]) -> list[dict]:
    dual = [r for r in rows if r["cls"] == "dual"]
    rest = [r for r in rows if r["cls"] in ("A_only", "B_only")]
    corr_rows = []
    for arm in ("tpsa", "clogp", "heavy", "mw", "vina_mean", "rtm_min_z", "gnina_cnn_min"):
        pos = [r[arm] for r in dual if arm in r]
        neg = [r[arm] for r in rest if arm in r]
        if len(pos) < 5 or len(neg) < 5:
            continue
        corr_rows.append(
            {
                "pair": "AChE/BChE",
                "feature": arm,
                "auroc_dual_vs_hardneg_pooled": round(auroc(pos, neg), 4),
                "mean_dual": round(st.mean(pos), 4),
                "mean_hardneg": round(st.mean(neg), 4),
                "delta_mean": round(st.mean(pos) - st.mean(neg), 4),
            }
        )
    return corr_rows


def continuous_spearman(packs: dict[str, list[dict]]) -> list[dict]:
    if spearmanr is None:
        return []
    out = []
    for name, rows in packs.items():
        both = [r for r in rows if r.get("min_pchembl") is not None]
        for arm in ARMS:
            xs, ys = [], []
            for r in both:
                if arm in r and r[arm] is not None:
                    xs.append(r[arm])
                    ys.append(r["min_pchembl"])
            if len(xs) < 8:
                continue
            rho, pval = spearmanr(xs, ys)
            out.append(
                {
                    "pair": name,
                    "scope": "all_with_min_pchembl",
                    "arm": arm,
                    "n": len(xs),
                    "spearman_rho": round(float(rho), 4) if rho == rho else "",
                    "pvalue": float(pval) if pval == pval else "",
                }
            )
            for cls in ("dual", "A_only", "B_only"):
                xs2, ys2 = [], []
                for r in both:
                    if r["cls"] != cls or arm not in r or r[arm] is None:
                        continue
                    xs2.append(r[arm])
                    ys2.append(r["min_pchembl"])
                if len(xs2) < 5:
                    continue
                rho2, p2 = spearmanr(xs2, ys2)
                out.append(
                    {
                        "pair": name,
                        "scope": f"class_{cls}",
                        "arm": arm,
                        "n": len(xs2),
                        "spearman_rho": round(float(rho2), 4),
                        "pvalue": float(p2),
                    }
                )
    return out


def threshold_sensitivity(packs: dict[str, list[dict]]) -> list[dict]:
    out = []
    for name, rows in packs.items():
        usable = [r for r in rows if r.get("pA") is not None and r.get("pB") is not None]
        if len(usable) < 20:
            continue
        for cut in CUTOFFS:
            labeled = []
            for r in usable:
                lab = assign_fourclass(r["pA"], r["pB"], cut)
                if lab is None:
                    continue
                labeled.append({**r, "cls": lab})
            counts = {
                c: sum(1 for r in labeled if r["cls"] == c)
                for c in ("dual", "A_only", "B_only", "neither")
            }
            for arm in ("vina_mean", "rtm_min_z", "gnina_cnn_min", "heavy", "tpsa", "clogp"):
                sub = [
                    r
                    for r in labeled
                    if arm in r and r["cls"] in ("dual", "A_only", "B_only")
                ]
                if len(sub) < 8:
                    continue
                y = np.array([r["cls"] for r in sub])
                s = np.array([r[arm] for r in sub], dtype=float)
                m = directional_point(s, y)
                if m is None:
                    continue
                out.append(
                    {
                        "pair": name,
                        "cutoff": cut,
                        "arm": arm,
                        "n_dual": counts["dual"],
                        "n_A_only": counts["A_only"],
                        "n_B_only": counts["B_only"],
                        "n_neither": counts["neither"],
                        "underpowered": int(min(counts["dual"], counts["A_only"], counts["B_only"]) < 8),
                        "auroc_D_vs_A": round(m[0], 4),
                        "auroc_D_vs_B": round(m[1], 4),
                        "summary_min": round(m[2], 4),
                        "auroc_pooled": round(pooled_auroc(s, y), 4),
                    }
                )
    return out


def prep_sensitivity_table() -> list[dict]:
    path = ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/directional_by_prep.csv"
    if not path.exists():
        return []
    return load_csv(path)


def round_floats(rows: list[dict], nd: int = 4):
    for r in rows:
        for k, v in list(r.items()):
            if isinstance(v, float):
                r[k] = round(v, nd)


def main():
    packs = {
        "EGFR/HER2": assemble_egfr(),
        "AChE/BChE": assemble_generic(
            "AChE/BChE",
            ROOT / "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
            ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict_with_smiles.csv",
            ROOT / "data/ache_bche_panel_v0/tables/scores_gnina_best.csv",
            "strict_6.5_5.5",
            "pchembl_ACHE",
            "pchembl_BCHE",
        ),
        "PIK3CA/PIK3CB": assemble_generic(
            "PIK3CA/PIK3CB",
            ROOT / "data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
            ROOT / "data/pik3ca_pik3cb_panel_v0/tables/panel_v0_strict_with_smiles.csv",
            ROOT / "data/pik3ca_pik3cb_panel_v0/tables/scores_gnina_best.csv",
            "strict_6.5_5.5",
            "pchembl_PIK3CA",
            "pchembl_PIK3CB",
        ),
        "PIK3CA/mTOR": assemble_generic(
            "PIK3CA/mTOR",
            ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
            ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
            ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/scores_gnina_best.csv",
            "theta_6_panel48",
            "pchembl_PIK3CA",
            "pchembl_MTOR",
        ),
    }

    long_fields = [
        "pair",
        "ligand",
        "cls",
        "smiles",
        "label_rule",
        "prep",
        "pA",
        "pB",
        "min_pchembl",
        "heavy",
        "mw",
        "clogp",
        "tpsa",
        "vina_mean",
        "vina_min",
        "rtm_min",
        "rtm_mean",
        "rtm_min_z",
        "gnina_cnn_min",
    ]
    long_rows = []
    for name, rows in packs.items():
        for r in rows:
            long_rows.append(r)
        write_csv(TAB / f"assembled_{name.replace('/', '_')}.csv", rows, long_fields)
    write_csv(TAB / "assembled_all_pairs_long.csv", long_rows, long_fields)

    inv = []
    for name, rows in packs.items():
        inv.append(
            {
                "pair": name,
                "n": len(rows),
                "n_dual": sum(r["cls"] == "dual" for r in rows),
                "n_A_only": sum(r["cls"] == "A_only" for r in rows),
                "n_B_only": sum(r["cls"] == "B_only" for r in rows),
                "n_neither": sum(r["cls"] == "neither" for r in rows),
                "n_with_gnina": sum("gnina_cnn_min" in r for r in rows),
                "n_with_pchembl_both": sum(r.get("pA") is not None and r.get("pB") is not None for r in rows),
                "label_rule": rows[0]["label_rule"] if rows else "",
            }
        )
    write_csv(TAB / "inventory_v1.csv", inv)

    boot_rows = []
    for name, rows in packs.items():
        for arm in ARMS:
            print(f"boot {name} {arm} ...", flush=True)
            m = bootstrap_directional(rows, arm, seed=SEED + stable_offset(name, arm))
            if not m:
                continue
            boot_rows.append(
                {"pair": name, "arm": arm, "family": "docking" if arm in DOCK else "baseline", **m}
            )
    round_floats(boot_rows)
    write_csv(TAB / "bootstrap_directional_ci_v1.csv", boot_rows)

    # asymmetry + pooled contrast
    asym = []
    for r in boot_rows:
        asym.append(
            {
                "pair": r["pair"],
                "arm": r["arm"],
                "family": r["family"],
                "auroc_D_vs_A": r["auroc_D_vs_A"],
                "auroc_D_vs_B": r["auroc_D_vs_B"],
                "delta_DA_minus_DB": round(r["auroc_D_vs_A"] - r["auroc_D_vs_B"], 4),
                "summary_min": r["summary_min"],
                "auroc_pooled": r["auroc_pooled"],
                "pooled_minus_summary_min": round(r["auroc_pooled"] - r["summary_min"], 4),
                "ci_excludes_0p5_on_weak_arm": bool(
                    (r["auroc_D_vs_A"] < r["auroc_D_vs_B"] and r["auroc_D_vs_A_ci_hi"] < 0.5)
                    or (r["auroc_D_vs_B"] <= r["auroc_D_vs_A"] and r["auroc_D_vs_B_ci_hi"] < 0.5)
                ),
            }
        )
    write_csv(TAB / "asymmetry_pooled_vs_directional_v1.csv", asym)

    # baseline gate with joint bootstrap Δ
    gate = []
    for name, rows in packs.items():
        sub = [r for r in boot_rows if r["pair"] == name]
        bases = [r for r in sub if r["family"] == "baseline"]
        docks = [r for r in sub if r["family"] == "docking"]
        if not bases or not docks:
            continue
        bestb = max(bases, key=lambda r: r["summary_min"])
        for darm in ("vina_mean", "rtm_min_z", "gnina_cnn_min"):
            drow = next((x for x in docks if x["arm"] == darm), None)
            if not drow:
                continue
            pairs = [
                (r["cls"], r[darm], r[bestb["arm"]])
                for r in rows
                if darm in r and bestb["arm"] in r and r["cls"] in ("dual", "A_only", "B_only")
            ]
            if len(pairs) < 8:
                continue
            y = np.array([p[0] for p in pairs])
            sd = np.array([p[1] for p in pairs], dtype=float)
            sb = np.array([p[2] for p in pairs], dtype=float)
            rng = np.random.default_rng(SEED + stable_offset(name, darm, bestb["arm"]))
            idx = np.arange(len(pairs))
            deltas = []
            for _ in range(N_BOOT):
                ii = rng.choice(idx, size=len(idx), replace=True)
                y_b = y[ii]
                md = directional_point(sd[ii], y_b)
                mb = directional_point(sb[ii], y_b)
                if md is None or mb is None:
                    continue
                deltas.append(md[2] - mb[2])
            if len(deltas) < N_BOOT // 2:
                continue
            lo, hi = np.percentile(deltas, [2.5, 97.5])
            point = drow["summary_min"] - bestb["summary_min"]
            gate.append(
                {
                    "pair": name,
                    "dock_arm": darm,
                    "best_baseline_arm": bestb["arm"],
                    "dock_summary_min": drow["summary_min"],
                    "dock_summary_min_ci_lo": drow["summary_min_ci_lo"],
                    "dock_summary_min_ci_hi": drow["summary_min_ci_hi"],
                    "baseline_summary_min": bestb["summary_min"],
                    "baseline_summary_min_ci_lo": bestb["summary_min_ci_lo"],
                    "baseline_summary_min_ci_hi": bestb["summary_min_ci_hi"],
                    "delta_summary_min": round(point, 4),
                    "delta_boot_mean": round(float(np.mean(deltas)), 4),
                    "delta_ci_lo": round(float(lo), 4),
                    "delta_ci_hi": round(float(hi), 4),
                    "beats_baseline_point": point > 0,
                    "beats_baseline_ci_excl0": bool(lo > 0),
                    "loses_baseline_ci_excl0": bool(hi < 0),
                }
            )
    write_csv(TAB / "baseline_gate_bootstrap_v1.csv", gate)

    # Top10 hardneg
    top10_rows = []
    for name, rows in packs.items():
        for arm in ARMS:
            m = top10_hardneg(rows, arm, seed=SEED + stable_offset("t10", name, arm))
            if not m:
                continue
            top10_rows.append(
                {
                    "pair": name,
                    "arm": arm,
                    "family": "docking" if arm in DOCK else "baseline",
                    **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in m.items()},
                }
            )
    write_csv(TAB / "top10_hardneg_bootstrap_v1.csv", top10_rows)

    # continuous Spearman
    spear = continuous_spearman(packs)
    write_csv(TAB / "continuous_spearman_v1.csv", spear)

    # threshold sensitivity
    thresh = threshold_sensitivity(packs)
    write_csv(TAB / "threshold_sensitivity_v1.csv", thresh)

    # descriptor diagnostics
    desc_all = []
    for name, rows in packs.items():
        desc_all.extend(descriptor_by_class(name, rows))
    write_csv(TAB / "descriptor_by_class_v1.csv", desc_all)
    write_csv(TAB / "ache_descriptor_by_class_v1.csv", [r for r in desc_all if r["pair"] == "AChE/BChE"])
    write_csv(TAB / "ache_feature_dual_vs_hardneg_v1.csv", ache_feature_vs_hardneg(packs["AChE/BChE"]))

    prep = prep_sensitivity_table()
    write_csv(TAB / "pm48_directional_by_prep_v1.csv", prep)

    forest = []
    for r in boot_rows:
        forest.append(
            {
                "pair": r["pair"],
                "arm": r["arm"],
                "family": r["family"],
                "summary_min": r["summary_min"],
                "ci_lo": r["summary_min_ci_lo"],
                "ci_hi": r["summary_min_ci_hi"],
                "auroc_D_vs_A": r["auroc_D_vs_A"],
                "auroc_D_vs_A_ci_lo": r["auroc_D_vs_A_ci_lo"],
                "auroc_D_vs_A_ci_hi": r["auroc_D_vs_A_ci_hi"],
                "auroc_D_vs_B": r["auroc_D_vs_B"],
                "auroc_D_vs_B_ci_lo": r["auroc_D_vs_B_ci_lo"],
                "auroc_D_vs_B_ci_hi": r["auroc_D_vs_B_ci_hi"],
                "auroc_pooled": r["auroc_pooled"],
                "n": r["n"],
            }
        )
    write_csv(TAB / "forest_summary_min_ci_v1.csv", forest)

    # failure-mode one-pager JSON for manuscript
    failures = []
    for name in packs:
        g_vina = next((g for g in gate if g["pair"] == name and g["dock_arm"] == "vina_mean"), None)
        boot_v = next((b for b in boot_rows if b["pair"] == name and b["arm"] == "vina_mean"), None)
        if not g_vina or not boot_v:
            continue
        mode = []
        if boot_v["auroc_D_vs_A"] - boot_v["auroc_D_vs_B"] > 0.15:
            mode.append("directional_asymmetry_weak_B")
        elif boot_v["auroc_D_vs_B"] - boot_v["auroc_D_vs_A"] > 0.15:
            mode.append("directional_asymmetry_weak_A")
        if g_vina["loses_baseline_ci_excl0"]:
            mode.append("significantly_below_best_baseline")
        elif not g_vina["beats_baseline_point"]:
            mode.append("point_estimate_below_best_baseline")
        elif g_vina["beats_baseline_point"] and not g_vina["beats_baseline_ci_excl0"]:
            mode.append("beats_baseline_point_but_CI_includes_0")
        elif g_vina["beats_baseline_ci_excl0"]:
            mode.append("significantly_above_best_baseline")
        if boot_v["auroc_pooled"] - boot_v["summary_min"] > 0.1:
            mode.append("pooled_masks_weak_arm")
        failures.append(
            {
                "pair": name,
                "vina_summary_min": boot_v["summary_min"],
                "vina_summary_min_ci": [boot_v["summary_min_ci_lo"], boot_v["summary_min_ci_hi"]],
                "best_baseline": g_vina["best_baseline_arm"],
                "delta_vs_baseline": g_vina["delta_summary_min"],
                "delta_ci": [g_vina["delta_ci_lo"], g_vina["delta_ci_hi"]],
                "modes": mode,
            }
        )

    meta = {
        "n_boot": N_BOOT,
        "seed": SEED,
        "primary_summary": "min(AUROC dual-vs-A_only, AUROC dual-vs-B_only)",
        "pairs": list(packs.keys()),
        "inventory": inv,
        "baseline_gate_n_rows": len(gate),
        "n_arms_bootstrapped": len(boot_rows),
        "n_top10_rows": len(top10_rows),
        "n_spearman_rows": len(spear),
        "n_threshold_rows": len(thresh),
        "failure_modes": failures,
    }
    (TAB / "analysis_meta_v1.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))
    print("gate rows:")
    for g in gate:
        print(
            f"  {g['pair']:16s} {g['dock_arm']:14s} vs {g['best_baseline_arm']:6s} "
            f"Δ={g['delta_summary_min']:+.3f} CI[{g['delta_ci_lo']:+.3f},{g['delta_ci_hi']:+.3f}] "
            f"sig+={g['beats_baseline_ci_excl0']} sig-={g['loses_baseline_ci_excl0']}"
        )


if __name__ == "__main__":
    main()
