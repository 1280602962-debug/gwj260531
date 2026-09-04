#!/usr/bin/env python3
"""Holdout evaluation: pocket-matched directional AUROC + trivial baselines.

Uses the same definitions as jcim_bench_v0 pocket-matched diagnostics
(Methods 2.6 / PLAN_V2 §2.2):
  - D vs A_only scored on pocket B
  - D vs B_only scored on pocket A
  - summary_min = min(AUROC_D/A, AUROC_D/B)
  - ligand-layer bootstrap B=2000, seed=20260729
  - Vina affinities flipped so higher = better

Does NOT retune protocol or re-sample the holdout panel.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
HOLDOUT = ROOT / "data/jcim_holdout_v0"
TAB = HOLDOUT / "tables"
AN = HOLDOUT / "analysis"
TAB.mkdir(parents=True, exist_ok=True)
AN.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
SEED = 20260729

# prefix -> (pair label, receptor A name, receptor B name)
PAIRS = {
    "HOAB": ("AChE/BChE", "4EY7", "4BDS"),
    "HOAP": ("PIK3CA/PIK3CB", "4L23", "2WXF"),
    "HOPM": ("PIK3CA/mTOR", "4L23", "4JT6"),
}


def main_panel_summary_min(pair: str) -> float:
    """Table 2 pocket-matched summary_min from the canonical unified-threshold file."""
    path = ROOT / "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv"
    for row in csv.DictReader(path.open(encoding="utf-8", newline="")):
        if row["pair"] == pair and row["label_rule"] == "theta_6.0":
            return float(row["pocket_matched_summary_min"])
    raise KeyError(f"no Table 2 summary_min for {pair}")


def fnum(v):
    try:
        return float(v) if v not in ("", None) else None
    except (TypeError, ValueError):
        return None


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def directional(recs, key_da, key_db):
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    da = auroc([r[key_da] for r in D], [r[key_da] for r in A])
    db = auroc([r[key_db] for r in D], [r[key_db] for r in B])
    return da, db


def boot_ci(recs, key_da, key_db, n_boot=N_BOOT, seed=SEED):
    usable = [r for r in recs if r["cls"] in ("dual", "A_only", "B_only")]
    rng = np.random.default_rng(seed)
    idx = np.arange(len(usable))
    das, dbs, mns = [], [], []
    for _ in range(n_boot):
        sub = [usable[i] for i in rng.choice(idx, size=len(idx), replace=True)]
        da, db = directional(sub, key_da, key_db)
        if da != da or db != db:
            continue
        das.append(da)
        dbs.append(db)
        mns.append(min(da, db))
    if len(mns) < n_boot // 2:
        return None

    def ci(xs):
        lo, hi = np.percentile(xs, [2.5, 97.5])
        return float(lo), float(hi)

    return {
        "da_ci": ci(das),
        "db_ci": ci(dbs),
        "min_ci": ci(mns),
        "n_boot_ok": len(mns),
    }


def assemble(prefix: str) -> tuple[list[dict], dict]:
    pair, recA, recB = PAIRS[prefix]
    panel = {r["holdout_id"]: r for r in csv.DictReader((TAB / f"holdout_panel_{prefix}.csv").open())}
    scores = list(csv.DictReader((TAB / f"scores_vina_mode1_{prefix}.csv").open()))

    by_lig: dict[str, dict] = {}
    fail_notes = []
    for s in scores:
        lig = s["ligand"]
        row = by_lig.setdefault(lig, {"ligand": lig})
        if s["status"] != "success" or fnum(s["vina_mode1"]) is None:
            fail_notes.append(
                {"prefix": prefix, "receptor": s["receptor"], "ligand": lig, "reason": s.get("reason", "")}
            )
            continue
        if s["receptor"] == recA:
            row["vina_A_raw"] = fnum(s["vina_mode1"])
        elif s["receptor"] == recB:
            row["vina_B_raw"] = fnum(s["vina_mode1"])

    recs = []
    for lig, row in by_lig.items():
        p = panel.get(lig)
        if p is None:
            continue
        a, b = row.get("vina_A_raw"), row.get("vina_B_raw")
        if a is None or b is None:
            continue
        mol = Chem.MolFromSmiles(p["smiles"])
        if mol is None:
            continue
        ha = float(mol.GetNumHeavyAtoms())
        recs.append(
            {
                "prefix": prefix,
                "pair": pair,
                "ligand": lig,
                "chembl": p["molecule_chembl_id"],
                "cls": p["class"],
                "smiles": p["smiles"],
                "vina_A_raw": a,
                "vina_B_raw": b,
                "vina_A": -a,
                "vina_B": -b,
                "vina_mean": -(a + b) / 2.0,
                "heavy": ha,
                "mw": float(Descriptors.MolWt(mol)),
                "clogp": float(Descriptors.MolLogP(mol)),
                "tpsa": float(Descriptors.TPSA(mol)),
            }
        )
    meta = {
        "n_panel": len(panel),
        "n_score_rows": len(scores),
        "n_assembled": len(recs),
        "n_fail_rows": len(fail_notes),
        "fail_notes": fail_notes,
        "n_dual": sum(r["cls"] == "dual" for r in recs),
        "n_A_only": sum(r["cls"] == "A_only" for r in recs),
        "n_B_only": sum(r["cls"] == "B_only" for r in recs),
    }
    return recs, meta


def evaluate(prefix: str, recs: list[dict], main_summary: float) -> list[dict]:
    variants = [
        ("pocket_matched_vina", "vina_B", "vina_A", "primary: D/A via pocket B; D/B via pocket A"),
        ("wrong_pocket_control_vina", "vina_A", "vina_B", "control; ~0.5 if no molecule-level confound"),
        ("pooled_vina_mean", "vina_mean", "vina_mean", "pooled mean (legacy forest arm)"),
        ("heavy_atoms_baseline", "heavy", "heavy", "trivial size baseline"),
        ("mw_baseline", "mw", "mw", "trivial MW baseline"),
        ("clogp_baseline", "clogp", "clogp", "trivial cLogP baseline"),
        ("tpsa_baseline", "tpsa", "tpsa", "trivial TPSA baseline"),
    ]
    out = []
    for name, kda, kdb, note in variants:
        da, db = directional(recs, kda, kdb)
        if da != da or db != db:
            continue
        seed = SEED + abs(hash((prefix, name))) % 99991
        ci = boot_ci(recs, kda, kdb, seed=seed)
        sm = min(da, db)
        best_base = None
        row = {
            "prefix": prefix,
            "pair": PAIRS[prefix][0],
            "variant": name,
            "n_dual": sum(r["cls"] == "dual" for r in recs),
            "n_A_only": sum(r["cls"] == "A_only" for r in recs),
            "n_B_only": sum(r["cls"] == "B_only" for r in recs),
            "auroc_D_vs_A": round(da, 4),
            "auroc_D_vs_B": round(db, 4),
            "summary_min": round(sm, 4),
            "summary_min_ci_lo": round(ci["min_ci"][0], 4) if ci else "",
            "summary_min_ci_hi": round(ci["min_ci"][1], 4) if ci else "",
            "auroc_D_vs_A_ci_lo": round(ci["da_ci"][0], 4) if ci else "",
            "auroc_D_vs_A_ci_hi": round(ci["da_ci"][1], 4) if ci else "",
            "auroc_D_vs_B_ci_lo": round(ci["db_ci"][0], 4) if ci else "",
            "auroc_D_vs_B_ci_hi": round(ci["db_ci"][1], 4) if ci else "",
            "main_panel_pocket_matched_vina": main_summary if name == "pocket_matched_vina" else "",
            "delta_vs_main_panel": round(sm - main_summary, 4) if name == "pocket_matched_vina" else "",
            "note": note,
        }
        out.append(row)
    # gate vs best trivial baseline on same holdout
    dock = next(r for r in out if r["variant"] == "pocket_matched_vina")
    bases = [r for r in out if r["variant"].endswith("_baseline")]
    best = max(bases, key=lambda r: r["summary_min"])
    dock["best_trivial_baseline"] = best["variant"]
    dock["best_trivial_summary_min"] = best["summary_min"]
    dock["delta_vs_best_trivial"] = round(dock["summary_min"] - best["summary_min"], 4)
    return out


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    # union keys preserving order
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def write_verdict(metric_rows: list[dict], metas: dict, ligand_rows: list[dict]):
    lines = []
    lines.append("# Holdout Verdict v1 — unused-pool post-panel-freeze validation")
    lines.append("")
    lines.append("> Generated by `data/jcim_holdout_v0/scripts/analyze_holdout_v1.py`")
    lines.append("> Protocol: frozen receptors/boxes/exhaustiveness; RDKit ETKDGv3 seed=20260727 + meeko; Vina mode-1.")
    lines.append("> Metric: **pocket-matched directional AUROC** (D/A via pocket B; D/B via pocket A); `summary_min = min(arms)`; ligand bootstrap B=2000, seed=20260729.")
    lines.append("> Scope: ChEMBL unused-pool holdout (seed 20260731), **not** an independent external database.")
    lines.append("")
    lines.append("## Docking inventory")
    lines.append("")
    lines.append("| pair | prefix | panel | assembled (D/A/B) | score fails |")
    lines.append("|------|--------|-------|-------------------|-------------|")
    for pref, (pair, *_rest) in PAIRS.items():
        m = metas[pref]
        lines.append(
            f"| {pair} | {pref} | {m['n_panel']} | {m['n_assembled']} "
            f"({m['n_dual']}/{m['n_A_only']}/{m['n_B_only']}) | {m['n_fail_rows']} |"
        )
    lines.append("")
    if any(metas[p]["fail_notes"] for p in PAIRS):
        lines.append("### Recorded docking failures (not re-tuned)")
        lines.append("")
        for pref in PAIRS:
            for f in metas[pref]["fail_notes"]:
                reason = (f["reason"] or "").replace("\n", " ")[:160]
                lines.append(f"- `{pref}` `{f['receptor']}`/`{f['ligand']}`: {reason}")
        lines.append("")
        lines.append("HOAP_028 contains boron; AutoDock atom type `B` is unsupported — both pockets failed. Ligand excluded from assembled AUROC (needs both ends).")
        lines.append("")

    lines.append("## Primary metric vs main panel (pocket-matched Vina)")
    lines.append("")
    lines.append("| pair | holdout summary_min (95% CI) | D vs A | D vs B | main-panel summary_min | Δ(holdout−main) |")
    lines.append("|------|------------------------------|--------|--------|------------------------|-----------------|")
    for pref in ("HOPM", "HOAB", "HOAP"):
        r = next(x for x in metric_rows if x["prefix"] == pref and x["variant"] == "pocket_matched_vina")
        lines.append(
            f"| {r['pair']} | **{r['summary_min']:.3f}** "
            f"[{r['summary_min_ci_lo']}, {r['summary_min_ci_hi']}] | "
            f"{r['auroc_D_vs_A']:.3f} | {r['auroc_D_vs_B']:.3f} | "
            f"{r['main_panel_pocket_matched_vina']} | {r['delta_vs_main_panel']} |"
        )
    lines.append("")
    lines.append("## Baseline gate on holdout (same ligands)")
    lines.append("")
    lines.append("| pair | pocket_matched_vina | best trivial | best trivial AUROC | Δ(dock−baseline) |")
    lines.append("|------|---------------------|--------------|--------------------|------------------|")
    for pref in ("HOPM", "HOAB", "HOAP"):
        r = next(x for x in metric_rows if x["prefix"] == pref and x["variant"] == "pocket_matched_vina")
        lines.append(
            f"| {r['pair']} | {r['summary_min']:.3f} | {r.get('best_trivial_baseline','')} | "
            f"{r.get('best_trivial_summary_min','')} | {r.get('delta_vs_best_trivial','')} |"
        )
    lines.append("")
    lines.append("## Wrong-pocket control")
    lines.append("")
    lines.append("| pair | pocket_matched | wrong_pocket | gap (matched−wrong) |")
    lines.append("|------|----------------|--------------|---------------------|")
    for pref in ("HOPM", "HOAB", "HOAP"):
        m = next(x for x in metric_rows if x["prefix"] == pref and x["variant"] == "pocket_matched_vina")
        w = next(x for x in metric_rows if x["prefix"] == pref and x["variant"] == "wrong_pocket_control_vina")
        gap = round(m["summary_min"] - w["summary_min"], 4)
        lines.append(f"| {m['pair']} | {m['summary_min']:.3f} | {w['summary_min']:.3f} | {gap} |")
    lines.append("")
    lines.append("## Verdict (honest ceiling)")
    lines.append("")
    for pref in ("HOPM", "HOAB", "HOAP"):
        r = next(x for x in metric_rows if x["prefix"] == pref and x["variant"] == "pocket_matched_vina")
        delta_base = float(r.get("delta_vs_best_trivial") or 0)
        lo = float(r["summary_min_ci_lo"]) if r["summary_min_ci_lo"] != "" else float("nan")
        beats = delta_base > 0
        ci_above_half = lo > 0.5
        lines.append(f"### {r['pair']} (`{pref}`)")
        lines.append(
            f"- Holdout pocket-matched `summary_min` = **{r['summary_min']:.3f}** "
            f"(95% CI {r['summary_min_ci_lo']}–{r['summary_min_ci_hi']}); "
            f"main panel was {r['main_panel_pocket_matched_vina']} (Δ={r['delta_vs_main_panel']})."
        )
        lines.append(
            f"- vs best trivial (`{r.get('best_trivial_baseline')}`={r.get('best_trivial_summary_min')}): "
            f"Δ={r.get('delta_vs_best_trivial')} "
            f"({'beats baseline point estimate' if beats else 'does NOT beat baseline point estimate'})."
        )
        lines.append(
            f"- CI vs 0.5: lower bound {'>' if ci_above_half else '≤'} 0.5 "
            f"({'directionally supported under bootstrap' if ci_above_half else 'CI still includes/≤ chance'})."
        )
        lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `tables/holdout_ligand_scores_v1.csv` — assembled dual-end scores + physchem")
    lines.append("- `tables/holdout_pocket_matched_v1.csv` — all variants + CIs")
    lines.append("- `tables/scores_vina_mode1_{HOAB,HOAP,HOPM}.csv` — raw docking (already frozen)")
    lines.append("")
    lines.append("Do **not** re-sample holdout or change aggregation because a number looks better.")
    (AN / "HOLDOUT_VERDICT.md").write_text("\n".join(lines) + "\n")


def main():
    all_ligs = []
    all_metrics = []
    metas = {}
    for pref, (pair, _a, _b) in PAIRS.items():
        recs, meta = assemble(pref)
        metas[pref] = meta
        all_ligs.extend(recs)
        metrics = evaluate(pref, recs, main_panel_summary_min(pair))
        all_metrics.extend(metrics)
        print(
            f"{pref}: assembled={meta['n_assembled']} "
            f"D/A/B={meta['n_dual']}/{meta['n_A_only']}/{meta['n_B_only']} "
            f"fails={meta['n_fail_rows']}"
        )
        pm = next(r for r in metrics if r["variant"] == "pocket_matched_vina")
        print(
            f"  pocket_matched summary_min={pm['summary_min']} "
            f"CI=[{pm['summary_min_ci_lo']},{pm['summary_min_ci_hi']}] "
            f"Δmain={pm['delta_vs_main_panel']} Δbase={pm.get('delta_vs_best_trivial')}"
        )

    write_csv(TAB / "holdout_ligand_scores_v1.csv", all_ligs)
    write_csv(TAB / "holdout_pocket_matched_v1.csv", all_metrics)
    write_verdict(all_metrics, metas, all_ligs)
    print(f"wrote {TAB / 'holdout_ligand_scores_v1.csv'}")
    print(f"wrote {TAB / 'holdout_pocket_matched_v1.csv'}")
    print(f"wrote {AN / 'HOLDOUT_VERDICT.md'}")


if __name__ == "__main__":
    main()
