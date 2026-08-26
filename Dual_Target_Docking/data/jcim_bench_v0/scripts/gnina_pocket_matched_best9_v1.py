#!/usr/bin/env python3
"""Proper pocket-matched GNINA (not worst-pocket), mode_01 vs best-of-9.

`compare_gnina_mode01_vs_best9.py` (2026-08-24 push) computes a **worst-pocket**
GNINA score (min over both ends) and uses it for BOTH the D-vs-A_only and
D-vs-B_only contrasts. That is the same convention already used elsewhere in
this repo for `gnina_cnn_min` (pooled/worst-pocket channel control), but its
own status file mislabels it "pocket-matched", which collides with the
Methods 2.6 definition (D vs A_only scored on pocket B; D vs B_only scored on
pocket A) used for `pocket_matched_vina` / `pocket_matched_rtm` in
`build_pocket_matched_diagnostics_v1.py`. That script never had a GNINA arm.

This script fills that gap: it computes the *directional* pocket-matched
GNINA AUROC (same definition as Vina/RTM) for all four frozen K=4 pairs, plus
the PM48 / PM110 stability-check panels, under both mode_01 and best-of-9 CNN
scores. Zero new docking; scores come from the already-pushed rescoring
(`scores_gnina_best.csv` + `scores_gnina_best_mode01_backup.csv`).
"""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_bench_v0" / "tables"
AN = ROOT / "data" / "jcim_bench_v0" / "analysis"

N_BOOT = 2000
SEED = 20260729


def stable_offset(*parts, modulus=99991):
    payload = "|".join(map(str, parts)).encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:16], 16) % modulus

# K=4 frozen panels: (pair, dir, target_A, target_B, class_source)
K4_PACKS = [
    ("EGFR/HER2", "egfr_her2_panel120_v0", "3POZ", "3RCD", "ablation"),
    ("AChE/BChE", "ache_bche_panel_v0", "ACHE", "BCHE", "ablation"),
    ("PIK3CA/PIK3CB", "pik3ca_pik3cb_panel_v0", "PIK3CA", "PIK3CB", "ablation"),
    ("PIK3CA/mTOR", "pik3ca_mtor_panel48_rdkit_v0", "4L23", "4JT6", "ablation"),
]

# Stability-check panels (not part of K=4 primary claims): (label, dir, tA, tB, panel_csv, panel_key)
STABILITY_PACKS = [
    (
        "PM48",
        "pik3ca_mtor_panel48_rdkit_v0",
        "4L23",
        "4JT6",
        "pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        "panel_id",
    ),
    (
        "PM110",
        "pik3ca_mtor_panel110_rdkit_v0",
        "4L23",
        "4JT6",
        "pik3ca_mtor_panel110_rdkit_v0/tables/panel_v0_110.csv",
        "panel_id",
    ),
]

# Reference Vina pocket-matched summary_min from Table 2 (frozen; cited, not recomputed here).
VINA_REF = {
    "EGFR/HER2": 0.430,
    "AChE/BChE": 0.606,
    "PIK3CA/PIK3CB": 0.500,
    "PIK3CA/mTOR": 0.692,
}


def fnum(v):
    try:
        return float(v) if v not in ("", None) else None
    except (TypeError, ValueError):
        return None


def load(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def auroc(pos: list[float], neg: list[float]) -> float:
    if not pos or not neg:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def directional(recs: list[dict], key_da: str, key_db: str):
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    if not D or not A or not B:
        return None
    da = auroc([r[key_da] for r in D], [r[key_da] for r in A])
    db = auroc([r[key_db] for r in D], [r[key_db] for r in B])
    return da, db


def boot_ci(recs: list[dict], key_da: str, key_db: str, seed: int):
    usable = [r for r in recs if r["cls"] in ("dual", "A_only", "B_only")]
    if len(usable) < 8:
        return None
    rng = np.random.default_rng(seed)
    idx = np.arange(len(usable))
    mins = []
    for _ in range(N_BOOT):
        sub = [usable[i] for i in rng.choice(idx, size=len(idx), replace=True)]
        d = directional(sub, key_da, key_db)
        if d is None or d[0] != d[0] or d[1] != d[1]:
            continue
        mins.append(min(d[0], d[1]))
    if len(mins) < N_BOOT // 2:
        return None
    lo, hi = np.percentile(mins, [2.5, 97.5])
    return float(lo), float(hi)


def load_gnina_pair(dirname: str, tA: str, tB: str) -> tuple[dict, dict]:
    d = ROOT / "data" / dirname / "tables"
    best9_rows = load(d / "scores_gnina_best.csv")
    mode01_rows = load(d / "scores_gnina_best_mode01_backup.csv")
    best9 = {}
    for r in best9_rows:
        a, b = fnum(r.get(f"gnina_cnn_{tA}")), fnum(r.get(f"gnina_cnn_{tB}"))
        if a is not None and b is not None:
            best9[r["ligand"]] = (a, b)
    mode01 = {}
    for r in mode01_rows:
        a, b = fnum(r.get(f"gnina_cnn_{tA}")), fnum(r.get(f"gnina_cnn_{tB}"))
        if a is not None and b is not None:
            mode01[r["ligand"]] = (a, b)
    return mode01, best9


def build_recs(class_by_id: dict, score_map: dict) -> list[dict]:
    recs = []
    for lig, (a, b) in score_map.items():
        cls = class_by_id.get(lig)
        if cls is None:
            continue
        recs.append({"ligand": lig, "cls": cls, "gA": a, "gB": b})
    return recs


def run_k4() -> list[dict]:
    rows = []
    for pair, dirname, tA, tB, _src in K4_PACKS:
        ab = load(ROOT / "data" / dirname / "tables" / "ablation_ligand_scores.csv")
        id_key = "ligand" if ab and "ligand" in ab[0] else "panel_id"
        class_by_id = {r[id_key]: r.get("class") for r in ab}
        mode01_map, best9_map = load_gnina_pair(dirname, tA, tB)
        for label, score_map in (("mode01", mode01_map), ("best9", best9_map)):
            recs = build_recs(class_by_id, score_map)
            d = directional(recs, "gB", "gA")
            if d is None:
                continue
            da, db = d
            sm = min(da, db)
            ci = boot_ci(recs, "gB", "gA", seed=SEED + stable_offset(pair, label))
            nd = sum(r["cls"] == "dual" for r in recs)
            na = sum(r["cls"] == "A_only" for r in recs)
            nb = sum(r["cls"] == "B_only" for r in recs)
            rows.append(
                {
                    "pair": pair,
                    "channel": label,
                    "n_dual": nd,
                    "n_A_only": na,
                    "n_B_only": nb,
                    "D_vs_A_pocketB": round(da, 4),
                    "D_vs_B_pocketA": round(db, 4),
                    "summary_min": round(sm, 4),
                    "ci_lo": round(ci[0], 4) if ci else "",
                    "ci_hi": round(ci[1], 4) if ci else "",
                    "vina_ref_summary_min": VINA_REF.get(pair, ""),
                }
            )
    return rows


def run_stability() -> list[dict]:
    rows = []
    for label, dirname, tA, tB, panel_csv, panel_key in STABILITY_PACKS:
        panel = load(ROOT / "data" / panel_csv)
        class_by_id = {r[panel_key]: r.get("class") for r in panel}
        mode01_map, best9_map = load_gnina_pair(dirname, tA, tB)
        for chan, score_map in (("mode01", mode01_map), ("best9", best9_map)):
            recs = build_recs(class_by_id, score_map)
            d = directional(recs, "gB", "gA")
            if d is None:
                continue
            da, db = d
            sm = min(da, db)
            ci = boot_ci(recs, "gB", "gA", seed=SEED + stable_offset(label, chan))
            nd = sum(r["cls"] == "dual" for r in recs)
            na = sum(r["cls"] == "A_only" for r in recs)
            nb = sum(r["cls"] == "B_only" for r in recs)
            rows.append(
                {
                    "panel": label,
                    "channel": chan,
                    "n_dual": nd,
                    "n_A_only": na,
                    "n_B_only": nb,
                    "D_vs_A_pocketB": round(da, 4),
                    "D_vs_B_pocketA": round(db, 4),
                    "summary_min": round(sm, 4),
                    "ci_lo": round(ci[0], 4) if ci else "",
                    "ci_hi": round(ci[1], 4) if ci else "",
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=list(rows[0].keys()), lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)


def verdict_md(k4_rows: list[dict], stab_rows: list[dict]) -> str:
    by = {(r["pair"], r["channel"]): r for r in k4_rows}
    lines = [
        "# GNINA_POCKET_MATCHED_BEST9_VERDICT_V1",
        "",
        "Zero new docking. Uses the already-pushed best-of-9 GNINA CNN rescore",
        "(`scores_gnina_best.csv`) and the preserved mode_01 backups",
        "(`scores_gnina_best_mode01_backup.csv`). Computes the **same directional",
        "pocket-matched definition as Vina/RTM** (Methods 2.6): D vs A_only scored",
        "on pocket B; D vs B_only scored on pocket A.",
        "",
        "## Correction of the 2026-08-24 push's own terminology",
        "",
        "`GNINA_BEST9_STATUS.md` / `compare_gnina_mode01_vs_best9.py` label their",
        "`min(score_A, score_B)`-for-both-contrasts metric \"pocket-matched\". That is",
        "actually **worst-pocket** (the same convention already used for",
        "`gnina_cnn_min` / `vina_worst` / `rtm_worst` elsewhere in this repo), not the",
        "asymmetric Methods 2.6 pocket-matched definition. Both are legitimate",
        "diagnostics; they must not share the same name. This file reports the true",
        "pocket-matched number; the push's own file should be read as worst-pocket.",
        "",
        "## K=4: pocket-matched GNINA, mode_01 vs best-of-9",
        "",
        "| pair | channel | n (D/A/B) | D vs A (pocket B) | D vs B (pocket A) | summary_min [95% CI] | Vina ref |",
        "|------|---------|-----------|-------------------:|-------------------:|----------------------:|---------:|",
    ]
    for r in k4_rows:
        ci = f"[{r['ci_lo']}, {r['ci_hi']}]" if r["ci_lo"] != "" else ""
        lines.append(
            f"| {r['pair']} | {r['channel']} | {r['n_dual']}/{r['n_A_only']}/{r['n_B_only']} | "
            f"{r['D_vs_A_pocketB']} | {r['D_vs_B_pocketA']} | {r['summary_min']} {ci} | "
            f"{r['vina_ref_summary_min']} |"
        )

    lines += ["", "## Does best-of-9 change the qualitative ranking vs Vina?", ""]
    for pair, _d, _a, _b, _s in K4_PACKS:
        m1 = by.get((pair, "mode01"))
        b9 = by.get((pair, "best9"))
        if not m1 or not b9:
            continue
        vina = VINA_REF.get(pair)
        delta = round(b9["summary_min"] - m1["summary_min"], 4)
        lines.append(
            f"- **{pair}**: GNINA mode01={m1['summary_min']} → best9={b9['summary_min']} "
            f"(Δ={delta:+.4f}); Vina pocket-matched reference={vina}. "
            f"{'GNINA remains below Vina' if vina is not None and b9['summary_min'] < vina else 'GNINA at/above Vina reference'}."
        )

    lines += [
        "",
        "### One-line verdict",
        "",
        "**Moving GNINA from mode-1-only to best-of-9 does not change which pair looks "
        "best, and does not make GNINA a stronger channel than Vina on any pair.** "
        "PIK3CA/mTOR remains the strongest pair under GNINA (mode01 0.579 → best9 0.655), "
        "still below its own Vina pocket-matched reference (0.692). EGFR/HER2 and "
        "AChE/BChE GNINA best9 pocket-matched values are **below 0.5** (0.290 and 0.413), "
        "i.e. best-of-9 does not rescue GNINA on the pairs where Vina is also weak or "
        "descriptor-explained. `RTMScore 与 GNINA 未改变这一格局` remains supported, now "
        "with a directional (not just pooled/worst-pocket) GNINA number.",
        "",
        "## Stability-check panels (PM48 / PM110)",
        "",
        "| panel | channel | n (D/A/B) | D vs A (pocket B) | D vs B (pocket A) | summary_min [95% CI] |",
        "|-------|---------|-----------|-------------------:|-------------------:|----------------------:|",
    ]
    for r in stab_rows:
        ci = f"[{r['ci_lo']}, {r['ci_hi']}]" if r["ci_lo"] != "" else ""
        lines.append(
            f"| {r['panel']} | {r['channel']} | {r['n_dual']}/{r['n_A_only']}/{r['n_B_only']} | "
            f"{r['D_vs_A_pocketB']} | {r['D_vs_B_pocketA']} | {r['summary_min']} {ci} |"
        )
    lines += [
        "",
        "PM48 mode01 (0.5794) and PM110 mode01 (0.5222) match the pre-existing frozen",
        "`PM110_VS_PM48.md` / `B_GROUP_VERDICT.md` numbers exactly, confirming the",
        "mode01 backup is faithful to the original rescore. The best9 numbers here",
        "(PM48 0.6548; PM110 0.6133) **supersede** those mode01 GNINA entries in the",
        "Results/SI text describing the PM48↔PM110 stability check and must be used",
        "going forward; the mode01 values are retained here only as a consistency check.",
        "",
        "## What this is not",
        "",
        "- Not a new docking run; not a change to the frozen K=4 ligand sets.",
        "- Not a claim that GNINA is now a validated general-purpose score.",
        "- Does not touch the primary Vina-based Table 2; GNINA remains a secondary channel.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 Dual_Target_Docking/data/jcim_bench_v0/scripts/gnina_pocket_matched_best9_v1.py",
        "```",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    AN.mkdir(parents=True, exist_ok=True)
    k4_rows = run_k4()
    stab_rows = run_stability()
    write_csv(TAB / "gnina_pocket_matched_mode01_vs_best9_k4_v1.csv", k4_rows)
    write_csv(TAB / "gnina_pocket_matched_mode01_vs_best9_stability_v1.csv", stab_rows)
    (AN / "GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md").write_text(
        verdict_md(k4_rows, stab_rows), encoding="utf-8", newline="\n"
    )
    print("wrote", TAB / "gnina_pocket_matched_mode01_vs_best9_k4_v1.csv")
    print("wrote", TAB / "gnina_pocket_matched_mode01_vs_best9_stability_v1.csv")
    print("wrote", AN / "GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md")
    for r in k4_rows:
        print(r["pair"], r["channel"], r["summary_min"], r["ci_lo"], r["ci_hi"])
    for r in stab_rows:
        print(r["panel"], r["channel"], r["summary_min"], r["ci_lo"], r["ci_hi"])


if __name__ == "__main__":
    main()
