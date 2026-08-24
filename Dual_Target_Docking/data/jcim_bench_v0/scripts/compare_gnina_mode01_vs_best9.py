#!/usr/bin/env python3
"""Compare GNINA mode_01 vs best-of-9 CNN scores and directional AUROC.

Writes:
  tables/gnina_mode01_vs_best9_ligand.csv
  tables/gnina_mode01_vs_best9_auroc.csv
  analysis/GNINA_BEST9_STATUS.md
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SEED = 20260729
N_BOOT = 2000

PACKS = [
    {
        "pair": "AChE/BChE",
        "dir": "ache_bche_panel_v0",
        "targets": ("ACHE", "BCHE"),
        "panel_csv": "ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        "label_rule": "strict",
        "a_key": "pchembl_ACHE",
        "b_key": "pchembl_BCHE",
        "id_key": "ligand",
    },
    {
        "pair": "PIK3CA/PIK3CB",
        "dir": "pik3ca_pik3cb_panel_v0",
        "targets": ("PIK3CA", "PIK3CB"),
        "panel_csv": "pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        "label_rule": "strict",
        "a_key": "pchembl_PIK3CA",
        "b_key": "pchembl_PIK3CB",
        "id_key": "ligand",
    },
    {
        "pair": "PIK3CA/mTOR",
        "dir": "pik3ca_mtor_panel48_rdkit_v0",
        "targets": ("4L23", "4JT6"),
        "panel_csv": "pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        "label_rule": "theta6",
        "a_key": "pchembl_PIK3CA",
        "b_key": "pchembl_MTOR",
        "id_key": "ligand",
    },
    {
        "pair": "EGFR/HER2",
        "dir": "egfr_her2_panel120_v0",
        "targets": ("3POZ", "3RCD"),
        "panel_csv": "egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        "label_rule": "theta6",
        "a_key": "pchembl_EGFR",
        "b_key": "pchembl_HER2",
        "id_key": "ligand",
    },
]


def fnum(x):
    if x is None or x == "":
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if math.isnan(v):
        return None
    return v


def load_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def auroc(pos: list[float], neg: list[float]) -> float:
    # Mann–Whitney U / ROC AUC
    if not pos or not neg:
        return float("nan")
    scores = [(s, 1) for s in pos] + [(s, 0) for s in neg]
    scores.sort(key=lambda t: t[0])
    n_pos, n_neg = len(pos), len(neg)
    rank_sum = 0.0
    i = 0
    while i < len(scores):
        j = i
        while j < len(scores) and scores[j][0] == scores[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0  # 1-based average rank for ties
        for k in range(i, j):
            if scores[k][1] == 1:
                rank_sum += avg_rank
        i = j
    u = rank_sum - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)


def classify(a: float | None, b: float | None, rule: str) -> str | None:
    if a is None or b is None:
        return None
    if rule == "strict":
        if a >= 6.5 and b >= 6.5:
            return "dual"
        if a >= 6.5 and b <= 5.5:
            return "A_only"
        if b >= 6.5 and a <= 5.5:
            return "B_only"
        if a <= 5.5 and b <= 5.5:
            return "neither"
        return None
    # theta = 6.0
    if a >= 6.0 and b >= 6.0:
        return "dual"
    if a >= 6.0 and b < 6.0:
        return "A_only"
    if b >= 6.0 and a < 6.0:
        return "B_only"
    return "neither"


def gnina_min_from_best(row: dict, targets: tuple[str, str]) -> float | None:
    vals = []
    for t in targets:
        v = fnum(row.get(f"gnina_cnn_{t}"))
        if v is not None:
            vals.append(v)
    return min(vals) if len(vals) == 2 else None


def directional(rows: list[dict], key: str) -> tuple[float, float, float] | None:
    pos = [r[key] for r in rows if r["cls"] == "dual" and r.get(key) is not None]
    a = [r[key] for r in rows if r["cls"] == "A_only" and r.get(key) is not None]
    b = [r[key] for r in rows if r["cls"] == "B_only" and r.get(key) is not None]
    if not pos or not a or not b:
        return None
    da, db = auroc(pos, a), auroc(pos, b)
    return da, db, min(da, db)


def find_panel_csv(meta: dict) -> Path | None:
    candidates = [
        ROOT / "data" / meta["panel_csv"],
        ROOT / "data" / meta["dir"] / "tables" / "panel_ligands.csv",
        ROOT / "data" / meta["dir"] / "tables" / "ligands.csv",
        ROOT / "data" / meta["dir"] / "panel_ligands.csv",
    ]
    # broader search
    d = ROOT / "data" / meta["dir"] / "tables"
    if d.exists():
        for p in sorted(d.glob("*.csv")):
            candidates.append(p)
    for p in candidates:
        if not p.exists():
            continue
        rows = load_csv(p)
        if not rows:
            continue
        keys = set(rows[0].keys())
        # need ligand id + two pchembl-ish columns OR class column
        if any(k.lower() in ("ligand", "ligand_id", "id") for k in keys):
            return p
    return None


def resolve_id_key(rows: list[dict]) -> str:
    for k in ("ligand_id", "ligand", "id", "Ligand", "name"):
        if k in rows[0]:
            return k
    return list(rows[0].keys())[0]


def resolve_pchembl(rows: list[dict], prefer: str) -> str | None:
    if prefer in rows[0]:
        return prefer
    # fuzzy
    low = {k.lower(): k for k in rows[0]}
    for cand in (prefer.lower(), prefer.lower().replace("pchembl_", "pchembl "), "class", "label"):
        if cand in low:
            return low[cand]
    return None


def assemble_pack(meta: dict) -> tuple[list[dict], dict]:
    d = ROOT / "data" / meta["dir"] / "tables"
    best9 = {r["ligand"]: r for r in load_csv(d / "scores_gnina_best.csv")}
    mode01 = {r["ligand"]: r for r in load_csv(d / "scores_gnina_best_mode01_backup.csv")}
    long9 = load_csv(d / "scores_gnina_long.csv")
    # mode preference: which mode wins
    win_mode = defaultdict(list)
    for r in long9:
        if r.get("status") not in ("success", "exists"):
            continue
        lig, tgt = r["ligand"], r["target"]
        sc = fnum(r.get("cnn_score"))
        if sc is None:
            continue
        win_mode[(lig, tgt)].append((sc, r["mode"]))

    panel_path = find_panel_csv(meta)
    panel_rows = load_csv(panel_path) if panel_path else []
    info = {
        "panel_path": str(panel_path) if panel_path else "",
        "n_best9": len(best9),
        "n_mode01": len(mode01),
        "n_long": len(long9),
    }
    if not panel_rows:
        # fall back: ligands from score tables only, no AUROC
        out = []
        for lig in sorted(set(best9) | set(mode01)):
            b9 = best9.get(lig, {})
            m1 = mode01.get(lig, {})
            out.append(
                {
                    "pair": meta["pair"],
                    "ligand": lig,
                    "cls": "",
                    "gnina_min_mode01": gnina_min_from_best(m1, meta["targets"]),
                    "gnina_min_best9": gnina_min_from_best(b9, meta["targets"]),
                }
            )
        info["labels"] = False
        return out, info

    id_key = resolve_id_key(panel_rows)
    # class column?
    class_key = None
    for k in ("class", "label", "cls", "category", "four_class"):
        if k in panel_rows[0]:
            class_key = k
            break
    a_key = resolve_pchembl(panel_rows, meta["a_key"])
    b_key = resolve_pchembl(panel_rows, meta["b_key"])
    # also try common alternates
    if a_key is None or b_key is None:
        keys = list(panel_rows[0].keys())
        pcols = [k for k in keys if "pchembl" in k.lower() or "pChEMBL" in k]
        if len(pcols) >= 2:
            a_key, b_key = pcols[0], pcols[1]

    out = []
    for pr in panel_rows:
        lig = pr[id_key]
        if class_key:
            cls = pr[class_key]
        else:
            cls = classify(fnum(pr.get(a_key)), fnum(pr.get(b_key)), meta["label_rule"])
        b9 = best9.get(lig, {})
        m1 = mode01.get(lig, {})
        rec = {
            "pair": meta["pair"],
            "ligand": lig,
            "cls": cls or "",
            "gnina_min_mode01": gnina_min_from_best(m1, meta["targets"]),
            "gnina_min_best9": gnina_min_from_best(b9, meta["targets"]),
        }
        # winning modes
        for t in meta["targets"]:
            opts = win_mode.get((lig, t), [])
            if opts:
                opts.sort(key=lambda x: -x[0])
                rec[f"best_mode_{t}"] = opts[0][1]
                rec[f"best_cnn_{t}"] = opts[0][0]
                rec[f"mode01_is_best_{t}"] = int(opts[0][1] == "mode_01")
        out.append(rec)
    info["labels"] = True
    info["id_key"] = id_key
    info["class_key"] = class_key
    info["a_key"] = a_key
    info["b_key"] = b_key
    return out, info


def main():
    out_tab = ROOT / "data/jcim_bench_v0/tables"
    out_an = ROOT / "data/jcim_bench_v0/analysis"
    out_tab.mkdir(parents=True, exist_ok=True)
    out_an.mkdir(parents=True, exist_ok=True)

    all_lig = []
    auroc_rows = []
    notes = []
    for meta in PACKS:
        rows, info = assemble_pack(meta)
        all_lig.extend(rows)
        notes.append(f"- **{meta['pair']}**: panel=`{info.get('panel_path')}` best9={info['n_best9']} mode01={info['n_mode01']} long={info['n_long']}")
        d = directional(rows, "gnina_min_mode01")
        e = directional(rows, "gnina_min_best9")
        # mode01 win rate
        win_flags = []
        for r in rows:
            for t in meta["targets"]:
                k = f"mode01_is_best_{t}"
                if k in r:
                    win_flags.append(r[k])
        win_rate = (sum(win_flags) / len(win_flags)) if win_flags else float("nan")
        auroc_rows.append(
            {
                "pair": meta["pair"],
                "n_labeled": sum(1 for r in rows if r.get("cls") in ("dual", "A_only", "B_only")),
                "D_vs_A_mode01": "" if d is None else f"{d[0]:.4f}",
                "D_vs_B_mode01": "" if d is None else f"{d[1]:.4f}",
                "summary_min_mode01": "" if d is None else f"{d[2]:.4f}",
                "D_vs_A_best9": "" if e is None else f"{e[0]:.4f}",
                "D_vs_B_best9": "" if e is None else f"{e[1]:.4f}",
                "summary_min_best9": "" if e is None else f"{e[2]:.4f}",
                "delta_summary_min_best9_minus_mode01": ""
                if d is None or e is None
                else f"{e[2] - d[2]:+.4f}",
                "frac_mode01_is_best_pose": f"{win_rate:.3f}" if win_flags else "",
            }
        )

    lig_path = out_tab / "gnina_mode01_vs_best9_ligand.csv"
    with lig_path.open("w", newline="") as fh:
        fields = sorted({k for r in all_lig for k in r.keys()})
        # stable order
        prefer = [
            "pair",
            "ligand",
            "cls",
            "gnina_min_mode01",
            "gnina_min_best9",
        ]
        fields = prefer + [f for f in fields if f not in prefer]
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in all_lig:
            w.writerow(r)

    auroc_path = out_tab / "gnina_mode01_vs_best9_auroc.csv"
    with auroc_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(auroc_rows[0].keys()))
        w.writeheader()
        w.writerows(auroc_rows)

    md = out_an / "GNINA_BEST9_STATUS.md"
    lines = [
        "# GNINA rescore — STATUS: DONE (best-of-9 CNN minimize/rescore)",
        "",
        "**Updated:** after fair all-mode GNINA rescore",
        "",
        "## Binary / protocol",
        "`/mnt/d/CADD paper exercise/gnina/bin/gnina` (v1.3.2), CPU `--no_gpu`",
        "",
        "- Input: **all** Vina `mode_01`…`mode_09` PDBQT → SDF (Open Babel)",
        "- `gnina --cnn_scoring rescore --minimize --seed 20260727 --cpu 1`",
        "- Per ligand–target: take **max CNNscore** over up to 9 modes (ties → first max)",
        "- Pocket-matched arm still uses `min(score_A, score_B)` on the per-end best-of-9 CNN",
        "- mode_01-only tables retained as `scores_gnina_*_mode01_backup.csv`",
        "",
        "## Packs",
        *notes,
        "",
        "## mode_01 vs best-of-9 directional AUROC",
        "",
        "| pair | summary_min mode01 | summary_min best9 | Δ | frac mode01 wins |",
        "|------|-------------------:|------------------:|--:|-----------------:|",
    ]
    for r in auroc_rows:
        lines.append(
            f"| {r['pair']} | {r['summary_min_mode01']} | {r['summary_min_best9']} | "
            f"{r['delta_summary_min_best9_minus_mode01']} | {r['frac_mode01_is_best_pose']} |"
        )
    lines += [
        "",
        "## Artifacts",
        f"- `{lig_path.relative_to(ROOT)}`",
        f"- `{auroc_path.relative_to(ROOT)}`",
        "- Per-pack: `tables/scores_gnina_long.csv`, `tables/scores_gnina_best.csv`",
        "",
        "## Claim update",
        "RTMScore and GNINA now share the same pose coverage (best-of-9 over the same Vina modes).",
        "Three-engine contrast is pose-symmetric; still do **not** claim a universal docking decision rule.",
        "",
    ]
    md.write_text("\n".join(lines))
    print("wrote", lig_path)
    print("wrote", auroc_path)
    print("wrote", md)
    for r in auroc_rows:
        print(r)


if __name__ == "__main__":
    main()
