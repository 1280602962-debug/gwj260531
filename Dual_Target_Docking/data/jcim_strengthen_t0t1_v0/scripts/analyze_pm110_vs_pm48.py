#!/usr/bin/env python3
"""B3: Pocket-matched Vina summary_min for PM110 vs PM48."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
OUT_MD = REPO / "data/jcim_strengthen_t0t1_v0/analysis/PM110_VS_PM48.md"
OUT_CSV = REPO / "data/jcim_strengthen_t0t1_v0/tables/pm110_vs_pm48_pocket_matched_v1.csv"
N_BOOT = 2000
SEED = 20260729


def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    p, n = np.asarray(pos, float), np.asarray(neg, float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def boot_ci(D, A, B, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        d = [D[i] for i in rng.integers(0, len(D), len(D))]
        a = [A[i] for i in rng.integers(0, len(A), len(A))]
        b = [B[i] for i in rng.integers(0, len(B), len(B))]
        da = auroc([r["B"] for r in d], [r["B"] for r in a])
        db = auroc([r["A"] for r in d], [r["A"] for r in b])
        vals.append(min(da, db))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def load_scores(panel_csv: Path, score_map: dict[str, dict]):
    recs = []
    for r in csv.DictReader(panel_csv.open()):
        sc = score_map.get(r["panel_id"])
        if not sc:
            continue
        recs.append({"cls": r["class"], **sc})
    D = [x for x in recs if x["cls"] == "dual"]
    A = [x for x in recs if x["cls"] == "A_only"]
    B = [x for x in recs if x["cls"] == "B_only"]
    da = auroc([r["B"] for r in D], [r["B"] for r in A])
    db = auroc([r["A"] for r in D], [r["A"] for r in B])
    mn = min(da, db) if da == da and db == db else float("nan")
    lo, hi = (float("nan"), float("nan"))
    if D and A and B:
        lo, hi = boot_ci(D, A, B)
    return {
        "n_dual": len(D),
        "n_A": len(A),
        "n_B": len(B),
        "auroc_DA_pocketB": da,
        "auroc_DB_pocketA": db,
        "summary_min": mn,
        "ci_lo": lo,
        "ci_hi": hi,
    }


def scores_from_ablation(path: Path):
    out = {}
    for r in csv.DictReader(path.open()):
        out[r["ligand"]] = {
            "A": -float(r["4L23_affinity"]),
            "B": -float(r["4JT6_affinity"]),
        }
    return out


def scores_from_long(path: Path):
    tmp: dict[str, dict] = {}
    for r in csv.DictReader(path.open()):
        lig = r.get("ligand") or r.get("panel_id")
        t = r["target"]
        aff = r.get("vina_best") or r.get("vina_score")
        if aff in ("", None):
            continue
        tmp.setdefault(lig, {})
        key = "A" if t == "4L23" else ("B" if t == "4JT6" else None)
        if key:
            tmp[lig][key] = -float(aff)
    return {k: v for k, v in tmp.items() if "A" in v and "B" in v}


def main():
    pm48_panel = REPO / "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv"
    pm110_panel = REPO / "data/pik3ca_mtor_panel110_rdkit_v0/tables/panel_v0_110.csv"
    abl = REPO / "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv"
    best110 = REPO / "data/pik3ca_mtor_panel110_rdkit_v0/tables/scores_vina_best.csv"
    if not best110.exists():
        best110 = Path(
            "/mnt/d/CADD paper exercise/dual target docking/results/"
            "pik3ca_mtor_panel110_rdkit_v0/tables/scores_vina_best.csv"
        )

    m48 = load_scores(pm48_panel, scores_from_ablation(abl))
    m110 = load_scores(pm110_panel, scores_from_long(best110))

    def round_row(panel, m):
        return {
            "panel": panel,
            "metric": "vina_pocket_matched",
            **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in m.items()},
        }

    rows = [round_row("PM48", m48), round_row("PM110", m110)]
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    delta = m110["summary_min"] - m48["summary_min"]
    lines = [
        "# PM110 vs PM48 — pocket-matched Vina",
        "",
        "> 主指标：`summary_min = min(AUROC D vs A on pocket B, AUROC D vs B on pocket A)`；分数取 −vina。",
        "> PM48 carryover 复用既有 E=16 姿态；新增 PM110_* 本机重对接。",
        "",
        "| Panel | n(D/A/B) | AUROC D↔A (B) | AUROC D↔B (A) | **summary_min** | 95% CI |",
        "|-------|----------|---------------|---------------|-----------------|--------|",
        f"| PM48 | {m48['n_dual']}/{m48['n_A']}/{m48['n_B']} | {m48['auroc_DA_pocketB']:.4f} | {m48['auroc_DB_pocketA']:.4f} | **{m48['summary_min']:.4f}** | [{m48['ci_lo']:.2f}, {m48['ci_hi']:.2f}] |",
        f"| PM110 | {m110['n_dual']}/{m110['n_A']}/{m110['n_B']} | {m110['auroc_DA_pocketB']:.4f} | {m110['auroc_DB_pocketA']:.4f} | **{m110['summary_min']:.4f}** | [{m110['ci_lo']:.2f}, {m110['ci_hi']:.2f}] |",
        "",
        f"**Δ summary_min (PM110 − PM48)**: {delta:+.4f}",
        "",
        "## 结论",
        "",
    ]
    if abs(delta) < 0.05:
        lines.append(
            "扩面后 summary_min 与 PM48 接近（|Δ|<0.05）；CI 应收窄。"
            "声称仍受 CLAIM_CEILING 约束（非通用决策臂）。"
        )
    elif delta > 0:
        lines.append("扩面后 summary_min 上升；仍仅作 PM 对证据加强，不外推为四对通用。")
    else:
        lines.append(
            "扩面后 summary_min 下降 → PM48 信号可能受小样本波动；"
            "文中应以 PM110 CI 为准并降调措辞。"
        )

    OUT_MD.write_text("\n".join(lines) + "\n")
    print(OUT_MD.read_text())


if __name__ == "__main__":
    main()
