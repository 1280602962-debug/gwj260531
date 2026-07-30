#!/usr/bin/env python3
"""Compare PM48 E=8 vs E=16 pocket-matched summary_min."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
PM48 = REPO / "data/pik3ca_mtor_panel48_rdkit_v0"
PANEL = PM48 / "tables/panel_v0_48.csv"
E16 = PM48 / "tables/ablation_ligand_scores.csv"
E8 = PM48 / "tables/scores_vina_E8_best.csv"
OUT = REPO / "data/jcim_strengthen_t0t1_v0/analysis/EXHAUSTIVENESS_E8_VS_E16.md"


def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    p, n = np.asarray(pos, float), np.asarray(neg, float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def load_panel():
    return {r["panel_id"]: r for r in csv.DictReader(PANEL.open())}


def main():
    panel = load_panel()
    e16 = {r["ligand"]: r for r in csv.DictReader(E16.open())}
    e8 = {r["ligand"]: r for r in csv.DictReader(E8.open())} if E8.exists() else {}

    def metrics(score_fn):
        recs = []
        for lig, meta in panel.items():
            sc = score_fn(lig)
            if sc is None:
                continue
            recs.append({"cls": meta["class"], **sc})
        D = [r for r in recs if r["cls"] == "dual"]
        A = [r for r in recs if r["cls"] == "A_only"]
        B = [r for r in recs if r["cls"] == "B_only"]
        da = auroc([r["B"] for r in D], [r["B"] for r in A])
        db = auroc([r["A"] for r in D], [r["A"] for r in B])
        return da, db, min(da, db), len(D), len(A), len(B)

    def e16_fn(lig):
        r = e16.get(lig)
        if not r:
            return None
        return {"A": -float(r["4L23_affinity"]), "B": -float(r["4JT6_affinity"])}

    def e8_fn(lig):
        r = e8.get(lig)
        if not r:
            return None
        a, b = r.get("4L23_affinity_E8"), r.get("4JT6_affinity_E8")
        if not a or not b:
            return None
        try:
            return {"A": -float(a), "B": -float(b)}
        except ValueError:
            return None

    m16 = metrics(e16_fn)
    m8 = metrics(e8_fn)

    lines = [
        "# EXHAUSTIVENESS E8 vs E16 — PIK3CA/mTOR PM48",
        "",
        "> 协议：同一 RDKit/meeko 配体；仅 exhaustiveness 8 vs 16；口袋匹配 summary_min。",
        "",
        "| E | AUROC D vs A (pocket B) | AUROC D vs B (pocket A) | **summary_min** | n(D/A/B) |",
        "|---|-------------------------|-------------------------|-----------------|----------|",
        f"| 16 | {m16[0]:.4f} | {m16[1]:.4f} | **{m16[2]:.4f}** | {m16[3]}/{m16[4]}/{m16[5]} |",
        f"| 8 | {m8[0]:.4f} | {m8[1]:.4f} | **{m8[2]:.4f}** | {m8[3]}/{m8[4]}/{m8[5]} |",
        "",
        f"**Δ summary_min (E16 − E8)**: {m16[2]-m8[2]:+.4f}",
        "",
        "## 结论",
        "",
    ]
    delta = m16[2] - m8[2]
    if abs(delta) < 0.05:
        lines.append("E=8 与 E=16 口袋匹配 summary_min 差异 <0.05 → **exhaustiveness 不足以单独解释 PM 优势**；需在文中并列报告。")
    elif delta > 0.05:
        lines.append("E=16 明显高于 E=8 → 部分 PM 信号可能来自更充分采样；主文应报告 E8 对照并谨慎措辞。")
    else:
        lines.append("E=8 不低于 E=16 → PM 优势**不能**归因于 exhaustiveness 差异。")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")
    print(OUT.read_text())


if __name__ == "__main__":
    main()
