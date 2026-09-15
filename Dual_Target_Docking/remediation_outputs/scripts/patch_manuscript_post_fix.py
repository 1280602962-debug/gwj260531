#!/usr/bin/env python3
"""Patch EN/ZH manuscripts from post-fix machine-readable CSVs.

No handwritten result numbers. EGFR/HER2 matched-pocket advantage is not claimed
when the corrected 95% CI includes 0.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
CMP = ROOT / "remediation_outputs/canonical_tables/egfr_her2_corrected_comparison.csv"
EGFR_MET = ROOT / "remediation_outputs/egfr_her2_pre_vs_post_box_metrics.csv"
GNINA = ROOT / "remediation_outputs/phase_gnina_independent/independent_dock_formulation_EGFR_HER2.csv"
COG = ROOT / "remediation_outputs/phase1_cognate_redock/cognate_rank_rmsd_corrected_box.csv"


def f3(x) -> str:
    return f"{float(x):.3f}"


def ci_br(lo, hi) -> str:
    return f"[{float(lo):.3f}, {float(hi):.3f}]"


def load_cmp():
    rows = {r["pair"]: r for r in csv.DictReader(CMP.open())}
    return rows


def load_egfr_ci():
    out = {}
    with EGFR_MET.open() as fh:
        for r in csv.DictReader(fh):
            ci = r.get("CI_new") or ""
            if "," in ci:
                lo, hi = ci.split(",", 1)
                out[r["metric"]] = (float(r["new"]) if r["new"][0].isdigit() or r["new"].startswith("-") else r["new"], float(lo), float(hi))
            else:
                out[r["metric"]] = (r["new"], None, None)
    return out


def fmt_auroc_ci(val, lo, hi) -> str:
    return f"{f3(val)} {ci_br(lo, hi)}"


def patch_text(text: str, egfr, ache, egfr_ci, gnina, cog) -> str:
    da, da_lo, da_hi = egfr_ci["AUROC_D_vs_A_pocketB"]
    db, db_lo, db_hi = egfr_ci["AUROC_D_vs_B_pocketA"]
    sm, sm_lo, sm_hi = egfr_ci["summary_min"]
    dlt, dlt_lo, dlt_hi = egfr_ci["fixed_score_delta_pocketA"]
    nei_pA = float(egfr["fixed_score_D_vs_neither_pocketA"])
    nei_mean, nei_lo, nei_hi = egfr_ci["AUROC_mean_D_vs_neither"]
    mm, mm_lo, mm_hi = egfr_ci["delta_summary_min_matched_minus_mismatched"]
    ache_da = float(ache["AUROC_D_vs_A_pocketB"])
    ache_db = float(ache["AUROC_D_vs_B_pocketA"])
    ache_mm = float(ache["matched_minus_mismatched"])
    ache_mm_lo, ache_mm_hi = ache["matched_minus_mismatched_ci"].split(",")

    # Abstract / results EGFR fixed-score sentence
    text = text.replace(
        "the EGFR-pocket AUROC for dual-versus-B-only was 0.430 and rose to 0.808 against neither (difference 0.378 [0.205, 0.547])",
        f"the EGFR-pocket AUROC for dual-versus-B-only was {f3(db)} and rose to {f3(nei_pA)} against neither (difference {f3(dlt)} {ci_br(dlt_lo, dlt_hi)})",
    )
    text = text.replace(
        "EGFR 口袋评分对 dual–B-only 的 AUROC 为 0.430，对 neither 升至 0.808（差值 0.378 [0.205, 0.547]）",
        f"EGFR 口袋评分对 dual–B-only 的 AUROC 为 {f3(db)}，对 neither 升至 {f3(nei_pA)}（差值 {f3(dlt)} {ci_br(dlt_lo, dlt_hi)}）",
    )
    text = text.replace(
        "the EGFR/HER2 dual-versus-B-only AUROC was 0.430 and rose to 0.808 against neither (difference 0.378 [0.205, 0.547])",
        f"the EGFR/HER2 dual-versus-B-only AUROC was {f3(db)} and rose to {f3(nei_pA)} against neither (difference {f3(dlt)} {ci_br(dlt_lo, dlt_hi)})",
    )
    text = text.replace(
        "dual–B-only 的 AUROC 为 0.430，将对照换成 neither 后升至 0.808，差值 0.378 [0.205, 0.547]",
        f"dual–B-only 的 AUROC 为 {f3(db)}，将对照换成 neither 后升至 {f3(nei_pA)}，差值 {f3(dlt)} {ci_br(dlt_lo, dlt_hi)}",
    )

    old_t2 = "| EGFR/HER2 | 28 / 38 / 32 | 0.666 [0.524, 0.793] | 0.430 [0.282, 0.579] | 0.430 [0.282, 0.578] |"
    new_t2 = f"| EGFR/HER2 | 28 / 38 / 32 | {fmt_auroc_ci(da, da_lo, da_hi)} | {fmt_auroc_ci(db, db_lo, db_hi)} | {fmt_auroc_ci(sm, sm_lo, sm_hi)} |"
    text = text.replace(old_t2, new_t2)

    old_ache_t2 = "| AChE/BChE | 27 / 25 / 28 | 0.650 [0.483, 0.801] | 0.606 [0.442, 0.751] | 0.606 [0.437, 0.730] |"
    new_ache_t2 = f"| AChE/BChE | 27 / 26 / 28 | {f3(ache_da)} | {f3(ache_db)} | {f3(ache['summary_min'])} |"
    text = text.replace(old_ache_t2, new_ache_t2)

    text = text.replace(
        "EGFR/HER2 reached 0.756 [0.562, 0.920]",
        f"EGFR/HER2 reached {fmt_auroc_ci(nei_mean, nei_lo, nei_hi)}",
    )
    text = text.replace("EGFR/HER2 为 0.756 [0.562, 0.920]", f"EGFR/HER2 为 {fmt_auroc_ci(nei_mean, nei_lo, nei_hi)}")
    text = text.replace(
        "whereas their directional $\\mathrm{summary}_{\\min}$ values were 0.430 and 0.365",
        f"whereas their directional $\\mathrm{{summary}}_{{\\min}}$ values were {f3(sm)} and 0.365",
    )
    text = text.replace(
        "而其方向性 $\\mathrm{summary}_{\\min}$ 分别为 0.430 和 0.365",
        f"而其方向性 $\\mathrm{{summary}}_{{\\min}}$ 分别为 {f3(sm)} 和 0.365",
    )

    old_t3 = "| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 |"
    new_t3 = f"| EGFR/HER2 | {fmt_auroc_ci(sm, sm_lo, sm_hi)} | {fmt_auroc_ci(nei_mean, nei_lo, nei_hi)} | 12 |"
    text = text.replace(old_t3, new_t3)

    # Matched-pocket: do not claim EGFR advantage
    en_old = (
        "In the main panels, only EGFR/HER2 and AChE/BChE had matched-minus-mismatched "
        "$\\mathrm{summary}_{\\min}$ 95% intervals that excluded zero, with differences of "
        "0.170 [0.060, 0.280] and 0.161 [0.037, 0.269]. The other six main-panel intervals included zero."
    )
    en_new = (
        "In the main panels, only AChE/BChE had a matched-minus-mismatched "
        f"$\\mathrm{{summary}}_{{\\min}}$ 95% interval that excluded zero "
        f"({f3(ache_mm)} {ci_br(ache_mm_lo, ache_mm_hi)}). "
        f"EGFR/HER2 was {f3(mm)} {ci_br(mm_lo, mm_hi)}, which includes 0, so this pair is not interpreted as having a matched-pocket advantage. "
        "The other six main-panel intervals included zero."
    )
    text = text.replace(en_old, en_new)

    zh_old = (
        "主评价集中仅 EGFR/HER2 和 AChE/BChE 的 matched−mismatched $\\mathrm{summary}_{\\min}$ 95% 区间排除 0，"
        "差值分别为 0.170 [0.060, 0.280] 和 0.161 [0.037, 0.269]。其余六对主集区间包含 0。"
    )
    zh_new = (
        f"主评价集中仅 AChE/BChE 的 matched−mismatched $\\mathrm{{summary}}_{{\\min}}$ 95% 区间排除 0"
        f"（{f3(ache_mm)} {ci_br(ache_mm_lo, ache_mm_hi)}）。"
        f"EGFR/HER2 为 {f3(mm)} {ci_br(mm_lo, mm_hi)}，区间包含 0，因此不将该对解释为具有 matched-pocket 优势。"
        "其余六对主集区间包含 0。"
    )
    text = text.replace(zh_old, zh_new)

    text = text.replace(
        "On the main panels, only EGFR/HER2 and AChE/BChE had positive matched-minus-mismatched intervals that excluded 0.",
        "On the main panels, only AChE/BChE had a positive matched-minus-mismatched interval that excluded 0; the EGFR/HER2 interval included 0.",
    )
    text = text.replace(
        "八个主评价集中仅两个 matched−mismatched 的 95% 区间排除 0",
        "八个主评价集中仅 AChE/BChE 的 matched−mismatched 95% 区间排除 0",
    )

    text = text.replace("1/95 AChE/BChE classes", "1/96 AChE/BChE classes")
    text = text.replace("AChE/BChE 有 1/95 个翻转", "AChE/BChE 有 1/96 个翻转")
    text = text.replace("($\\mathrm{summary}_{\\min}$ 0.430 to 0.424)", f"($\\mathrm{{summary}}_{{\\min}}$ {f3(sm)} to 0.424)")
    text = text.replace("（$\\mathrm{summary}_{\\min}$ 0.430 变为 0.424）", f"（$\\mathrm{{summary}}_{{\\min}}$ {f3(sm)} 变为 0.424）")

    if cog:
        text = text.replace(
            "EGFR 3POZ top-1 was 9.505 Å, with lowest saved-pose RMSD 0.760 Å.",
            f"EGFR 3POZ top-1 was {cog['3POZ_top1']:.3f} Å, with lowest saved-pose RMSD {cog['3POZ_best']:.3f} Å.",
        )
    if gnina:
        text = re.sub(
            r"For EGFR/HER2, the dual-versus-neither AUROC was 0\.783 \[0\.610, 0\.922\] \(\\?n_\{\\mathrm\{neither\}\}=11\), whereas dual-versus-B-only was 0\.220 \[0\.109, 0\.343\]\.",
            f"For EGFR/HER2, the dual-versus-neither AUROC was {gnina['nei']} ({gnina['n_nei']}), whereas dual-versus-B-only was {gnina['db']}.",
            text,
        )
    return text


def read_gnina():
    if not GNINA.exists():
        return None
    rows = {r["contrast"]: r for r in csv.DictReader(GNINA.open())}
    db = rows["D_vs_B_pocketA"]
    nei = rows["D_vs_neither_mean"]
    return {
        "db": fmt_auroc_ci(db["auroc"], db["ci_lo"], db["ci_hi"]),
        "nei": fmt_auroc_ci(nei["auroc"], nei["ci_lo"], nei["ci_hi"]),
        "n_nei": f"$n_{{\\mathrm{{neither}}}}={int(float(nei['n_neg']))}$",
    }


def read_cog():
    if not COG.exists():
        return None
    by = {}
    with COG.open() as fh:
        for r in csv.DictReader(fh):
            if int(r["pose_rank"]) == 1:
                by[f"{r['pdb']}_top1"] = float(r["rmsd_A"])
                by[f"{r['pdb']}_best"] = float(r["best_all_deposited_A"])
    return by


def main() -> int:
    rows = load_cmp()
    egfr, ache = rows["EGFR/HER2"], rows["AChE/BChE"]
    egfr_ci = load_egfr_ci()
    gnina = read_gnina()
    cog = read_cog()
    for rel in ("docs/MANUSCRIPT_JCIM_EN.md", "docs/MANUSCRIPT_JCIM_ZH.md",
                "docs/RESULTS_SECTION_JCIM_EN_V1.md", "docs/DISCUSSION_SECTION_JCIM_EN_V1.md",
                "docs/TITLE_AND_ABSTRACT_JCIM_EN_V1.md", "docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md"):
        path = ROOT / rel
        if not path.exists():
            continue
        old = path.read_text(encoding="utf-8")
        new = patch_text(old, egfr, ache, egfr_ci, gnina, cog)
        if new != old:
            path.write_text(new, encoding="utf-8")
            print("patched", rel)
        else:
            print("unchanged", rel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
