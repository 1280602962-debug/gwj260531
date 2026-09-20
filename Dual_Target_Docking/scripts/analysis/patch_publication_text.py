#!/usr/bin/env python3
"""Rewrite publication-facing numbers from results/canonical (scheme B)."""
from __future__ import annotations

import csv
import math
import re
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CANON = ROOT / "results" / "canonical"
DOCS = ROOT / "docs"
PAIRS = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def r3(x) -> str:
    d = Decimal(str(float(x))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)
    return f"{d:.3f}"


def r3m(x) -> str:
    s = r3(x)
    return "−" + s[1:] if s.startswith("-") else s


def ci(lo, hi) -> str:
    if lo in ("", None) or hi in ("", None):
        return "not recomputed"
    return f"[{r3m(lo)}, {r3m(hi)}]"


def _ensure_once(text: str, needle: str, addition: str) -> str:
    """Keep exactly one copy of addition immediately after needle."""
    if needle not in text:
        return text
    i = text.index(needle) + len(needle)
    rest = text[i:]
    while rest.startswith(addition):
        rest = rest[len(addition):]
    return text[:i] + addition + rest


def signed(x) -> str:
    v = float(x)
    if abs(v) < 5e-4:
        return r3(0)
    return ("+" if v > 0 else "−") + r3(abs(v))


def _maxmed_lookup(d):
    by = {}
    for r in d.get("maxmed") or []:
        by.setdefault(r["pair"], {})[r["aggregation"]] = r
    return by


def _cluster_lookup(d, pair):
    return {r["estimator"]: r for r in d.get("cluster") or [] if r["pair"] == pair}


def apply_boundary_patches(text: str, d: dict, zh: bool) -> str:
    """Scientific-boundary wording that must not be restored by number-only patches."""
    s_egfr = d["smin"]["EGFR/HER2"]
    mx = _maxmed_lookup(d)
    egfr_md = mx.get("EGFR/HER2", {}).get("median", {})
    ache_md = mx.get("AChE/BChE", {}).get("median", {})
    ache_mx = mx.get("AChE/BChE", {}).get("max", {})
    ppar_md = mx.get("PPARA/PPARD", {}).get("median", {})
    cl_egfr = _cluster_lookup(d, "EGFR/HER2")
    cl_jak = _cluster_lookup(d, "JAK1/TYK2")
    egfr_flips = egfr_md.get("class_flips_vs_max", "")
    egfr_n = egfr_md.get("n_ligands", "")
    ache_flips = ache_md.get("class_flips_vs_max", "")
    ache_n = ache_md.get("n_ligands", "")
    ppar_flips = ppar_md.get("class_flips_vs_max", "")
    ppar_n = ppar_md.get("n_ligands", "")
    jak_doc = cl_jak.get("document_cluster") or {}
    if jak_doc.get("delta_ci_lo") not in ("", None):
        jak_doc_en = (
            f"document-cluster was {ci(jak_doc['delta_ci_lo'], jak_doc['delta_ci_hi'])} "
            "(includes 0; current recompute)"
        )
        jak_doc_zh = (
            f"文献簇区间为 {ci(jak_doc['delta_ci_lo'], jak_doc['delta_ci_hi'])}（包含 0；当前重算）"
        )
    else:
        jak_doc_en = (
            "document-cluster could not be recomputed because the ligand–document map and ChEMBL 37 sqlite "
            "are unavailable; a previously deposited interval included 0 and is not treated as a current calculation"
        )
        jak_doc_zh = (
            "文献簇无法在本冻结中重算（配体–文献分组映射与 ChEMBL 37 sqlite 均不可用）；"
            "此前存档区间包含 0，不作为当前计算结果"
        )
    if cl_egfr.get("document_cluster") and cl_egfr["document_cluster"].get("delta_ci_lo"):
        egfr_doc = ci(cl_egfr["document_cluster"]["delta_ci_lo"], cl_egfr["document_cluster"]["delta_ci_hi"])
        egfr_scaf = ci(cl_egfr["scaffold_cluster"]["delta_ci_lo"], cl_egfr["scaffold_cluster"]["delta_ci_hi"])
        jak_scaf = ci(cl_jak["scaffold_cluster"]["delta_ci_lo"], cl_jak["scaffold_cluster"]["delta_ci_hi"])
        if zh:
            text = text.replace(
                "EGFR/HER2 簇区间为 [0.235, 0.665]（骨架）和 [0.125, 0.644]（文献），均排除 0。"
                "JAK1/TYK2 骨架簇区间为 [0.220, 0.631]，排除 0，文献簇区间为 [−0.034, 0.682]，包含 0",
                f"EGFR/HER2 簇区间为 {egfr_scaf}（骨架）和 {egfr_doc}（文献），均排除 0。"
                f"JAK1/TYK2 骨架簇区间为 {jak_scaf}，排除 0；{jak_doc_zh}",
            )
        else:
            text = text.replace(
                "the EGFR/HER2 cluster intervals were [0.235, 0.665] (scaffold) and [0.125, 0.644] (document); both exclude 0. "
                "JAK1/TYK2 scaffold-cluster was [0.220, 0.631] (excludes 0) and document-cluster was [−0.034, 0.682] (includes 0)",
                f"the EGFR/HER2 cluster intervals were {egfr_scaf} (scaffold) and {egfr_doc} (document); both exclude 0. "
                f"JAK1/TYK2 scaffold-cluster was {jak_scaf} (excludes 0) and {jak_doc_en}",
            )
    if zh:
        text = text.replace(
            "在同一套 ChEMBL 37 记录上将最大 pChEMBL 改为中位数后，EGFR/HER2 有 6/110 个类别翻转（主分析 \(\mathrm{summary}_{\min}\) 0.324），AChE/BChE 有 1/96 个翻转（CHEMBL659；0.606 变为 0.629），PPARA/PPARD 有 1/110 个翻转（CHEMBL121；\(\mathrm{summary}_{\min}\) 仍为 0.446）。其余五对保持类别组成和 \(\mathrm{summary}_{\min}\) 点估计（Table S3）。",
            f"最大与中位数聚合使用同一套合格原始记录、同一分子交集和当前对接分数，而不是把全部已评分配体与 dump 缺失记录混为一谈。"
            f"EGFR/HER2 有 {egfr_flips}/{egfr_n} 个类别翻转（主分析 \(\mathrm{{summary}}_{{\min}}\) {r3(s_egfr['summary_min'])}）；"
            f"AChE/BChE 有 {ache_flips}/{ache_n} 个翻转（CHEMBL659；max {r3(ache_mx.get('summary_min', 0.6058))} 变为 {r3(ache_md.get('summary_min', 0.6291))}）；"
            f"PPARA/PPARD 有 {ppar_flips}/{ppar_n} 个翻转（CHEMBL121；\(\mathrm{{summary}}_{{\min}}\) 仍为 {r3(mx.get('PPARA/PPARD', {}).get('max', {}).get('summary_min', 0.4463))}）。"
            "其余靶对保持类别组成和 \(\mathrm{summary}_{\min}\) 点估计（Table S3）。",
        )
        text = text.replace(
            "最大–中位数聚合在全部八对上对同一套 ChEMBL 37 记录重复进行，从而只改变实验标签。",
            "最大–中位数聚合按靶对使用同一套合格记录和当前分数；EGFR/HER2、PIK3CA/mTOR 与 AChE/BChE 使用裁决后的高置信记录，其余靶对使用 dump-gated 记录。缺少 dump 记录不等于从主分析删除该分子。",
        )
        text = text.replace(
            "制备脚本将受体链、altloc 和残基模板选择记录在 Table S2。",
            "Table S2 记录 PDB、共晶配体、分辨率和对接盒坐标。链、altLoc 和残基模板见已提交的受体 PDBQT 与制备记录，不在 Table S2。从提交 PDBQT 复现不同于从原始 PDB 重新制备。",
        )
        text = text.replace("| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 38 / 32 | 8 |",
                            "| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 37 / 32 | 8 |")
        text = text.replace("| EGFR/HER2 | \(\theta=6.0\) | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 38 / 32 | 8 |",
                            "| EGFR/HER2 | \(\theta=6.0\) | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 37 / 32 | 8 |")
        text = text.replace("| EGFR/HER2 | \\(\\theta=6.0\\) | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 38 / 32 | 8 |",
                            "| EGFR/HER2 | \\(\\theta=6.0\\) | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 37 / 32 | 8 |")
        if "实际对接的 mTOR 结构 4JT6 分辨率为 3.60 Å" not in text:
            text = text.replace(
                "两端各至少 5 个人源全配体结构（\\(\\leq 3.5\\) Å，至少一个非聚合物配体）后剩 19 对。",
                "两端各至少 5 个人源全配体结构（\\(\\leq 3.5\\) Å，至少一个非聚合物配体）后剩 19 对。"
                "该 3.5 Å 门槛是靶对供给筛选。实际对接的 mTOR 结构 4JT6 分辨率为 3.60 Å，不改写为 3.5 Å 结构。",
            )
    else:
        text = text.replace(
            "Relabeling the same ChEMBL 37 records by median rather than maximum pChEMBL flipped 6/110 EGFR/HER2 classes (primary \(\mathrm{summary}_{\min}\) 0.324), 1/96 AChE/BChE classes (CHEMBL659; 0.606 to 0.629), and 1/110 PPARA/PPARD classes (CHEMBL121; \(\mathrm{summary}_{\min}\) remained 0.446). The other five pairs kept class composition and \(\mathrm{summary}_{\min}\) point estimates (Table S3).",
            f"Maximum and median aggregation used the same qualified source records, the same ligand intersection, and the current docking scores; dump-missing ligands were not mixed into that denominator. "
            f"Class flips: EGFR/HER2 {egfr_flips}/{egfr_n} (primary \(\mathrm{{summary}}_{{\min}}\) {r3(s_egfr['summary_min'])}); "
            f"AChE/BChE {ache_flips}/{ache_n} (CHEMBL659; max {r3(ache_mx.get('summary_min', 0.6058))} to {r3(ache_md.get('summary_min', 0.6291))}); "
            f"PPARA/PPARD {ppar_flips}/{ppar_n} (CHEMBL121; \(\mathrm{{summary}}_{{\min}}\) remained {r3(mx.get('PPARA/PPARD', {}).get('max', {}).get('summary_min', 0.4463))}). "
            "The other pairs kept class composition and \(\mathrm{summary}_{\min}\) point estimates (Table S3).",
        )
        text = text.replace(
            "Maximum-versus-median aggregation was repeated on the same ChEMBL 37 records for all eight pairs, so that only the experimental labels changed. That check does not replace Table 2.",
            "Maximum-versus-median aggregation used one qualified record set per pair and the current scores, so that only the experimental labels changed on that intersection. EGFR/HER2, PIK3CA/mTOR, and AChE/BChE used adjudicated high-confidence rows; the remaining pairs used dump-gated rows. A missing dump row is not automatic removal from the primary analysis. That check does not replace Table 2.",
        )
        text = text.replace(
            "The preparation scripts record receptor-specific chain selection, alternate-location handling, and residue-template choices in Table S2. The deposited preparation records do not specify a common pH-dependent protonation or missing-loop reconstruction procedure.",
            "Table S2 records PDB identifiers, cognate ligands, resolutions, and docking-box coordinates. Chain, altLoc, and residue-template choices are in the deposited receptor PDBQT files and receptor-prep records, not in Table S2. Reproducing scores from those PDBQT files is not the same as re-preparing receptors from the original PDBs. A common pH-dependent protonation or missing-loop reconstruction procedure was not uniformly recorded and is not reconstructed here.",
        )
        text = text.replace("| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 38 / 32 | 8 |",
                            "| EGFR/HER2 | θ = 6.0 | 28 / 38 / 32 / 12 | 3POZ / 3RCD | 28 / 37 / 32 | 8 |")
        if "The selected mTOR structure 4JT6 is 3.60 Å" not in text:
            text = text.replace(
                "Requiring at least five human holo structures per end (\(\leq 3.5\) Å, at least one non-polymer ligand) left 19 pairs.",
                "Requiring at least five human holo structures per end (\(\leq 3.5\) Å, at least one non-polymer ligand) left 19 pairs. That 3.5 Å cutoff was a pair-supply screen. The selected mTOR structure 4JT6 is 3.60 Å and is the receptor actually docked; it is not rewritten as a 3.5 Å structure.",
            )
        text = text.replace(
            "lines join the two points for one direction and are not confidence intervals.",
            "horizontal segments connect each \(\Delta\)AUROC to 0 and are not confidence intervals.",
        )
    text = text.replace(
        "all eight pairs used a fully identical complete record set",
        "the eight pairs used one analysis protocol with pair-specific candidate pools and sampling rules",
    )
    return text


def load():
    smin = {r["pair"]: r for r in read_csv(CANON / "primary_summary_min.csv")}
    direc = {(r["pair"], r["estimand"]): r for r in read_csv(CANON / "primary_directional_auroc.csv")}
    two = {r["pair"]: r for r in read_csv(CANON / "two_pocket_mean_ranking.csv")}
    rank = {r["pair"]: r for r in read_csv(CANON / "top10_operating_points.csv")}
    andf = {r["pair"]: r for r in read_csv(ROOT / "data/jcim_novelty_v0/tables/operating_point_examples_review_v1.csv")}
    fixed = {(r["pair"], r["contrast"]): r for r in read_csv(CANON / "fixed_score_negative_class_delta.csv")}
    mm = {r["pair"]: r for r in read_csv(CANON / "matched_mismatched_pocket.csv")}
    hold = {r["pair"]: r for r in read_csv(CANON / "holdout_metrics.csv")}
    wp = read_csv(ROOT / "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv")
    hold_wp = {r["pair"]: r for r in wp if r["set"] == "unused_pool_holdout"}
    desc = {r["pair"]: r for r in read_csv(CANON / "descriptor_baselines.csv")}
    inc = {(r["pair"], r["contrast"]): r for r in read_csv(CANON / "ecfp4_incremental_information.csv")}
    labels = read_csv(CANON / "label_aggregation_sensitivity.csv")
    maxmed = read_csv(CANON / "max_vs_median_sensitivity.csv")
    cluster = read_csv(CANON / "cluster_bootstrap_sensitivity.csv")
    gnina = {r["pair"]: r for r in read_csv(CANON / "computational_robustness.csv") if r["engine"] == "gnina_dock_mode1"}
    rec = {r["replacement"]: r for r in read_csv(CANON / "receptor_substitution.csv")}
    scaler = read_csv(CANON / "ecfp4_scaler_sensitivity.csv")
    det_path = CANON / "detectable_effect_simulation.csv"
    det = read_csv(det_path) if det_path.is_file() else []
    return locals()


def table2(d, zh=False):
    h = (
        "| 靶对 | n_scored (dual / A-only / B-only) | Dual vs A-only（口袋 B）[95% CI] | Dual vs B-only（口袋 A）[95% CI] | summary_min [95% CI] |"
        if zh
        else "| Pair | n_scored (dual / A-only / B-only) | dual vs A-only (pocket B) [95% CI] | dual vs B-only (pocket A) [95% CI] | summary_min [95% CI] |"
    )
    rows = [h, "|------|---------------------------:|-------------------------:|-------------------------:|----------------------|"]
    for p in PAIRS:
        s = d["smin"][p]
        da = d["direc"][(p, "AUROC_D_vs_A_pocketB")]
        db = d["direc"][(p, "AUROC_D_vs_B_pocketA")]
        rows.append(
            f"| {p} | {s['n_dual']} / {s['n_A_only']} / {s['n_B_only']} | "
            f"{r3(da['point'])} {ci(da['ci_lo'], da['ci_hi'])} | "
            f"{r3(db['point'])} {ci(db['ci_lo'], db['ci_hi'])} | "
            f"{r3(s['summary_min'])} {ci(s['ci_lo'], s['ci_hi'])} |"
        )
    return "\n".join(rows)


def table3(d, zh=False):
    h = (
        "| 靶对 | 双口袋平均 D-vs-neither AUROC [95% CI] | n_neither | Top 10% D/A/B/N | \(\mathrm{EF}_{\mathrm{dual},10\%}\) |"
        if zh
        else "| Pair | two-pocket mean D-vs-neither AUROC [95% CI] | n_neither | Top 10% D/A/B/N | \(\mathrm{EF}_{\mathrm{dual},10\%}\) |"
    )
    rows = [h, "|------|--------------------------------------------:|----------:|----------------:|-------------------------------------:|"]
    for p in PAIRS:
        t = d["two"][p]
        rk = d["rank"][p]
        rows.append(
            f"| {p} | {r3(t['two_pocket_mean_D_vs_neither'])} {ci(t['ci_lo'], t['ci_hi'])} | "
            f"{rk['n_neither']} | {rk['top_dual']} / {rk['top_A_only']} / {rk['top_B_only']} / {rk['top_neither']} | "
            f"{r3(rk['ef_dual_10pct'])} |"
        )
    return "\n".join(rows)


METHODS_EN = """Table 2 reports pointwise 95% confidence intervals for both directional AUROCs and their descriptive minimum. Primary intervals used a class-stratified nonparametric percentile bootstrap, with B = 2000 replicates and seed 20260729. For each directional AUROC, ligands were resampled with replacement within the two experimental-state classes, preserving the original class sizes. For \(\\mathrm{summary}_{\\min}\), dual ligands were resampled once per replicate and that same dual draw was applied to both pocket scores; A-only and B-only were resampled independently within class; the two directional AUROCs were recomputed and their minimum was taken inside the replicate. Fixed-score \(\\Delta\)AUROC used the same shared dual resample for dual-versus-selective and dual-versus-neither, with the two negative classes resampled independently. Matched-versus-mismatched pocket comparisons used paired ligand resampling so that a ligand's matched and mismatched scores were drawn together. Point estimates used the full analysis sample. A ligand-level non-stratified bootstrap was retained as a sensitivity analysis and did not replace the class-stratified intervals. Cluster bootstrap by Bemis–Murcko scaffold and literature-connected document groups was a source-dependence sensitivity for the two largest fixed-score differences (Table S4) and did not replace the ligand-level primary intervals. A binormal detectable-effect simulation reused the same class-stratified shared-dual bootstrap, the current eight-pair class sizes, B = 2000, and seed 20260729; it estimates the probability that a CI excludes 0.5 under a specified true AUROC and is not observed power. Intervals were not adjusted for multiple comparisons."""

METHODS_ZH = """Table 2 同时报告两条方向性 AUROC 及其描述性最小值的逐项 95% 置信区间。主分析采用类别分层非参数百分位 bootstrap，B = 2000，种子 20260729。方向性 AUROC 在对应的两个实验状态类别内分别有放回抽样，并保持原类别样本量。对于 \(\\mathrm{summary}_{\\min}\)，每次重采样只抽取一次 dual 配体，并将同一次 dual 抽取同时用于两个口袋评分；A-only 与 B-only 在各自类别内独立重采样；在该次重采样内重新计算两个方向 AUROC 后再取较小值。固定评分 \(\\Delta\)AUROC 中，dual–selective 与 dual–neither 共用同一次 dual 抽取，两类阴性对照各自独立重采样。matched 与 mismatched 口袋比较采用配体水平配对重采样，使同一配体的两条口袋评分同时被抽中或同时不被抽中。点估计由完整分析样本直接计算。配体水平非分层 bootstrap 仅作为敏感性分析，不替换类别分层主区间。按 Bemis–Murcko 骨架簇和文献连通簇进行的簇 bootstrap 用于两个最大固定评分差的来源依赖性敏感性（Table S4），不替换配体水平主区间。二项正态 detectable-effect 仿真复用同一类别分层、共享 dual 的 bootstrap、当前八对类别样本量和 B = 2000、种子 20260729；它估计在指定真 AUROC 下 CI 排除 0.5 的概率，不是观测功效。区间未做多重比较校正。"""


def replace_between(text, start, end, new, count=1):
    i = text.find(start)
    if i < 0:
        raise SystemExit(f"missing start {start[:60]!r}")
    j = text.find(end, i)
    if j < 0:
        raise SystemExit(f"missing end {end[:60]!r}")
    return text[:i] + new + text[j:]


def patch_table_block(text, header_line, new_table):
    i = text.find(header_line)
    if i < 0:
        return text
    # from this header through the last consecutive table row
    k = i
    lines = text[i:].splitlines()
    n = 0
    for line in lines:
        if line.startswith("|") or line.startswith("|-"):
            n += 1
        else:
            break
    old = "\n".join(lines[:n])
    return text[:i] + new_table + text[i + len(old) :]


def smin_cross(d):
    above, below, cross = [], [], []
    for p in PAIRS:
        s = d["smin"][p]
        lo, hi = float(s["ci_lo"]), float(s["ci_hi"])
        if lo > 0.5:
            above.append(p)
        elif hi < 0.5:
            below.append(p)
        else:
            cross.append(p)
    return above, below, cross


def max_inc(d):
    best = max(d["inc"].values(), key=lambda r: abs(float(r["delta_ECFP4_plus_docking_minus_ECFP4"])))
    return best, abs(float(best["delta_ECFP4_plus_docking_minus_ECFP4"]))


def desc_excl(d):
    excl, incl = [], []
    for p in PAIRS:
        r = d["desc"][p]
        lo, hi = float(r["delta_ci_lo"]), float(r["delta_ci_hi"])
        if lo > 0 or hi < 0:
            excl.append(p)
        else:
            incl.append(p)
    return excl, incl


def patch_en(text, d):
    egfr_fix = d["fixed"][("EGFR/HER2", "D_vs_B_or_neither_pocketA")]
    jak_fix = d["fixed"][("JAK1/TYK2", "D_vs_B_or_neither_pocketA")]
    s_egfr = d["smin"]["EGFR/HER2"]
    s_jak = d["smin"]["JAK1/TYK2"]
    s_pparg = d["smin"]["PPARG/PPARA"]
    s_f2 = d["smin"]["F2/F10"]
    two_egfr = d["two"]["EGFR/HER2"]
    two_jak = d["two"]["JAK1/TYK2"]
    mm_egfr = d["mm"]["EGFR/HER2"]
    mm_ache = d["mm"]["AChE/BChE"]
    ache_t = d["desc"]["AChE/BChE"]
    g_egfr = d["gnina"]["EGFR/HER2"]
    g_jak = d["gnina"]["JAK1/TYK2"]
    g_pm = d["gnina"]["PIK3CA/mTOR"]
    above, below, cross = smin_cross(d)
    inc_row, inc_abs = max_inc(d)
    excl, incl = desc_excl(d)

    text = text.replace(
        "The original ligand-level non-stratified percentile bootstrap was retained, with 2000 replicates and the recorded deterministic seeds. Each replicate sampled the original number of ligands with replacement from the pooled dual, A-only, and B-only set, recomputed both directional AUROCs using the same resampled dual ligands, and took their minimum. Replicates missing a required class were omitted. The 2.5th and 97.5th percentiles were calculated separately for each directional AUROC and for the replicate-wise minimum. Point estimates used the full analysis sample. A supplementary class-stratified bootstrap preserved the three class sizes and reused the same dual draw in both directions. This post hoc sensitivity analysis is archived with the score tables and did not replace the original intervals. Dual-versus-neither and dual-versus-all-nonduals intervals used class-stratified percentile bootstrap. Intervals were not adjusted for multiple comparisons.",
        METHODS_EN.split("Table 2 reports pointwise 95% confidence intervals for both directional AUROCs and their descriptive minimum. ", 1)[1],
    )
    # If METHODS already starts with Table 2... the replace above only hits the old tail.
    text = text.replace(
        "Table 2 reports pointwise 95% confidence intervals for both directional AUROCs and their descriptive minimum. Primary intervals used a class-stratified nonparametric percentile bootstrap",
        "Table 2 reports pointwise 95% confidence intervals for both directional AUROCs and their descriptive minimum. Primary intervals used a class-stratified nonparametric percentile bootstrap",
    )
    if "non-stratified percentile bootstrap was retained" in text:
            text = re.sub(
                r"Table 2 reports pointwise 95% confidence intervals for both directional AUROCs and their descriptive minimum\..*?Intervals were not adjusted for multiple comparisons\.",
                lambda _: METHODS_EN,
                text,
                count=1,
                flags=re.S,
            )

    text = patch_table_block(
        text,
        "| Pair | n_scored (dual / A-only / B-only) | dual vs A-only (pocket B) [95% CI] |",
        table2(d, zh=False),
    )
    text = patch_table_block(
        text,
        "| Pair | two-pocket mean D-vs-neither AUROC [95% CI] | n_neither |",
        table3(d, zh=False),
    )

    text = text.replace(
        "the EGFR/HER2 dual-versus-B-only AUROC was 0.324 and rose to 0.786 against neither (difference 0.462 [0.262, 0.651]). Using the JAK1-pocket score, replacing B-only with neither increased the AUROC by 0.444 [0.263, 0.620].",
        f"the EGFR/HER2 dual-versus-B-only AUROC was {r3(egfr_fix['auroc_dual_vs_selective'])} and rose to {r3(egfr_fix['auroc_dual_vs_neither'])} against neither (difference {r3(egfr_fix['delta_neither_minus_selective'])} {ci(egfr_fix['delta_ci_lo'], egfr_fix['delta_ci_hi'])}). Using the JAK1-pocket score, replacing B-only with neither increased the AUROC by {r3(jak_fix['delta_neither_minus_selective'])} {ci(jak_fix['delta_ci_lo'], jak_fix['delta_ci_hi'])}.",
    )
    text = text.replace(
        "The PPARG/PPARA \(\mathrm{summary}_{\min}\) interval lay entirely above 0.5 (0.649 [0.504, 0.751]), the F2/F10 interval lay entirely below 0.5 (0.345 [0.211, 0.477]), and the remaining six pairs crossed 0.5",
        f"The PPARG/PPARA \(\mathrm{{summary}}_{{\min}}\) interval lay entirely above 0.5 ({r3(s_pparg['summary_min'])} {ci(s_pparg['ci_lo'], s_pparg['ci_hi'])}), the EGFR/HER2 and F2/F10 intervals lay entirely below 0.5 ({r3(s_egfr['summary_min'])} {ci(s_egfr['ci_lo'], s_egfr['ci_hi'])} and {r3(s_f2['summary_min'])} {ci(s_f2['ci_lo'], s_f2['ci_hi'])}), and the remaining {len(cross)} pairs crossed 0.5",
    )
    text = text.replace(
        f"EGFR/HER2 reached 0.759 [0.557, 0.923] and JAK1/TYK2 reached 0.770 [0.597, 0.906]",
        f"EGFR/HER2 reached {r3(two_egfr['two_pocket_mean_D_vs_neither'])} {ci(two_egfr['ci_lo'], two_egfr['ci_hi'])} and JAK1/TYK2 reached {r3(two_jak['two_pocket_mean_D_vs_neither'])} {ci(two_jak['ci_lo'], two_jak['ci_hi'])}",
    )
    text = text.replace(
        "TPSA alone gave dual-versus-A-only and dual-versus-B-only AUROCs of 0.733 and 0.801",
        f"TPSA alone gave dual-versus-A-only and dual-versus-B-only AUROCs of {r3(ache_t['tpsa_D_vs_A'])} and {r3(ache_t['tpsa_D_vs_B'])}",
    )
    text = text.replace(
        "six of eight 95% intervals included 0, whereas F2/F10 and JAK1/TYK2 excluded 0",
        f"{len(incl)} of eight 95% intervals included 0, whereas {', '.join(excl)} excluded 0",
    )
    text = text.replace(
        "changed AUROC by at most 0.023 across the 16 directional comparisons",
        f"changed AUROC by at most {r3(inc_abs)} across the 16 directional comparisons",
    )
    text = text.replace(
        "changed AUROC by at most 0.011 across the 16 directional comparisons",
        f"changed AUROC by at most {r3(inc_abs)} across the 16 directional comparisons",
    )
    text = text.replace(
        "changed AUROC by at most 0.023 across 16 directions",
        f"changed AUROC by at most {r3(inc_abs)} across 16 directions",
    )
    text = text.replace(
        "changed AUROC by at most 0.011 across 16 directions",
        f"changed AUROC by at most {r3(inc_abs)} across 16 directions",
    )
    text = text.replace(
        "only AChE/BChE had a matched-minus-mismatched \(\mathrm{summary}_{\min}\) 95% interval that excluded zero (0.177 [0.050, 0.297]). EGFR/HER2 was 0.056 [−0.044, 0.160]",
        f"only AChE/BChE had a matched-minus-mismatched \(\mathrm{{summary}}_{{\min}}\) 95% interval that excluded zero ({r3(mm_ache['delta'])} {ci(mm_ache['delta_ci_lo'], mm_ache['delta_ci_hi'])}). EGFR/HER2 was {r3(mm_egfr['delta'])} {ci(mm_egfr['delta_ci_lo'], mm_egfr['delta_ci_hi'])}",
    )
    text = text.replace(
        "only AChE/BChE had a matched-minus-mismatched \(\mathrm{summary}_{\min}\) 95% interval that excluded zero (0.177 [0.053, 0.291]). EGFR/HER2 was 0.056 [−0.029, 0.157]",
        f"only AChE/BChE had a matched-minus-mismatched \(\mathrm{{summary}}_{{\min}}\) 95% interval that excluded zero ({r3(mm_ache['delta'])} {ci(mm_ache['delta_ci_lo'], mm_ache['delta_ci_hi'])}). EGFR/HER2 was {r3(mm_egfr['delta'])} {ci(mm_egfr['delta_ci_lo'], mm_egfr['delta_ci_hi'])}",
    )
    gnina_en = (
        f"the dual-versus-neither AUROC was {r3(g_egfr['auroc_D_vs_neither_mean'])} {ci(g_egfr['d_vs_neither_ci_lo'], g_egfr['d_vs_neither_ci_hi'])} ($n_{{\mathrm{{neither}}}}={g_egfr['n_neither']}$), whereas dual-versus-B-only was {r3(g_egfr['auroc_D_vs_B_pocketA'])} {ci(g_egfr['auroc_D_vs_B_ci_lo'], g_egfr['auroc_D_vs_B_ci_hi'])}. For JAK1/TYK2, dual-versus-neither was {r3(g_jak['auroc_D_vs_neither_mean'])} {ci(g_jak['d_vs_neither_ci_lo'], g_jak['d_vs_neither_ci_hi'])} and directional \(\mathrm{{summary}}_{{\min}}\) was {r3(g_jak['summary_min'])} {ci(g_jak['summary_min_ci_lo'], g_jak['summary_min_ci_hi'])}. For PIK3CA/mTOR, \(\mathrm{{summary}}_{{\min}}\) was {r3(g_pm['summary_min'])} {ci(g_pm['summary_min_ci_lo'], g_pm['summary_min_ci_hi'])}, the weaker arm was dual-versus-A-only {r3(g_pm['auroc_D_vs_A_pocketB'])} {ci(g_pm['auroc_D_vs_A_ci_lo'], g_pm['auroc_D_vs_A_ci_hi'])}, and dual-versus-neither was {r3(g_pm['auroc_D_vs_neither_mean'])} {ci(g_pm['d_vs_neither_ci_lo'], g_pm['d_vs_neither_ci_hi'])}"
    )
    for old in (
        "the dual-versus-neither AUROC was 0.737 [0.536, 0.903] ($n_{\mathrm{neither}}=11$), whereas dual-versus-B-only was 0.265 [0.148, 0.394]. For JAK1/TYK2, dual-versus-neither was 0.705 [0.517, 0.876] and directional \(\mathrm{summary}_{\min}\) was 0.317 [0.183, 0.463]. For PIK3CA/mTOR, \(\mathrm{summary}_{\min}\) was 0.633, the weaker arm was dual-versus-A-only 0.633 [0.427, 0.825], and dual-versus-neither was 0.569 [0.222, 0.889]",
        "the dual-versus-neither AUROC was 0.737 [0.529, 0.919] ($n_{\mathrm{neither}}=11$), whereas dual-versus-B-only was 0.265 [0.144, 0.401]. For JAK1/TYK2, dual-versus-neither was 0.705 [0.524, 0.872] and directional \(\mathrm{summary}_{\min}\) was 0.317 [0.187, 0.455]. For PIK3CA/mTOR, \(\mathrm{summary}_{\min}\) was 0.633, the weaker arm was dual-versus-A-only 0.633 [0.410, 0.769], and dual-versus-neither was 0.569 [0.236, 0.889]",
        "the dual-versus-neither AUROC was 0.705 [0.465, 0.910] ($n_{\mathrm{neither}}=10$), whereas dual-versus-B-only was 0.227 [0.104, 0.373]. For JAK1/TYK2, dual-versus-neither was 0.705 [0.524, 0.872] and directional \(\mathrm{summary}_{\min}\) was 0.317 [0.187, 0.455]. For PIK3CA/mTOR, \(\mathrm{summary}_{\min}\) was 0.633, the weaker arm was dual-versus-A-only 0.633 [0.410, 0.769], and dual-versus-neither was 0.569 [0.236, 0.889]",
    ):
        text = text.replace(old, gnina_en)
    r4j, r5d, r4s = d["rec"]["4JPS"], d["rec"]["5DXT"], d["rec"]["4JSX"]
    s_pm = d["smin"]["PIK3CA/mTOR"]
    text = text.replace(
        f"replacing PIK3CA 4L23 with 4JPS lowered \(\mathrm{{summary}}_{{\min}}\) from 0.692 [0.470, 0.813] to 0.486 [0.259, 0.692]. Replacement with 5DXT gave 0.505 [0.292, 0.696]. Replacing mTOR 4JT6 with 4JSX gave 0.639 [0.418, 0.776]",
        f"replacing PIK3CA 4L23 with 4JPS lowered \(\mathrm{{summary}}_{{\min}}\) from {r3(s_pm['summary_min'])} {ci(s_pm['ci_lo'], s_pm['ci_hi'])} to {r3(r4j['summary_min'])} {ci(r4j['summary_min_ci_lo'], r4j['summary_min_ci_hi'])}. Replacement with 5DXT gave {r3(r5d['summary_min'])} {ci(r5d['summary_min_ci_lo'], r5d['summary_min_ci_hi'])}. Replacing mTOR 4JT6 with 4JSX gave {r3(r4s['summary_min'])} {ci(r4s['summary_min_ci_lo'], r4s['summary_min_ci_hi'])}",
    )
    text = text.replace(
        "PPARG/PPARA was the only pair whose primary Vina \(\mathrm{summary}_{\min}\) interval lay entirely above 0.5 (0.649 [0.504, 0.751])",
        f"PPARG/PPARA was the only pair whose primary Vina \(\mathrm{{summary}}_{{\min}}\) interval lay entirely above 0.5 ({r3(s_pparg['summary_min'])} {ci(s_pparg['ci_lo'], s_pparg['ci_hi'])})",
    )
    cl_egfr = {r["estimator"]: r for r in d["cluster"] if r["pair"] == "EGFR/HER2"}
    cl_jak = {r["estimator"]: r for r in d["cluster"] if r["pair"] == "JAK1/TYK2"}
    jak_doc = cl_jak.get("document_cluster") or {}
    if jak_doc.get("delta_ci_lo") not in ("", None):
        jak_doc_txt = f"document-cluster was {ci(jak_doc['delta_ci_lo'], jak_doc['delta_ci_hi'])} (includes 0)"
    else:
        jak_doc_txt = (
            "document-cluster could not be recomputed because the ligand–document map and ChEMBL 37 sqlite "
            "are unavailable; a previously deposited interval included 0 and is not treated as a current calculation"
        )
    text = text.replace(
        "The official ligand-level fixed-score difference on EGFR/HER2 pocket A is 0.462 [0.262, 0.651]. Cluster bootstrap was not recomputed after the canonical box correction. JAK1/TYK2 scaffold-cluster was [0.234, 0.633] (excludes 0) and document-cluster was [−0.034, 0.682] (includes 0)",
        f"The ligand-level fixed-score difference on EGFR/HER2 pocket A is {r3(egfr_fix['delta_neither_minus_selective'])} {ci(egfr_fix['delta_ci_lo'], egfr_fix['delta_ci_hi'])}. After substituting the corrected-box scores into the frozen scaffold and document groupings, the EGFR/HER2 cluster intervals were {ci(cl_egfr['scaffold_cluster']['delta_ci_lo'], cl_egfr['scaffold_cluster']['delta_ci_hi'])} (scaffold) and {ci(cl_egfr['document_cluster']['delta_ci_lo'], cl_egfr['document_cluster']['delta_ci_hi'])} (document); both exclude 0. JAK1/TYK2 scaffold-cluster was {ci(cl_jak['scaffold_cluster']['delta_ci_lo'], cl_jak['scaffold_cluster']['delta_ci_hi'])} (excludes 0) and {jak_doc_txt}",
    )
    h_pg = d["hold"]["PPARG/PPARA"]
    text = text.replace(
        "PPARG/PPARA fell from 0.649 to 0.535 [0.350, 0.717]",
        f"PPARG/PPARA fell from {r3(s_pparg['summary_min'])} to {r3(h_pg['summary_min'])} {ci(h_pg['ci_lo'], h_pg['ci_hi'])}",
    )
    text = text.replace(
        "Cluster bootstrap by Bemis–Murcko scaffold and literature-connected document groups was a source-dependence sensitivity for the two largest fixed-score differences (Table S4) and did not replace the ligand-level primary intervals. Intervals were not adjusted for multiple comparisons.",
        "Cluster bootstrap by Bemis–Murcko scaffold and literature-connected document groups was a source-dependence sensitivity for the two largest fixed-score differences (Table S4) and did not replace the ligand-level primary intervals. A binormal detectable-effect simulation reused the same class-stratified shared-dual bootstrap, the current eight-pair class sizes, B = 2000, and seed 20260729; it estimates the probability that a CI excludes 0.5 under a specified true AUROC and is not observed power. Intervals were not adjusted for multiple comparisons.",
    )
    text = _ensure_once(
        text,
        "Complete-case counts and fixed-membership intersections are archived in the repository.",
        " A binormal detectable-effect simulation used the current eight-pair class sizes and the same class-stratified shared-dual bootstrap as Table 2.",
    )
    # Abstract
    text = text.replace(
        "On EGFR/HER2, the EGFR-pocket AUROC for dual-versus-B-only was 0.324 and rose to 0.786 against neither (difference 0.462 [0.262, 0.651]); JAK1/TYK2 showed a similar difference (0.444 [0.263, 0.620]).",
        f"On EGFR/HER2, the EGFR-pocket AUROC for dual-versus-B-only was {r3(egfr_fix['auroc_dual_vs_selective'])} and rose to {r3(egfr_fix['auroc_dual_vs_neither'])} against neither (difference {r3(egfr_fix['delta_neither_minus_selective'])} {ci(egfr_fix['delta_ci_lo'], egfr_fix['delta_ci_hi'])}); JAK1/TYK2 showed a similar difference ({r3(jak_fix['delta_neither_minus_selective'])} {ci(jak_fix['delta_ci_lo'], jak_fix['delta_ci_hi'])}).",
    )
    det_en = _detectable_sentence(d, zh=False)
    if det_en:
        text = re.sub(
            r"Under a binormal simulation that reused these class sizes and the Table 2 bootstrap, "
            r"the probability that a `summary_min` CI excludes 0\.5 was at most 0\.\d+ when the true weaker-arm AUROC was 0\.60, "
            r"and at least 0\.\d+ at 0\.75 except for PIK3CA/mTOR \(0\.\d+\)\. "
            r"The simulation is not observed power \(Table S4\)\.",
            det_en.rstrip(),
            text,
        )
        if "Under a binormal simulation that reused these class sizes" not in text:
            text = text.replace(
                "These values are pair-specific and are not an eight-pair ranking.",
                "These values are pair-specific and are not an eight-pair ranking. " + det_en,
            )
    text = text.replace("adding the matched-pocket docking score changed AUROC by at most 0.023.", f"adding the matched-pocket docking score changed AUROC by at most {r3(inc_abs)}.")
    rk_e = d["rank"]["EGFR/HER2"]
    rk_a = d["rank"]["AChE/BChE"]
    text = text.replace("and 0.357 (EGFR/HER2)", f"and {r3(rk_e['ef_dual_10pct'])} (EGFR/HER2)")
    text = text.replace(
        "\\(\\mathrm{EF}_{\\mathrm{dual},10\\%}=0.357\\)",
        "\\(\\mathrm{EF}_{\\mathrm{dual},10\\%}=" + r3(rk_e["ef_dual_10pct"]) + "\\)",
    )
    text = text.replace(
        "\\(\\mathrm{EF}_{\\mathrm{dual},10\\%}=1.778\\)",
        "\\(\\mathrm{EF}_{\\mathrm{dual},10\\%}=" + r3(rk_a["ef_dual_10pct"]) + "\\)",
    )
    return apply_boundary_patches(text, d, zh=False)


def patch_zh(text, d):
    egfr_fix = d["fixed"][("EGFR/HER2", "D_vs_B_or_neither_pocketA")]
    jak_fix = d["fixed"][("JAK1/TYK2", "D_vs_B_or_neither_pocketA")]
    s_egfr = d["smin"]["EGFR/HER2"]
    s_pparg = d["smin"]["PPARG/PPARA"]
    s_f2 = d["smin"]["F2/F10"]
    two_egfr = d["two"]["EGFR/HER2"]
    two_jak = d["two"]["JAK1/TYK2"]
    mm_egfr = d["mm"]["EGFR/HER2"]
    mm_ache = d["mm"]["AChE/BChE"]
    ache_t = d["desc"]["AChE/BChE"]
    g_egfr = d["gnina"]["EGFR/HER2"]
    g_jak = d["gnina"]["JAK1/TYK2"]
    g_pm = d["gnina"]["PIK3CA/mTOR"]
    _, _, cross = smin_cross(d)
    _, inc_abs = max_inc(d)
    excl, incl = desc_excl(d)
    if "主分析保留 2000 次配体水平非分层百分位 bootstrap" in text:
        text = re.sub(
            r"Table 2 同时报告两条方向性 AUROC 及其描述性最小值的逐项 95% 置信区间。.*?区间未做多重比较校正。",
            lambda _: METHODS_ZH,
            text,
            count=1,
            flags=re.S,
        )
    text = patch_table_block(
        text,
        "| 靶对 | n_scored (dual / A-only / B-only) | Dual vs A-only（口袋 B）[95% CI] |",
        table2(d, zh=True),
    )
    text = patch_table_block(
        text,
        "| 靶对 | 双口袋平均 D-vs-neither AUROC [95% CI] | n_neither |",
        table3(d, zh=True),
    )
    text = text.replace(
        "EGFR/HER2 使用 EGFR 口袋评分时，dual–B-only 的 AUROC 为 0.324，将对照换成 neither 后升至 0.786，差值 0.462 [0.262, 0.651]。JAK1/TYK2 使用 JAK1 口袋评分时，对应差值为 0.444 [0.263, 0.620]。",
        f"EGFR/HER2 使用 EGFR 口袋评分时，dual–B-only 的 AUROC 为 {r3(egfr_fix['auroc_dual_vs_selective'])}，将对照换成 neither 后升至 {r3(egfr_fix['auroc_dual_vs_neither'])}，差值 {r3(egfr_fix['delta_neither_minus_selective'])} {ci(egfr_fix['delta_ci_lo'], egfr_fix['delta_ci_hi'])}。JAK1/TYK2 使用 JAK1 口袋评分时，对应差值为 {r3(jak_fix['delta_neither_minus_selective'])} {ci(jak_fix['delta_ci_lo'], jak_fix['delta_ci_hi'])}。",
    )
    text = text.replace(
        "PPARG/PPARA 的 \(\mathrm{summary}_{\min}\) 区间完全高于 0.5（0.649 [0.504, 0.751]），F2/F10 完全低于 0.5（0.345 [0.211, 0.477]），其余六对穿过 0.5",
        f"PPARG/PPARA 的 \(\mathrm{{summary}}_{{\min}}\) 区间完全高于 0.5（{r3(s_pparg['summary_min'])} {ci(s_pparg['ci_lo'], s_pparg['ci_hi'])}），EGFR/HER2 与 F2/F10 完全低于 0.5（{r3(s_egfr['summary_min'])} {ci(s_egfr['ci_lo'], s_egfr['ci_hi'])} 和 {r3(s_f2['summary_min'])} {ci(s_f2['ci_lo'], s_f2['ci_hi'])}），其余 {len(cross)} 对穿过 0.5",
    )
    text = text.replace(
        "EGFR/HER2 为 0.759 [0.557, 0.923]，JAK1/TYK2 为 0.770 [0.597, 0.906]",
        f"EGFR/HER2 为 {r3(two_egfr['two_pocket_mean_D_vs_neither'])} {ci(two_egfr['ci_lo'], two_egfr['ci_hi'])}，JAK1/TYK2 为 {r3(two_jak['two_pocket_mean_D_vs_neither'])} {ci(two_jak['ci_lo'], two_jak['ci_hi'])}",
    )
    text = text.replace(
        "AChE/BChE 仅用 TPSA 时，dual–A-only 和 dual–B-only 的 AUROC 分别为 0.733 和 0.801",
        f"AChE/BChE 仅用 TPSA 时，dual–A-only 和 dual–B-only 的 AUROC 分别为 {r3(ache_t['tpsa_D_vs_A'])} 和 {r3(ache_t['tpsa_D_vs_B'])}",
    )
    text = text.replace(
        "八对中六对的 95% 区间包含 0，F2/F10 与 JAK1/TYK2 不包含 0",
        f"八对中 {len(incl)} 对的 95% 区间包含 0，{ '、'.join(excl) } 不包含 0",
    )
    text = text.replace("16 个方向的 AUROC 最多变化 0.023", f"16 个方向的 AUROC 最多变化 {r3(inc_abs)}")
    text = text.replace("16 个方向的 AUROC 最多变化 0.011", f"16 个方向的 AUROC 最多变化 {r3(inc_abs)}")
    text = text.replace("AUROC 最大绝对变化仅为 0.023", f"AUROC 最大绝对变化为 {r3(inc_abs)}")
    text = text.replace("AUROC 最大绝对变化为 0.023", f"AUROC 最大绝对变化为 {r3(inc_abs)}")
    text = text.replace("AUROC 最大绝对变化为 0.011", f"AUROC 最大绝对变化为 {r3(inc_abs)}")
    text = text.replace(
        "仅 AChE/BChE 的 matched−mismatched \(\mathrm{summary}_{\min}\) 95% 区间排除 0（0.177 [0.050, 0.297]）。EGFR/HER2 为 0.056 [−0.044, 0.160]",
        f"仅 AChE/BChE 的 matched−mismatched \(\mathrm{{summary}}_{{\min}}\) 95% 区间排除 0（{r3(mm_ache['delta'])} {ci(mm_ache['delta_ci_lo'], mm_ache['delta_ci_hi'])}）。EGFR/HER2 为 {r3(mm_egfr['delta'])} {ci(mm_egfr['delta_ci_lo'], mm_egfr['delta_ci_hi'])}",
    )
    text = text.replace(
        "仅 AChE/BChE 的 matched−mismatched \(\mathrm{summary}_{\min}\) 95% 区间排除 0（0.177 [0.053, 0.291]）。EGFR/HER2 为 0.056 [−0.029, 0.157]",
        f"仅 AChE/BChE 的 matched−mismatched \(\mathrm{{summary}}_{{\min}}\) 95% 区间排除 0（{r3(mm_ache['delta'])} {ci(mm_ache['delta_ci_lo'], mm_ache['delta_ci_hi'])}）。EGFR/HER2 为 {r3(mm_egfr['delta'])} {ci(mm_egfr['delta_ci_lo'], mm_egfr['delta_ci_hi'])}",
    )
    gnina_zh = (
        f"EGFR/HER2 的 dual–neither AUROC 为 {r3(g_egfr['auroc_D_vs_neither_mean'])} {ci(g_egfr['d_vs_neither_ci_lo'], g_egfr['d_vs_neither_ci_hi'])}（\(n_{{\mathrm{{neither}}}}={g_egfr['n_neither']}\)），dual–B-only 为 {r3(g_egfr['auroc_D_vs_B_pocketA'])} {ci(g_egfr['auroc_D_vs_B_ci_lo'], g_egfr['auroc_D_vs_B_ci_hi'])}。JAK1/TYK2 的 dual–neither 为 {r3(g_jak['auroc_D_vs_neither_mean'])} {ci(g_jak['d_vs_neither_ci_lo'], g_jak['d_vs_neither_ci_hi'])}，方向性 \(\mathrm{{summary}}_{{\min}}\) 为 {r3(g_jak['summary_min'])} {ci(g_jak['summary_min_ci_lo'], g_jak['summary_min_ci_hi'])}。PIK3CA/mTOR 的 \(\mathrm{{summary}}_{{\min}}\) 为 {r3(g_pm['summary_min'])} {ci(g_pm['summary_min_ci_lo'], g_pm['summary_min_ci_hi'])}，较弱臂 dual–A-only 为 {r3(g_pm['auroc_D_vs_A_pocketB'])} {ci(g_pm['auroc_D_vs_A_ci_lo'], g_pm['auroc_D_vs_A_ci_hi'])}，dual–neither 为 {r3(g_pm['auroc_D_vs_neither_mean'])} {ci(g_pm['d_vs_neither_ci_lo'], g_pm['d_vs_neither_ci_hi'])}"
    )
    for old in (
        "EGFR/HER2 的 dual–neither AUROC 为 0.737 [0.536, 0.903]（\(n_{\mathrm{neither}}=11\)），dual–B-only 为 0.265 [0.148, 0.394]。JAK1/TYK2 的 dual–neither 为 0.705 [0.517, 0.876]，方向性 \(\mathrm{summary}_{\min}\) 为 0.317 [0.183, 0.463]。PIK3CA/mTOR 的 \(\mathrm{summary}_{\min}\) 为 0.633，较弱臂 dual–A-only 为 0.633 [0.427, 0.825]，dual–neither 为 0.569 [0.222, 0.889]",
        "EGFR/HER2 的 dual–neither AUROC 为 0.737 [0.529, 0.919]（\(n_{\mathrm{neither}}=11\)），dual–B-only 为 0.265 [0.144, 0.401]。JAK1/TYK2 的 dual–neither 为 0.705 [0.524, 0.872]，方向性 \(\mathrm{summary}_{\min}\) 为 0.317 [0.187, 0.455]。PIK3CA/mTOR 的 \(\mathrm{summary}_{\min}\) 为 0.633，较弱臂 dual–A-only 为 0.633 [0.410, 0.769]，dual–neither 为 0.569 [0.236, 0.889]",
        "EGFR/HER2 的 dual–neither AUROC 为 0.705 [0.465, 0.910]（\(n_{\mathrm{neither}}=10\)），dual–B-only 为 0.227 [0.104, 0.373]。JAK1/TYK2 的 dual–neither 为 0.705 [0.524, 0.872]，方向性 \(\mathrm{summary}_{\min}\) 为 0.317 [0.187, 0.455]。PIK3CA/mTOR 的 \(\mathrm{summary}_{\min}\) 为 0.633，较弱臂 dual–A-only 为 0.633 [0.410, 0.769]，dual–neither 为 0.569 [0.236, 0.889]",
    ):
        text = text.replace(old, gnina_zh)
    r4j, r5d, r4s = d["rec"]["4JPS"], d["rec"]["5DXT"], d["rec"]["4JSX"]
    s_pm = d["smin"]["PIK3CA/mTOR"]
    text = text.replace(
        "将 PIK3CA 受体由 4L23 替换为 4JPS 后，\(\mathrm{summary}_{\min}\) 从 0.692 [0.470, 0.813] 降至 0.486 [0.259, 0.692]；替换为 5DXT 后为 0.505 [0.292, 0.696]；将 mTOR 4JT6 替换为 4JSX 后为 0.639 [0.418, 0.776]",
        f"将 PIK3CA 受体由 4L23 替换为 4JPS 后，\(\mathrm{{summary}}_{{\min}}\) 从 {r3(s_pm['summary_min'])} {ci(s_pm['ci_lo'], s_pm['ci_hi'])} 降至 {r3(r4j['summary_min'])} {ci(r4j['summary_min_ci_lo'], r4j['summary_min_ci_hi'])}；替换为 5DXT 后为 {r3(r5d['summary_min'])} {ci(r5d['summary_min_ci_lo'], r5d['summary_min_ci_hi'])}；将 mTOR 4JT6 替换为 4JSX 后为 {r3(r4s['summary_min'])} {ci(r4s['summary_min_ci_lo'], r4s['summary_min_ci_hi'])}",
    )
    text = text.replace(
        "PPARG/PPARA 是唯一主分析 Vina \(\mathrm{summary}_{\min}\) 区间完全高于 0.5 的靶对（0.649 [0.504, 0.751]）",
        f"PPARG/PPARA 是唯一主分析 Vina \(\mathrm{{summary}}_{{\min}}\) 区间完全高于 0.5 的靶对（{r3(s_pparg['summary_min'])} {ci(s_pparg['ci_lo'], s_pparg['ci_hi'])}）",
    )
    cl_egfr = {r["estimator"]: r for r in d["cluster"] if r["pair"] == "EGFR/HER2"}
    cl_jak = {r["estimator"]: r for r in d["cluster"] if r["pair"] == "JAK1/TYK2"}
    jak_doc = cl_jak.get("document_cluster") or {}
    if jak_doc.get("delta_ci_lo") not in ("", None):
        jak_doc_txt = f"文献簇区间为 {ci(jak_doc['delta_ci_lo'], jak_doc['delta_ci_hi'])}，包含 0"
    else:
        jak_doc_txt = (
            "文献簇无法在本冻结中重算（配体–文献分组映射与 ChEMBL 37 sqlite 均不可用）；"
            "此前存档区间包含 0，不作为当前计算结果"
        )
    text = text.replace(
        "EGFR/HER2 口袋 A 的正式配体层固定评分差值为 0.462 [0.262, 0.651]。规范盒校正后未重算簇 bootstrap。JAK1/TYK2 骨架簇区间为 [0.234, 0.633]，排除 0，文献簇区间为 [−0.034, 0.682]，包含 0",
        f"EGFR/HER2 口袋 A 的配体层固定评分差值为 {r3(egfr_fix['delta_neither_minus_selective'])} {ci(egfr_fix['delta_ci_lo'], egfr_fix['delta_ci_hi'])}。将校正盒评分代入原冻结的骨架簇与文献簇后，EGFR/HER2 簇区间为 {ci(cl_egfr['scaffold_cluster']['delta_ci_lo'], cl_egfr['scaffold_cluster']['delta_ci_hi'])}（骨架）和 {ci(cl_egfr['document_cluster']['delta_ci_lo'], cl_egfr['document_cluster']['delta_ci_hi'])}（文献），均排除 0。JAK1/TYK2 骨架簇区间为 {ci(cl_jak['scaffold_cluster']['delta_ci_lo'], cl_jak['scaffold_cluster']['delta_ci_hi'])}，排除 0，{jak_doc_txt}",
    )
    h_pg = d["hold"]["PPARG/PPARA"]
    text = text.replace(
        "PPARG/PPARA 则由 0.649 降至 0.535 [0.350, 0.717]",
        f"PPARG/PPARA 则由 {r3(s_pparg['summary_min'])} 降至 {r3(h_pg['summary_min'])} {ci(h_pg['ci_lo'], h_pg['ci_hi'])}",
    )
    text = text.replace(
        "EGFR/HER2 中 EGFR 口袋评分对 dual–B-only 的 AUROC 为 0.324，对 neither 升至 0.786（差值 0.462 [0.262, 0.651]）；JAK1/TYK2 出现类似差值（0.444 [0.263, 0.620]）。",
        f"EGFR/HER2 中 EGFR 口袋评分对 dual–B-only 的 AUROC 为 {r3(egfr_fix['auroc_dual_vs_selective'])}，对 neither 升至 {r3(egfr_fix['auroc_dual_vs_neither'])}（差值 {r3(egfr_fix['delta_neither_minus_selective'])} {ci(egfr_fix['delta_ci_lo'], egfr_fix['delta_ci_hi'])}）；JAK1/TYK2 出现类似差值（{r3(jak_fix['delta_neither_minus_selective'])} {ci(jak_fix['delta_ci_lo'], jak_fix['delta_ci_hi'])}）。",
    )
    text = text.replace("AUROC 最多变化 0.023。", f"AUROC 最多变化 {r3(inc_abs)}。")
    rk_e = d["rank"]["EGFR/HER2"]
    rk_a = d["rank"]["AChE/BChE"]
    text = text.replace(
        "按 Bemis–Murcko 骨架簇和文献连通簇进行的簇 bootstrap 用于两个最大固定评分差的来源依赖性敏感性（Table S4），不替换配体水平主区间。区间未做多重比较校正。",
        "按 Bemis–Murcko 骨架簇和文献连通簇进行的簇 bootstrap 用于两个最大固定评分差的来源依赖性敏感性（Table S4），不替换配体水平主区间。二项正态 detectable-effect 仿真复用同一类别分层、共享 dual 的 bootstrap、当前八对类别样本量和 B = 2000、种子 20260729；它估计在指定真 AUROC 下 CI 排除 0.5 的概率，不是观测功效。区间未做多重比较校正。",
    )
    text = _ensure_once(
        text,
        "完整病例计数和固定成员交集留在仓库。",
        "二项正态 detectable-effect 仿真使用当前八对类别样本量和与 Table 2 相同的类别分层、共享 dual bootstrap。",
    )
    det_zh = _detectable_sentence(d, zh=True)
    if det_zh:
        text = re.sub(
            r"在复用当前类别样本量和 Table 2 同类 bootstrap 的二项正态仿真中，当真较弱臂 AUROC 为 0\.60 时，"
            r"`summary_min` CI 排除 0\.5 的概率最高为 0\.\d+；为 0\.75 时，除 PIK3CA/mTOR（0\.\d+）外均不低于 0\.\d+。"
            r"该仿真不是观测功效（Table S4）。",
            det_zh.rstrip(),
            text,
        )
    if det_zh and "在复用当前类别样本量和 Table 2 同类 bootstrap" not in text:
        text = text.replace(
            "这些数值是靶对特异的，不是八对排行。",
            "这些数值是靶对特异的，不是八对排行。" + det_zh,
            1,
        )
    rk_e = d["rank"]["EGFR/HER2"]
    rk_a = d["rank"]["AChE/BChE"]
    text = text.replace("和 0.357（EGFR/HER2）", f"和 {r3(rk_e['ef_dual_10pct'])}（EGFR/HER2）")
    text = text.replace(
        "\\(\\mathrm{EF}_{\\mathrm{dual},10\\%}=0.357\\)",
        "\\(\\mathrm{EF}_{\\mathrm{dual},10\\%}=" + r3(rk_e["ef_dual_10pct"]) + "\\)",
    )
    text = text.replace(
        "\\(\\mathrm{EF}_{\\mathrm{dual},10\\%}=1.778\\)",
        "\\(\\mathrm{EF}_{\\mathrm{dual},10\\%}=" + r3(rk_a["ef_dual_10pct"]) + "\\)",
    )
    return apply_boundary_patches(text, d, zh=True)


def si_s3(d):
    rows = ["| Pair | Label rule | n (D / A / B) | summary_min | 95% CI |", "|------|----------|--------------:|------------:|--------|"]
    show = {
        "EGFR/HER2": ["theta_5.5", "theta_6.0", "theta_6.5", "strict_6.5_5.5"],
        "PIK3CA/mTOR": ["theta_5.5", "theta_6.0", "theta_6.5", "strict_6.5_5.5"],
    }
    labels = {(r["pair"], r["label_rule"]): r for r in d["labels"]}
    name = {"theta_5.5": "θ = 5.5", "theta_6.0": "θ = 6.0", "theta_6.5": "θ = 6.5", "strict_6.5_5.5": "strict 6.5/5.5"}
    for p in PAIRS:
        rules = show.get(p, ["theta_6.0"])
        for rule in rules:
            r = labels[(p, rule)]
            rows.append(
                f"| {p} | {name[rule]} | {r['n_dual']} / {r['n_A_only']} / {r['n_B_only']} | {r3(r['summary_min'])} | {ci(r['ci_lo'], r['ci_hi'])} |"
            )
    return "\n".join(rows)


def si_s4_fixed(d):
    rows = [
        "| Pair | Score channel | dual vs selective | dual vs neither | Δ | 95% CI | neither underpowered |",
        "|------|----------|---------------:|----------------:|--:|--------|:----------------:|",
    ]
    order = [
        ("EGFR/HER2", "D_vs_B_or_neither_pocketA", "pocket A (vs B-only)"),
        ("EGFR/HER2", "D_vs_A_or_neither_pocketB", "pocket B (vs A-only)"),
        ("AChE/BChE", "D_vs_B_or_neither_pocketA", "pocket A"),
        ("AChE/BChE", "D_vs_A_or_neither_pocketB", "pocket B"),
        ("PIK3CA/mTOR", "D_vs_B_or_neither_pocketA", "pocket A"),
        ("PIK3CA/mTOR", "D_vs_A_or_neither_pocketB", "pocket B"),
        ("F2/F10", "D_vs_B_or_neither_pocketA", "pocket A"),
        ("F2/F10", "D_vs_A_or_neither_pocketB", "pocket B"),
        ("JAK1/TYK2", "D_vs_B_or_neither_pocketA", "pocket A"),
        ("JAK1/TYK2", "D_vs_A_or_neither_pocketB", "pocket B"),
        ("JAK1/JAK2", "D_vs_B_or_neither_pocketA", "pocket A"),
        ("JAK1/JAK2", "D_vs_A_or_neither_pocketB", "pocket B"),
        ("PPARG/PPARA", "D_vs_B_or_neither_pocketA", "pocket A"),
        ("PPARG/PPARA", "D_vs_A_or_neither_pocketB", "pocket B"),
        ("PPARA/PPARD", "D_vs_B_or_neither_pocketA", "pocket A"),
        ("PPARA/PPARD", "D_vs_A_or_neither_pocketB", "pocket B"),
    ]
    for p, c, lab in order:
        r = d["fixed"][(p, c)]
        und = "yes" if int(r["n_neither"]) < 10 else "no"
        rows.append(
            f"| {p} | {lab} | {r3(r['auroc_dual_vs_selective'])} | {r3(r['auroc_dual_vs_neither'])} | "
            f"{r3m(r['delta_neither_minus_selective'])} | {ci(r['delta_ci_lo'], r['delta_ci_hi'])} | {und} |"
        )
    return "\n".join(rows)


def _detectable_sentence(d, zh=False) -> str:
    rows = d.get("det") or []
    if not rows:
        return ""
    vals_60 = []
    vals_75 = []
    pik_75 = None
    for r in rows:
        if r["contrast"] != "summary_min":
            continue
        p = float(r["p_ci_excludes_0p5"])
        if r["true_auroc"] == "0.60":
            vals_60.append(p)
        elif r["true_auroc"] == "0.75":
            if r["pair"] == "PIK3CA/mTOR":
                pik_75 = p
            else:
                vals_75.append(p)
    if not vals_60 or not vals_75 or pik_75 is None:
        return ""
    if zh:
        return (
            f"在复用当前类别样本量和 Table 2 同类 bootstrap 的二项正态仿真中，当真较弱臂 AUROC 为 0.60 时，"
            f"`summary_min` CI 排除 0.5 的概率最高为 {r3(max(vals_60))}；为 0.75 时，除 PIK3CA/mTOR（{r3(pik_75)}）外均不低于 {r3(min(vals_75))}。"
            "该仿真不是观测功效（Table S4）。"
        )
    return (
        f"Under a binormal simulation that reused these class sizes and the Table 2 bootstrap, "
        f"the probability that a `summary_min` CI excludes 0.5 was at most {r3(max(vals_60))} when the true weaker-arm AUROC was 0.60, "
        f"and at least {r3(min(vals_75))} at 0.75 except for PIK3CA/mTOR ({r3(pik_75)}). "
        "The simulation is not observed power (Table S4)."
    )


def si_detectable(d, zh=False):
    rows = d.get("det") or []
    if not rows:
        return ""
    by = {(r["pair"], r["contrast"], r["true_auroc"]): r for r in rows}
    smin = d["smin"]
    if zh:
        lines = [
            "**Detectable-effect 仿真（二项正态；不是观测功效）。** 类别样本量取自当前 `current_score_master.csv` 的 complete-case 主集。内部 95% CI 与 Table 2 相同：类别分层、共享 dual 的百分位 bootstrap（B = 2000，种子 20260729）。N_MC = 1000。表中为 `summary_min` 区间排除 0.5 的蒙特卡洛概率。该仿真不替代 Table 2。",
            "",
            "| 靶对 | n (D / A / B) | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 |",
            "|------|--------------:|-----:|-----:|-----:|-----:|-----:|",
        ]
        src = "源：`results/canonical/detectable_effect_simulation.csv`。"
    else:
        lines = [
            "**Detectable-effect simulation (binormal; not observed power).** Class sizes are the current complete-case main-panel counts in `current_score_master.csv`. Inner 95% CIs use the same class-stratified shared-dual percentile bootstrap as Table 2 (B = 2000, seed 20260729). N_MC = 1000. Cells are the Monte Carlo probability that the `summary_min` CI excludes 0.5. This simulation does not replace Table 2.",
            "",
            "| Pair | n (D / A / B) | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 |",
            "|------|--------------:|-----:|-----:|-----:|-----:|-----:|",
        ]
        src = "Source: `results/canonical/detectable_effect_simulation.csv`."
    for pair in PAIRS:
        rec = smin[pair]
        ns = f"{rec['n_dual']} / {rec['n_A_only']} / {rec['n_B_only']}"
        cells = []
        for auc in ("0.55", "0.60", "0.65", "0.70", "0.75"):
            hit = by.get((pair, "summary_min", auc))
            if hit is None:
                return ""
            cells.append(r3(hit["p_ci_excludes_0p5"]))
        lines.append(f"| {pair} | {ns} | {' | '.join(cells)} |")
    lines += ["", src]
    return "\n".join(lines)


def si_s4_cluster(d):
    rows = [
        "| Pair | Resampling unit | Δ point | 95% CI | CI excludes 0 |",
        "|------|-----------------|--------:|--------|:-------------:|",
    ]
    labels = {
        "ligand_stratified": "ligand-level (primary)",
        "scaffold_cluster": "scaffold cluster",
        "document_cluster": "document cluster",
    }
    for r in d["cluster"]:
        if r.get("status") == "unresolved_mapping_unavailable" or r.get("delta_ci_lo") in ("", None):
            rows.append(
                f"| {r['pair']} | {labels.get(r['estimator'], r['estimator'])} | {r3(r['delta_point'])} | "
                f"not recomputed | — |"
            )
            continue
        rows.append(
            f"| {r['pair']} | {labels.get(r['estimator'], r['estimator'])} | {r3(r['delta_point'])} | "
            f"{ci(r['delta_ci_lo'], r['delta_ci_hi'])} | {'yes' if str(r['excludes_zero']) in {'1', 'True'} else 'no'} |"
        )
    return "\n".join(rows)


def si_s5_ecfp(d):
    rows = [
        "| Pair | Arm | ECFP4 | ECFP4+docking | Δ | Vina ranking AUROC (Table 2) |",
        "|------|------|------:|--------------:|--:|--------------------------:|",
    ]
    for p in PAIRS:
        for arm, lab in (("D_vs_A", "D vs A"), ("D_vs_B", "D vs B")):
            r = d["inc"][(p, arm)]
            rows.append(
                f"| {p} | {lab} | {r3(r['cv_auroc_ECFP4'])} | {r3(r['cv_auroc_ECFP4_docking'])} | "
                f"{signed(r['delta_ECFP4_plus_docking_minus_ECFP4'])} | {r3(r['rank_auroc_docking'])} |"
            )
    return "\n".join(rows)


def si_s5_desc(d):
    rows = [
        "| Pair | Best descriptor | Descriptor summary_min | Δ | 95% CI | CI excludes 0 |",
        "|------|-----------------|-----------------------:|--:|--------|:-------------:|",
    ]
    names = {"clogp": "cLogP", "tpsa": "TPSA", "heavy": "heavy", "mw": "MW"}
    for p in PAIRS:
        r = d["desc"][p]
        excl = "yes" if (float(r["delta_ci_lo"]) > 0 or float(r["delta_ci_hi"]) < 0) else "no"
        rows.append(
            f"| {p} | {names.get(r['best_descriptor'], r['best_descriptor'])} | {r3(r['best_descriptor_summary_min'])} | "
            f"{r3m(r['vina_minus_best_descriptor'])} | {ci(r['delta_ci_lo'], r['delta_ci_hi'])} | {excl} |"
        )
    return "\n".join(rows)


def si_s5_nested(d, zh=False):
    if zh:
        rows = [
            "| 靶对 | D vs A OOF AUROC [95% CI] | D vs B OOF AUROC [95% CI] | summary_min [95% CI] | 各折所选描述符 (A / B) |",
            "|------|---------------------------:|---------------------------:|----------------------|------------------------|",
        ]
    else:
        rows = [
            "| Pair | D vs A OOF AUROC [95% CI] | D vs B OOF AUROC [95% CI] | summary_min [95% CI] | Descriptors by fold (A / B) |",
            "|------|---------------------------:|---------------------------:|----------------------|-----------------------------|",
        ]
    names = {"clogp": "cLogP", "tpsa": "TPSA", "heavy": "heavy", "mw": "MW"}
    for p in PAIRS:
        r = d["desc"][p]
        a_mode = names.get(r.get("nested_selected_D_vs_A_mode", ""), r.get("nested_selected_D_vs_A_mode", ""))
        b_mode = names.get(r.get("nested_selected_D_vs_B_mode", ""), r.get("nested_selected_D_vs_B_mode", ""))
        rows.append(
            f"| {p} | {r3(r['nested_scaffold_cv_oof_D_vs_A'])} {ci(r['nested_scaffold_cv_oof_D_vs_A_ci_lo'], r['nested_scaffold_cv_oof_D_vs_A_ci_hi'])} | "
            f"{r3(r['nested_scaffold_cv_oof_D_vs_B'])} {ci(r['nested_scaffold_cv_oof_D_vs_B_ci_lo'], r['nested_scaffold_cv_oof_D_vs_B_ci_hi'])} | "
            f"{r3(r['nested_scaffold_cv_oof_summary_min'])} {ci(r['nested_scaffold_cv_oof_summary_min_ci_lo'], r['nested_scaffold_cv_oof_summary_min_ci_hi'])} | "
            f"{a_mode} / {b_mode} |"
        )
    return "\n".join(rows)


def _gnina_counts(gp) -> str:
    return f"{gp['n_dual']} / {gp['n_A_only']} / {gp['n_B_only']} / {gp['n_neither']}"


def _gnina_weaker(gp, zh=False):
    da = float(gp["auroc_D_vs_A_pocketB"])
    db = float(gp["auroc_D_vs_B_pocketA"])
    if db <= da:
        lab = "dual–B-only（口袋 A）" if zh else "dual–B-only (pocket A)"
        return lab, gp["auroc_D_vs_B_pocketA"], gp["auroc_D_vs_B_ci_lo"], gp["auroc_D_vs_B_ci_hi"]
    lab = "dual–A-only（口袋 B）" if zh else "dual–A-only (pocket B)"
    return lab, gp["auroc_D_vs_A_pocketB"], gp["auroc_D_vs_A_ci_lo"], gp["auroc_D_vs_A_ci_hi"]


def si_s6(d):
    rows = [
        "| Pair | Set | Δ | 95% CI | CI excludes 0 | Weaker arm switched |",
        "|------|------|--:|--------|:---------:|:-------------------:|",
    ]
    for p in PAIRS:
        r = d["mm"][p]
        rows.append(
            f"| {p} | main | {r3m(r['delta'])} | {ci(r['delta_ci_lo'], r['delta_ci_hi'])} | "
            f"{'yes' if str(r['ci_excludes_zero']) in {'1','True'} else 'no'} | "
            f"{'yes' if str(r['weaker_arm_switched']) in {'1','True'} else 'no'} |"
        )
    for p in PAIRS:
        if p not in d["hold_wp"]:
            continue
        r = d["hold_wp"][p]
        rows.append(
            f"| {p} | holdout | {r3m(r['delta_matched_minus_wrong'])} | {ci(r['delta_ci_lo'], r['delta_ci_hi'])} | "
            f"{'yes' if r['ci_excludes_zero']=='True' else 'no'} | "
            f"{'yes' if str(r['weaker_arm_switched']) in {'1','True'} else 'no'} |"
        )
    return "\n".join(rows)


def si_s6_hold(d):
    rows = [
        "| Pair | Main summary_min [95% CI] | holdout n (D / A / B) | holdout summary_min [95% CI] |",
        "|------|----------------------------:|----------------------:|------------------------------|",
    ]
    for p in PAIRS:
        if p not in d["hold"]:
            continue
        s = d["smin"][p]
        h = d["hold"][p]
        rows.append(
            f"| {p} | {r3(s['summary_min'])} {ci(s['ci_lo'], s['ci_hi'])} | "
            f"{h['n_dual']} / {h['n_A_only']} / {h['n_B_only']} | "
            f"{r3(h['summary_min'])} {ci(h['ci_lo'], h['ci_hi'])} |"
        )
    return "\n".join(rows)


def si_s10(d):
    rows = [
        "| Pair | n_ranked (D / A / B / N) | k | Top 10% D / A / B / N | dual / k | n_dual / n | EF_dual,10% | AND input → pass (D / A / B) | AND dual precision |",
        "|------|-------------------------:|--:|----------------------:|---------:|-----------:|------------:|------------------------------:|-------------------:|",
    ]
    for p in PAIRS:
        rk = d["rank"][p]
        a = d["andf"][p]
        n = int(rk["n_ranked"])
        nd = int(rk["n_dual"])
        rows.append(
            f"| {p} | {n} ({rk['n_dual']} / {rk['n_A_only']} / {rk['n_B_only']} / {rk['n_neither']}) | {rk['top_k']} | "
            f"{rk['top_dual']} / {rk['top_A_only']} / {rk['top_B_only']} / {rk['top_neither']} | "
            f"{r3(int(rk['top_dual'])/int(rk['top_k']))} | {r3(nd/n)} | {r3(rk['ef_dual_10pct'])} | "
            f"{a['n_filter_input']} → {int(a['retained_dual'])+int(a['retained_A_only'])+int(a['retained_B_only'])} "
            f"({a['retained_dual']} / {a['retained_A_only']} / {a['retained_B_only']}) | {r3(a['dual_precision'])} |"
        )
    return "\n".join(rows)


def patch_si(text, d):
    text = text.replace(
        "class-stratified bootstrap alternatives, sample-size simulations, raw BindingDB/PubChem supply counts",
        "non-stratified bootstrap sensitivity, the eight-pair detectable-effect simulation summarized after Table S4 (`results/canonical/detectable_effect_simulation.csv`), raw BindingDB/PubChem supply counts",
    )
    text = patch_table_block(text, "| Pair | Label rule | n (D / A / B) | summary_min | 95% CI |", si_s3(d))
    text = patch_table_block(
        text,
        "| Pair | Score channel | dual vs selective | dual vs neither | Δ | 95% CI | neither underpowered |",
        si_s4_fixed(d),
    )
    text = text.replace(
        "EGFR/HER2 cluster bootstrap was not recomputed after the canonical heavy-atom box; the official ligand-level difference is 0.462 [0.262, 0.651]. JAK1/TYK2 cluster rows are unchanged because that pair was not remade.",
        "EGFR/HER2 cluster bootstrap was recomputed with the corrected-box scores and the previously frozen scaffold and document groupings. Cluster intervals remain a sensitivity analysis and do not replace the ligand-level primary interval.",
    )
    text = patch_table_block(text, "| Pair | Resampling unit | Δ point | 95% CI | CI excludes 0 |", si_s4_cluster(d))
    det_block = si_detectable(d, zh=False)
    if det_block:
        cluster_src = (
            "Source: `results/canonical/fixed_score_negative_class_delta.csv`; "
            "`results/canonical/cluster_bootstrap_sensitivity.csv`. Dual-versus-neither with two-pocket mean scores is main-text Table 3."
        )
        text = re.sub(
            r"\*\*Detectable-effect simulation \(binormal; not observed power\)\.\*\*.*?results/canonical/detectable_effect_simulation\.csv`\.",
            det_block,
            text,
            count=1,
            flags=re.S,
        )
        if "Detectable-effect simulation (binormal; not observed power)." not in text:
            text = text.replace(cluster_src, cluster_src + "\n\n" + det_block, 1)
    _, inc_abs = max_inc(d)
    text = text.replace("Across 16 arms the largest |Δ| is 0.023.", f"Across 16 arms the largest |Δ| is {r3(inc_abs)}.")
    text = patch_table_block(
        text,
        "| Pair | Arm | ECFP4 | ECFP4+docking | Δ | Vina ranking AUROC (Table 2) |",
        si_s5_ecfp(d),
    )
    ache = d["desc"]["AChE/BChE"]
    excl, incl = desc_excl(d)
    text = text.replace(
        f"AChE/BChE TPSA directional AUROCs are 0.733 / 0.801. Six of eight 95% CIs include 0; F2/F10 and JAK1/TYK2 exclude 0.",
        f"AChE/BChE TPSA directional AUROCs are {r3(ache['tpsa_D_vs_A'])} / {r3(ache['tpsa_D_vs_B'])} on the full panel. {len(incl)} of eight 95% CIs for Vina minus full-panel-best include 0; {', '.join(excl)} exclude 0.",
    )
    text = text.replace(
        "The chemistry-control predictive number is nested scaffold-GroupKFold train-only descriptor selection OOF AUROC in the same CSV (`nested_scaffold_cv_oof_summary_min`; fold table `descriptor_nested_scaffold_cv.csv`).",
        "The chemistry-control predictive number is nested scaffold-GroupKFold inner-CV descriptor selection followed by train-only StandardScaler + univariate logistic regression; pooled AUROC uses held-out P(dual) (`nested_scaffold_cv_oof_summary_min`; OOF table `descriptor_nested_oof_predictions.csv`). The `summary_min` interval is a ligand-level shared-dual bootstrap of those fixed OOF probabilities.",
    )
    text = patch_table_block(
        text,
        "| Pair | Best descriptor | Descriptor summary_min | Δ | 95% CI | CI excludes 0 |",
        si_s5_desc(d),
    )
    nested_en = (
        "**Nested scaffold-GroupKFold single-descriptor baseline.** "
        "Each pair × arm uses independent outer scaffold GroupKFold. "
        "The descriptor is selected by inner scaffold CV on outer-training data only, "
        "then a train-only StandardScaler + univariate logistic regression emits held-out P(dual). "
        "Pooled AUROC uses those probabilities, not raw TPSA/cLogP/MW/heavy values. "
        "`summary_min` CI is a ligand-level shared-dual bootstrap of the min of the two directional OOF AUROCs "
        "on fixed OOF probabilities; it is not a re-run of model selection."
    )
    if "| Pair | D vs A OOF AUROC [95% CI] |" in text:
        text = patch_table_block(text, "| Pair | D vs A OOF AUROC [95% CI] |", si_s5_nested(d, zh=False))
    else:
        text = text.replace(
            "Source: `results/canonical/descriptor_baselines.csv`; `results/canonical/descriptor_nested_scaffold_cv.csv`; `results/canonical/ecfp4_incremental_information.csv`.",
            nested_en + "\n\n" + si_s5_nested(d, zh=False) + "\n\n"
            "Source: `results/canonical/descriptor_baselines.csv`; `results/canonical/descriptor_nested_scaffold_cv.csv`; "
            "`results/canonical/descriptor_nested_oof_predictions.csv`; `results/canonical/ecfp4_incremental_information.csv`.",
        )
    text = text.replace("does not replace the unscaled 0.023 primary result.", f"does not replace the unscaled {r3(inc_abs)} primary result.")
    text = patch_table_block(
        text,
        "| Pair | Set | Δ | 95% CI | CI excludes 0 | Weaker arm switched |",
        si_s6(d),
    )
    text = patch_table_block(
        text,
        "| Pair | Main summary_min [95% CI] | holdout n (D / A / B) | holdout summary_min [95% CI] |",
        si_s6_hold(d),
    )
    text = patch_table_block(
        text,
        "| Pair | n_ranked (D / A / B / N) | k | Top 10% D / A / B / N | dual / k | n_dual / n | EF_dual,10% | AND input → pass (D / A / B) | AND dual precision |",
        si_s10(d),
    )
    # S7 GNINA / receptor
    s = d["smin"]
    two = d["two"]
    g = d["gnina"]
    db = d["direc"]
    s7a = [
        "| Pair | Engine | n_dual / n_A / n_B / n_neither | summary_min | Weaker-arm AUROC [95% CI] | Dual vs neither |",
        "|------|------|------|------------:|---------------------------|----------------:|",
    ]
    for p in ("EGFR/HER2", "PIK3CA/mTOR", "JAK1/TYK2"):
        sp = s[p]
        tp = two[p]
        gp = g[p]
        weak = "dual–B-only (pocket A)" if float(sp["auroc_D_vs_B_pocketA"]) <= float(sp["auroc_D_vs_A_pocketB"]) else "dual–A-only (pocket B)"
        gweak, garm, glo, ghi = _gnina_weaker(gp, zh=False)
        s7a.append(
            f"| {p} | Vina primary | {sp['n_dual']} / {sp['n_A_only']} / {sp['n_B_only']} / {d['rank'][p]['n_neither']} | "
            f"{r3(sp['summary_min'])} {ci(sp['ci_lo'], sp['ci_hi'])} | {weak} {r3(sp['auroc_D_vs_B_pocketA'] if 'B-only' in weak else sp['auroc_D_vs_A_pocketB'])} "
            f"{ci(db[(p,'AUROC_D_vs_B_pocketA' if 'B-only' in weak else 'AUROC_D_vs_A_pocketB')]['ci_lo'], db[(p,'AUROC_D_vs_B_pocketA' if 'B-only' in weak else 'AUROC_D_vs_A_pocketB')]['ci_hi'])} | "
            f"{r3(tp['two_pocket_mean_D_vs_neither'])} {ci(tp['ci_lo'], tp['ci_hi'])} |"
        )
        s7a.append(
            f"| {p} | GNINA independent | {_gnina_counts(gp)} | {r3(gp['summary_min'])} {ci(gp['summary_min_ci_lo'], gp['summary_min_ci_hi'])} | "
            f"{gweak} {r3(garm)} {ci(glo, ghi)} | "
            f"{r3(gp['auroc_D_vs_neither_mean'])} {ci(gp['d_vs_neither_ci_lo'], gp['d_vs_neither_ci_hi'])} |"
        )
    text = patch_table_block(text, "| Pair | Engine | n_dual / n_A / n_B / n_neither | summary_min | Weaker-arm AUROC [95% CI] | Dual vs neither |", "\n".join(s7a))
    r4j, r5d, r4s = d["rec"]["4JPS"], d["rec"]["5DXT"], d["rec"]["4JSX"]
    sp = s["PIK3CA/mTOR"]
    s7c = [
        "| Replacement | Pocket replaced | D vs A | D vs B | summary_min [95% CI] |",
        "|------|------------|-------:|-------:|----------------------|",
        f"| Main 4L23 / 4JT6 | — | {r3(sp['auroc_D_vs_A_pocketB'])} | {r3(sp['auroc_D_vs_B_pocketA'])} | {r3(sp['summary_min'])} {ci(sp['ci_lo'], sp['ci_hi'])} |",
        f"| PIK3CA → 4JPS | A | {r3(r4j['auroc_D_vs_A'])} | {r3(r4j['auroc_D_vs_B'])} | {r3(r4j['summary_min'])} {ci(r4j['summary_min_ci_lo'], r4j['summary_min_ci_hi'])} |",
        f"| PIK3CA → 5DXT | A | {r3(r5d['auroc_D_vs_A'])} | {r3(r5d['auroc_D_vs_B'])} | {r3(r5d['summary_min'])} {ci(r5d['summary_min_ci_lo'], r5d['summary_min_ci_hi'])} |",
        f"| mTOR → 4JSX | B | {r3(r4s['auroc_D_vs_A'])} | {r3(r4s['auroc_D_vs_B'])} | {r3(r4s['summary_min'])} {ci(r4s['summary_min_ci_lo'], r4s['summary_min_ci_hi'])} |",
    ]
    text = patch_table_block(text, "| Replacement | Pocket replaced | D vs A | D vs B | summary_min [95% CI] |", "\n".join(s7c))
    text = text.replace(
        "EGFR/HER2 and PIK3CA/mTOR independent-GNINA `summary_min` rows have empty CI columns in the source file; intervals below are labeled on the corresponding single arm and are not min-of-two bootstrap intervals.",
        "Independent-GNINA `summary_min` intervals below use the same class-stratified shared-dual protocol as Table 2.",
    )
    text = text.replace(
        "Source: `pocket_matched_vs_best_descriptor_delta_v1.csv`; `descriptor_paired_delta_s19_v1.csv`; `incremental_information_v1.csv`; `ecfp4_incremental_s20s24_v1.csv`; `ligand_ml_baseline_scaffold_cv_v1.csv`.",
        "Source: `results/canonical/descriptor_baselines.csv`; `results/canonical/ecfp4_incremental_information.csv`.",
    )
    text = text.replace(
        "Source: `wrong_pocket_paired_delta_bootstrap_v1.csv` (`set=main_panel` / `unused_pool_holdout`); `wrong_pocket_by_channel_v1.csv`; `pocket_unidirectional_delta_v1.csv`; `holdout_pocket_matched_v1.csv`; `table2_comparable_by_channel_v1.csv` (`holdout_vina_20260727`).",
        "Source: `results/canonical/matched_mismatched_pocket.csv`; `results/canonical/holdout_metrics.csv`.",
    )
    text = text.replace(
        "Source: `eight_pair_ranking_operating_point_v1.csv`; `review_scored_membership_v1.csv`.",
        "Source: `results/canonical/top10_operating_points.csv`; `results/canonical/current_score_master.csv`.",
    )
    pg = d["smin"]["PPARG/PPARA"]
    tpg = d["two"]["PPARG/PPARA"]
    hpg = d["hold"]["PPARG/PPARA"]
    text = re.sub(
        r"\| Vina primary \| 0\.649 \[0\.\d+, 0\.\d+\] \| 0\.685 \[0\.\d+, 0\.\d+\] \|",
        f"| Vina primary | {r3(pg['summary_min'])} {ci(pg['ci_lo'], pg['ci_hi'])} | {r3(tpg['two_pocket_mean_D_vs_neither'])} {ci(tpg['ci_lo'], tpg['ci_hi'])} |",
        text,
    )
    text = re.sub(
        r"\| unused-pool holdout \| 0\.535 \[0\.\d+, 0\.\d+\] \| — \|",
        f"| unused-pool holdout | {r3(hpg['summary_min'])} {ci(hpg['ci_lo'], hpg['ci_hi'])} | — |",
        text,
    )
    mx = _maxmed_lookup(d)
    egfr_md = mx.get("EGFR/HER2", {}).get("median", {})
    ache_md = mx.get("AChE/BChE", {}).get("median", {})
    ache_mx = mx.get("AChE/BChE", {}).get("max", {})
    ppar_md = mx.get("PPARA/PPARD", {}).get("median", {})
    text = text.replace(
        "Max-to-median class flips: EGFR/HER2 6/110 (primary `summary_min` 0.324); AChE/BChE 1/96 (CHEMBL659; 0.606 → 0.629); PPARA/PPARD 1/110 (CHEMBL121; A-only 32→31; dual-versus-A-only 0.646 → 0.636; `summary_min` remained 0.446).",
        f"Max-to-median class flips on the qualified-record intersection: EGFR/HER2 {egfr_md.get('class_flips_vs_max','')}/{egfr_md.get('n_ligands','')} "
        f"(primary `summary_min` {r3(d['smin']['EGFR/HER2']['summary_min'])}); "
        f"AChE/BChE {ache_md.get('class_flips_vs_max','')}/{ache_md.get('n_ligands','')} "
        f"(CHEMBL659; {r3(ache_mx.get('summary_min', 0.6058))} → {r3(ache_md.get('summary_min', 0.6291))}); "
        f"PPARA/PPARD {ppar_md.get('class_flips_vs_max','')}/{ppar_md.get('n_ligands','')} "
        f"(CHEMBL121; `summary_min` remained {r3(mx.get('PPARA/PPARD', {}).get('max', {}).get('summary_min', 0.4463))}).",
    )
    text = text.replace(
        "Source: `eight_pair_dump_gated_v1/max_vs_median_auroc_v1.csv`; `eight_pair_dump_gated_v1/parity_v1.csv`.",
        "Source: `results/canonical/max_vs_median_sensitivity.csv`.",
    )
    return apply_boundary_patches(text, d, zh=False)


def yesno(flag, zh=False) -> str:
    truth = str(flag) in {"1", "True", "true", "yes"}
    if zh:
        return "是" if truth else "否"
    return "yes" if truth else "no"


def si_s3_zh(d):
    rows = ["| 靶对 | 标签规则 | n (D / A / B) | summary_min | 95% CI |", "|------|----------|--------------:|------------:|--------|"]
    show = {
        "EGFR/HER2": ["theta_5.5", "theta_6.0", "theta_6.5", "strict_6.5_5.5"],
        "PIK3CA/mTOR": ["theta_5.5", "theta_6.0", "theta_6.5", "strict_6.5_5.5"],
    }
    labels = {(r["pair"], r["label_rule"]): r for r in d["labels"]}
    name = {"theta_5.5": "θ = 5.5", "theta_6.0": "θ = 6.0", "theta_6.5": "θ = 6.5", "strict_6.5_5.5": "严格 6.5/5.5"}
    for p in PAIRS:
        for rule in show.get(p, ["theta_6.0"]):
            r = labels[(p, rule)]
            rows.append(
                f"| {p} | {name[rule]} | {r['n_dual']} / {r['n_A_only']} / {r['n_B_only']} | {r3(r['summary_min'])} | {ci(r['ci_lo'], r['ci_hi'])} |"
            )
    return "\n".join(rows)


def si_s4_fixed_zh(d):
    rows = [
        "| 靶对 | 评分通道 | dual vs 选择性 | dual vs neither | Δ | 95% CI | neither 效能不足 |",
        "|------|----------|---------------:|----------------:|--:|--------|:----------------:|",
    ]
    order = [
        ("EGFR/HER2", "D_vs_B_or_neither_pocketA", "口袋 A（对 B-only）"),
        ("EGFR/HER2", "D_vs_A_or_neither_pocketB", "口袋 B（对 A-only）"),
        ("AChE/BChE", "D_vs_B_or_neither_pocketA", "口袋 A"),
        ("AChE/BChE", "D_vs_A_or_neither_pocketB", "口袋 B"),
        ("PIK3CA/mTOR", "D_vs_B_or_neither_pocketA", "口袋 A"),
        ("PIK3CA/mTOR", "D_vs_A_or_neither_pocketB", "口袋 B"),
        ("F2/F10", "D_vs_B_or_neither_pocketA", "口袋 A"),
        ("F2/F10", "D_vs_A_or_neither_pocketB", "口袋 B"),
        ("JAK1/TYK2", "D_vs_B_or_neither_pocketA", "口袋 A"),
        ("JAK1/TYK2", "D_vs_A_or_neither_pocketB", "口袋 B"),
        ("JAK1/JAK2", "D_vs_B_or_neither_pocketA", "口袋 A"),
        ("JAK1/JAK2", "D_vs_A_or_neither_pocketB", "口袋 B"),
        ("PPARG/PPARA", "D_vs_B_or_neither_pocketA", "口袋 A"),
        ("PPARG/PPARA", "D_vs_A_or_neither_pocketB", "口袋 B"),
        ("PPARA/PPARD", "D_vs_B_or_neither_pocketA", "口袋 A"),
        ("PPARA/PPARD", "D_vs_A_or_neither_pocketB", "口袋 B"),
    ]
    for p, c, lab in order:
        r = d["fixed"][(p, c)]
        und = "是" if int(r["n_neither"]) < 10 else "否"
        rows.append(
            f"| {p} | {lab} | {r3(r['auroc_dual_vs_selective'])} | {r3(r['auroc_dual_vs_neither'])} | "
            f"{r3m(r['delta_neither_minus_selective'])} | {ci(r['delta_ci_lo'], r['delta_ci_hi'])} | {und} |"
        )
    return "\n".join(rows)


def si_s4_cluster_zh(d):
    rows = [
        "| 靶对 | 重采样单位 | Δ 点估计 | 95% CI | CI 排除 0 |",
        "|------|-----------------|--------:|--------|:-------------:|",
    ]
    labels = {
        "ligand_stratified": "配体层（正式）",
        "scaffold_cluster": "骨架簇",
        "document_cluster": "文献簇",
    }
    for r in d["cluster"]:
        if r.get("status") == "unresolved_mapping_unavailable" or r.get("delta_ci_lo") in ("", None):
            rows.append(
                f"| {r['pair']} | {labels.get(r['estimator'], r['estimator'])} | {r3(r['delta_point'])} | "
                f"未重算 | — |"
            )
            continue
        rows.append(
            f"| {r['pair']} | {labels.get(r['estimator'], r['estimator'])} | {r3(r['delta_point'])} | "
            f"{ci(r['delta_ci_lo'], r['delta_ci_hi'])} | {yesno(r['excludes_zero'], zh=True)} |"
        )
    return "\n".join(rows)


def si_s5_ecfp_zh(d):
    rows = [
        "| 靶对 | 方向 | ECFP4 | ECFP4+对接 | Δ | Vina 排序 AUROC（Table 2） |",
        "|------|------|------:|--------------:|--:|--------------------------:|",
    ]
    for p in PAIRS:
        for arm, lab in (("D_vs_A", "D vs A"), ("D_vs_B", "D vs B")):
            r = d["inc"][(p, arm)]
            rows.append(
                f"| {p} | {lab} | {r3(r['cv_auroc_ECFP4'])} | {r3(r['cv_auroc_ECFP4_docking'])} | "
                f"{signed(r['delta_ECFP4_plus_docking_minus_ECFP4'])} | {r3(r['rank_auroc_docking'])} |"
            )
    return "\n".join(rows)


def si_s5_desc_zh(d):
    rows = [
        "| 靶对 | 最佳描述符 | 描述符 summary_min | Δ | 95% CI | CI 排除 0 |",
        "|------|-----------------|-----------------------:|--:|--------|:-------------:|",
    ]
    names = {"clogp": "cLogP", "tpsa": "TPSA", "heavy": "重原子数", "mw": "MW"}
    for p in PAIRS:
        r = d["desc"][p]
        excl = "是" if (float(r["delta_ci_lo"]) > 0 or float(r["delta_ci_hi"]) < 0) else "否"
        rows.append(
            f"| {p} | {names.get(r['best_descriptor'], r['best_descriptor'])} | {r3(r['best_descriptor_summary_min'])} | "
            f"{r3m(r['vina_minus_best_descriptor'])} | {ci(r['delta_ci_lo'], r['delta_ci_hi'])} | {excl} |"
        )
    return "\n".join(rows)


def si_s6_zh(d):
    rows = [
        "| 靶对 | 集合 | Δ | 95% CI | CI 排除 0 | 较弱方向切换 |",
        "|------|------|--:|--------|:---------:|:-------------------:|",
    ]
    for p in PAIRS:
        r = d["mm"][p]
        rows.append(
            f"| {p} | 主集 | {r3m(r['delta'])} | {ci(r['delta_ci_lo'], r['delta_ci_hi'])} | "
            f"{yesno(r['ci_excludes_zero'], zh=True)} | {yesno(r['weaker_arm_switched'], zh=True)} |"
        )
    for p in PAIRS:
        if p not in d["hold_wp"]:
            continue
        r = d["hold_wp"][p]
        rows.append(
            f"| {p} | 留出 | {r3m(r['delta_matched_minus_wrong'])} | {ci(r['delta_ci_lo'], r['delta_ci_hi'])} | "
            f"{yesno(r['ci_excludes_zero'], zh=True)} | {yesno(r['weaker_arm_switched'], zh=True)} |"
        )
    return "\n".join(rows)


def si_s6_hold_zh(d):
    rows = [
        "| 靶对 | 主集 summary_min [95% CI] | 留出 n (D / A / B) | 留出 summary_min [95% CI] |",
        "|------|----------------------------:|----------------------:|------------------------------|",
    ]
    for p in PAIRS:
        if p not in d["hold"]:
            continue
        s = d["smin"][p]
        h = d["hold"][p]
        rows.append(
            f"| {p} | {r3(s['summary_min'])} {ci(s['ci_lo'], s['ci_hi'])} | "
            f"{h['n_dual']} / {h['n_A_only']} / {h['n_B_only']} | "
            f"{r3(h['summary_min'])} {ci(h['ci_lo'], h['ci_hi'])} |"
        )
    return "\n".join(rows)


def patch_si_zh(text, d):
    text = text.replace(
        "Table 2 报告锁定的配体层非分层百分位 bootstrap 下两条方向性区间及逐次 `summary_min` 区间。dual–neither 为两类分层百分位区间；对应/非对应口袋差值为配对 bootstrap。类别分层 `summary_min` 敏感性见归档表 `summary_min_stratified_sensitivity_review_v1.csv`，不替换 Table 2",
        "Table 2 报告类别分层非参数百分位 bootstrap（B = 2000，种子 20260729）下两条方向性区间及逐次 `summary_min` 区间。每次重采样内 dual 配体只抽一次并同时用于两个方向。dual–neither 使用同一类别分层协议；对应/非对应口袋差值为配体水平配对 bootstrap。非分层配体 bootstrap 仅作敏感性（`non_stratified_bootstrap_sensitivity.csv`），不替换 Table 2",
    )
    text = patch_table_block(text, "| 靶对 | 标签规则 | n (D / A / B) | summary_min | 95% CI |", si_s3_zh(d))
    text = patch_table_block(
        text,
        "| 靶对 | 评分通道 | dual vs 选择性 | dual vs neither | Δ | 95% CI | neither 效能不足 |",
        si_s4_fixed_zh(d),
    )
    text = text.replace(
        "EGFR/HER2 的簇 bootstrap 未在规范重原子盒上重算；正式配体层差值为 0.462 [0.262, 0.651]。JAK1/TYK2 簇行未改，因为该对未重做。",
        "EGFR/HER2 的簇 bootstrap 已用校正盒评分和原冻结的骨架/文献分组重算。簇区间仍为敏感性分析，不替换配体层主区间。",
    )
    text = patch_table_block(text, "| 靶对 | 重采样单位 | Δ 点估计 | 95% CI | CI 排除 0 |", si_s4_cluster_zh(d))
    det_zh = si_detectable(d, zh=True)
    if det_zh:
        cluster_src_zh = (
            "源：`results/canonical/fixed_score_negative_class_delta.csv`；"
            "`results/canonical/cluster_bootstrap_sensitivity.csv`。双口袋平均分的 dual–neither 比较见正文 Table 3。"
        )
        text = re.sub(
            r"\*\*Detectable-effect 仿真（二项正态；不是观测功效）。\*\*.*?results/canonical/detectable_effect_simulation\.csv`。",
            det_zh,
            text,
            count=1,
            flags=re.S,
        )
        if "Detectable-effect 仿真（二项正态；不是观测功效）。" not in text:
            text = text.replace(cluster_src_zh, cluster_src_zh + "\n\n" + det_zh, 1)
    text = text.replace(
        "源：`formulation_equal_score_negative_v1.csv`；`equal_score_negative_s34_v1.csv`；`equal_score_cluster_bootstrap_v1.csv`（仅 JAK1/TYK2）。",
        "源：`results/canonical/fixed_score_negative_class_delta.csv`；`results/canonical/cluster_bootstrap_sensitivity.csv`。",
    )
    text = text.replace(
        "源：`pocket_matched_vs_best_descriptor_delta_v1.csv`；`descriptor_paired_delta_s19_v1.csv`；`incremental_information_v1.csv`；`ecfp4_incremental_s20s24_v1.csv`；`ligand_ml_baseline_scaffold_cv_v1.csv`。",
        "源：`results/canonical/descriptor_baselines.csv`；`results/canonical/ecfp4_incremental_information.csv`。",
    )
    text = text.replace(
        "源：`wrong_pocket_paired_delta_bootstrap_v1.csv`（`set=main_panel` / `unused_pool_holdout`）；`wrong_pocket_by_channel_v1.csv`；`pocket_unidirectional_delta_v1.csv`；`holdout_pocket_matched_v1.csv`；`table2_comparable_by_channel_v1.csv`（`holdout_vina_20260727`）。",
        "源：`results/canonical/matched_mismatched_pocket.csv`；`results/canonical/holdout_metrics.csv`。",
    )
    text = text.replace(
        "源：`eight_pair_ranking_operating_point_v1.csv`；`review_scored_membership_v1.csv`。",
        "源：`results/canonical/top10_operating_points.csv`；`results/canonical/current_score_master.csv`。",
    )
    _, inc_abs = max_inc(d)
    text = text.replace("16 个方向上最大 |Δ| 为 0.023。", f"16 个方向上最大 |Δ| 为 {r3(inc_abs)}。")
    text = patch_table_block(
        text,
        "| 靶对 | 方向 | ECFP4 | ECFP4+对接 | Δ | Vina 排序 AUROC（Table 2） |",
        si_s5_ecfp_zh(d),
    )
    ache = d["desc"]["AChE/BChE"]
    excl, incl = desc_excl(d)
    text = text.replace(
        "AChE/BChE 的 TPSA 方向性 AUROC 为 0.733 / 0.801。八对中六对 95% CI 包含 0；F2/F10 与 JAK1/TYK2 不包含 0。",
        f"AChE/BChE 的 TPSA 方向性 AUROC 为 {r3(ache['tpsa_D_vs_A'])} / {r3(ache['tpsa_D_vs_B'])}。八对中 {len(incl)} 对 95% CI 包含 0；{'、'.join(excl)} 不包含 0。",
    )
    text = patch_table_block(
        text,
        "| 靶对 | 最佳描述符 | 描述符 summary_min | Δ | 95% CI | CI 排除 0 |",
        si_s5_desc_zh(d),
    )
    nested_zh = (
        "**嵌套骨架 GroupKFold 单描述符基线。** "
        "每个靶对 × 方向独立做外层骨架 GroupKFold。"
        "描述符仅在外层训练折上通过内层骨架 CV 选择，再在外层训练集上拟合 StandardScaler + 单变量逻辑回归，输出折外 P(dual)。"
        "合并 AUROC 使用该概率，不再拼接原始 TPSA/cLogP/MW/重原子数。"
        "`summary_min` CI 是在固定折外概率上的配体层共享 dual bootstrap，不是重新执行模型选择的不确定度。"
    )
    if "| 靶对 | D vs A OOF AUROC [95% CI] |" in text:
        text = patch_table_block(text, "| 靶对 | D vs A OOF AUROC [95% CI] |", si_s5_nested(d, zh=True))
    else:
        text = text.replace(
            "源：`results/canonical/descriptor_baselines.csv`；`results/canonical/descriptor_nested_scaffold_cv.csv`；`results/canonical/ecfp4_incremental_information.csv`。",
            nested_zh + "\n\n" + si_s5_nested(d, zh=True) + "\n\n"
            "源：`results/canonical/descriptor_baselines.csv`；`results/canonical/descriptor_nested_scaffold_cv.csv`；"
            "`results/canonical/descriptor_nested_oof_predictions.csv`；`results/canonical/ecfp4_incremental_information.csv`。",
        )
    text = text.replace("缩放不替代未缩放的 0.023 主结果。", f"缩放不替代未缩放的 {r3(inc_abs)} 主结果。")
    text = patch_table_block(text, "| 靶对 | 集合 | Δ | 95% CI | CI 排除 0 | 较弱方向切换 |", si_s6_zh(d))
    text = patch_table_block(
        text,
        "| 靶对 | 主集 summary_min [95% CI] | 留出 n (D / A / B) | 留出 summary_min [95% CI] |",
        si_s6_hold_zh(d),
    )
    s = d["smin"]
    two = d["two"]
    g = d["gnina"]
    db = d["direc"]
    s7a = [
        "| 靶对 | 引擎 | n_dual / n_A / n_B / n_neither | summary_min | 较弱臂 AUROC [95% CI] | Dual vs neither |",
        "|------|------|------|------------:|---------------------------|----------------:|",
    ]
    for p in ("EGFR/HER2", "PIK3CA/mTOR", "JAK1/TYK2"):
        sp = s[p]
        tp = two[p]
        gp = g[p]
        weak = "dual–B-only（口袋 A）" if float(sp["auroc_D_vs_B_pocketA"]) <= float(sp["auroc_D_vs_A_pocketB"]) else "dual–A-only（口袋 B）"
        gweak, garm, glo, ghi = _gnina_weaker(gp, zh=True)
        s7a.append(
            f"| {p} | Vina 主分析 | {sp['n_dual']} / {sp['n_A_only']} / {sp['n_B_only']} / {d['rank'][p]['n_neither']} | "
            f"{r3(sp['summary_min'])} {ci(sp['ci_lo'], sp['ci_hi'])} | {weak} {r3(sp['auroc_D_vs_B_pocketA'] if 'B-only' in weak else sp['auroc_D_vs_A_pocketB'])} "
            f"{ci(db[(p,'AUROC_D_vs_B_pocketA' if 'B-only' in weak else 'AUROC_D_vs_A_pocketB')]['ci_lo'], db[(p,'AUROC_D_vs_B_pocketA' if 'B-only' in weak else 'AUROC_D_vs_A_pocketB')]['ci_hi'])} | "
            f"{r3(tp['two_pocket_mean_D_vs_neither'])} {ci(tp['ci_lo'], tp['ci_hi'])} |"
        )
        s7a.append(
            f"| {p} | GNINA 独立 | {_gnina_counts(gp)} | {r3(gp['summary_min'])} {ci(gp['summary_min_ci_lo'], gp['summary_min_ci_hi'])} | "
            f"{gweak} {r3(garm)} {ci(glo, ghi)} | "
            f"{r3(gp['auroc_D_vs_neither_mean'])} {ci(gp['d_vs_neither_ci_lo'], gp['d_vs_neither_ci_hi'])} |"
        )
    text = patch_table_block(
        text,
        "| 靶对 | 引擎 | n_dual / n_A / n_B / n_neither | summary_min | 较弱臂 AUROC [95% CI] | Dual vs neither |",
        "\n".join(s7a),
    )
    r4j, r5d, r4s = d["rec"]["4JPS"], d["rec"]["5DXT"], d["rec"]["4JSX"]
    sp = s["PIK3CA/mTOR"]
    s7c = [
        "| 替换 | 被替换口袋 | D vs A | D vs B | summary_min [95% CI] |",
        "|------|------------|-------:|-------:|----------------------|",
        f"| 主集 4L23 / 4JT6 | — | {r3(sp['auroc_D_vs_A_pocketB'])} | {r3(sp['auroc_D_vs_B_pocketA'])} | {r3(sp['summary_min'])} {ci(sp['ci_lo'], sp['ci_hi'])} |",
        f"| PIK3CA → 4JPS | A | {r3(r4j['auroc_D_vs_A'])} | {r3(r4j['auroc_D_vs_B'])} | {r3(r4j['summary_min'])} {ci(r4j['summary_min_ci_lo'], r4j['summary_min_ci_hi'])} |",
        f"| PIK3CA → 5DXT | A | {r3(r5d['auroc_D_vs_A'])} | {r3(r5d['auroc_D_vs_B'])} | {r3(r5d['summary_min'])} {ci(r5d['summary_min_ci_lo'], r5d['summary_min_ci_hi'])} |",
        f"| mTOR → 4JSX | B | {r3(r4s['auroc_D_vs_A'])} | {r3(r4s['auroc_D_vs_B'])} | {r3(r4s['summary_min'])} {ci(r4s['summary_min_ci_lo'], r4s['summary_min_ci_hi'])} |",
    ]
    text = patch_table_block(text, "| 替换 | 被替换口袋 | D vs A | D vs B | summary_min [95% CI] |", "\n".join(s7c))
    text = text.replace(
        "EGFR/HER2 与 PIK3CA/mTOR 独立 GNINA 的 `summary_min` 行在源文件中 CI 列为空；下表区间标在对应单臂上，不是对两臂逐次取最小值的区间。",
        "独立 GNINA 的 `summary_min` 区间使用与 Table 2 相同的类别分层、共享 dual 协议。",
    )
    pg = d["smin"]["PPARG/PPARA"]
    tpg = d["two"]["PPARG/PPARA"]
    hpg = d["hold"]["PPARG/PPARA"]
    text = re.sub(
        r"\| Vina 主分析 \| 0\.649 \[0\.\d+, 0\.\d+\] \| 0\.685 \[0\.\d+, 0\.\d+\] \|",
        f"| Vina 主分析 | {r3(pg['summary_min'])} {ci(pg['ci_lo'], pg['ci_hi'])} | {r3(tpg['two_pocket_mean_D_vs_neither'])} {ci(tpg['ci_lo'], tpg['ci_hi'])} |",
        text,
    )
    text = re.sub(
        r"\| 未使用池留出集 \| 0\.535 \[0\.\d+, 0\.\d+\] \| — \|",
        f"| 未使用池留出集 | {r3(hpg['summary_min'])} {ci(hpg['ci_lo'], hpg['ci_hi'])} | — |",
        text,
    )
    text = patch_table_block(
        text,
        "| 靶对 | n_ranked (D / A / B / N) | k | Top 10% D / A / B / N | dual / k | n_dual / n | EF_dual,10% | AND 输入 → 通过 (D / A / B) | AND dual precision |",
        si_s10(d).replace("| Pair |", "| 靶对 |").replace("AND input → pass", "AND 输入 → 通过"),
    )
    text = text.replace(
        "EGFR/HER2 簇 bootstrap 在规范盒校正后未重算；正式配体层差值为 0.462 [0.262, 0.651]。该图不重复 Figure 4 的留出集口袋对调。",
        "EGFR/HER2 簇 bootstrap 已用校正盒评分和冻结的骨架/文献分组重算。该图不重复 Figure 4 的留出集口袋对调。",
    )
    mx = _maxmed_lookup(d)
    egfr_md = mx.get("EGFR/HER2", {}).get("median", {})
    ache_md = mx.get("AChE/BChE", {}).get("median", {})
    ache_mx = mx.get("AChE/BChE", {}).get("max", {})
    ppar_md = mx.get("PPARA/PPARD", {}).get("median", {})
    text = text.replace(
        "最大值到中位数的类别翻转：EGFR/HER2 6/110（主分析 `summary_min` 0.324）；AChE/BChE 1/96（CHEMBL659；0.606 → 0.629）；PPARA/PPARD 1/110（CHEMBL121；A-only 32→31；dual–A-only 0.646 → 0.636；`summary_min` 仍为 0.446）。",
        f"合格记录交集上最大值到中位数的类别翻转：EGFR/HER2 {egfr_md.get('class_flips_vs_max','')}/{egfr_md.get('n_ligands','')}（主分析 `summary_min` {r3(d['smin']['EGFR/HER2']['summary_min'])}）；"
        f"AChE/BChE {ache_md.get('class_flips_vs_max','')}/{ache_md.get('n_ligands','')}（CHEMBL659；{r3(ache_mx.get('summary_min', 0.6058))} → {r3(ache_md.get('summary_min', 0.6291))}）；"
        f"PPARA/PPARD {ppar_md.get('class_flips_vs_max','')}/{ppar_md.get('n_ligands','')}（CHEMBL121；`summary_min` 仍为 {r3(mx.get('PPARA/PPARD', {}).get('max', {}).get('summary_min', 0.4463))}）。",
    )
    text = text.replace(
        "Table 2 所用的同一套 ChEMBL 37 记录在 θ = 6.0 下按中位数重新汇总。最大 pChEMBL 与全部已评分配体一致（缺失端 0；不一致 0）。",
        "最大与中位数使用同一套合格记录、同一分子交集和当前分数。缺少 dump 记录不等于从主分析删除该分子。",
    )
    return apply_boundary_patches(text, d, zh=True)


def patch_lock(text, d):
    egfr_fix = d["fixed"][("EGFR/HER2", "D_vs_B_or_neither_pocketA")]
    s = d["smin"]["EGFR/HER2"]
    two = d["two"]["EGFR/HER2"]
    mm = d["mm"]["EGFR/HER2"]
    mm_a = d["mm"]["AChE/BChE"]
    rk = d["rank"]["EGFR/HER2"]
    rk_a = d["rank"]["AChE/BChE"]
    g = d["gnina"]["EGFR/HER2"]
    text = re.sub(
        r"EGFR/HER2: D vs A-only pocket B = 0\.661; D vs B-only pocket A = 0\.324; summary_min = 0\.324\.\s*"
        r"Same pocket-A score: D vs B-only = 0\.324; D vs neither = 0\.786; fixed-score ΔAUROC = 0\.462 \[[^\]]+\]\.\s*"
        r"Two-pocket mean D-vs-neither = 0\.759\.\s*"
        r"Matched-minus-mismatched = 0\.056 \[[^\]]+\]; CI includes zero; no matched-pocket advantage\.\s*"
        r"AChE/BChE: n_scored = 27 / 26 / 28; matched-minus-mismatched = 0\.177; CI excludes 0 \(only main-panel interval excluding zero\)\.\s*"
        r"EGFR/HER2 Top 10% = 1 / 5 / 5 / 0; EFdual,10% = 0\.357\.\s*"
        r"AChE/BChE Top 10% = 5 / 3 / 1 / 1; EFdual,10% ≈ 1\.778\.\s*"
        r"GNINA EGFR: D-vs-B-only ≈ 0\.265; D-vs-neither ≈ 0\.737\.",
        f"EGFR/HER2: D vs A-only pocket B = {r3(s['auroc_D_vs_A_pocketB'])}; D vs B-only pocket A = {r3(s['auroc_D_vs_B_pocketA'])}; summary_min = {r3(s['summary_min'])} {ci(s['ci_lo'], s['ci_hi'])}.  \n"
        f"Same pocket-A score: D vs B-only = {r3(egfr_fix['auroc_dual_vs_selective'])}; D vs neither = {r3(egfr_fix['auroc_dual_vs_neither'])}; fixed-score ΔAUROC = {r3(egfr_fix['delta_neither_minus_selective'])} {ci(egfr_fix['delta_ci_lo'], egfr_fix['delta_ci_hi'])}.  \n"
        f"Two-pocket mean D-vs-neither = {r3(two['two_pocket_mean_D_vs_neither'])} {ci(two['ci_lo'], two['ci_hi'])}.  \n"
        f"Matched-minus-mismatched = {r3(mm['delta'])} {ci(mm['delta_ci_lo'], mm['delta_ci_hi'])}; CI includes zero; no matched-pocket advantage.  \n"
        f"AChE/BChE: n_scored = 27 / 26 / 28; matched-minus-mismatched = {r3(mm_a['delta'])} {ci(mm_a['delta_ci_lo'], mm_a['delta_ci_hi'])}; CI excludes 0 (only main-panel interval excluding zero).  \n"
        f"EGFR/HER2 Top 10% = {rk['top_dual']} / {rk['top_A_only']} / {rk['top_B_only']} / {rk['top_neither']}; EFdual,10% = {r3(rk['ef_dual_10pct'])}.  \n"
        f"AChE/BChE Top 10% = {rk_a['top_dual']} / {rk_a['top_A_only']} / {rk_a['top_B_only']} / {rk_a['top_neither']}; EFdual,10% = {r3(rk_a['ef_dual_10pct'])}.  \n"
        f"GNINA EGFR: D-vs-B-only = {r3(g['auroc_D_vs_B_pocketA'])}; D-vs-neither = {r3(g['auroc_D_vs_neither_mean'])}.",
        text,
    )
    text = text.replace(
        "EGFR cluster was not recomputed after the box correction",
        "EGFR cluster recomputed with corrected-box scores and frozen groupings",
    )
    text = text.replace("Post-fix EGFR pocket A = 0.462", f"EGFR pocket A Δ = {r3(egfr_fix['delta_neither_minus_selective'])}")
    return text


def main() -> int:
    d = load()
    en = DOCS / "MANUSCRIPT_JCIM_EN.md"
    zh = DOCS / "MANUSCRIPT_JCIM_ZH.md"
    si = DOCS / "SUPPORTING_INFORMATION_JCIM_EN_V1.md"
    sizh = DOCS / "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md"
    lock = DOCS / "FIGURE_TABLE_LOCK_POSTFIX_V4.md"
    en.write_text(patch_en(en.read_text(encoding="utf-8"), d), encoding="utf-8")
    zh.write_text(patch_zh(zh.read_text(encoding="utf-8"), d), encoding="utf-8")
    si.write_text(patch_si(si.read_text(encoding="utf-8"), d), encoding="utf-8")
    if sizh.is_file():
        t = sizh.read_text(encoding="utf-8")
        t = patch_zh(t, d)
        t = patch_si_zh(t, d)
        sizh.write_text(t, encoding="utf-8")
    lock.write_text(patch_lock(lock.read_text(encoding="utf-8"), d), encoding="utf-8")
    extra_zh = [
        DOCS / "METHODS_DRAFT_ZH_JCIM_V1.md",
        DOCS / "RESULTS_DRAFT_ZH_JCIM_V1.md",
        DOCS / "DISCUSSION_DRAFT_ZH_JCIM_V1.md",
        DOCS / "TITLE_AND_ABSTRACT_JCIM_ZH_V1.md",
    ]
    extra_en = [
        DOCS / "METHODS_SECTION_JCIM_EN_V1.md",
        DOCS / "RESULTS_SECTION_JCIM_EN_V1.md",
        DOCS / "DISCUSSION_SECTION_JCIM_EN_V1.md",
        DOCS / "TITLE_AND_ABSTRACT_JCIM_EN_V1.md",
    ]
    for extra in extra_en:
        if extra.is_file():
            extra.write_text(patch_en(extra.read_text(encoding="utf-8"), d), encoding="utf-8")
    for extra in extra_zh:
        if extra.is_file():
            extra.write_text(patch_zh(extra.read_text(encoding="utf-8"), d), encoding="utf-8")
    print("patched manuscripts and SI from canonical tables")
    print("EGFR summary_min", r3(d["smin"]["EGFR/HER2"]["summary_min"]), ci(d["smin"]["EGFR/HER2"]["ci_lo"], d["smin"]["EGFR/HER2"]["ci_hi"]))
    print("max |ECFP Δ|", r3(max_inc(d)[1]), max_inc(d)[0]["pair"], max_inc(d)[0]["contrast"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
