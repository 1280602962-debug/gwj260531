"""Independent display audit. Reads original CSVs and manuscript text.

Does not import the figure or integrate scripts. Does not rewrite concordance counts.
"""

from __future__ import annotations

import csv
import hashlib
import subprocess
import zipfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parent
MAN = ROOT / "manuscript"
ORDER = [
    "EGFR/HER2", "JAK1/JAK2", "JAK1/TYK2", "PIK3CA/mTOR",
    "AChE/BChE", "F2/F10", "PPARG/PPARA", "PPARA/PPARD",
]
CHECKS = []


def add(name: str, ok: bool, detail: str = "") -> None:
    CHECKS.append((name, "PASS" if ok else "FAIL", detail))


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def f4(value: str) -> str:
    return f"{float(value):.4f}"


def activity_source_status(si: str) -> tuple[bool, str]:
    """Check the restored ChEMBL 37 attachment against the current source wording."""
    attach = ROOT / "activity_label_source_audit_20261008"
    problems = []
    wide = read_csv(attach / "chembl37_dual_endpoint_crosswalk_805.csv")
    long_rows = read_csv(attach / "chembl37_activity_rows_long.csv")
    cover = read_csv(attach / "document_identifier_coverage_805.csv")
    diff = read_csv(attach / "historical_adjudication_differences.csv")
    keys = [(row["pair"], row["panel_id"]) for row in wide]
    if len(wide) != 805 or len(set(keys)) != 805:
        problems.append(f"wide {len(wide)} unique {len(set(keys))}")
    if any(int(row["n_rows_A"]) < 1 or int(row["n_rows_B"]) < 1 for row in wide):
        problems.append("missing endpoint")
    if len(long_rows) != 3555:
        problems.append(f"long {len(long_rows)}")
    if any(row["class_matches_frozen"] != "YES" for row in wide):
        problems.append("class replay")
    if any(float(row["delta_vs_panel_A"]) != 0 or float(row["delta_vs_panel_B"]) != 0 for row in wide):
        problems.append("panel delta")
    both = sum(row["tie_literature_A"] == "present" and row["tie_literature_B"] == "present" for row in wide)
    if both != 794:
        problems.append(f"wide literature {both}")
    if len(cover) != 805:
        problems.append(f"cover {len(cover)}")
    cover_yes = [row for row in cover if row["all_tied_maxima_both_ends"] == "YES"]
    cover_no = [row for row in cover if row["all_tied_maxima_both_ends"] == "NO"]
    if len(cover_yes) != 794 or len(cover_no) != 11:
        problems.append(f"cover split {len(cover_yes)}/{len(cover_no)}")
    if any("全部并列" not in row["coverage_rule"] for row in cover):
        problems.append("coverage rule")
    if any(row["raw_activity_present_both_ends"] != "YES" or "原始活性存在" not in row["judgment"] or not row["dataset_titles"] or not row["document_urls"] for row in cover_no):
        problems.append("dataset gap rows")
    if any("缺少原始活性" in row["judgment"] for row in cover_no):
        problems.append("dataset gap misstated")
    expected = {
        "PM48_04": ("8.85", "9.35", "8.4", "9.35", "dual"),
        "PM48_05": ("10", "8.52", "9.12", "7.32", "dual"),
        "PM48_22": ("8.77", "5.83", "8.77", "5.52", "A_only"),
    }
    for ligand, values in expected.items():
        rows = [row for row in diff if row["ligand"] == ligand]
        if not rows or any((row["dump_panel_max_A"], row["dump_panel_max_B"], row["audit_aggregate_max_A"], row["audit_aggregate_max_B"], row["frozen_class"]) != values for row in rows):
            problems.append(f"diff {ligand}")
        elif any(row["substantive_reason"] != "排除步骤可追溯；实质审核理由待核实；冻结类别一致。" or row["raw_record_called_unreliable"] != "NO" for row in rows):
            problems.append(f"diff reason {ligand}")
    phrases = [
        "805 条合格观察均已与恢复的 ChEMBL 37 原始活性记录建立对应",
        "阈值 6.0 回放类别与冻结类别 805/805 一致",
        "该核对不等同于逐篇人工核查原始论文",
        "也不表示所有记录来自统一测定平台",
        "11 条数据集文档缺少指定文献标识",
        "三条后续审核差异的实质理由尚待核实",
        "Dual_Target_Docking/activity_label_source_audit_20261008/",
    ]
    for phrase in phrases:
        if phrase not in si:
            problems.append("SI missing " + phrase)
    if "不能写成 805 条标签都已逐条核对" in si:
        problems.append("old SI sentence")
    register = (MAN / "Source_Register.csv").read_text()
    if "ACTIVITY-CHEMBL37" not in register or "activity_label_source_audit_20261008" not in register:
        problems.append("register attachment")
    if "不能写成 805 条标签均已逐条回溯" in register:
        problems.append("old register sentence")
    return (not problems, "; ".join(problems))


def main() -> int:
    results = (MAN / "Results_CN.md").read_text()
    captions = (MAN / "Figure_Captions_CN.md").read_text()
    si = (MAN / "Supporting_Information_CN.md").read_text()
    joined = results + "\n" + captions + "\n" + si

    pop = read_csv(ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv")
    mapping = read_csv(ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/pair_ligand_mapping.csv")
    eligible = [r for r in pop if r["activity_eligible"] == "1" and r["class"] in {"dual", "A_only", "B_only", "neither"}]
    counts = Counter(r["class"] for r in eligible)
    add("counts_805", len(eligible) == 805 and counts["dual"] == 234 and counts["A_only"] == 238 and counts["B_only"] == 232 and counts["neither"] == 101, str(dict(counts)))
    add("population_807_mapping_808", len(pop) == 807 and len(mapping) == 808, f"{len(pop)}/{len(mapping)}")
    add("938_identity", 2 * counts["dual"] + counts["A_only"] + counts["B_only"] == 938, "")

    layers = read_csv(ROOT / "analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv")
    add("layers", sum(int(r["n_labeled_population"]) for r in layers) == 938 and sum(int(r["n_m0_directional_population"]) for r in layers) == 934 and sum(int(r["n_final_chemistry_oof"]) for r in layers) == 928, "")

    primary = read_csv(ROOT / "results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv")
    deltas = read_csv(ROOT / "results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv")
    add("primary_rows", len(primary) == 120 and len(deltas) == 32, "")
    egfr_b = next(r for r in primary if r["pair"] == "EGFR/HER2" and r["method"] == "M0" and r["contrast"] == "dual_vs_B_only")
    egfr_s = next(r for r in primary if r["pair"] == "EGFR/HER2" and r["method"] == "M0" and r["contrast"] == "summary_min")
    jak_s = next(r for r in primary if r["pair"] == "JAK1/JAK2" and r["method"] == "M0" and r["contrast"] == "summary_min")
    jak_a = next(r for r in primary if r["pair"] == "JAK1/JAK2" and r["method"] == "M0" and r["contrast"] == "dual_vs_A_only")
    add("summary_min_ci_not_copied", jak_s["ci_hi"] != jak_a["ci_hi"], f"{jak_s['ci_hi']} vs {jak_a['ci_hi']}")
    for token in (f4(egfr_b["auroc"]), f4(jak_s["ci_hi"]), f4(jak_a["ci_hi"])):
        add(f"text_has_{token}", token in results, token)

    negative = read_csv(ROOT / "results/jcim_stage2/module_B_negative_class.csv")
    pocket = read_csv(ROOT / "results/jcim_stage2/module_C_wrong_pocket.csv")
    signs_b = Counter("+" if float(r["delta_negative"]) > 0 else "-" for r in negative)
    signs_c = Counter("0" if float(r["delta_pocket"]) == 0 else "+" if float(r["delta_pocket"]) > 0 else "-" for r in pocket)
    add("B_11_5", signs_b["+"] == 11 and signs_b["-"] == 5, str(signs_b))
    add("C_11_4_1", signs_c["+"] == 11 and signs_c["-"] == 4 and signs_c["0"] == 1, str(signs_c))
    add("extra_drop_zero", all(r["extra_drop_vs_directional_m0"] == "0" for r in pocket), "")

    chem = read_csv(ROOT / "results/jcim_stage2/module_D_chemistry_units.csv")
    ecfp = [float(r["auroc_ecfp4"]) for r in chem]
    inc = [float(r["delta_ecfp4_plus_m0"]) for r in chem]
    add("ecfp_above_matched", all(float(r["auroc_ecfp4"]) > float(r["auroc_SUPPORTING_CHEMISTRY_MATCHED_M0"]) for r in chem), "")
    add("text_ecfp_range", f4(min(ecfp)) in results and f4(max(ecfp)) in results, f"{f4(min(ecfp))} {f4(max(ecfp))}")
    add("text_increment_range", f4(min(inc)) in results and f4(max(inc)) in results, f"{f4(min(inc))} {f4(max(inc))}")

    oof = read_csv(ROOT / "results/jcim_stage2/module_D_chemistry_oof_ligands.csv")
    add("oof_928_includes_fold0", len(oof) == 928 and any(r["fold_id"] == "0" for r in oof), "")

    compare = read_csv(ROOT / "results/jcim_stage2/module_I_method_compare.csv")
    counted = Counter((r["method"], r["concordance"]) for r in compare)
    stored = read_csv(ROOT / "results/jcim_stage2_figures/module_I_concordance_counts.csv")
    mismatch = []
    for row in stored:
        for state in ("CONCORDANT_IMPROVE", "CONCORDANT_WORSEN", "DISCORDANT", "NO_CHANGE"):
            if int(row[state]) != counted[row["method"], state]:
                mismatch.append((row["method"], state, row[state], counted[row["method"], state]))
    add("concordance_independent_match", not mismatch, str(mismatch))
    add("m2_no_change_zero", counted["M2", "NO_CHANGE"] == 0, "")
    add("text_mentions_m2_zero", "M2 为 4/3/1/0" in results and "M2 没有无变化" in results, "")

    required_prefixes = [
        "Table S2", "Figures S1–S2", "Table S3", "Table S4", "Figure 2", "Table S5",
        "Figure 4", "Table S6", "Figure S5", "Tables S10–S11", "Figure S6", "Table S12", "Table S17",
    ]
    for token in required_prefixes:
        add(f"prefix_{token}", token in joined, token)
    add("bare_s6_not_used_for_pocket", "及S6" not in joined and "在S6" not in joined, "")
    add("activity_source_attachment", *activity_source_status(si))
    add("figure5a_scope", "Figure 5a 只显示 summary_min" in results and "两个方向的完整结果在 Figure S5" in results, "")

    stems = [
        "Scheme01_Study_Design",
        "Fig01_Primary_Directional_and_Deltas",
        "Fig02_Negative_Class_Contrast",
        "Fig03_Chemistry_and_Docking_Increment",
        "Fig04_Pocket_Contrast",
        "Fig05_Seed_and_Receptor_Sensitivity",
        "Fig06_Directional_and_Topk_Changes",
        "FigS01_Descriptors_D_vs_A",
        "FigS02_Descriptors_D_vs_B",
        "FigS03_Cross_Class_Similarity",
        "FigS04_All_Method_AUROC",
        "FigS05_Five_Seeds",
        "FigS06_Alternative_Receptors",
        "FigS07_Missingness",
        "FigS08_Topk_and_Concordance",
    ]
    missing_fig = [stem for stem in stems if not all((ROOT / "results/jcim_stage2_figures" / f"{stem}{ext}").exists() for ext in (".pdf", ".svg", ".png"))]
    add("figure_files", not missing_fig, ",".join(missing_fig))

    zip_path = MAN / "Supplementary_Source_Data.zip"
    if zip_path.exists():
        with zipfile.ZipFile(zip_path) as archive:
            bad = []
            for info in archive.infolist():
                if info.filename == "SHA256SUMS.csv":
                    continue
                source = REPO / info.filename
                packed = hashlib.sha256(archive.read(info.filename)).hexdigest()
                if not source.exists() or hashlib.sha256(source.read_bytes()).hexdigest() != packed:
                    bad.append(info.filename)
        add("zip_bytes", not bad, str(bad[:3]))
    else:
        add("zip_bytes", False, "missing")

    protected = [
        "Dual_Target_Docking/scripts/rerun_v4_2/jcim_stage2_lib.py",
        "Dual_Target_Docking/scripts/rerun_v4_2/jcim_stage2_compute.py",
        "Dual_Target_Docking/scripts/rerun_v4_2/jcim_stage2_audit.py",
        "Dual_Target_Docking/scripts/rerun_v4_2/jcim_stage2_impl_check.py",
        "Dual_Target_Docking/scripts/rerun_v4_2/jcim_stage2_result_audit.py",
        "Dual_Target_Docking/results/jcim_stage2/RUN_MANIFEST.json",
        "Dual_Target_Docking/results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv",
        "Dual_Target_Docking/results/jcim_stage2_figures/module_I_concordance_counts.csv",
    ]
    dirty = []
    for rel in protected:
        diff = subprocess.check_output(["git", "diff", "--name-only", "HEAD", "--", rel], cwd=REPO, text=True).strip()
        if diff:
            dirty.append(rel)
    add("protected_unchanged_vs_head", not dirty, ",".join(dirty))

    failed = [item for item in CHECKS if item[1] != "PASS"]
    report = MAN / "source_views" / "display_audit_report.txt"
    report.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{status}\t{name}\t{detail}" for name, status, detail in CHECKS]
    lines.append(f"AUDIT {len(CHECKS) - len(failed)}/{len(CHECKS)} PASS")
    report.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
