#!/usr/bin/env python3
"""Build popular-science Word document and data_tables folder from project report."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Cm, Pt
from docx.oxml.ns import qn

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "popular_science"
DATA_DIR = OUT / "data_tables"
FIG_DIR = OUT / "figures"
DOC_PATH = OUT / "JNK1筛选项目通俗解读报告.docx"

plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial Unicode MS", "SimHei", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

# --- Data tables to copy / synthesize -----------------------------------------

COPY_MAP = {
    "01_文献benchmark化合物.csv": "data/benchmarks/literature_benchmarks.csv",
    "02_ML模型性能comparison.json": "results/model_comparison/comparison.json",
    "03_ML外部decoy验证指标.json": "results/ml_external_validation/ml_external_validation_metrics_9bd8.json",
    "04_再对接验证redocking.csv": "results/docking_validation/redocking_summary_7725.csv",
    "05_benchmark对接差值.csv": "results/docking_validation/benchmark_deltas_51c1.csv",
    "06_benchmark方向混淆矩阵.csv": "results/docking_validation/direction_confusion_27c3.csv",
    "07_benchmark_MMGBSA标定.csv": "results/docking_validation/benchmark_mmgbsa_calibration.csv",
    "08_Gly87自检.csv": "results/docking_validation/gly87_selfcheck_16be.csv",
    "09_各亚型排序相关.csv": "results/docking_validation/isoform_rank_correlations_299a.csv",
    "10_选择性标签统计.csv": "results/similarity/sel_class_counts.csv",
    "11_采购清单purchase_after_md.csv": "data/purchase/purchase_after_md.csv",
    "12_2231延伸MD_RMSD分位数.csv": "results/md_2231_200ns/tables/09_production_rmsd_percentiles.csv",
    "13_2231延伸MD_MMGBSA分量.csv": "results/md_2231_200ns/tables/14_mm_gbsa_components.csv",
    "14_ML阈值校准threshold.json": "results/calibration/threshold_recommendation.json",
    "15_ML虚拟筛选demo_screening.json": "results/screening_v2/screening_report.json",
}


def synthesize_tables() -> None:
  """Create summary tables documented in JNK1_PROJECT_REPORT.md."""
  pd.DataFrame([
    {"阶段": "ML初筛后对接库(F0)", "数量": 4983, "数据来源": "md_shortlist_report_23c8.md"},
    {"阶段": "Glide XP VSW有效记录", "数量": 4979, "数据来源": "JNK1_SELECTIVITY_FINAL_REPORT_41d9.md"},
    {"阶段": "VSW pass_selectivity(遗留探索)", "数量": 233, "数据来源": "同上;未用于MD短名单"},
    {"阶段": "VSW Tier 1′", "数量": 57, "数据来源": "同上"},
    {"阶段": "MD短名单 F1∧F2(ADMET前)", "数量": 157, "数据来源": "md_shortlist_report_23c8.md"},
    {"阶段": "MD shortlist(ADMET后)", "数量": 25, "数据来源": "同上"},
    {"阶段": "MD pose QC输入", "数量": 16, "数据来源": "MD_QC_report_cf26.md"},
    {"阶段": "最终采购推荐", "数量": 10, "数据来源": "data/purchase/purchase_after_md.csv"},
  ]).to_csv(DATA_DIR / "00_端到端漏斗汇总.csv", index=False, encoding="utf-8-sig")

  pd.DataFrame([
    {"阶段": "ML F1后", "数量": 4983, "条件": "p_family ≥ 6.0"},
    {"阶段": "VSW有效", "数量": 4979, "条件": "三isoform均有XP分"},
    {"阶段": "pass_pose", "数量": 3234, "条件": "Glide pose质量门"},
    {"阶段": "pass_potency", "数量": 1681, "条件": "score_JNK1 ≤ -7.43"},
    {"阶段": "pass_selectivity(遗留)", "数量": 233, "条件": "Δsel_dock>0 且 Δsel_MMGBSA≥2"},
    {"阶段": "has_selectivity_contact", "数量": 63, "条件": "铰链H-bond代理+Δsel启发式"},
  ]).to_csv(DATA_DIR / "16_VSW分支A各阶段数量.csv", index=False, encoding="utf-8-sig")

  pd.DataFrame([
    {"Tier": "Tier 1′", "数量": 57, "条件": "pose+potency+pass_selectivity+contact"},
    {"Tier": "Tier 2", "数量": 92, "条件": "pose+potency+pass_selectivity"},
    {"Tier": "Tier 3", "数量": 1191, "条件": "pose+potency"},
    {"Tier": "Tier 0", "数量": 3639, "条件": "未达Tier3"},
  ]).to_csv(DATA_DIR / "17_Tier分布.csv", index=False, encoding="utf-8-sig")

  pd.DataFrame([
    {"阶段": "输入(F0后)", "数量": 4983},
    {"阶段": "F1 pose QC通过", "数量": 3125},
    {"阶段": "F2活性+配体效率通过", "数量": 182},
    {"阶段": "F1∧F2通过", "数量": 157},
    {"阶段": "F7 ADMET剔除", "数量": 9},
    {"阶段": "ADMET后shortlist", "数量": 25},
    {"阶段": "进入MD", "数量": 16},
    {"阶段": "采购推荐", "数量": 10},
  ]).to_csv(DATA_DIR / "18_MD分支B漏斗.csv", index=False, encoding="utf-8-sig")

  pd.DataFrame([
    {"指标": "Spearman(Δsel_dock_vsw, -ΔpIC50_sel)", "数值": 0.750, "n": 7},
    {"指标": "Spearman(Δsel_mmgbsa_vsw, -ΔpIC50_sel)", "数值": 0.786, "n": 7},
    {"指标": "方向准确率(docking,VSW PDB)", "数值": "43%(3/7)", "n": 7},
    {"指标": "方向准确率(MM-GBSA,VSW PDB)", "数值": "43%(3/7)", "n": 7},
    {"指标": "方向准确率(docking,ensemble归档)", "数值": "29%(2/7)", "n": 7},
  ]).to_csv(DATA_DIR / "19_benchmark定量结果.csv", index=False, encoding="utf-8-sig")

  pd.DataFrame([
    {"ID": 2231, "组": "G2", "JNK1": 0.91, "JNK2": 0.00, "JNK3": 0.10},
    {"ID": 2157, "组": "G1", "JNK1": 0.85, "JNK2": 0.46, "JNK3": 0.02},
    {"ID": 2232, "组": "G1", "JNK1": 1.00, "JNK2": 1.00, "JNK3": 0.04},
    {"ID": 690, "组": "G1", "JNK1": 1.00, "JNK2": 0.51, "JNK3": 0.77},
    {"ID": 1280, "组": "G2", "JNK1": 0.00, "JNK2": 0.58, "JNK3": 0.43},
    {"ID": 4795, "组": "G2", "JNK1": 0.04, "JNK2": 0.04, "JNK3": 0.49},
  ]).to_csv(DATA_DIR / "20_采购分子铰链占有率.csv", index=False, encoding="utf-8-sig")

  pd.DataFrame([
    {"ID": 2231, "JNK1": 0.48, "JNK2": 1.17, "JNK3": 0.66},
    {"ID": 2157, "JNK1": 0.49, "JNK2": 1.13, "JNK3": 0.35},
    {"ID": 690, "JNK1": 0.72, "JNK2": 1.98, "JNK3": 0.61},
    {"ID": 2232, "JNK1": 0.57, "JNK2": 0.29, "JNK3": 1.31},
  ]).to_csv(DATA_DIR / "21_采购分子配体RMSD中位数_Angstrom.csv", index=False, encoding="utf-8-sig")

  pd.DataFrame([
    {"亚型": "JNK1", "化合物数": 444, "Holdout_R2": 0.697, "Holdout_Spearman": 0.858, "Holdout_n": 31},
    {"亚型": "JNK2", "化合物数": 610, "Holdout_R2": 0.574, "Holdout_Spearman": 0.780, "Holdout_n": 67},
    {"亚型": "JNK3", "化合物数": 1147, "Holdout_R2": 0.774, "Holdout_Spearman": 0.869, "Holdout_n": 98},
  ]).to_csv(DATA_DIR / "22_ML模型性能汇总.csv", index=False, encoding="utf-8-sig")

  pd.DataFrame([
    {"组分": "Decoys(Taosu)", "n": 10000, "标签": "假定无活性"},
    {"组分": "Benchmarks", "n": 9, "标签": "已知活性"},
    {"组分": "ChEMBL actives", "n": 1210, "标签": "已知活性"},
  ]).to_csv(DATA_DIR / "23_外部decoy验证设计.csv", index=False, encoding="utf-8-sig")

  pd.DataFrame([
    {"指标": "Sensitivity(recall)", "数值": "99.3%"},
    {"指标": "Specificity", "数值": "4.7%"},
    {"指标": "Decoy FPR", "数值": "95.3%"},
    {"指标": "Precision", "数值": "11.3%"},
    {"指标": "ROC-AUC", "数值": 0.876},
    {"指标": "EF1%", "数值": 9.20},
  ]).to_csv(DATA_DIR / "24_外部decoy验证指标汇总.csv", index=False, encoding="utf-8-sig")


def setup_data_folder() -> None:
  DATA_DIR.mkdir(parents=True, exist_ok=True)
  FIG_DIR.mkdir(parents=True, exist_ok=True)
  for dst, src in COPY_MAP.items():
    shutil.copy2(ROOT / src, DATA_DIR / dst)
  synthesize_tables()
  readme = """# 数据表格文件夹

本文件夹收录 JNK1/2/3 计算筛选项目报告中引用的全部可复现数据表。
所有数值均可回溯至仓库内原始文件；合成汇总表（00、16–24 号）的字段与
`docs/JNK1_PROJECT_REPORT.md`（v2.7）正文一致。

## 文件索引

| 文件名 | 说明 | 原始路径 |
|--------|------|----------|
"""
  for dst, src in sorted(COPY_MAP.items()):
    readme += f"| {dst} | 见文件名 | `{src}` |\n"
  readme += """
| 00_端到端漏斗汇总.csv | 摘要漏斗各阶段数量 | 报告§摘要 |
| 16–24 号 CSV | 报告正文汇总表 | 报告对应章节 |

生成脚本：`scripts/build_popular_science_doc.py`
"""
  (DATA_DIR / "README.md").write_text(readme, encoding="utf-8")


def make_figures() -> dict[str, Path]:
  paths = {}

  funnel = pd.read_csv(DATA_DIR / "00_端到端漏斗汇总.csv")
  fig, ax = plt.subplots(figsize=(10, 5))
  y = funnel["数量"][::-1]
  labels = funnel["阶段"][::-1]
  ax.barh(range(len(y)), y, color="#4C72B0")
  ax.set_yticks(range(len(labels)))
  ax.set_yticklabels(labels, fontsize=9)
  ax.set_xlabel("化合物数量")
  ax.set_title("端到端筛选漏斗（有数据支撑的阶段）")
  for i, v in enumerate(y):
    ax.text(v + 30, i, str(int(v)), va="center", fontsize=9)
  fig.tight_layout()
  p = FIG_DIR / "fig01_funnel.png"
  fig.savefig(p, dpi=150)
  plt.close(fig)
  paths["funnel"] = p

  redock = pd.read_csv(DATA_DIR / "04_再对接验证redocking.csv")
  fig, ax = plt.subplots(figsize=(7, 4))
  ax.bar(redock["PDB_ID"], redock["rmsd_A"], color=["#55A868" if x else "#C44E52" for x in redock["pass(<2A)"]])
  ax.axhline(2.0, color="red", linestyle="--", label="阈值 2.0 Å")
  ax.set_ylabel("RMSD (Å)")
  ax.set_title("共晶再对接验证（5/5 通过）")
  ax.legend()
  fig.tight_layout()
  p = FIG_DIR / "fig02_redocking.png"
  fig.savefig(p, dpi=150)
  plt.close(fig)
  paths["redock"] = p

  bench = pd.read_csv(DATA_DIR / "06_benchmark方向混淆矩阵.csv")
  fig, ax = plt.subplots(figsize=(6, 4))
  colors = ["#55A868" if m else "#C44E52" for m in bench["direction_match"]]
  ax.bar(bench["name"], bench["delta_sel_dock"], color=colors)
  ax.axhline(0, color="black", linewidth=0.8)
  ax.set_ylabel("Δsel_dock (kcal/mol 量级)")
  ax.set_title("文献对照：对接选择性方向预测（绿=正确，红=错误）")
  plt.xticks(rotation=30, ha="right")
  fig.tight_layout()
  p = FIG_DIR / "fig03_benchmark_direction.png"
  fig.savefig(p, dpi=150)
  plt.close(fig)
  paths["bench_dir"] = p

  rmsd2231 = pd.read_csv(DATA_DIR / "12_2231延伸MD_RMSD分位数.csv")
  lig = rmsd2231[rmsd2231["Metric"] == "Ligand"]
  fig, ax = plt.subplots(figsize=(6, 4))
  systems = lig["System"].tolist()
  med = lig["median"].tolist()
  desmond = [0.48, 1.17, 0.66]
  x = np.arange(len(systems))
  w = 0.35
  ax.bar(x - w/2, med, w, label="Amber 200ns 中位数")
  ax.bar(x + w/2, desmond, w, label="Desmond 短MD 中位数(报告§6.4.2)")
  ax.set_xticks(x)
  ax.set_xticklabels(systems)
  ax.set_ylabel("配体 RMSD (Å)")
  ax.set_title("化合物 2231：三 isoform 配体稳定性对比")
  ax.legend(fontsize=8)
  fig.tight_layout()
  p = FIG_DIR / "fig04_2231_rmsd.png"
  fig.savefig(p, dpi=150)
  plt.close(fig)
  paths["2231_rmsd"] = p

  hinge = pd.read_csv(DATA_DIR / "20_采购分子铰链占有率.csv")
  fig, ax = plt.subplots(figsize=(8, 4))
  ids = hinge["ID"].astype(str)
  x = np.arange(len(ids))
  w = 0.25
  ax.bar(x - w, hinge["JNK1"], w, label="JNK1")
  ax.bar(x, hinge["JNK2"], w, label="JNK2")
  ax.bar(x + w, hinge["JNK3"], w, label="JNK3")
  ax.axhline(0.3, color="gray", linestyle="--", label="pass 阈值 0.3")
  ax.set_xticks(x)
  ax.set_xticklabels(ids)
  ax.set_ylabel("铰链氢键占有率 (0–1)")
  ax.set_title("采购分子：三亚型铰链氢键占有率")
  ax.legend(fontsize=8)
  fig.tight_layout()
  p = FIG_DIR / "fig05_hinge_occupancy.png"
  fig.savefig(p, dpi=150)
  plt.close(fig)
  paths["hinge"] = p

  metrics = pd.read_csv(DATA_DIR / "19_benchmark定量结果.csv")
  fig, ax = plt.subplots(figsize=(7, 4))
  dirs = ["VSW对接", "VSW MM-GBSA", "ensemble对接"]
  acc = [43, 43, 29]
  ax.bar(dirs, acc, color="#4C72B0")
  ax.axhline(55, color="red", linestyle="--", label="期望阈值 55%")
  ax.set_ylabel("方向准确率 (%)")
  ax.set_title("计算选择性方向预测准确率（均不达标）")
  ax.set_ylim(0, 70)
  ax.legend()
  fig.tight_layout()
  p = FIG_DIR / "fig06_direction_accuracy.png"
  fig.savefig(p, dpi=150)
  plt.close(fig)
  paths["dir_acc"] = p

  return paths


def set_cell_font(cell, size=9):
  for p in cell.paragraphs:
    for r in p.runs:
      r.font.size = Pt(size)


def add_table(doc: Document, df: pd.DataFrame, title: str, max_rows: int | None = None):
  doc.add_heading(title, level=3)
  if max_rows:
    df = df.head(max_rows)
  table = doc.add_table(rows=1, cols=len(df.columns))
  table.style = "Table Grid"
  hdr = table.rows[0].cells
  for i, c in enumerate(df.columns):
    hdr[i].text = str(c)
    set_cell_font(hdr[i], 9)
  for _, row in df.iterrows():
    cells = table.add_row().cells
    for i, v in enumerate(row):
      cells[i].text = str(v) if pd.notna(v) else ""
      set_cell_font(cells[i], 8)
  doc.add_paragraph()


def add_image(doc: Document, path: Path, caption: str, width_cm=15):
  doc.add_picture(str(path), width=Cm(width_cm))
  p = doc.add_paragraph(caption)
  p.alignment = WD_ALIGN_PARAGRAPH.CENTER
  if p.runs:
    p.runs[0].italic = True
  doc.add_paragraph()


def build_word(fig_paths: dict[str, Path]) -> None:
  doc = Document()
  style = doc.styles["Normal"]
  style.font.name = "宋体"
  style._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
  style.font.size = Pt(11)

  title = doc.add_heading("JNK1/2/3 抑制剂计算机辅助筛选项目\n通俗解读报告", 0)
  title.alignment = WD_ALIGN_PARAGRAPH.CENTER
  doc.add_paragraph("版本：2.7（面向非 CADD 专业读者）")
  doc.add_paragraph("数据来源：项目仓库 JNK1_Selectivity_Project，所有数值可回溯至 data_tables/ 文件夹")
  doc.add_paragraph("技术报告原文：docs/JNK1_PROJECT_REPORT.md")
  doc.add_page_break()

  sections = [
    ("写给完全不懂计算机辅助药物设计（CADD）的读者", """
本报告用通俗语言说明：我们如何从约五千个候选分子出发，经过机器学习、分子对接、成药性评估和分子动力学模拟，最终选出 10 个分子去做湿实验（试管酶学测定）。

需要事先明确的三件事：
第一，这是一套「计算预筛选」流程，目的是在花钱合成/购买化合物之前，用计算机尽量排除明显不靠谱的分子。
第二，本项目最初希望找到「主要抑制 JNK1、少抑制 JNK2/JNK3」的选择性分子，但计算验证表明：现有方法无法可靠判断亚型方向，因此策略已调整为寻找「可能与 JNK 家族结合的候选分子」，选择性只能交给实验。
第三，文中每一个数字都有原始数据文件支撑，汇总表见本报告配套文件夹 data_tables/。
"""),
    ("核心术语解释（阅读前请先浏览）", """
• JNK（c-Jun N-terminal kinase）：一类与炎症、纤维化等相关的蛋白激酶，有 JNK1、JNK2、JNK3 三个亚型（同工酶）。
• 亚型 / isoform：同一基因家族下结构相似但功能略有差别的蛋白版本。
• 抑制剂 / 配体 / 小分子：有望结合蛋白并调节其活性的化合物。
• IC50：使酶活性下降一半所需化合物浓度，越小通常活性越强。
• pIC50 = −log10(IC50/M)：活性常用对数刻度，数值越大活性越强。
• ChEMBL：公开的药物化学数据库，收录化合物生物活性（参考文献 [1]）。
• 虚拟筛选：用计算机从大量分子中快速缩小候选范围，而非逐个做实验。
• 机器学习（ML）：让计算机从已知「分子结构—活性」数据中学习规律，再预测新分子活性。
• XGBoost：一种常用的机器学习算法（参考文献 [9]）。
• 分子对接（docking）：把小分子「摆放」到蛋白结合口袋，用打分函数估计结合强弱（Glide，参考文献 [7]）。
• PDB：蛋白三维结构数据库编号；本项目用 3ELJ(JNK1)、3E7O(JNK2)、3TTI(JNK3) 等共晶结构。
• RMSD：两个分子结构相差多少埃（Å），越小越相似；再对接 RMSD < 2 Å 表示计算能复现实验结构。
• MM-GBSA：在对接 pose 基础上估算结合自由能（本项目用 Schrödinger Prime 实现）。
• ADMET：吸收、分布、代谢、排泄、毒性等成药性指标；本项目用 QikProp 预测。
• 分子动力学（MD）：模拟分子在溶液中随时间运动，观察结合是否稳定。
• Decoy（诱饵分子）：假定无活性的对照分子，用于检验筛选是否会把「噪音」当成活性。
• pan-JNK：对三个亚型都有一定抑制活性的「广谱」分子，而非只偏向 JNK1。
"""),
    ("一、项目背景：我们在找什么药？", """
JNK 三个亚型由不同基因编码（JNK1/MAPK8、JNK2/MAPK9、JNK3/MAPK10）。JNK1 与特发性肺纤维化（IPF）、非酒精性脂肪性肝炎（NASH）等疾病相关。

临床上已有 JNK 相关药物研究，例如 CC-90001（参考文献 [2]）在细胞层面偏向 JNK1 功能抑制；SP600125（参考文献 [3]）是经典的 pan-JNK 研究工具药。本项目还纳入文献报道的强 JNK1 抑制剂 E1（参考文献 [4]）作为对照。

项目初期目标：通过计算筛选找到 preferentially 抑制 JNK1 的候选分子。
项目后期结论：计算无法可靠区分亚型方向，因此把候选定位为「JNK 家族结合剂」，用同批次 JNK1/2/3 酶学 IC50 实验做最终判断（参考文献 [8] 综述了 JNK 成药挑战）。
"""),
    ("二、整体筛选流程（双分支设计）", """
从 Enamine 商业库约 5000 个分子出发，经过以下主链路：

步骤 1 — 机器学习粗筛：用三个亚型活性模型判断分子是否「可能有 JNK 家族活性」，保留 4983 个。
步骤 2 — Glide XP 对接：把每个分子分别对接到 JNK1/JNK2/JNK3 三个蛋白结构，得到 4979 个有效结果。
步骤 3 — 分叉：
  • 分支 A（遗留探索标签）：根据对接分数打 Tier、pass_selectivity 等标签，得到 233/57 等统计数字，仅用于回顾分析，不决定买哪些分子。
  • 分支 B（真正采购决策链）：按 pose 质量、JNK1 活性门槛、成药性、ADMET 筛选 → 157 → 25 → 16 个做 MD → 最终 10 个采购。

为什么要有两个分支？因为早期流水线曾尝试用「对接能量差 Δsel」判断 JNK1 选择性，但文献 benchmark 证明该指标方向准确率仅 43%（VSW 单结构）或 29%（双结构平均），不能用于采购决策。MD 短名单因此完全不用 Δsel 标签。
"""),
    ("三、机器学习：如何粗筛「可能有活性」的分子？", """
3.1 训练数据与模型
从 ChEMBL（参考文献 [1]）收集 JNK1/2/3 活性数据，清洗后：JNK1 444 个、JNK2 610 个、JNK3 1147 个化合物。
采用 XGBoost（参考文献 [9]）分别训练三个回归模型，预测 pActivity（活性对数刻度）。
数据划分采用 Murcko 骨架分组：化学结构相似的分子不会同时出现在训练集和测试集，避免「背答案」。

独立 holdout 测试集性能（数据表 22_ML模型性能汇总.csv）：
JNK1 R²=0.697，JNK2 R²=0.574，JNK3 R²=0.774。
说明：模型对「这个分子强不强」有一定预测力，但对「偏向哪个亚型」不可靠。

3.2 为什么用 F1 门槛 p_family ≥ 6.0？
九个文献对照分子在 6.0 阈值下全部通过（9/9），目的是高召回——不让已知活性分子被漏掉。
外部 decoy 验证（10,000 个 Taosu 随机分子 + 1,219 个已知活性）显示：decoy 假阳性率 95.3%，特异性仅 4.7%，ROC-AUC=0.876，EF1%=9.20（数据表 23–24）。
解读：F1 是「宽网」，负责不漏真活性；真正去噪音靠综合排序分 final_score（含 QED、合成可及性 SA 等）和后续对接。

3.3 为什么 ML 不能判断亚型选择性？
ChEMBL 中明确标注 JNK1-selective 的分子仅 8 个（数据表 10），数据太少无法训练可靠分类器（测试 F1=0）。
对 E1、TCS JNK 6O 等对照，ML 预测的最高活性亚型与实验不符（数据表 05）。
"""),
    ("四、分子对接：如何把分子「放进」蛋白口袋？", """
4.1 软件与结构
使用 Schrödinger Glide XP 模式（参考文献 [7]）。每个亚型选一个主 PDB：JNK1→3ELJ，JNK2→3E7O，JNK3→3TTI。
打分采用 XP 终分 gscore，越负表示预测结合越强。

4.2 为什么要做共晶再对接？
在筛选大量分子之前，先拿共晶配体重新对接回自己的蛋白结构。若 RMSD < 2 Å，说明蛋白结构处理和对接参数可信。
结果 5/5 全部通过（数据表 04、图 2）。注意：这只能证明「结构准备正确」，不能证明「能预测选择性」。

4.3 选择性指标 Δsel 的含义与失败
Δsel_dock = min(score_JNK2, score_JNK3) − score_JNK1。若 > 0，计算上「更偏向 JNK1」。
对 9 个有实验 IC50 的文献分子做回顾：Spearman 秩相关尚可，但离散方向准确率仅 43%（3/7），低于 55% 可用线（数据表 19、图 6）。
因此：233 个 pass_selectivity 和 57 个 Tier 1′ 只是探索性统计，不用于 MD 进门。
"""),
    ("五、曾尝试的选择性策略与失败原因", """
5.1 Gly87 占据假说
JNK1 在铰链区有甘氨酸 Gly87，JNK2/JNK3 分别为 Ser/Met，体积更大。曾假设：配体占据 JNK1 特有的小空间可得选择性（机制类比参考文献 [10]）。
回顾测试显示：所有 benchmark 配体距 Gly87 仅 0.59–1.18 Å，E1 与 CC-930 均显示占据，无法区分（数据表 08）。策略放弃。

5.2 MM-GBSA 选择性门槛
曾要求 Δsel_MMGBSA ≥ 2 kcal/mol。标定发现：九个对照的 |Δsel| 中位数噪声达 8.13 kcal/mol，2.0 门槛远低于噪声，方向准确率仍仅 43%（数据表 07）。
注意：MMGBSA_JNK1 ≤ −51.6 作为「JNK1 单点活性」门槛仍用于 MD 短名单，与选择性差值无关。
"""),
    ("六、MD 短名单与成药性：如何选出 25 个、再做 16 个 MD？", """
6.1 分支 B 漏斗（数据表 18）
3125 个通过 pose 质量 → 182 个通过 JNK1 活性双门槛（Glide ≤ −7.43 且 MMGBSA_JNK1 ≤ −51.6）→ 157 个 F1∧F2 → ADMET 剔除 9 个 → 25 个 shortlist。
25 个按化学策略分四组：G1 文献 chemotype 邻近组(9)、G2 新骨架(10)、G3 已知活性对照(4)、G4 阴性锚点(2)。
按组配额取 16 个做 Desmond 分子动力学（48 个任务 = 16 分子 × 3 蛋白）。

6.2 为什么 2231 能进 MD 但 pass_selectivity=No？
2231 的 JNK1 对接分极强（−11.22），满足活性门槛；未过 pass_selectivity 是因为遗留双门槛中的 Δsel_MMGBSA 探索标签，MD 短名单从不读取该标签（详见技术报告 §6.3.1）。
"""),
    ("七、分子动力学 QC：结合姿势稳不稳？", """
对每个分子在三个亚型各跑一条 MD 轨迹（约 20–50 ns），用两项指标质检：
• 配体 RMSD ≤ 3 Å：分子在口袋里没有「飞出去」
• 铰链氢键占有率 ≥ 30%：与铰链区有稳定极性接触

pass_md_overall = JNK1 通过 且 (JNK2 或 JNK3 至少一个通过)。
结果：G1 组 3/4 通过，G2 组 0/6 通过——但 G2 的 2231 在 JNK1 上 hinge=0.91、RMSD 最低，仍有探索价值（图 5、数据表 20–21）。

重要：SP600125、E1 等已知活性药在 MD 指标上可 fail，说明 MD 不能代替活性实验，只能辅助判断 pose 可信度。
"""),
    ("八、化合物 2231 的 200 ns 延伸 MD", """
对 MD 相对排名第一的 G2 分子 2231，补做 200 ns Amber 模拟（三 isoform 各一条轨迹）。
与 Desmond 短 MD 对比：JNK1 配体 RMSD 中位数最低（Amber 0.57 Å vs Desmond 0.48 Å），JNK2 最高（图 4，数据表 12）。
JNK1 中 Asn108–配体氧原子氢键占有率约 68%（技术报告 §6.5.4）。

限制：仅单条轨迹、生产期对配体有位置约束、不能与 Desmond 数值直接对比绝对值。
结论：支持「2231 在 JNK1 中 pose 相对更稳」的方向性假说，但不能确认选择性，也不能改写 §6.2 中 overall MD fail 的记录。
"""),
    ("九、最终采购的 10 个分子与实验计划", """
采购清单见数据表 11，分为：
• G3 对照 4 个：SP600125、CC-90001、CC-930、E1——用于校准实验体系
• G1 主力 3 个：690、2232、2157——MD 总体通过、接近文献 chemotype
• G2 探索 3 个：2231、1280、4795——新骨架与 off-target pose 假说

花钱买的不是「已经算出的选择性 hit」，而是：
1) 验证计算管线能否富集有活性分子；2) 比较 G1 与 G2 哪类骨架更易出活性；3) 用 2231/2157 等检验 JNK1 偏好假说。
必做实验：同批次 JNK1 + JNK2 + JNK3 重组酶 IC50。
"""),
    ("十、诚实结论（给决策者与合作者）", """
1. 受体准备与对接流程可信（5/5 再对接）。
2. 机器学习能粗筛家族活性，不能判断亚型方向。
3. Glide Δsel 与 MM-GBSA Δsel 不能用于 isoform 分型（方向准确率 43% 或更低）。
4. MD 可辅助看结合姿势，不能单独确认选择性。
5. 没有任何分子被计算「确认」为 JNK1 选择性抑制剂；选择性只能由湿实验回答。

若只有经费买 2 个分子优先赌 JNK1 选择性：MD 可视化支持 2231 + 2157 的组合；690 更适合作为 pan-JNK 活性验证。
"""),
    ("参考文献与方法的准确对应", """
以下编号与项目技术报告 §11 一致；未列入者为本项目使用的商业软件/通用规则，不强行编造文献。

[1] ChEMBL 数据库 — 训练数据来源。Zdrazil B, et al. Nucleic Acids Res. 2024;52(D1):D1180-D1192. doi:10.1093/nar/gkad1004

[2] CC-90001 临床候选与 JNK1 功能偏向。Bennett BL, et al. J. Med. Chem. 2021;64(3):1776-1795. doi:10.1021/acs.jmedchem.0c01843. PMID:33404223

[3] SP600125 经典 pan-JNK 工具药。Bennett BL, et al. Proc. Natl. Acad. Sci. USA 2001;98(24):13681-13686. doi:10.1073/pnas.251194298. PMID:11717429

[4] E1 与 JNK1 抑制剂系列。Pan X, et al. J. Med. Chem. 2024. doi:10.1021/acs.jmedchem.4c01764

[5] TCS JNK 6O 氨基吡唑 JNK 抑制剂。Szczepankiewicz BG, et al. J. Med. Chem. 2006;49(14):3563-3566. doi:10.1021/jm060150w

[6] CC-930 共晶结构 3TTI。Plantevin-Krenitsky V, et al. Bioorg. Med. Chem. Lett. 2012;22(3):1433-1438. doi:10.1016/j.bmcl.2011.12.111

[7] Glide 分子对接。Friesner RA, et al. J. Med. Chem. 2004;47(7):1739-1749. doi:10.1021/jm0306430

[8] JNK 成药靶点综述。Manning BD, Davis RJ. Nat. Rev. Drug Discov. 2003;2(7):554-565. doi:10.1038/nrd1132

[9] XGBoost 机器学习。Chen T, Guestrin C. Proc. 22nd ACM SIGKDD 2016. doi:10.1145/2939672.2939785

[10] JNK2/3 选择性共价抑制剂与 Gly87/Leu106 机制讨论。Bennett BL, et al. J. Med. Chem. 2022. doi:10.1021/acs.jmedchem.2c01834

方法—工具补充说明（无额外编造文献）：
• MM-GBSA、Desmond MD、QikProp ADMET：Schrödinger 商业软件套件（Prime / Desmond / QikProp），流程参数见项目配置与工作区归档报告。
• 2231 延伸 MD：AMBER 分子动力学 + cpptraj 分析 + MMPBSA.py；力场与输入见 results/md_2231_200ns/README.md。
• Lipinski 五规则、Murcko 骨架划分、Butina 聚类、RDKit：药物化学与化学信息学常用开源/工业界标准流程；本项目按脚本实现，详见仓库 scripts/ 与 data/processed/。
• JNK-IN-8 共价抑制剂对照来源：Cell Chem Biol 2012（见 data/benchmarks/literature_benchmarks.csv 的 source 字段）。
"""),
  ]

  for heading, body in sections:
    doc.add_heading(heading, level=1)
    for para in body.strip().split("\n"):
      if para.strip():
        doc.add_paragraph(para.strip())
    doc.add_paragraph()

  doc.add_heading("附图与附表", level=1)
  add_image(doc, fig_paths["funnel"], "图1 端到端筛选漏斗各阶段化合物数量（数据表 00）")
  add_table(doc, pd.read_csv(DATA_DIR / "00_端到端漏斗汇总.csv"), "表1 端到端漏斗汇总")
  add_image(doc, fig_paths["redock"], "图2 五个蛋白结构的共晶再对接 RMSD（数据表 04）")
  add_table(doc, pd.read_csv(DATA_DIR / "04_再对接验证redocking.csv"), "表2 再对接验证明细")
  add_image(doc, fig_paths["dir_acc"], "图3 Benchmark 方向准确率对比（数据表 19）")
  add_image(doc, fig_paths["bench_dir"], "图4 四个关键对照的 Δsel_dock 与方向匹配（数据表 06）")
  add_table(doc, pd.read_csv(DATA_DIR / "06_benchmark方向混淆矩阵.csv"), "表3 Benchmark 方向混淆")
  add_table(doc, pd.read_csv(DATA_DIR / "22_ML模型性能汇总.csv"), "表4 ML 模型 holdout 性能")
  add_table(doc, pd.read_csv(DATA_DIR / "24_外部decoy验证指标汇总.csv"), "表5 外部 decoy 验证指标")
  add_table(doc, pd.read_csv(DATA_DIR / "18_MD分支B漏斗.csv"), "表6 MD 短名单漏斗（分支 B）")
  add_table(doc, pd.read_csv(DATA_DIR / "17_Tier分布.csv"), "表7 VSW Tier 分布（分支 A，探索用）")
  add_image(doc, fig_paths["hinge"], "图5 采购分子铰链氢键占有率（数据表 20）")
  add_table(doc, pd.read_csv(DATA_DIR / "11_采购清单purchase_after_md.csv"), "表8 最终采购清单（完整）")
  add_image(doc, fig_paths["2231_rmsd"], "图6 化合物 2231 延伸 MD 配体 RMSD 对比（数据表 12）")
  add_table(doc, pd.read_csv(DATA_DIR / "12_2231延伸MD_RMSD分位数.csv"), "表9 2231 延伸 MD RMSD 统计")
  add_table(doc, pd.read_csv(DATA_DIR / "13_2231延伸MD_MMGBSA分量.csv"), "表10 2231 MM-GBSA 分量（辅助记录，不作选择性裁决）")

  add_table(doc, pd.read_csv(DATA_DIR / "01_文献benchmark化合物.csv"), "表11 九个文献 benchmark 化合物与实验活性")
  add_table(doc, pd.read_csv(DATA_DIR / "08_Gly87自检.csv"), "表12 Gly87 占据策略回顾性测试")
  add_table(doc, pd.read_csv(DATA_DIR / "16_VSW分支A各阶段数量.csv"), "表13 VSW 分支 A 各过滤阶段数量")
  add_table(doc, pd.read_csv(DATA_DIR / "20_采购分子铰链占有率.csv"), "表14 采购分子铰链氢键占有率")
  add_table(doc, pd.read_csv(DATA_DIR / "21_采购分子配体RMSD中位数_Angstrom.csv"), "表15 采购分子配体 RMSD 中位数 (Å)")

  doc.add_heading("数据文件索引", level=1)
  doc.add_paragraph(
    f"全部 {len(list(DATA_DIR.glob('*.csv')))} 个 CSV 与 JSON 数据表位于：docs/popular_science/data_tables/\n"
    "详见该文件夹 README.md。"
  )

  doc.save(DOC_PATH)
  print(f"Wrote {DOC_PATH}")


def main():
  setup_data_folder()
  figs = make_figures()
  build_word(figs)
  print(f"Data tables: {DATA_DIR}")
  print(f"Figures: {FIG_DIR}")


if __name__ == "__main__":
  main()
