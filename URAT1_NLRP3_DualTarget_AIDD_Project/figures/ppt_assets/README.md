# PPT 汇报素材（`figures/ppt_assets/`）

本目录收录课题汇报 PPT 可直接插入的图片：PDB 结构图、开放获取（OA）文献图页、项目结果图与自制示意图。完整清单见 [`manifest.json`](manifest.json)。

重新生成：

```bash
python3 scripts/build_ppt_assets.py
```

（需 `matplotlib`；从 Nature Communications PDF 提取页面时需 `pymupdf`。源 PDF 在本地缓存，默认不纳入 git。）

---

## 建议幻灯片映射

| 汇报环节 | 推荐文件 | 说明 |
|----------|----------|------|
| **背景：痛风双轴** | `schematics/schematic_gout_dual_pathway.png` | 代谢轴（URAT1）+ 炎症轴（NLRP3） |
| **背景：NLRP3 与痛风** | `literature/frontiers2023_liu_nlrp3_gout_g001–g003.jpg` | Liu et al. *Front Immunol* 2023 |
| **背景：MSU–NLRP3 机制** | `literature/chen2023_jir_page1–3.png` | Chen et al. *J Inflamm Res* 2023（CC BY） |
| **动机：Eurycoma 对照** | `schematics/schematic_inspiration_vs_this_project.png` | 与 Zhang et al. *Nat Commun* 2025 互补定位 |
| **动机/先例** | `literature/eurycoma2025_pdf_page2–5.png` | 表型筛选 → 化合物 32 双靶点验证 |
| **靶点结构：URAT1** | `structures/9dkb_assembly-1.jpeg` | Fedor et al. 2025 共晶（PDB 9DKB） |
| **靶点结构：NLRP3** | `structures/7alv_assembly-1.jpeg`, `8etr_assembly-1.jpeg` | 7ALV（MCC950 类）、8ETR（GDC-2394 类） |
| **URAT1 药理（文献）** | `literature/fedor2025_urat1_pdf_page1–4.png` | lesinurad 等抑制剂结合模式 |
| **方法流程** | `schematics/schematic_project_workflow.png` | 8319 → NLRP3 ML 缩库 → **gnina P2** 双靶对接 → Pareto 审计 → 化学提名 |
| **数据不对称（SI）** | `project_results/si_data_asymmetry.png` | ChEMBL 临床库 vs 8973 |
| **NLRP3 ML（Fig 2）** | `project_results/fig02_nlrp3_screening_composite.png` | OOF 与筛选漏斗 |
| **8973 回顾（Fig 3）** | `project_results/fig03_urat1_retrospective_composite.png` | URAT1 8973 SI（非生产 P2） |
| **Pareto 对接（Fig 4）** | `project_results/fig04_pareto_dual_docking_9dkb_7alv.png` | 历史双靶表的裸 Pareto **审计**（EGCG 标为 PAINS 降级；不是 P2 提名） |
| **模型验证（SI）** | `project_results/si_nlrp3_oof_roc_pr.png` | NLRP3 ROC/PR |

---

## 目录结构

```
figures/ppt_assets/
├── manifest.json          # 机器可读清单（来源、许可、用途）
├── structures/            # RCSB PDB 组装图（JPEG）
├── literature/            # OA 文献图（JPG/PNG；PDF 本地缓存）
├── schematics/            # 项目自制示意图（PNG）
└── project_results/       # 从 figures/generated/ 复制的项目图
```

---

## 版权与引用（插入 PPT 时请标注）

| 来源 | 许可 | 引用 |
|------|------|------|
| RCSB PDB 结构图 | CC0（结构图）；需引用 PDB ID | 如 9DKB: Fedor et al. *Nat Commun* 2025 |
| Liu et al. Front Immunol 2023 | CC BY 4.0 | doi:10.3389/fimmu.2023.1137822 |
| Chen et al. J Inflamm Res 2023 | CC BY 4.0 | doi:10.2147/JIR.S413477 |
| Zhang et al. Nat Commun 2025 (Eurycoma) | CC BY 4.0 | doi:10.1038/s41467-025-62645-6 |
| Fedor et al. Nat Commun 2025 (URAT1) | CC BY 4.0 | doi:10.1038/s41467-025-60480-3 |
| `schematics/`、`project_results/` | 本项目（MIT） | URAT1/NLRP3 dual-target AIDD project |

PDF 页面提取自 OA 全文，仅用于学术汇报；正式发表请使用出版社原图。

---

## 使用提示

1. **优先用 PNG/JPG**：`literature/*_pdf_page*.png` 为 2.2× 渲染，适合投影；不必在 PPT 里嵌整份 PDF。
2. **三张示意图**可改字后重跑 `build_ppt_assets.py`，或直接在 PPT 里叠加文字框。
3. **项目主图**仍以 `figures/generated/main/` 为准；此处为汇报副本，便于一次性打包。
