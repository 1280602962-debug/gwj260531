# P0 完成指南 — Zenodo 成品 + 英文稿（投 JCIM 前必做）

> 目的：补齐相对 Vu et al. JCIM 2025 的「成品感」缺口。  
> **不需要新对接。**  
> 仓库已清理为 DualFourClass 投稿核心；入口见 [`../README.md`](../README.md)。  
> **写作前逐项核对：** [`JCIM_PREWRITING_CHECKLIST_V1.md`](JCIM_PREWRITING_CHECKLIST_V1.md)

---

## P0 是什么（两件事）

| 编号 | 内容 | 谁做 | 完成标志 |
|------|------|------|----------|
| **P0-A** | Zenodo（或同等）公开数据包 + DOI + 一键复现 | 本地打包上传 + 云端可写说明 | 有 DOI；外人按 README 能重算主表 |
| **P0-B** | 英文 Article 初稿 + 主图定稿 | **正文与 Fig1–7 已完成；待期刊排版** | `MANUSCRIPT_JCIM_EN.md` + `figures/jcim_article/` |

---

## P0-A — Zenodo 包（具体怎么做）

### A1. 打包目录（建议本地组好再上传）

在本地（有姿态的机器）建：

```text
DualFourClass-Bench_v1/
├── README.md                 # 英文：任务定义、主指标、复现命令
├── LICENSE                   # 建议 CC-BY 4.0（数据）+ MIT（脚本）
├── CITATION.cff
├── CLAIM_CEILING.md          # 从仓库复制
├── PROTOCOL_BENCH.md         # 口袋匹配定义 + 标签规则 + prep
├── environment.yml           # 或 ENV_PIN.md 内容
├── panels/                   # 每对：panel CSV + SMILES
├── receptors_boxes/          # PDBQT + box JSON（无大姿态也要）
├── scores/                   # ablation / vina_long / rtm / gnina 表
├── analysis/                 # jcim_bench_v0 + jcim_strengthen 主表与图
├── scripts/                  # 复现主表的最少脚本
└── poses_top1/               # 可选但强烈建议：每配体每受体 top1
```

姿态来源清单见：  
`data/jcim_strengthen_t0t1_v0/POSE_UPLOAD_CHECKLIST.md`

### A2. 最少必须进包（没有姿态也能先发「分数版」）

1. `assembled_all_pairs_long.csv`（或四对 assembled）  
2. `PRIMARY_METRIC_V2` 相关表 + `pocket_matched_*`  
3. `jcim_strengthen_t0t1_v0/tables/*` 主表  
4. K=4 `ablation_ligand_scores.csv` + panel CSV  
5. `protocols/PAIR_ROLES_APPROVED_JCIM.yaml` + `CLAIM_CEILING.md`  
6. `build_pocket_matched_diagnostics_v1.py` + `build_t0_strengthen_v1.py`  
7. `ENV_PIN.md`

**有姿态再补 top1**（审稿加分）；无姿态也可先上传分数包拿 DOI，再发 v1.1 补姿态。

### A3. 上传步骤（白话）

1. 注册 [Zenodo](https://zenodo.org)，登录后 New Upload  
2. 上传 zip；Title 例：`DualFourClass-Bench: panels, scores, and analysis tables for dual-target docking evaluation`  
3. Upload type: Dataset；License: CC-BY-4.0  
4. 关联 GitHub 本仓库 URL  
5. Publish → 得到 **DOI**  
6. 把 DOI 写回仓库 `Dual_Target_Docking/README.md` 与英文稿 Data Availability

### A4. 复现命令（写进 Zenodo README）

```bash
# 最小：重算口袋匹配主表 + 混淆包
python3 data/jcim_bench_v0/scripts/build_pocket_matched_diagnostics_v1.py
python3 data/jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py
python3 data/jcim_bench_v0/scripts/plot_forest_ci_v1.py
```

---

## P0-B — 英文稿 + 主图（已完成科学组装）

### B1. 当前主稿

- 英文投稿主稿：[`MANUSCRIPT_JCIM_EN.md`](MANUSCRIPT_JCIM_EN.md)
- 中文核对稿：[`MANUSCRIPT_JCIM_ZH.md`](MANUSCRIPT_JCIM_ZH.md)
- 重建命令：`python3 docs/assemble_manuscript_en.py` / `assemble_manuscript_zh.py`
- 题名：*Benchmark Formulation and Chemical Confounding in Docking-Based Dual-Target Recognition*

### B2. 当前主线

- formulation：Dual-versus-neither 与方向性 selective-hard-negative 任务并不等价；EGFR/HER2 为清晰案例，不是四对定律；
- confounder-aware evaluation：docking 必须相对 ligand-only、描述符与 wrong-pocket 对照解释；
- evaluation-condition sensitivity：max/median 对主结论影响小，unused-pool 与 receptor realization 显示条件依赖；
- K=4 是供给受限案例面板，不是 comprehensive suite。

### B3. 当前主图

| Fig | 内容 |
|-----|------|
| 1 | 四状态任务 + 两条口袋匹配方向主终点 |
| 2 | 49→4 供给审计 |
| 3 | Dual-versus-neither vs directional `summary_min` |
| 4 | 弱臂与物化混淆 |
| 5 | 两靶对 receptor realization（PM 含 4JSX） |
| 6 | wrong-pocket 在 holdout 上的未解决反转 |
| 7 | ECFP4 / 全描述符 / 协变量 / 匹配子集 |

原主森林图与 unused-pool 对照分别为 Figure S4 / S5。所有图由 `plot_jcim_article_figures_v1.py` 从冻结 CSV 重绘并校验。

### B4. 尚未完成

- 按 ACS 模板排版 Word/PDF；
- Zenodo DOI 写回 Data and Software Availability；
- cover letter、作者信息、利益冲突与投稿元数据。

---

## P0 完成检查清单

- [ ] Zenodo 数据包已 publish，DOI 写入 README  
- [ ] 最少复现脚本在干净环境跑通主表  
- [x] 英文稿完整（Abstract→Conclusions + References）
- [x] Fig1–7 + SI 图由冻结 CSV 生成并通过脚本校验
- [x] Limitations 写清：K=4、ChEMBL/assay、receptor、无 prospective validation
- [x] 全文无「通吃 scorer / 显著通吃四对」表述

**做完 P0 = 可以按 JCIM Evaluation Article 投稿形态组装；不必再为 P0 开对接。**
