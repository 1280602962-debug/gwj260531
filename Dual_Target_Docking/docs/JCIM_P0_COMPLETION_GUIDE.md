# P0 完成指南 — Zenodo 成品 + 英文稿（投 JCIM 前必做）

> 目的：补齐相对 Vu et al. JCIM 2025 的「成品感」缺口。  
> **不需要新对接。**  
> 仓库已清理为 DualFourClass 投稿核心；入口见 [`../README.md`](../README.md)。

---

## P0 是什么（两件事）

| 编号 | 内容 | 谁做 | 完成标志 |
|------|------|------|----------|
| **P0-A** | Zenodo（或同等）公开数据包 + DOI + 一键复现 | 本地打包上传 + 云端可写说明 | 有 DOI；外人按 README 能重算主表 |
| **P0-B** | 英文 Article 初稿 + 主图定稿 | 写作（可人写 / agent 起稿） | 完整 IMRaD + Fig1–5 可投稿级 |

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

## P0-B — 英文稿 + 主图（具体怎么做）

### B1. 主张（开写前钉死，勿漂移）

- 评测/基准文，不是新 scorer  
- 主指标 = **pocket-matched directional AUROC**  
- K=4 冻结集；EGFR = 供给案例  
- 不写「通用决策臂已验证」

### B2. 建议章节 ↔ 已有材料

| 章节 | 用这些文件 |
|------|------------|
| Intro | 任务缺口；可引 Vu 2025 / VSDS 作单靶对照 |
| Methods 2.1 Dataset | `jcim_j0j1_v0` 供给；`public_pair_selection`；四对 panel CSV |
| Methods 2.2 Docking | 各 panel `protocol.yaml`；PM E=8→16 cognate |
| Methods 2.3 Metrics | `PRIMARY_METRIC_V2.md`；CLAIM_CEILING |
| Results 3.1–3.7（英文 JCIM 稿） | **[`RESULTS_SECTION_JCIM_EN_V1.md`](RESULTS_SECTION_JCIM_EN_V1.md)** |
| Results 3.1–3.7（中文对齐稿） | [`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) |
| Results 3.1 Supply | J0 表 |
| Results 3.2 Pooling vs directional | asymmetry / PRIMARY_METRIC_V2 |
| Results 3.3 Main forest + baselines | forest 图；baseline gate |
| Results 3.4 Confounds | wrong-pocket；matched subset；LE；covariate；scaffold ML |
| Results 3.5 Robustness | E8 vs E16；PM110；prep；θ；enrichment |
| Discussion / Conclusions | CLAIM_CEILING + Limitations |
| Data Availability | Zenodo DOI |

### B3. 主图最少 5 张（定稿级）

| Fig | 内容 | 现成素材 |
|-----|------|----------|
| 1 | 四类任务 + 口袋匹配定义 | 新画 schematic |
| 2 | 供给审计（49 对硬负） | `jcim_j0j1_v0` |
| 3 | K=4 森林图（口袋匹配 ± CI + 基线） | `figures/forest_*` 需改标注为 pocket-matched |
| 4 | 混淆解剖（错口袋 / LE / 匹配） | `pocket_specificity_gap` + matched subset |
| 5 | 稳健性（E8/E16、PM48/110、单靶 enrichment） | B 组 md/表 |

### B4. 写作执行方式（任选）

- **人写：** 按上表填空，先写 Results 数字段再写 Intro  
- **Agent 起稿：** 另开任务「按 `JCIM_P0_COMPLETION_GUIDE.md` 的 B2/B3 生成英文 Draft v0」  
- Cover letter 一句：*evaluation/benchmark article; not a new scoring function*

---

## P0 完成检查清单

- [ ] Zenodo 数据包已 publish，DOI 写入 README  
- [ ] 最少复现脚本在干净环境跑通主表  
- [ ] 英文稿完整（Abstract→Conclusions + Data Availability）  
- [ ] Fig1–5 定稿（矢量 PDF/PNG）  
- [ ] Limitations 写清：K 小、singleton 支架、ChEMBL max 聚合、无湿实验  
- [ ] 全文无「通吃 scorer / 显著通吃四对」表述  

**做完 P0 = 可以按 JCIM Evaluation Article 投稿形态组装；不必再为 P0 开对接。**
