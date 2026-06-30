# 项目准备清单（TAPE-GATE v2.0）

> 按优先级排列。打勾项为发表 **最低要求**；加星项为 **冲高期刊建议**。

**v2.0 变更**：主路径为不对称双证据 + 库筛/生成双路径；MTL 降为消融；须准备 PLK1-style baseline 对照。

---

## 一、数据准备

### 1.1 生物活性数据 ✅

| 任务 | 来源 | 操作 | 状态 |
|------|------|------|------|
| URAT1 ChEMBL | CHEMBL6120 / 用户 CSV | 清洗后 **822** SMILES | ☐ |
| NLRP3 ChEMBL | CHEMBL1741208 / 用户 CSV | IL-1β + Assay B，**513** SMILES（609 records） | ☐ |
| NLRP3 assay 元数据 | ChEMBL 导出列 | 保留 assay_id, cell_line, assay_type | ☐ |
| 重叠检查 | — | 确认 0 shared SMILES → 独立双模型 | ☐ |
| SLC22 辅助 | OAT1 `CHEMBL1641347`、OAT3 `CHEMBL1641348`（**主迁移**）；OCT1/2（**脱靶**） | 见 `SLC22_AUXILIARY_RATIONALE.md` | ☐ |
| 骨架划分 | Murcko + GroupKFold | train/val/test CSV | ☐ |

**预期输出**：`data/processed/urat1_curated.csv`, `nlrp3_records.csv`（含 assay 列）, `splits/`, `data_summary.json`

> 所有规模与 ChEMBL ID 以 [`docs/DATA_FACT_CHECK.md`](DATA_FACT_CHECK.md) 为准；改清洗规则后须重跑 `scripts/00_prepare_data.py`。
### 1.2 结构数据 ✅

| 任务 | PDB ID | 操作 |
|------|--------|------|
| URAT1 inward（主 grid） | **9DKB**（lesinurad, 2.55 Å） | RCSB 下载 → Prep Wizard → Grid |
| URAT1 occluded | **9B1K**（urate, occluded） | 独立 PDB；勿用 9JDZ |
| URAT1 outward | **9B1L**（urate, outward-facing） | 独立 PDB；勿用 9JDZ |
| 备用 inward / benchmark | 9B1H, 9DKA, 9JDY, 9JE1 | redock 与敏感性分析 |
| NLRP3 NACHT + inhibitor | 7ALV | 提取 NACHT 域或全链 |
| NLRP3 + GDC-2394 | 8ETR | cryo-EM 低分辨率，作补充 |
| 膜蛋白准备 | — | CHARMM-GUI 或 OPM 取向，嵌入 POPC 脂双层（MD 用）|

**工具**：PDBFixer, Modeller（若有缺失 loop）, PyMOL/ChimeraX 可视化

### 1.3 筛选库

| 库 | 规模 | 用途 |
|----|------|------|
| Enamine REAL | ~10^6 | 主筛选 |
| ChEMBL 全库 | ~2×10^6 | 回顾性/补充 |
| 自建 dual-benchmark | 见 `data/benchmarks/` | 强制回收测试 |

---

## 二、软件环境

### 2.1 核心依赖 ✅

```bash
# Python 3.10+
pip install rdkit-pypi pandas numpy scikit-learn xgboost pyyaml
pip install chemprop  # v2.x
pip install shap matplotlib seaborn umap-learn
```

### 2.2 分子基础模型

| 包 | 安装 | 用途 |
|----|------|------|
| graphium / MiniMol | `pip install graphium` 或 GitHub | 预训练指纹 |
| 备选：chemprop pretrained | Chemprop checkpoint | baseline |

### 2.3 对接与 MD

| 软件 | 许可 | 用途 |
|------|------|------|
| AutoDock Vina | 免费 | 快速系综对接 |
| GNINA | 免费 | CNN 重打分（可选）|
| Schrödinger Glide | 商业（学校许可）| SP/XP + MM-GBSA |
| GROMACS 2024 | 免费 | MD 模拟 |
| AMBER/OpenMM | 免费 | 备选 MD 引擎 |

### 2.4 生成式模块 ⭐

| 工具 | 用途 |
|------|------|
| HuggingFace Transformers | 化学 LLM |
| MolGPT / ChemGPT 权重 | CLM 预训练基座 |
| REINFORCE 自实现 | RL 微调 |

### 2.5 ADMET 预测

| 工具 | 网址 |
|------|------|
| SwissADME | http://www.swissadme.ch |
| pkCSM | https://biosig.lab.uq.edu.au/pkcsm |
| ADMETlab 2.0 | https://admetmesh.scbdd.com |

---

## 三、文献与 Benchmark 准备 ✅

### 3.1 必须精读的论文（Methods 引用）

1. Dai et al., 2024 — URAT1 cryo-EM (*Cell Res*)
2. Fedor et al., 2025 — URAT1 药物结合 (*Nat Commun*)
3. Dekker A et al., 2021 — NLRP3 NACHT + inhibitor analog（PDB 7ALV；*J Mol Biol*）
4. Beaini et al., 2024 — Graphium/MiniMol (ICLR)
5. Zhao et al., 2024 — NLRP3 ML 筛选流程 (*BMC Chemistry*)
6. Ferreira et al., 2024 — POLYGON (*Nat Commun*)
7. Schneider et al., 2024 — CLM dual-target (*Nat Commun*)

### 3.2 Benchmark 化合物（必须准备 SMILES）

从 ChEMBL 或 PubChem 获取精确结构：

- [ ] lesinurad (CHEMBL3301572)
- [ ] benzbromarone (CHEMBL892)
- [ ] verinurad
- [ ] dotinurad
- [ ] MCC950
- [ ] GDC-2394
- [ ] 阴性对照：allopurinol, colchicine

---

## 四、算力与存储

详见 [`COMPUTE_REQUIREMENTS.md`](COMPUTE_REQUIREMENTS.md)。

**最低配置**：1× RTX 3090/4090 (24GB) + 32GB RAM + 200GB SSD

---

## 五、论文素材准备

### 5.1 图表清单（见 `MANUSCRIPT_OUTLINE.md`）

- [ ] Fig 1: STAD-AIDD 流程图
- [ ] Fig 2: 数据集化学空间 UMAP
- [ ] Fig 3: MTL 模型 CV 性能 + 消融
- [ ] Fig 4: URAT1 构象系综对接模式（与共晶对比）
- [ ] Fig 5: NLRP3 对接 + MD
- [ ] Fig 6: 筛选漏斗 + benchmark 回收
- [ ] Fig 7: 生成候选结构与评分分布
- [ ] Table 1: 与相关工作对比
- [ ] Table 2: Top 20 候选化合物性质
- [ ] SI: 完整候选列表、对接参数、超参

### 5.2 开源准备

- [ ] GitHub 仓库（本项目）
- [ ] Zenodo 归档（数据 + 结果快照）
- [ ] README 含复现命令

---

## 六、时间线建议（按技术难度，非日历）

| 阶段 | 核心交付物 |
|------|-----------|
| 数据 + 结构准备 | 清洗 CSV、PDB 就绪 |
| 数据集分析 | UMAP、骨架统计、Fig 2 |
| MTL 训练 + CV | 模型权重、消融、Fig 3 |
| 系综对接 + 筛选 | Top 500、Fig 4–6 |
| 生成式优化（可选）| 新候选、Fig 7 |
| 论文撰写 | 初稿 |

---

## 七、合作建议（弥补无湿实验）

若需增强论文说服力，可考虑：

| 合作方向 | 可提供 |
|---------|--------|
| 药学院 URAT1 细胞摄取 lab | 验证 Top 5–10 化合物 |
| 免疫学 lab | NLRP3 IL-1β 抑制测定 |
| CRO | 合成 3–5 个生成候选 |

**即使无合作**：纯计算论文在 JCIM / J Cheminformatics 仍可发表，但需在 Limitations 明确说明。
