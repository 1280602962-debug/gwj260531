# 开源分子对接流程（AutoDock Vina）

替代 Schrödinger Glide XP / Maestro / Canvas，保证无商业许可亦可完整复现。

---

## 工具栈

| 步骤 | 工具 | 说明 |
|------|------|------|
| 受体 | gemmi + Open Babel | 去配体/水，加氢，导出 PDBQT |
| 配体 | RDKit + Meeko | 3D 嵌入 + PDBQT |
| 对接 | **AutoDock Vina 1.2.5** | `exhaustiveness=32`（约等于 Glide SP 级搜索） |
| 可选精修 | smina | 同一搜索盒内局部优化 |
| 后处理 | `normalize_docking_export.py` | 统一 `dock_score` 列 |

配置文件：`config/docking_open_source.yaml`

---

## 一键流程（P≥0.5 池，1588 化合物）

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project

# 0. 依赖（conda 或 pip）
pip install -r requirements.txt
# Vina 二进制：见 https://github.com/ccsb-scripps/AutoDock-Vina/releases

# 1. 下载结构（若未包含）
mkdir -p data/structures/pdb
curl -fsSL -o data/structures/pdb/9DKB.cif https://files.rcsb.org/download/9DKB.cif
curl -fsSL -o data/structures/pdb/7ALV.pdb https://files.rcsb.org/download/7ALV.pdb

# 2. 受体 PDBQT
python3 scripts/prepare_receptor_vina.py --target urat1_9dkb
python3 scripts/prepare_receptor_vina.py --target nlrp3_7alv

# 3. 配体库
python3 scripts/prepare_ligands_vina.py \
  --input data/repurposing/screening/docking_pool_p05.csv \
  --output-dir results/repurposing/ligands_p05

# 4. 批量对接（并行）
python3 scripts/run_vina_batch.py \
  --target urat1_9dkb \
  --manifest results/repurposing/ligands_p05/ligand_manifest.csv \
  --output-dir results/repurposing/docking_vina/9dkb \
  --jobs 8

python3 scripts/run_vina_batch.py \
  --target nlrp3_7alv \
  --manifest results/repurposing/ligands_p05/ligand_manifest.csv \
  --output-dir results/repurposing/docking_vina/7alv \
  --jobs 8

# 5. 规范化 → Pareto
python3 scripts/normalize_docking_export.py \
  --input results/repurposing/docking_vina/9dkb/docking_9dkb_vina.csv \
  --pdb 9DKB \
  --output results/repurposing/docking_raw/urat1_9dkb_p05.csv \
  --engine vina

python3 scripts/normalize_docking_export.py \
  --input results/repurposing/docking_vina/7alv/docking_7alv_vina.csv \
  --pdb 7ALV \
  --output results/repurposing/docking_raw/nlrp3_7alv_p05.csv \
  --engine vina

python3 scripts/merge_docking_pareto.py \
  --ml-scores data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv \
  --urat1-dock results/repurposing/docking_raw/urat1_9dkb_p05.csv \
  --nlrp3-dock results/repurposing/docking_raw/nlrp3_7alv_p05.csv \
  --pool data/repurposing/screening/docking_pool_p05.csv
```

---

## 8973 URAT1 回顾轨（8928 化合物）

对 `data/distill/distill_manifest.csv` 重复步骤 3–4（配体目录可单独建 `results/docking/ligands_8973`），再：

```bash
python3 scripts/merge_8973_docking_results.py \
  --glide-csv results/docking/vina_8973/docking_9dkb_vina.csv \
  --pdb 9DKB
python3 scripts/analyze_urat1_docking_vs_ml.py
```

> CLI 参数名 `--glide-csv` 保留向后兼容，接受 Vina 导出。

---

## 搜索盒（与 Glide 网格对齐）

| 靶点 | PDB | 中心 (Å) | 盒边长 (Å) | 参考配体 |
|------|-----|----------|------------|----------|
| URAT1 | 9DKB | 99.97, 102.97, 105.70 | 22³ | lesinurad (A1AIL) |
| NLRP3 | 7ALV | 16.76, 35.45, 125.71 | 20³ | RM5 (MCC950-class) |

---

## 与 Glide 的等价性说明（Methods 用语）

- **搜索深度**：Vina `exhaustiveness=32` 为文献常用的「高精度」设置，计算量与 Glide SP 同量级；**不等同** Glide XP 的额外经验势与 WaterMap。
- **评分**：使用 Vina `minimizedAffinity`（kcal/mol）；**仅用于池内百分位**，不声称与实验亲和力或 Glide 分数可换算。
- **可重复性**：全开源栈 + GitHub 脚本；满足 JMM 对计算可复现的要求。

---

## 迁移状态

| 数据集 | 历史引擎 | 目标引擎 | 状态 |
|--------|----------|----------|------|
| P05 双靶 (1588) | Glide XP | Vina 1.2.5 | **待重跑**（流水线已就绪） |
| 8973 URAT1 | Glide XP | Vina 1.2.5 | **待重跑** |
| Pareto / 图表 | Glide 分数 | Vina 百分位 | 重跑后更新 |

`data/repurposing/pareto/` 中现有 CSV 为 **Glide 时代结果**，投稿前须用 Vina 重算并替换。

---

## 烟雾测试（3 个 benchmark 配体）

```bash
python3 scripts/prepare_ligands_vina.py \
  --input config/docking_open_source.yaml \
  ...
# 或 --limit 3 于 run_vina_batch.py
python3 scripts/run_vina_batch.py --target urat1_9dkb --manifest ... --limit 3
```
