# SLC22 辅助库（OAT 迁移 + OCT 脱靶）

完整科学逻辑见 [`docs/SLC22_AUXILIARY_RATIONALE.md`](../../docs/SLC22_AUXILIARY_RATIONALE.md)。

**要点**：URAT1（SLC22A12）是 **OAT 亚家族**有机阴离子/尿酸交换体，**不是 OCT**。  
- **Tier A（主）**：OAT1 + OAT3 → ML 迁移预训练  
- **Tier B（辅）**：OCT1 + OCT2 → 阳离子转运 **脱靶** 讨论（Tier 3 对接假说）

配置：`config/targets.yaml` → `auxiliary_targets`

---

## Tier A — OAT 迁移预训练（优先导出）

| 靶点 | 基因 | UniProt | ChEMBL ID | **当前清洗后规模** | 建议规模 | 活性类型 |
|------|------|---------|-----------|-------------------|----------|----------|
| **OAT1** | SLC22A6 | O95742 | **CHEMBL1641347** | **48** SMILES（2026-06 用户导出） | 500–2000 条 | 摄取/抑制 IC50、Ki |
| **OAT3** | SLC22A8 | O95816 | **CHEMBL1641348** | **23** SMILES | 500–2000 条 | 同上 |
| **OAT 合并** | — | — | — | **55** unique（`oat_combined_transfer.csv`） | — | — |

> **重要**：当前仅为 **IC50 子集**，清洗后远小于论文建议规模；迁移学习前请 **重新导出**（含 Ki/EC50、Assay Type B 优先）。  
> 运行：`python3 scripts/00b_prepare_auxiliary_data.py --copy-raw` → 见 `auxiliary_data_summary.json`

- OAT1: https://www.ebi.ac.uk/chembl/explore/target/CHEMBL1641347  
- OAT3: https://www.ebi.ac.uk/chembl/explore/target/CHEMBL1641348  

输出：

- `data/auxiliary/oat1_chembl_curated.csv`
- `data/auxiliary/oat3_chembl_curated.csv`
- `data/auxiliary/oat_combined_transfer.csv`

---

## Tier B — OCT 脱靶（非主迁移）

| 靶点 | 基因 | UniProt | ChEMBL ID（配置） | **导出常见 ID** | **当前清洗后** | 用途 |
|------|------|---------|-------------------|-----------------|---------------|------|
| OCT1 | SLC22A1 | O15245 | CHEMBL2073664 | **CHEMBL5685**（同基因别名） | **108** SMILES | 肝摄取脱靶 |
| OCT2 | SLC22A2 | O15244 | CHEMBL1770032 | **CHEMBL1743122**（同基因别名） | **105** SMILES | 肾阳离子分泌脱靶 |

输出：

- `data/auxiliary/oct1_chembl_curated.csv`
- `data/auxiliary/oct2_chembl_curated.csv`

**勿用错误 ID**：`CHEMBL242`、`CHEMBL3968`、`CHEMBL1906`、`CHEMBL210`（见 `DATA_FACT_CHECK.md`）。

---

## 清洗规则

与 URAT1 主任务相同（`config/targets.yaml` → `data_curation`）：

1. `standard_type` ∈ {IC50, Ki, EC50}，`standard_relation` = `=`
2. 关键词：`uptake`, `transport`, `inhibition`
3. `pactivity_range` 4.0–10.0；冲突按 median / discard 规则

---

## 与 URAT1 主任务关系

| 成员 | ChEMBL | 关系 |
|------|--------|------|
| URAT1 | CHEMBL6120 | 主任务（822 curated SMILES） |
| OAT1/3 | 见上 | **主迁移**（阴离子化学空间） |
| OCT1/2 | 见上 | **脱靶**（阳离子；抑制剂 ≠ URAT1 抑制剂） |
