# SLC22 家族辅助训练库（OCT1 / OCT2）

**用途**：TC-Encoder 迁移学习预训练；URAT1 选择性讨论（URAT1 vs OCT 脱靶）；对接比值 $R_{\text{sel}}$（Tier 3 计算假说，非实验选择性）。

配置见 `config/targets.yaml` → `auxiliary_targets.slc22_family`。

---

## 5. 靶点与 ChEMBL ID（已校正）

| 靶点 | 基因 | UniProt | ChEMBL ID | 建议规模 | 活性类型 |
|------|------|---------|-----------|----------|----------|
| OCT1 | SLC22A1 | O15245 | **CHEMBL2073664** | 500–2000 条 | 细胞摄取抑制 / 转运 IC50、Ki |
| OCT2 | SLC22A2 | O15244 | **CHEMBL1770032** | 500–2000 条 | 同上 |

ChEMBL 页面：

- OCT1: https://www.ebi.ac.uk/chembl/explore/target/CHEMBL2073664
- OCT2: https://www.ebi.ac.uk/chembl/explore/target/CHEMBL1770032

---

## 常见错误 ID（勿用）

| 错误 ID | 实际靶点 | 说明 |
|---------|----------|------|
| CHEMBL242 | Estrogen receptor β (ESR2) | 与 OCT1 无关 |
| CHEMBL3968 | Acidic phospholipase A2 (PLA2G1B) | 与 OCT2 无关 |
| CHEMBL1906 | RAF1 kinase | 旧版项目配置误用 |
| CHEMBL210 | Beta-2 adrenergic receptor (ADRB2) | 旧版项目配置误用 |

---

## 数据导出与清洗

1. 自 ChEMBL 按 `target_chembl_id` 导出 bioactivity（与 URAT1 相同 `cf12` 列格式）。
2. 保留 `standard_type` ∈ {IC50, Ki, EC50}，`standard_relation` = `=`。
3. Assay 关键词优先：`uptake`, `transport`, `inhibition`, `substrate`；排除纯结合（Kd）若无细胞功能上下文。
4. `pactivity_range`: 4.0–10.0；同一 SMILES 冲突按 `config/targets.yaml` → `data_curation` 规则聚合。
5. 输出路径（待生成）：
   - `data/auxiliary/oct1_chembl_curated.csv`
   - `data/auxiliary/oct2_chembl_curated.csv`
   - `data/auxiliary/slc22_combined_transfer.csv`（合并后用于序贯微调）

---

## 与 URAT1 的关系

| 成员 | ChEMBL | 迁移 / 讨论角色 |
|------|--------|-----------------|
| URAT1 | CHEMBL6120 | 主任务（822 条 curated） |
| OCT1 | CHEMBL2073664 | 家族表示预训练；肝摄取脱靶 |
| OCT2 | CHEMBL1770032 | 家族表示预训练；肾排泄脱靶 |

**注意**：SLC22 迁移是表示学习层面；OCT 抑制剂 ≠ URAT1 抑制剂。论文中须区分 ML 迁移与 $R_{\text{sel}}$ 计算假说。
