# 数据与文献事实核验清单

> 本文档记录 **可复现、已核对** 的项目事实。所有数字应来自 `data/processed/data_summary.json` 或 `scripts/00_prepare_data.py` 输出；文献 ID 以 ChEMBL / PubMed / RCSB 为准。  
> 更新：`python3 scripts/00_prepare_data.py` 后同步本节 NLRP3 统计。

---

## 1. 训练集规模（`data/processed/data_summary.json`）

| 靶点 | 化合物 (unique SMILES) | 记录数 | 其他 |
|------|------------------------|--------|------|
| URAT1 | **822** | 822（每 SMILES 一条，median 聚合） | 218 Murcko scaffolds |
| NLRP3 | **513** | **609**（assay-conditioned 记录） | **39** assays；THP-1 子集 **302** SMILES / **313** records |
| 重叠 | **0** shared SMILES | — | 独立双模型依据 |

### NLRP3 assay 异质性（实测，非估计）

| 指标 | 值 | 定义 |
|------|-----|------|
| 多 assay 化合物 | **80 / 513** (15.6%) | 同一 canonical SMILES 出现在 ≥2 个 assay |
| 跨 assay >1 log 离散 | **37 / 513** (**7.2%**) | 同一 SMILES 在不同 assay 的 pActivity 极差 >1.0 |

**勿用** 旧文档中的「503 化合物」「92 assays」「47% 冲突」——与当前 curated 数据不符。

---

## 2. ChEMBL Target ID（`config/targets.yaml`）

| 靶点 | 基因 | UniProt | ChEMBL ID |
|------|------|---------|-----------|
| URAT1 | SLC22A12 | Q96S37 | CHEMBL6120 |
| NLRP3 | NLRP3 | Q96P20 | CHEMBL1741208 |
| OCT1 | SLC22A1 | O15245 | CHEMBL2073664 |
| OCT2 | SLC22A2 | O15244 | CHEMBL1770032 |
| OAT1 | SLC22A6 | O95742 | CHEMBL1641347 |
| OAT3 | SLC22A8 | O95816 | CHEMBL1641348 |

### OAT 辅助库清洗后规模（`auxiliary_data_summary.json`，全量 ChEMBL 导出）

| 库 | 原始行数 | 清洗后 SMILES | 与 URAT1 重叠 |
|----|----------|---------------|---------------|
| OAT1 | 280 | **63** | 13 |
| OAT3 | 254 | **41** | 5 |
| OAT 合并 | — | **73** | 13 |

清洗规则：仅 IC50/Ki/EC50；`=`；pActivity 4–10；冲突丢弃；缺 pChEMBL 时由 nM/µM 换算。

### ChEMBL 靶点 ID 别名（OCT 导出常见）

| 配置 ID | 导出中常见 ID | 基因 | 说明 |
|---------|---------------|------|------|
| CHEMBL2073664 | **CHEMBL5685** | SLC22A1 (OCT1) | 同一蛋白，ChEMBL 历史条目 |
| CHEMBL1770032 | **CHEMBL1743122** | SLC22A2 (OCT2) | 同一蛋白，ChEMBL 历史条目 |

### 常见错误 ID（禁止用于导出）

| 错误 ID | 实际靶点 |
|---------|----------|
| CHEMBL242 | ESR2 |
| CHEMBL3968 | PLA2G1B |
| CHEMBL1906 | RAF1 |
| CHEMBL210 | ADRB2 |
| CHEMBL3989876 | 非 MCC950（MCC950 用 **CHEMBL230208**） |
| CHEMBL1777665 | **大鼠** OAT1（人源 OAT1 用 **CHEMBL1641347**） |

辅助库分层见 `docs/SLC22_AUXILIARY_RATIONALE.md`。

## 3. Benchmark 化合物 ChEMBL

| 化合物 | ChEMBL ID |
|--------|-----------|
| lesinurad | CHEMBL3301572 |
| benzbromarone | CHEMBL892 |
| verinurad | CHEMBL3989871 |
| dotinurad | CHEMBL4594374 |
| MCC950 | CHEMBL230208 |
| allopurinol | CHEMBL1467 |

---

## 4. PDB 与配体对应（结构考试）

| PDB | 配体 / 说明 |
|-----|-------------|
| 9B1H | lesinurad inward（Dai 2024, *Cell Res*） |
| 9DKB | lesinurad inward（Fedor/Suo 2025, *Nat Commun*）— **三态对接 inward 主 grid** |
| 9B1K | urate **occluded**（Dai 2024）— **三态对接 occluded** |
| 9B1L | urate **outward-facing**（Dai 2024）— **三态对接 outward** |
| 9B1J | urate inward-facing（Dai 2024） |
| 9DKA | **benzbromarone**（勿与 9DKB 混淆） |
| 9JDZ | lesinurad inward（Wu 2025, *Cell Discov*）— **非** occluded/outward |
| 9JDY / 9JE1 | verinurad / dotinurad |
| 7ALV | **MCC950 类类似物 NP3-146**（非 MCC950 共晶；药效团模板） |
| 8ETR | GDC-2394（McBride 2022, *J Med Chem*） |

---

## 5. 关键文献 PMID（benchmark 主来源）

| 主题 | PMID | DOI |
|------|------|-----|
| Dai URAT1 2024 | **39245778** | 10.1038/s41422-024-01023-1 |
| Fedor/Suo URAT1 2025 | 40467597 | 10.1038/s41467-025-60480-3 |
| Burns URAT1 IC50 2016 | 27716403 | 10.1186/s13075-016-1107-x |
| verinurad 25 nM | 28386072 | 10.1038/s41598-017-00706-7（Tan PK et al.） |
| Nakamura dotinurad 2019 | 31371478 | 10.1124/jpet.119.262741 |
| MCC950 IL-1β | 25686105 | 10.1038/nm.3806 |
| GDC-2394 | 36279149 | 10.1021/acs.jmedchem.2c01250 |
| 7ALV 结构 | 34687713 | 10.1016/j.jmb.2021.167189 |

---

## 6. 维护规则

1. 改 ChEMBL 导出或清洗规则后 **必须** 重跑 `00_prepare_data.py` 并更新本节表 1。  
2. 新增 benchmark 行须填 **真实 PMID**，不得用 ChEMBL compound ID / IUPHAR ligand ID 冒充。  
3. 论文 Methods 引用规模数字时写：`data_summary.json` 生成日期 + 清洗过滤器字符串。  
4. 结构 redock 使用 7ALV 时须写明 **analog-based template**，不得称 MCC950 共晶。

---

## 7. 实现状态与审计记录（2026-06-29）

| 类别 | 状态 |
|------|------|
| 数据清洗 822/513/39 assays/7.2% 冲突 | ✅ `00_prepare_data.py` 可复现 |
| NLRP3 THP-1 子集 | **302** unique SMILES，**313** records | `data_summary.json` → `nlrp3.thp1_*`（`assay_cell_type` 含 THP） |
| 双模型训练 + benchmark | ✅ `run_model_build_and_validate.py` |
| ChEMBL/PDB/PMID 黑名单 | ✅ 见 §2–§5 |
| $S_{\text{trap}}$、Path A/B、PLK1 消融 | ☐ 骨架脚本，**不可写进 Results** |
| OAT 辅助库 CSV | ☐ 待 ChEMBL 导出 |
| MASFL v3.1 全管线 | ☐ 设计稿 |

**维护**：每次改清洗规则或 benchmark 行后更新本节日期并重跑 `00_prepare_data.py`。
