# 双靶重定位药物库选择与获取指南

> **用途**：路线 2 — 在 **合理共享库** 上做 URAT1（9DKB 对接）+ NLRP3（ML）双证据重定位  
> **不是**：8973 URAT1 蒸馏集；不是百万随机库  
> **目标规模**：去盐、去肽后 **约 1,500–2,500** 个小分子

---

## 1. 库应满足什么条件

| 条件 | 原因 |
|------|------|
| **已上市 / 临床阶段小分子** | 重定位叙事；类药性合理 |
| **有可靠 SMILES** | Glide + RDKit 去重 |
| **与训练化学空间部分重叠** | NLRP3 ML 不在适用域外乱预测 |
| **规模 1.5k–3k** | 全库 9DKB XP 可算；NLRP3 ML 秒级 |
| **可引用、可复现** | Methods 写清来源与过滤 |

---

## 2. 推荐数据源（按优先级）

### 首选：ChEMBL（免费、项目已用）

| 项目 | 说明 |
|------|------|
| 入口 | https://www.ebi.ac.uk/chembl/ → Downloads；FTP `ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/` |
| 筛选 | `max_phase` 3=临床，4=上市；`molecule_type = Small molecule` |
| 规模 | 过滤后约 **2,000–2,800** |
| Methods 表述 | *ChEMBL clinical-phase and approved small molecules (max_phase ≥ 3)* |

### 备选：DrugBank（需学术注册）

| 项目 | 说明 |
|------|------|
| 入口 | https://go.drugbank.com/releases/latest |
| 文件 | `approved` + 可选 `investigational` 的 Structure External Links CSV |
| 规模 | Approved ~1.5–2k；+ investigational ~3k |
| 注意 | 不可再分发；Methods 引用 DrugBank 版本 |

### 不推荐

| 来源 | 原因 |
|------|------|
| 8973 distill | URAT1 偏置 |
| Enamine 百万库 | 算力 + ML 适用域 + 无实验 |
| 随机 PubChem | 非药物，重定位不成立 |

---

## 3. 构建四步（ChEMBL）

1. **下载** ChEMBL molecule 表（SQLite 或 CSV）  
2. **过滤**：max_phase≥3，small molecule，有 SMILES；MW 150–800；去盐；InChIKey 去重  
3. **追加 benchmark**：lesinurad 等六药 + allopurinol/colchicine（见 `literature_benchmarks_summary.csv`）  
4. **输出**：`data/repurposing/repurposing_manifest.csv`（目标 ~2k）

---

## 4. 与 8000 XP 的分工

| 数据 | 用途 |
|------|------|
| 822+8000 @ 9DKB | URAT1 对接 **回顾性验证**（ML vs 结构） |
| ~2k 药物库 | **双靶重定位**（URAT1 XP + NLRP3 ML，同一批药） |

---

## 5. 待建脚本

`scripts/build_repurposing_library.py` — 从 ChEMBL SQLite 筛 max_phase≥3 并输出 manifest。
