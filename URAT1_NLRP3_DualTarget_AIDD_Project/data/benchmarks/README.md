# Literature Benchmark 化合物库

> **文件**：`literature_benchmarks.csv`（多行/多文献来源）  
> **用途**：回顾性回收检验、结构 redocking、CPDL 方向测试、阴性对照  
> **非用途**：融合网络训练、等渗校准拟合（样本过少）

---

## 文件结构说明

`literature_benchmarks.csv` 采用 **长表格式**：同一化合物可有 **多行**，每行对应一条 **文献活性测定** 或 **参考结构**。

| 列名 | 含义 |
|------|------|
| `compound_id` | 项目内部 ID（如 `URAT1_POS_01`） |
| `validation_tier` | `Tier1a` 外推 / `Tier1b` 训练集内 / `Tier_neg` 阴性 |
| `validation_types` | 考试类型（分号分隔） |
| `ic50_nm` / `pactivity` | 该文献行的活性 |
| `ref_pmid` / `ref_doi` | 可追溯文献 |
| `pdb_id_primary` | 首选共晶结构 |
| `project_curation_note` | 与项目训练集/ChEMBL 的关系 |

---

## 化合物一览（主推荐标签）

### URAT1 阳性（必回收）

| 化合物 | 主 IC50（推荐） | 主文献 | PDB | Tier |
|--------|----------------|--------|-----|------|
| **lesinurad** | 3.53 µM (HEK-URAT1) 或 39 µM (Dai EM) | Burns 2016 / Dai 2024 | 9B1H, 9DKB | 1a |
| **benzbromarone** | 0.29 µM | Burns 2016 | 9DKA | 1a |
| **verinurad** | 25 nM | Tan et al. 2017 (*Sci Rep*) | 9JDY | 1b（在训练集） |
| **dotinurad** | 37.2 nM | Nakamura 2019 | 9JE1, 9B1G | 1a |

### NLRP3 阳性（必回收）

| 化合物 | 主 IC50（推荐） | 主文献 | PDB | Tier |
|--------|----------------|--------|-----|------|
| **MCC950** | 7.5 nM IL-1β (BMDM) | Coll 2015 | 7ALV（MCC950 类类似物共晶，非 MCC950 本身） | 1b |
| **GDC-2394** | 16 nM IL-1β (HMDM) | McBride 2022 | 8ETR | 1b |

### 阴性对照

| 化合物 | 为何选它 | 期望 |
|--------|---------|------|
| **allopurinol** | XO 抑制剂，降尿酸但非 URAT1 | URAT1 排名后 20% |
| **colchicine** | 微管抑制剂，间接影响 NLRP3 | NLRP3 排名后 20% |

---

## 关键文献（Methods 必引）

### URAT1 结构药理学

1. Dai Y, Lee CH. *Cell Res* 2024;34:776-787. doi:10.1038/s41422-024-01023-1 — PDB 9B1H/9B1G 等  
2. Suo Y et al. *Nat Commun* 2025;16:5178. doi:10.1038/s41467-025-60480-3 — PDB 9DKA/9DKB  
3. Wu C et al. *Cell Discov* 2025 — native URAT1 drug inward structures (9JDZ/9JDY/9JE1); urate occ/out **not** separate PDB — use Dai 9B1K/9B1L
4. Burns RL et al. *Arthritis Res Ther* 2016;18:214. doi:10.1186/s13075-016-1107-x — 细胞 IC50 面板  
5. Nakamura M et al. *J Pharmacol Exp Ther* 2019;371:162-177. doi:10.1124/jpet.119.262741 — dotinurad 选择性  
6. Tan PK, Liu S, Gunic E, et al. *Sci Rep* 2017;7:665. doi:10.1038/s41598-017-00706-7 — verinurad 25 nM  

### NLRP3

1. Coll RC et al. *Nat Med* 2015;21:248-255. doi:10.1038/nm.3806 — MCC950  
2. McBride C et al. *J Med Chem* 2022;65:14721-14739. doi:10.1021/acs.jmedchem.2c01250 — GDC-2394  
3. Dekker A et al. *J Mol Biol* 2021 — 7ALV（NP3-146 类抑制剂，MCC950 药效团模板）  

---

## 使用规则（避免误用）

1. **ML 回归**：URAT1 优先用 Burns 2016 HEK-URAT1 行；结构考试用 Dai 2024  
2. **不要混用** 不同细胞系/construct 的 IC50 训练单一模型而不条件化  
3. **verinurad / MCC950 / GDC-2394** 在训练集 → 只报 sanity，不宣称外推  
4. **lesinurad / dotinurad** 是 scaffold-novel 硬考试  
5. 阴性 pass = **排名低**，不是「活性数字高」

---

## 脚本读取示例

```python
import pandas as pd

df = pd.read_csv("data/benchmarks/literature_benchmarks.csv")
# 每个化合物主推荐行（每 compound_id 取第一行）
primary = df.drop_duplicates("compound_id", keep="first")
# 仅 scaffold-novel Tier1a
novel = df[(df.validation_tier == "Tier1a") & (df.scaffold_novel_expected == True)].drop_duplicates("compound_name")
```
