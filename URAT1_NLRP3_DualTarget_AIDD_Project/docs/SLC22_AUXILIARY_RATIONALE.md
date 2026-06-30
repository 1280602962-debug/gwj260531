# SLC22 辅助库：科学逻辑、历史错误与可信边界

## 1. 为什么会出现 ChEMBL / 文献 / 数字错误？

| 错误类型 | 典型表现 | 根因 |
|----------|----------|------|
| **ID 张冠李戴** | OCT1 写成 `CHEMBL242`（ESR2）、`CHEMBL1906`（RAF1） | 手工抄 ChEMBL 编号、未用 target 页面核对 UniProt |
| **字段混用** | `ref_pmid=1467`（ChEMBL 化合物 ID）、`12280`（IUPHAR ligand ID） | 把非 PMID 填进 PMID 列 |
| **PDB–配体错配** | benzbromarone 标成 `9DKB`（实为 lesinurad） | 同一论文多结构，未查 RCSB 配体名 |
| **统计捏造** | NLRP3「503 / 92 assays / 47% 冲突」 | 未跑 `00_prepare_data.py` 就从旧笔记抄写 |
| **子家族混淆** | 把 **OCT（阳离子）** 当作 URAT1 迁移主库 | 仅见「SLC22 家族」统称，忽略 URAT1 的 **OAT（阴离子）** 分类 |

**结论**：此前对话与旧版文档中，**仅下列内容可默认可信**——其余须对照 `DATA_FACT_CHECK.md` 或重跑脚本：

- 已写入 `config/`、`data/benchmarks/` 且经本轮审计的 ID  
- `data/processed/data_summary.json`（本地运行 `00_prepare_data.py` 生成）  
- `docs/MODEL_QUALITY_REPORT.md` 等 **已执行** 的训练/回测输出  
- **不可信**：未实现的 MASFL v3.1 模块、URAT1 单模型大规模筛选结论（`URAT1_NO_GO`）、任何未标注来源的百分比

---

## 2. URAT1 是 OAT，不是 OCT

| 项目 | URAT1 | OCT1/OCT2 |
|------|-------|-----------|
| 基因 | **SLC22A12** | SLC22A1 / SLC22A2 |
| 亚家族 | **OAT（有机阴离子转运体）** | **OCT（有机阳离子转运体）** |
| 底物电荷 | 尿酸盐、有机**阴离子**交换 | 有机**阳离子** |
| 肾定位 | 近端小管 **顶膜**（重吸收尿酸） | OCT2 顶膜分泌阳离子；OCT1 肝摄取 |
| 与促尿酸排泄药 | **直接靶点** | **脱靶**（肝/肾阳离子转运） |

文献依据：NCBI Gene / OMIM 将 SLC22A12 定义为 **OAT 家族成员**；URAT1 与 OAT4 同属进化枝 F，OAT1/OAT3 属另一 OAT 枝（Inui, *Pharmacol Ther* 2012 综述 Fig.1）。

**因此**：用 OCT1/OCT2 作 URAT1 **主迁移库** 在化学空间与机制上 **弱于** OAT1/OAT3；原设计把 OCT 放在「辅助训练库」首位是 **科学表述不当**（算法上可用 OCT 大数据预训练，但论文必须写清 **OAT 优先 + OCT 脱靶**）。

---

## 3. 修正后的辅助库分层（项目标准）

### Tier A — **OAT 迁移预训练**（主）

| 靶点 | 基因 | UniProt | ChEMBL | 角色 |
|------|------|---------|--------|------|
| OAT1 | SLC22A6 | O95742 | CHEMBL1641347 | 阴离子底物/抑制剂化学空间；与 probenecid、苯溴马隆等 uricosuric 重叠 |
| OAT3 | SLC22A8 | O95816 | CHEMBL1641348 | 肾阴离子分泌；URAT1 通路下游/并行肾排泄 |

**用途**：TC-Encoder / MLP head 序贯微调；**须做消融** `Abl-7a`（无 OAT 迁移）vs `Abl-7b`（无 OCT 脱靶特征）。

### Tier B — **OCT 脱靶与选择性**（辅，非主迁移）

| 靶点 | 基因 | UniProt | ChEMBL | 角色 |
|------|------|---------|--------|------|
| OCT1 | SLC22A1 | O15245 | CHEMBL2073664 | 肝摄取脱靶；阳离子药物 DDI |
| OCT2 | SLC22A2 | O15244 | CHEMBL1770032 | 肾分泌脱靶 |

**用途**：对接比值 $R_{\text{sel}}$（Tier 3 **计算假说**）；**不**声称 OCT 抑制剂 = URAT1 抑制剂。

### 可选 Tier C — OAT4（讨论）

- **SLC22A11 / OAT4**：lesinurad 亦抑制 OAT4（Burns 2016）；与 URAT1 同 apical OAT 枝，数据量小于 OAT1/3，可作敏感性分析，**非当前必需导出**。

---

## 4. 与 verinurad / dotinurad 选择性表述一致

Benchmark 中 verinurad「相对 OAT1/OAT4 ~200-fold 选择性」说明：**临床 URAT1 选择性讨论本就包含 OAT 成员**，而非 OCT。OCT 脱靶是 **另一类风险**（阳离子），应在 Discussion 分句表述。

---

## 5. 论文 Methods 推荐写法（英文模板）

> *URAT1 (SLC22A12) is an organic anion transporter (OAT) family member. Auxiliary ChEMBL bioactivity for **human OAT1 and OAT3** was used for representation pre-training before URAT1 fine-tuning. **OCT1/OCT2** data were used only for in silico off-target ranking (hepatic/renal organic cation transport), not as the primary transfer source. Transfer learning does not imply that OAT or OCT inhibitors are URAT1 inhibitors.*

---

## 6. 配置与文件

- `config/targets.yaml` → `auxiliary_targets.slc22_oat_transfer` / `slc22_oct_detargeting`  
- `data/auxiliary/README.md` — 导出路径与 ChEMBL ID  
- `docs/DATA_FACT_CHECK.md` — ID 与规模数字  
