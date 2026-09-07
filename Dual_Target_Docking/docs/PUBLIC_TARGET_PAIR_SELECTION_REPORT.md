# 公开双靶靶点对硬门槛评估报告（v1）

**日期：** 2026-07-23  
**数据源：** ChEMBL Web API（`pchembl_value` 非空；类型 ∈ IC50/Ki/Kd/EC50 等）+ RCSB Search API（分辨率 ≤ 3.5 Å 且含非聚合物配体的 holo 条目粗计数）  
**标签规则：** 同一分子在该靶上取 **最大 pChEMBL**；**pChEMBL ≥ 6 = 活**；**已测且 pChEMBL < 6 = 弱/阴**；**未测 ≠ 阴**。  
**机器表：** [`../data/public_pair_selection/chembl_pair_fourclass.csv`](../data/public_pair_selection/chembl_pair_fourclass.csv) · [`pdb_holo_counts.csv`](../data/public_pair_selection/pdb_holo_counts.csv) · 冻结清单 [`FROZEN_PUBLIC_PAIRS.yaml`](../data/public_pair_selection/FROZEN_PUBLIC_PAIRS.yaml)  
**复现脚本：** [`../scripts/audit_public_target_pairs.py`](../scripts/audit_public_target_pairs.py)

---

## 1. 评估范围

共审计 **12 个候选靶点对**（覆盖文献高频、共晶金标准、经典双抑制剂通路）：

| 候选对 | 动机 |
|--------|------|
| EGFR / HER2 | Tier A 双端共晶；经典双 TKI |
| Mcl-1 / Bcl-xL | Tier A 双端共晶；linked 金标准 |
| PIK3CA / mTOR | 经典双抑制剂；ChEMBL 预期规模大 |
| BRD4 / HDAC1、BRD4 / HDAC6 | 100 篇文献最热 HDAC–X |
| PARP1 / MET | 100 篇中有双靶文 |
| AChE / BChE | 双胆碱酯酶；公开数据厚 |
| JAK2 / HDAC1 | 文献双靶常见 |
| CDK6 / BRD4 | 100 篇有双靶文 |
| AKT1 / p70S6K | 通路双靶（M2698 类） |
| Mcl-1 / Bcl-2 | 与 Bcl-xL 并列的 BH3 对 |
| PIK3CA / PIK3CB | 同工酶对（对照「过近」） |

---

## 2. 硬门槛判定标准（本报告执行版）

| # | 门槛 | 本审计操作化 |
|---|------|----------------|
| H1 | 两端定量活性分子足够 | 每端 unique 分子 **≥ 200**（有 pChEMBL 的 IC50/Ki/Kd/EC50） |
| H2 | 配对后能分出四类（至少三类） | 两端**都测过**的分子中：dual / A-only / B-only 各自 **≥ 10** |
| H3 | 两端有可用 holo PDB | RCSB：reso ≤ 3.5 Å 且含配体实例；建议每端 **≥ 5** 条 |
| H4 | 可独立二元对接 | 常规小分子口袋优先；**Zn-HDAC / 纯 PPI 免疫检查点**不作**唯一**主对 |

未测分子不进入 A-only/B-only/双弱计数。

---

## 3. ChEMBL 四类结果总表

| 靶点对 | n_A | n_B | 两端都测 | dual | A-only | B-only | 双弱 | H1 | 三类≥10 | 硬门槛活性 |
|--------|-----|-----|----------|------|--------|--------|------|----|---------|------------|
| **PIK3CA / mTOR** | 7732 | 5209 | **2713** | 2002 | 266 | 236 | 209 | ✓ | ✓ | **通过** |
| **EGFR / HER2** | 11198 | 2619 | **1751** | 1182 | 207 | 46 | 316 | ✓ | ✓ | **通过** |
| **AChE / BChE** | 6197 | 3798 | **2537** | 986 | 483 | 225 | 843 | ✓ | ✓ | **通过** |
| **PIK3CA / PIK3CB** | 7732 | 2786 | 1990 | 988 | 299 | 213 | 490 | ✓ | ✓ | 通过* |
| **Mcl-1 / Bcl-2** | 3412 | 3171 | 371 | 160 | 87 | 34 | 90 | ✓ | ✓ | **通过** |
| **Mcl-1 / Bcl-xL** | 3412 | 1612 | **305** | 82 | 77 | 24 | 122 | ✓ | ✓ | **通过** |
| **AKT1 / p70S6K** | 3699 | 1918 | 601 | 453 | **11** | 100 | 37 | ✓ | ✓† | 通过但脆 |
| BRD4 / HDAC1 | 9338 | 8348 | 82 | 66 | **0** | 10 | 6 | ✓ | ✗ | **否** |
| JAK2 / HDAC1 | 12711 | 8348 | 76 | 52 | 23 | **1** | 0 | ✓ | ✗ | **否** |
| BRD4 / HDAC6 | 9338 | 6554 | 39 | 27 | 3 | 9 | 0 | ✓ | ✗ | **否** |
| PARP1 / MET | 4576 | 4658 | 13 | 10 | 1 | 0 | 2 | ✓ | ✗ | **否** |
| CDK6 / BRD4 | 838 | 9338 | 6 | 5 | 1 | 0 | 0 | ✓ | ✗ | **否** |

\*同工酶对：活性门槛过，但叙事上「过近」，不作三主对之一。  
†A-only 仅 11，统计脆弱。

---

## 4. 结构 / 对接 / 加分项

| 靶点对 | holo PDB 粗计 (A/B, ≤3.5Å) | 独立二元对接 | 双端同配体共晶 | 文献双靶 |
|--------|----------------------------|--------------|----------------|----------|
| PIK3CA / mTOR | 106 / 41 | ✓ 常规 | 本 catalog 无 | 大量经典双抑制（另有 PROTAC，需分表） |
| EGFR / HER2 | 343 / 28 | ✓ 激酶 ATP | **有** TAK-285 (3POZ/3RCD) | 大量双 TKI |
| Mcl-1 / Bcl-xL | 122 / 72 | ✓ 但 PPI 沟槽更难 | **有** LC6 (3WIY/3WIZ) | Tanaka JMC 等 |
| AChE / BChE | 74 / 108 | ✓ | 非本 catalog Tier A | 大量双胆碱酯酶 |
| BRD4 / HDAC* | 584 / 6(HDAC1) | HDAC=Zn **金属酶** | 多为仅 BRD4 端 | **100 篇最热** |
| PARP1 / MET | 90 / 113 | ✓ | 无（配对活性也无） | 有个例 |

---

## 5. 逐对结论（过 / 备 / 否）

### 锁定为公开主对（3）

#### DTPAIR-01 — **PIK3CA / mTOR**（主规模对）

- **硬门槛：** 全过；配对 2713，四类齐全且 A-only/B-only 都很厚 → **最适合打假朴素融合**。  
- **强烈建议：** 文献双靶充分；化学型以经典 dual 为主，**PROTAC 另表**。  
- **加分（2026-07-27 纠偏）：** 有双端同配体共晶 — **PI-103**（chem_comp **X6K**）在 **4L23**（PIK3Cα, P42336, 2.5 Å）与 **4JT6**（mTOR, P42345, 3.6 Å）。对接主结构冻结为 **4L23 + 4JT6**；勿用 7L1C（HLA+肽）、9CMK（RBD glue）、4DRI（FKBP+rapamycin/FRB）。  
- **风险：** 通路相关，存在「本来就会双活」的化学型 → 必须保留硬负（A-only/B-only）并做 scaffold 泄漏控制。

#### DTPAIR-02 — **EGFR / HER2**（激酶 + 姿态金标准）

- **硬门槛：** 全过；B-only=46 可用。  
- **加分：** Tier A 双端共晶（TAK-285）。  
- **风险：** 口袋同源，交叉活性高；dual 比例高属预期，仍以 A-only/B-only 污染率做主证据。

#### DTPAIR-03 — **Mcl-1 / Bcl-xL**（异质折叠 + 姿态金标准）

- **硬门槛：** 全过；配对 305，规模小于前两对但仍可做 dual-vs-single。  
- **加分：** Tier A linked 共晶（3WIY/3WIZ）。  
- **风险：** PPI 对接更难 → 单靶 QC（共晶 RMSD + EF）必须先过关；不过关则该对只作姿态案例、不作融合主表。

### 强烈建议满足，但未进「三主对」

| 对 | 理由 |
|----|------|
| **AChE / BChE** | 硬门槛极好，可作 **第 4 对 / leave-pair**；口袋过近，不宜占三个主名额中的多样性位。 |
| **AKT1 / p70S6K** | 三类刚过线（A-only=11），仅 backup。 |
| **Mcl-1 / Bcl-2** | 与 DTPAIR-03 化学空间重叠；姿态金标准更认 Bcl-xL。 |

### 硬门槛或一票否决 — 不作公开主对

| 对 | 主要否决点 |
|----|------------|
| BRD4 / HDAC1、BRD4 / HDAC6 | **拼不出够用的 A-only**；HDAC Zn；文献热 ≠ 配对标签可用 |
| JAK2 / HDAC1 | B-only≈0 |
| PARP1 / MET | 两端都测仅 13 分子 |
| CDK6 / BRD4 | 配对仅 6 |
| PIK3CA / PIK3CB | 同工酶过近，不作主叙事对 |
| 纯 PROTAC 三元 | **另表**，不进普通双抑制 D2 |

---

## 6. 最终冻结（请按此建 D2 / 冻 YAML / 做单靶 QC）

| 优先级 | pair_id | 靶点对 | UniProt | ChEMBL | 角色 |
|--------|---------|--------|---------|--------|------|
| 1 | DTPAIR-01 | **PIK3CA / mTOR** | P42336 / P42345 | CHEMBL4005 / CHEMBL2842 | 主评测规模 |
| 2 | DTPAIR-02 | **EGFR / HER2** | P00533 / P04626 | CHEMBL203 / CHEMBL1824 | 激酶 + 双端共晶 |
| 3 | DTPAIR-03 | **Mcl-1 / Bcl-xL** | Q07820 / Q07817 | CHEMBL4361 / CHEMBL4625 | 异质口袋 + linked 共晶 |

**可选第 4 对（未冻结）：** AChE / BChE，用于 leave-pair-out 或敏感性分析。

---

## 7. 与「100 篇文献」的关系（避免误解）

- 100 篇里 **HDAC 最热**，但 ChEMBL 上 BRD4–HDAC / JAK–HDAC **几乎没有 A-only** → **不能**当公开主对。  
- 冻结的三对里，**EGFR/HER2、Mcl-1/Bcl-xL 几乎不在那 100 篇**，靠的是 **公开库规模 +（后两对）双端共晶**，不是文献热度。  
- 文献 curated 仍可从 100 篇抽 HDAC–X 等进 **D2-curated**，但主打假表必须跑在上述冻结公开对上。

---

## 8. 局限（写进 Methods 的诚实句）

1. ChEMBL 配对是「同分子两端都有记录」，**不是**同一实验条件下的平行测定。  
2. max-pChEMBL 会偏乐观；正式入库需再滤 assay type、去盐立体、InChIKey 去重。  
3. RCSB holo 计数含溶剂/添加剂配体噪声，冻对接 YAML 时需人工挑 **真正抑制剂 holo**。  
4. 本审计的面板建造以 ChEMBL pChEMBL 为准。冻结 K=4 的 BindingDB/PubChem **计数核对**（零对接）见 `data/jcim_supply_crossdb_v0/` 与 Supporting Information Table S12：等式测定下厚面板门槛不翻转；EGFR/HER2 升至薄面板（min HN ≈ 30）仍达不到 ≥50。

---

## 9. 下一步（计算侧，可立刻开干）

1. 按 `FROZEN_PUBLIC_PAIRS.yaml` 写三对的对接 YAML（受体 PDB、box、种子）。  
2. 每靶单靶 QC：共晶自对接 RMSD + 已知活性 vs decoy EF。  
3. 脚本化导出三对的 D2-public 四分类表（带 scaffold split）。  
4. 合成侧文献金标准：优先挂到这三对上的分子；其余对进 curated 旁支。
