# JNK2 选择性共价抑制剂筛选方案

**项目类型：** 计算—化学—生化一体化先导发现  
**靶点：** 人源 JNK2（MAPK9, c-Jun N-terminal kinase 2）  
**共价位点：** **Cys116**（ATP 结合区近端，JNK1–3 中 JNK1/2 保守；JNK3 为 Cys154）  
**选择性目标：** **JNK2 > JNK1**（JNK2/JNK3 可评估；JNK3 脑组织表达需单独考虑）  
**编制日期：** 2026 年 7 月  
**版本：** v1.0

---

## 摘要

JNK1 与 JNK2 在 ATP 口袋序列同源性高达 **98%**，仅 **Met77/Leu77** 与 **Ile106/Leu106** 两处差异，传统可逆抑制剂极难实现亚型选择性。Gray/Westover 组通过 **JNK-IN-8（pan-JNK 共价）** 的结构优化，获得 **YL5084**——首个对 JNK2 显示 **kinact/KI 较 JNK1 高约 21 倍** 的共价抑制剂，并经 **8ELC 共晶** 确证 Cys116 共价键 [1]。

本方案以 **YL5084 化学系列为锚点**，采用 **MedChem 类似物库 + 可选 AF3/共价对接缩库** 的策略，以 **kinact/KI、C116S 突变拯救、细胞 biotin 探针竞争** 为选择性核心 readout，**明确不使用** JNK1 非共价项目中已证伪的 **Δsel 对接排名** 作为采购决策依据 [2]。

---

## 一、科学背景与立项依据

### 1.1 生物学 rationale

| 亚型 | 主要功能倾向 | 选择性抑制剂价值 |
|------|-------------|-----------------|
| JNK1 | 促凋亡、TNFα/UV 应激 | JNK1 选择性探针 |
| **JNK2** | 促存活（情境依赖）、MM 中组成型激活 | **JNK2 选择性探针** |
| JNK3 | 神经退行（脑表达为主） | CNS 靶点 |

JNK1 与 JNK2 在多种细胞情境下功能 **拮抗**（如 UV 应激、MM 细胞存活）[1,3]。亚型选择性化学探针是解析信号通路金标准，但 ATP 口袋高度保守使非共价选择性极难实现 [4]。

### 1.2 共价策略的优势

- **Cys116** 位于 JNK1/2 催化域 ATP 结合 motif 近端，JNK-IN-8 已证明为 **可药化共价位点** [5]。
- 共价 **kinact/KI** 可同时编码非共价预定位（KI）与成键速率（kinact），为亚型选择性提供 **独立于 IC50 单点读数** 的维度 [1,6]。
- YL5084 通过 **(R)-3-氨基吡咯烷** 替换苯环、**flag methyl** 等修饰，在保持 Cys116 共价前提下，利用 **Leu106（JNK2）vs Ile106（JNK1）** 及 P-loop 动态差异实现选择性 [1]。

### 1.3 与仓库内 JNK1 非共价项目的衔接

`JNK1_Selectivity_Project` 结论：**Δsel 对接/MM-GBSA 无法区分 benchmark 选择性方向**（标定通过率 43%，ML 选择性 F1=0）[2]。  
本课题 **不复用** 该选择性 gate，但可复用：

- **3E7O**（JNK2 DFG-in 非共价共晶）作非共价姿态交叉验证；
- Schrödinger IFP/MM-GBSA 批处理脚本（`jnk_docking_export/`，非共价 benchmark 专用）作 **辅助**，不作选择性主排序。

---

## 二、靶点、结构与选择性机制

### 2.1 共价位点

| 项目 | 内容 |
|------|------|
| 靶点蛋白 | JNK2（MAPK9），催化域残基 4–364（与 8ELC 一致） |
| 共价残基 | **Cys116**（SG 与丙烯酰胺 β-碳 Michael 加成） |
| 弹头 | **丙烯酰胺**（与 JNK-IN-8 / YL5084 / BTK 共价药一致） |
| JNK1 对应 | Cys116（序列保守；**反应性相当**，选择性来自非共价差异） |
| JNK3 对应 | **Cys154**（需单独评估 cross-reactivity） |

### 2.2 结构模板选用（重要）

| PDB | 配体 | DFG | 用途 | 备注 |
|-----|------|-----|------|------|
| **8ELC** | YL2056（共价） | **In** | **主受体** — 共价对接 / AF3 校准 | 2.0 Å；Cys116 共价键清晰 [1] |
| 7N8T | AMP | In | 活化环 / 构象参考 | 1.6 Å [1] |
| 3E7O | 35F（indazole） | In | 非共价交叉验证 | 勿与 8ELC 混为共价主模板 [2] |
| 3ELJ / 4L7F | 多种 | In | JNK1 同源建模 / 选择性对接 | JNK1 选择性参照 |
| **3NPC** | 抑制剂 | **Out** | **排除** | 8ELC 结构解算 MR 模型；**非**共价 DFG-in pose [1] |

> **原则：** 共价筛选 **仅使用 DFG-in 共晶 8ELC**；DFG-out（3NPC）与非共价 indazole（3E7O）不得作为共价排名主受体。

### 2.3 选择性结构假说（来自 YL5084 论文）

1. **Leu106（JNK2/3）vs Ile106（JNK1）**：后口袋芳环/ scaffold 占据；JNK1 Ile106 导致 **back pocket clash** [1]。
2. **Val54 vs Ile54**：主链柔性差异，影响配体诱导 fit [1]。
3. **Arg50–Glu109 盐桥（JNK1）**：稳定 P-loop，不利于该类配体结合；JNK2 Ile50 无此桥 [1]。
4. **Activation loop 动态**：YL2056 结合诱导活化环大幅移动；共价抑制剂对 **DFG-in 构象** 有构象锁定需求 [1]。

---

## 三、总体筛选架构

```
┌─────────────────────────────────────────────────────────────┐
│  Phase 0：回顾性校准（必过 gate）                             │
│  YL5084 / JNK-IN-8 / YL5084R @ Cys116                        │
│  kinact/KI、共晶 RMSD、C116S 拯救                             │
└──────────────────────────┬──────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  Phase 1：化合物库构建                                        │
│  Tier 0 对照 + Tier 1 YL5084 类似物 + Tier 2 丙烯酰胺 SMARTS │
└──────────────────────────┬──────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  Phase 2：计算缩库（可选并行）                                 │
│  共价 Glide SP/XP 或 AF3 @ Cys116 → 几何 QC → 短 MD           │
│  交叉对接 JNK1/JNK3 同源模型（非共价 ΔΔG 仅作参考）            │
└──────────────────────────┬──────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  Phase 3：湿实验确证（选择性主 readout）                       │
│  GSH t½ → LC-MS 位点 → kinact/KI → C116S → 细胞探针 → 表型   │
└─────────────────────────────────────────────────────────────┘
```

### 3.1 设计原则

1. **MedChem 优先，计算辅助** — 沿用 YL5084 论文路径 [1]，非 AF3 大库盲筛；
2. **选择性以 kinact/KI 为准**，不以对接 Δsel 为准 [1,2]；
3. **每个 hit 必须过 C116S 突变拯救** [5]；
4. **共价与非共价活性分离**：配置 **YL5084R（弹头饱和）** 阴性对照 [1]；
5. **AF3 使用需过 COValid 式校准**（见 §5.3）[7]。

---

## 四、化合物库构建

### 4.1 化学 scaffold 与弹头

**母核：** 2-芳基-吡唑并[1,5-a]吡啶 + 2-氨基嘧啶 + **(R)-3-氨基吡咯烷** 连接子 + **丙烯酰胺** 弹头（YL5084 系列）[1]。

**起点化合物：** JNK-IN-8 → THZ-3-60-1 / YL2056 → YL5084 演化路径 [1,5]。

### 4.2 分层建库

| 层级 | 内容 | 规模建议 | 来源 |
|------|------|----------|------|
| **Tier 0** | YL5084、YL2056、JNK-IN-8、YL5084R、CC-930 | 5–8 | 文献 / 供应商 |
| **Tier 1** | YL5084 类似物（吡咯烷取代、嘧啶/吡唑修饰） | 50–200 | 文献化合物 + 理性设计 |
| **Tier 2** | **丙烯酰胺** SMARTS 商业子库 | 2,000–20,000 | Enamine/ZINC 等 [7] |
| **Tier 3** | 以 8ELC pose 为种子的 R-group 枚举 | 200–1,000 | 内部生成 |

**Tier 2 SMARTS（丙烯酰胺核心）：**

```
[CH2]=[CH]C(=O)N
```

**预过滤（统一）：**

| 步骤 | 标准 |
|------|------|
| 分子量 | 350–650 Da（YL5084 MW ~601） |
| cLogP | 1–5 |
| TPSA | 80–150 Å² |
| 排除 | 卤乙酰胺、硝基烯、多弹头、PAINS |
| GSH 反应性 | t½ > 30 min（参考 YL5084 ~46 min）[1] |
| 去冗余 | 对 YL5084 Tc > 0.35 标注为类似物 |

### 4.3 不宜照搬的策略

| 策略 | 原因 |
|------|------|
| 906K 盲筛无校准（BTK 式） | 需先证明 JNK2 Cys116 上 AF3/对接能富集 YL5084 [7] |
| 非共价 Δsel 门槛 | JNK1 项目 benchmark 已否定 [2] |
| 3NPC DFG-out 对接 | 与 8ELC 共价 DFG-in 几何不一致 [1] |

---

## 五、计算筛选流程

### 5.1 Phase 0：回顾性校准 Gate（必做）

**未通过则停止 Tier 2 大库计算。**

| 对照 | 预期 |
|------|------|
| YL5084 @ JNK2 Cys116 | 共价 pose RMSD vs 8ELC < 2 Å；共价 S–C ~1.8 Å |
| YL5084 vs JNK1 | 计算结合模式显示 **Ile106 clash** 或 MM-GBSA(JNK1) > MM-GBSA(JNK2) + 2 kcal/mol |
| JNK-IN-8 | JNK1/JNK2 均高置信（pan 阳性） |
| YL5084R | 无法形成共价键 / 评分显著差于 YL5084 |

**AF3 专用 gate（若启用 AF3）：**

- 输入：JNK2 序列（4–364）+ 配体 3D + **bondedAtomPairs: Cys116-SG ↔ acrylamide β-C**
- 排序：mPAE（参考 Shamir COValid [7]）
- 通过标准：YL5084 mPAE **优于** 性质匹配 decoy 中位数；EF@1% ≥ 2

### 5.2 Phase 1：共价对接（Glide Covalent / Schrödinger）

| 参数 | 建议 |
|------|------|
| 受体 | 8ELC prep（ProtPrep，Cys116 质子化状态按 pH 7.4 优化） |
| 共价约束 | Cys116-SG — 丙烯酰胺 Cβ |
| 采样 | SP → XP 重打分 Top 500 |
| 排序 | GlideScore + 共价几何（SG–Cβ 距离、Cβ 方位角） |
| 交叉对接 | 同源 JNK1（3ELJ 为模板）/ JNK3 仅 **pose 可视化**，不作采购排序 |

### 5.3 Phase 2：AF3 共价共折叠（可选）

**适用条件：** Phase 0 AF3 gate 通过；靶点具备 **单 Cys + 深 ATP 口袋 + 统一丙烯酰胺** — 与 BTK 成功条件一致 [7]。

| 步骤 | 操作 |
|------|------|
| 1 | Tier 1 + Tier 2 子集（≤5,000）送 AF3 @ Cys116 |
| 2 | mPAE 排序 → 阈值参照 YL5084 校准 |
| 3 | 去已知相似 Tc > 0.35 |
| 4 | Top 50 → 50–100 ns Desmond MD（共价键稳定性） |
| 5 | 输出 Top 10–15 进入合成/采购 |

### 5.4 计算输出物

- `docking_results/covalent_scores.csv`
- `af3_results/mpae_ranked.csv`（若启用）
- `md_results/ligand_rmsd_cys116_distance.csv`
- **禁止输出：** 单独的 `delta_sel_dock` 采购清单

---

## 六、实验确证流程

### 6.1 化学与反应性

| 实验 | 方法 | 通过标准 |
|------|------|----------|
| 结构确证 | HRMS, NMR | 结构正确 |
| GSH 稳定性 | PBS + GSH 1–5 mM, 37°C, LC-MS | t½ ≥ 30 min（排除过高反应性） |
| 共价标记 | 重组 JNK2 + 化合物, LC-MS | 质量位移 +1 分子量 [1,5] |

### 6.2 生化选择性（核心）

| 实验 | 方法 | 通过标准 |
|------|------|----------|
| **kinact/KI** | Sox 荧光肽激酶 assay [1,8] | **JNK1/JNK2 kinact/KI ≥ 10×**（YL5084 ~21×） |
| 固定时间 IC50 | 1 h 预孵育后测 IC50 [1] | JNK1/JNK2 IC50 ratio ≥ 20×（YL5084 ~31×） |
| 位点确证 | 胰蛋白酶消化 + MS/MS | **仅 Cys116 修饰** [1] |
| **C116S 突变** | 重组 JNK2 C116S | IC50 或 kinact/KI **≥100× 右移**（JNK-IN-8 类 ≥100×）[5] |
| 激酶组 | KINOMEscan @ 1 μM | 仅 JNK 家族为主要 hit（YL5084 已优化掉 PIKFYVE）[1] |

**YL5084 参考数值 [1]：**

| 参数 | JNK1 | JNK2 |
|------|------|------|
| kinact/KI (M⁻¹ s⁻¹) | 335 | 7166 |
| 固定时间 IC50 (1 h) | 2173 nM | 70 nM |
| 细胞 JNK2 占据 | — | ~500 nM |
| 细胞 JNK1 占据 | 2 μM 仍弱 | — |

### 6.3 细胞选择性

| 实验 | 方法 | 预期 |
|------|------|------|
| Biotin-JNK-IN-7 竞争 | MM.1S 或 HEK293T, 6 h 预孵育 [1,5] | YL5084 类：JNK2 pull-down ↓，JNK1 需更高浓度 |
| NanoBRET | JNK1 vs JNK2 荧光标记靶标 [1] | EC50(JNK2) << EC50(JNK1) |
| c-Jun pSer63 | Western blot | 需 C116S 或 YL5084R 证明依赖共价 JNK2 |

### 6.4 表型与 on-target 验证

| 实验 | 说明 |
|------|------|
| JNK2 KO / C116S 过表达 | 抗增殖等表型是否仍依赖 JNK2（YL5084 在 MM 中 **部分 off-target**）[1] |
| CC-930 共处理 | 可逆 JNK2 抑制剂不能复制共价化合物全部表型 → 区分机制 |

> **注意：** YL5084 在 MM 细胞抗增殖 **不完全依赖 JNK2** [1]；本方案以 **生化/细胞占据选择性** 为 hit 标准，表型作为 secondary。

---

## 七、Hit 标准与 Go/No-Go

### 7.1 最小 hit 标准（MVP）

| 层级 | 标准 |
|------|------|
| 化学 | 结构确证；GSH t½ 可接受 |
| 共价 | MS 定位 **Cys116**；C116S **≥100× 拯救** |
| 选择性 | **kinact/KI (JNK1/JNK2) ≥ 10×**；细胞 JNK2 占据 EC50 ≤ 1 μM 且 JNK1 弱 |
| 选择性（理想） | kinact/KI ≥ 20×；接近 YL5084 |
| 安全性 | KINOMEscan @ 1 μM：无 PIKFYVE 级 off-target（参考 YL5084 优化路径）[1] |

### 7.2 阶段性 Gate

| 阶段 | Gate | 失败动作 |
|------|------|----------|
| Phase 0 | YL5084 计算/AF3 回顾性通过 | 不启动 Tier 2 大库 |
| Phase 1 | ≥2/5 Tier 1 类似物 kinact/KI ≥ 10× | 缩小 scaffold，回到 MedChem |
| Phase 2 | Top 10 中 ≥1 通过 C116S + MS | 调整共价对接参数或 AF3 输入 |
| Phase 3 | 1 化合物满足 MVP + 共晶尝试 | 论文/专利级先导 |

---

## 八、里程碑与资源

| 阶段 | 工作包 | 产出 |
|------|--------|------|
| W1–W2 | 8ELC 受体准备；Tier 0 采购；Phase 0 校准 | 校准报告 |
| W3–W6 | Tier 1 共价对接 / 可选 AF3；Top 20 清单 | 计算 hit list |
| W7–W12 | 10 化合物 GSH + MS + kinact/KI | 2–3 生化 hit |
| W13–W16 | 细胞探针 + C116S；激酶组 | 1 选择性 lead |
| W17–W20 | 共晶（可选）；初步 SAR | 论文数据包 |

**首轮合成/采购建议：** 8–12 个（Tier 0 全部 + Tier 1 Top 5 + 2 阴性对照）。

---

## 九、风险与应对

| 风险 | 说明 | 应对 |
|------|------|------|
| 选择性窗口极窄 | JNK1/JNK2 仅 2 个口袋差异氨基酸 [1] | 以 kinact/KI 非 IC50 单点；固定 (R)-吡咯烷立体化学 |
| 共价 off-target | 丙烯酰胺 全蛋白质组反应 [6] | GSH 预筛 + KINOMEscan + ABPP（可选） |
| 计算误导 | JNK1 项目 Δsel 失败 [2] | 禁止 Δsel 采购；MS + 突变为准 |
| 表型 off-target | YL5084 MM 抗增殖非 JNK2 依赖 [1] | 早期强调 C116S；不以外推表型为 hit 唯一标准 |
| AF3 过拟合 | NLRP3 Cys409 已有失败先例 | JNK2 必须先过 YL5084 AF3 gate |
| DFG 构象 | 活化环高度柔性 [1] | 8ELC 为主；短 MD 验证共价键 |

---

## 十、预期成果

1. **1 类** Cys116 共价、**JNK2 选择性优于 JNK1 ≥10×（kinact/KI）** 的化学探针；
2. 建立 **JNK2 共价选择性筛选 SOP**（含 Phase 0 gate 与禁用 Δsel 的文档化决策）；
3. 可选：JNK2–YL 系列共晶（延伸 8ELC 系列）；
4. 方法学：评估 AF3 @ JNK2 Cys116 相对 Glide covalent 的富集性能（COValid 式 decoy）。

---

## 参考文献

完整书目与链接见 [REFERENCES.md](./REFERENCES.md)。

| 编号 | 引用 |
|------|------|
| [1] | Lu W, Liu Y, Gao Y, et al. Development of a Covalent Inhibitor of JNK 2/3 with Selectivity over JNK1. *J Med Chem.* 2023;66(5):3356-3371. doi:[10.1021/acs.jmedchem.2c01834](https://doi.org/10.1021/acs.jmedchem.2c01834) |
| [2] | 本项目 JNK1 非共价筛选报告 v2.8（`JNK1_Selectivity_Project/docs/JNK1_PROJECT_REPORT.md`，分支 `cursor/purchase-list-md-05df`） |
| [3] | Liu J, Lin A. Role of JNK activation in apoptosis. *J Biol Chem.* 2005;280(22):21453-21456. doi:[10.1074/jbc.R500018200](https://doi.org/10.1074/jbc.R500018200) |
| [4] | Wang Y, et al. Unraveling the Design and Discovery of JNK Inhibitors. *J Med Chem.* 2022 Perspective. doi:[10.1021/acs.jmedchem.1c01947](https://doi.org/10.1021/acs.jmedchem.1c01947) |
| [5] | Zhang T, et al. Discovery of Potent and Selective Covalent Inhibitors of JNK. *Chem Biol.* 2012;19(1):140-154. doi:[10.1016/j.chembiol.2011.11.010](https://doi.org/10.1016/j.chembiol.2011.11.010) |
| [6] | Boike L, Nomura DK, Cravatt BF. Advances in covalent drug discovery. *Nat Rev Drug Discov.* 2022;21(12):881-898. doi:[10.1038/s41573-022-00516-z](https://doi.org/10.1038/s41573-022-00516-z) |
| [7] | Shamir Y, et al. Discovery of Covalent Ligands with AlphaFold3. *J Am Chem Soc.* 2026;148(12):13043-13054. doi:[10.1021/jacs.5c22222](https://doi.org/10.1021/jacs.5c22222) |
| [8] | Park H, et al. JNK2/3 isoform selective aminopyrazoles. *Sci Rep.* 2015;5:8047. doi:[10.1038/srep08047](https://doi.org/10.1038/srep08047) |

---

## 附录 A：Tier 0 对照化合物表

| 化合物 | 类型 | 用途 | 关键数据 |
|--------|------|------|----------|
| **YL5084** | 共价 JNK2 选择性 | 主阳性 / SAR 锚点 | kinact/KI JNK2/JNK1 ~21× [1] |
| **YL2056** | 共价 | 8ELC 共晶配体 | JNK1/JNK2 IC50 ratio ~33× [1] |
| **JNK-IN-8** | 共价 pan-JNK | 选择性阴性参照 | Cys116 必需 [5] |
| **YL5084R** | 非共价（饱和弹头） | 机制阴性对照 | 无共价标记 [1] |
| **CC-930** | 可逆 pan-JNK | 可逆对照 | Ki 5–61 nM 各亚型 [9] |

## 附录 B：PDB 结构索引

| PDB | 分辨率 (Å) | 内容 | DOI |
|-----|-----------|------|-----|
| [8ELC](https://www.rcsb.org/structure/8ELC) | 2.0 | JNK2–YL2056 共价 | [10.2210/pdb8elc/pdb](https://doi.org/10.2210/pdb8elc/pdb) |
| [7N8T](https://www.rcsb.org/structure/7N8T) | 1.6 | JNK2–AMP | [10.2210/pdb7n8t/pdb](https://doi.org/10.2210/pdb7n8t/pdb) |
| [3V6R](https://www.rcsb.org/structure/3V6R) | 2.97 | JNK3–JNK-IN-8 共价 | [10.2210/pdb3v6r/pdb](https://doi.org/10.2210/pdb3v6r/pdb) |
| [3E7O](https://www.rcsb.org/structure/3E7O) | 2.14 | JNK2–35F 非共价 | Shaw D, *J Mol Biol.* 2008;383(4):885-893 |
| [3ELJ](https://www.rcsb.org/structure/3ELJ) | — | JNK1 非共价 | JNK1 选择性对接参照 |

---

*文档版本 v1.0 | 2026-07*
