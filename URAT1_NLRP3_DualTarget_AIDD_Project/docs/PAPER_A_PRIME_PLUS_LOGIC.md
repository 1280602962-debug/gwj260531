# 论文 A′+ 完整逻辑 — 不对称双证据 + 共享库双靶优先（修订版）

> **⚠️ 已归档（2026-07）**：文中 8973 双靶 Pareto 部分已废弃；8973 **仅作 URAT1 回顾**。当前提纲：[`MANUSCRIPT_OUTLINE_CURRENT.md`](MANUSCRIPT_OUTLINE_CURRENT.md)

> **回应漏洞**：仅 6 药 MD 不够；不做 geometric mean 主图则「双靶抑制剂从哪来」；OAT 迁移不能当创新。  
> **首投**：*Journal of Molecular Modeling*（Subscription / 非 OA）  
> **前置文档**：`DUAL_TARGET_AND_FAST_JOURNALS.md`、`PAPER_A_PRIME_JMM_ACTION_PLAN.md`

---

## 0. 核心逻辑链（先读这一段）

```
痛风双节点（URAT1 代谢 + NLRP3 炎症）
        ↓
ChEMBL 训练集 0 SMILES 重叠 → 不能从「双靶监督学习」直接得到共抑制剂
        ↓
必须在【同一共享化合物库】上，对每条 SMILES 分别算 URAT1 证据 + NLRP3 证据
        ↓
两靶证据类型不对称（benchmark 证明：URAT1 靠结构，NLRP3 靠 ML）
        ↓
共享库上双目标优先（Pareto / 双百分位交集）→ 计算短名单（非实验验证的 hit）
        ↓
代表药 + 短名单 Top 做 MD → 机制解释
        ↓
回顾性检验：已知单靶药在双证据平面上的位置 + 富集指标
```

**本文交付物**：  
- **不是**「已发现双靶抑制剂」  
- **是**「在零训练重叠条件下，如何用不对称双证据在共享库上产生可实验跟进的双靶 **计算短名单**」

---

## 1. 原 A′+ 的三大漏洞与修补

| 漏洞 | 问题本质 | 修补 |
|------|----------|------|
| **仅 6 药 MD** | 无库筛、无定量双靶结果 | 共享库 **8973**（或 A+D+E 子集）双证据打分 + 主图 **双证据平面/Pareto** |
| **不做融合排序** | 双靶抑制剂无从定义 | 在 **同一 SMILES 列表** 上做 **双目标优先**（见 §3），不是 ChEMBL 训练集融合 |
| **OAT 迁移** | Δρ≈0.004，无方法学意义 | **全文删除**，不作创新、不作消融、不作 SI |

---

## 2. 双靶抑制剂在计算上怎么定义（无湿实验版）

### 2.1 不能用的路径（数据不支持）

| 路径 | 为何不行 |
|------|----------|
| 从 ChEMBL 822∩513 找双靶药 | **交集 = 0** |
| MTL / 双任务几何平均标签 | 无双标签样本 |
| OAT→URAT1 迁移提升筛选 | 无显著提升 |
| 三态 $S_\pi$ 全库 | B1K/B1L 刚性失败；留给论文 B |
| 宣称「虚拟筛选发现新双靶 hit」 | 无 IC50/IL-1β 实验 |

### 2.2 唯一自洽的定义：共享库上的双证据共优先

对集合 $\mathcal{L}$ 中每个分子 $x$：

| 臂 | 主证据（A′+ 默认） | 符号 | 依据 |
|----|-------------------|------|------|
| **URAT1** | 9DKB Glide **XP** 分数（越高越好） | $S_U(x)$ | ML benchmark **2/4**；转运体宜结构 |
| **NLRP3** | Assay-conditioned **$P_\text{active}(x)$** | $S_N(x)$ | AUROC≈0.89；数据够监督 |

可选（SI / Top 50 only）：NLRP3 8ETR XP、URAT1 ML $\hat{y}$ 作 **消融臂**。

**双靶计算候选**：

$$
\mathcal{C}_{\text{dual}} = \left\{ x \in \mathcal{L} : R_U(x) \ge \tau_U \;\land\; R_N(x) \ge \tau_N \right\}
$$

$R_U, R_N$ = 在 $\mathcal{L}$ 内的 **百分位排名**（或 decoy-adjusted percentile，见 §4）。  
$\tau_U, \tau_N$ 建议 **90th** 或 **95th**（主文做敏感性 85–95）。

**Pareto 非支配**（主图推荐，替代 geometric mean）：

- 目标：最大化 $(S_U, S_N)$  
- 报告 Pareto front 大小、前 20 结构、与已知药的相对位置  
- **消融对照**（SI）：固定 0.5/0.5 线性融合、几何平均 — 展示 **为何不用** 对称融合（URAT1 ML 拉低或结构分与 ML 分不可通约）

### 2.3 共享库 $\mathcal{L}$ 取什么

| 方案 | 规模 | 推荐 |
|------|------|------|
| **A′+ 主方案** | `distill_manifest.csv` **8973** 唯一 SMILES | ✅ 已建库，A 有标签、D 无标签 |
| 最小可发表 | A(822) + D(8000) | 你正在 dock D；**A 必须补 dock** |
| 过大 | 百万 Enamine | A′+ **不做**（无 NLRP3 实验跟进） |

你正在跑的 **8000 XP** = 子集 D；**还必须对子集 A（822 活性）补 9DKB XP**，否则无法做「活性 vs decoy」富集，也无法在统一标尺上定位四药。

---

## 3. 文章目的（修订，可写进 Abstract）

**中文**：  
在痛风 URAT1–NLRP3 双节点框架下，提出并检验一种 **不对称双证据优先策略**：对共享化合物库，URAT1 用 9DKB 结构对接、NLRP3 用 assay-conditioned 机器学习，通过 Pareto/双百分位产生双靶 **计算短名单**；以回顾性富集和已知药物定位验证；对代表分子与短名单 Top 做 MD。  

**不声称**：实验验证的双靶抑制剂、新融合算法、OAT 迁移。

**英文标题方向**：

*Asymmetric dual-evidence prioritization on a shared compound library for URAT1 and NLRP3 in gout: structure-based transporter scoring, assay-conditioned inflammasome modeling, and retrospective validation*

---

## 4. 创新点（可辩护的 4 条，无 OAT）

1. **问题设定**：系统说明 ChEMBL **0 重叠** 下双靶共抑制剂不能从训练集 MTL 得到，必须转向 **共享库双证据**（Gap + 方法学动机）。  
2. **不对称证据分配**（empirically motivated）：URAT1 以 **9DKB XP** 为主、NLRP3 以 **assay-conditioned ML** 为主；用 benchmark 与 **ML+ML 消融** 证明 URAT1 不宜单用指纹回归。  
3. **双目标优先短名单**：Pareto / 双百分位在 8973 库上生成 $\mathcal{C}_{\text{dual}}$，并用 **已知单靶药六元组** 在 $(S_U, S_N)$ 平面定位（lesinurad 是否只强 UR 臂等）。  
4. **回顾性 VS 指标**：A vs D（+ benchmark）的 EF@k、AUC；回答「结构对接能否富集 URAT1 活性」— 与双靶短名单 **同一套** $S_U$。

**不是创新**：OAT 迁移、Teacher 8973、三态 $\pi$、geometric mean 新公式、Path B 生成。

---

## 5. 主文结构（完整）

### Introduction
- 痛风双节点；0 重叠 Gap  
- 为何 PLK1 式对称融合/锚点相似性不适用（简表）  
- 本文：**共享库 + 不对称双证据 + 回顾验证**

### Methods
- **§2.1** 数据：URAT1 822、NLRP3 513、manifest 8973  
- **§2.2** URAT1 ML（XGBoost+conformal）— **用于消融**，非主排序  
- **§2.3** NLRP3 assay-conditioned ML — **主 $S_N$**  
- **§2.4** 9DKB Glide SP→XP — **主 $S_U$**（全 $\mathcal{L}$ 或 A+D）  
- **§2.5** 双证据优先：百分位、Pareto、$\mathcal{C}_{\text{dual}}$ 阈值  
- **§2.6** 回顾性指标：EF@5%、AUC（$S_U$ vs A/D 标签）  
- **§2.7** MD/MM-GBSA：6 benchmark + **Pareto Top 3–5**（见 §6）  
- **§2.8** 消融：ML+ML、线性 0.5/0.5、几何平均（**放 SI 或 Methods 一段**）

### Results
| 节 | 内容 | 主图 |
|----|------|------|
| 3.1 | URAT1 ML benchmark **2/4** | 表 S1 |
| 3.2 | NLRP3 ML AUROC + benchmark 2/2（训练集内） | ROC 图 |
| 3.3 | **$S_U$ 回顾性富集**（A vs D） | EF 曲线 |
| 3.4 | **双证据平面**：8973 散点 + 六药标注 + Pareto front | **Fig 4（主图）** |
| 3.5 | $\mathcal{C}_{\text{dual}}$ 数量、Top 20 化学空间、与 decoy 比 | 表 2 |
| 3.6 | URAT1 四药对接 + MD | Fig 5–6 |
| 3.7 | NLRP3 两药对接 + MD | Fig 7 |
| 3.8 | Pareto Top 3–5 MD（若 pose 合理） | Fig 8 或 SI |

### Discussion
- 短名单含义：**实验待验证**，非成药  
- 若 $\mathcal{C}_{\text{dual}}$ 富集 decoy/泛筛选命中 → **诚实写**（Sindt 2025 rescoring 语境）  
- 单靶药（lesinurad / MCC950）在平面上的 **单臂强** → 符合临床现实  
- 与 PLK1/NLRP3、对称融合的差异  
- 局限：无湿实验；URAT1 单态 9DKB；NLRP3 结构未全库对接  

### Conclusion
- 交付：**可复现的双证据优先流程 + 短名单 + 回顾指标**  
- 未来：论文 B 三态 URAT1；湿实验验证 Top hits

---

## 6. MD 做多少、对谁做（修订）

| 类别 | 数量 | 目的 |
|------|------|------|
| **Benchmark** | URAT1×4 + NLRP3×2 = **6** | 已知药机制、redock 验证 |
| **Pareto Top** | **3–5** | 证明短名单分子复合物 **动力学可模拟**（pose 来自对接或后续 NLRP3 dock） |
| **全库 8973** | **0** | 不做 |
| **$\mathcal{C}_{\text{dual}}$ 全做** | **否** | 只对 Top 做 |

若 Pareto Top 全是「类药 decoy」且 MD 不稳定 → **仍可作为结果**：计算短名单需实验过滤。

---

## 7. 计算任务清单（按依赖排序）

### Phase 1 — 共享库双分数（主文核心，与 8000 并行）

- [ ] **8000 D**：9DKB XP（进行中）→ 导出 `glide_xp_score`  
- [ ] **822 A**：补 9DKB XP（**必做**）  
- [ ] **E 等其余 manifest**：可选；至少保证 **六药 + Pareto top** 有分  
- [ ] **8973 全库**：NLRP3 `P_active` 批量预测（脚本快，**立即能做**）  
- [ ] 合并表：`results/paper_a_prime/dual_evidence_scores.csv`  

列：`smiles, subset, S_U, S_N, R_U, R_N, pareto_rank, in_C_dual`

### Phase 2 — 图表与回顾指标

- [ ] A vs D：$S_U$ 的 AUC、EF@5%  
- [ ] 双证据散点 + 六药标注 + Pareto  
- [ ] $\mathcal{C}_{\text{dual}}$ 统计（阈值敏感性）  
- [ ] 消融：若用 $\hat{y}_U$ 代替 $S_U$，四药/富集变化（SI）

### Phase 3 — 结构深化

- [ ] NLRP3：GDC-2394 @ 8ETR、MCC950 @ 7ALV 对接 + redock  
- [ ] Pareto Top 3–5：缺 NLRP3 pose 的补 dock 8ETR  
- [ ] MD：6 benchmark + 3–5 Top  

### Phase 4 — 写作与投稿

- [ ] 摘要：**asymmetric dual-evidence prioritization**，非 discovery  
- [ ] 投稿 JMM，Subscription  

### 明确删除

- [ ] ~~OAT 迁移~~  
- [ ] ~~geometric mean 主图~~  
- [ ] ~~Teacher 8973 三态~~  
- [ ] ~~宣称双靶新药~~  

---

## 8. 审稿人可能问什么 — 预设回答

| 问题 | 回答 |
|------|------|
| 双靶抑制剂在哪？ | $\mathcal{C}_{\text{dual}}$ 短名单 + Top MD；**待实验验证** |
| 为何不用 ML 筛 URAT1？ | Benchmark 2/4 + 富集对比；转运体构象依赖 |
| 为何不全库 NLRP3 对接？ | ML 已够；结构仅验证 benchmark + Top（资源限制） |
| 与 PLK1/NLRP3 何异？ | 不对称证据类型 + Pareto；非 5-anchor + 0.5 融合 |
| OAT 迁移？ | **未采用**（预实验无 practical gain） |
| 6 药 MD 够吗？ | 6 药 + Top 3–5 + 8973 层级的 **双证据图** 才是主贡献 |

---

## 9. 能否发表（诚实预期）

| 配置 | JMM 概率 | 说明 |
|------|----------|------|
| 仅 6 MD | 低 | 你已判断正确 |
| 6 MD + 8000 富集 | 中 | 仍缺双靶定义 |
| **A′+ 完整（§7）** | **中–中高** | 有方法+数据+短名单+局限 |
| + 湿实验 1–3 个 | 高 | 非当前 |

被拒可转 **Chemical Biology & Drug Design**（同叙事）。

---

## 10. `dual_evidence_scores.csv` 规格（待脚本）

```csv
canonical_smiles,subset,has_urat1_label,has_nlrp3_label,S_U_glide_xp,S_N_ml_pactive,R_U_pct,R_N_pct,pareto_layer,in_C_dual_90,benchmark_name
```

生成脚本建议：`scripts/build_dual_evidence_table.py`（NLRP3 推理 + 合并 Glide 导出）。

---

## 11. 文档关系

| 文件 | 角色 |
|------|------|
| **PAPER_A_PRIME_PLUS_LOGIC.md** | 本文件 — 完整逻辑 |
| `PAPER_A_PRIME_JMM_ACTION_PLAN.md` | 操作清单（需按 §7 更新） |
| `distill_manifest.csv` | 共享库 $\mathcal{L}$ |
| `DIFFERENTIATION_VS_PLK1_NLRP3.md` | 与对称融合差异化 |
