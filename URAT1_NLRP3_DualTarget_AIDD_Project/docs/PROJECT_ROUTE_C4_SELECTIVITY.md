# 路线 C4：以选择性为主轴的候选筛选（计算 SURI 画像）

> **现行路线。** 取代 C3（实验结构模板 + IFP 筛选）。机器可读锁：[`config/campaign_c4.yaml`](../config/campaign_c4.yaml)
> 冻结数据只读：`data/repurposing/p2/`、`data/campaigns/c1/`。

---

## 0. 为什么 C3 要废掉

C3 的内核是"换更好的受体结构 + 用 IFP 代替分数，再走一遍漏斗"。三个致命问题：

1. **体裁已被占满且对手更强。** Yang *RSC Adv.* 2023、Du 2024、荷叶碱文 2026 走的是同一条漏斗（对接 → 药效团/MM-GBSA → 候选），而且**都有 HEK 摄取实验**。我们没有湿实验，做同一体裁只能是更弱的版本。
2. **IFP 筛选不是新方法。** 相互作用指纹筛选是 2004 年 Deng/Rognan 就有的技术。C3 的"新颖性"被压缩成"9DK 系列 2025 年 6 月才公开、还没人拿它筛过库"——这是时间窗口，不是贡献。
3. **它没有解决真正让排序失效的原因。** C3 仍然在做**跨化学型的绝对强弱排序**：给不同分子在同一个口袋里打分再比大小。这正是 P2（True EF 2.587 → Random 0.215）和酸门（OR 0.970）失败的那个任务形态。换读出不换任务，失败机理照旧。

结论：不是换工具的问题，是**问的问题不对**。

---

## 1. 换掉的问题

| | 旧问题（C1–C3） | 新问题（C4） |
|---|---|---|
| 问什么 | 谁**结合** URAT1 更强 | 谁结合 URAT1 **而不**结合尿酸分泌转运体 |
| 计算形态 | 跨分子比较打分**幅值** | **同一分子在两个口袋之间的差** |
| 共模误差 | 分子大小/疏水体积直接进入排序 | 同一配体作差，共模成分抵消 |
| 验证数据 | 469 个 ChEMBL IC50（标签病理、27% 单骨架） | 4 药 × 4 转运体 IC50 矩阵，选择性比跨越约 1600 倍 |
| 临床含义 | "又一个虚筛" | 直接对应现有药物的失败原因 |

**选择性预测是同源蛋白间的同配体作差**，这是计算化学里成绩最好的任务形态之一；而跨化学型的绝对强弱排序是成绩最差的那一类。这次的任务形态从后者换成前者。

---

## 2. 为什么这是领域的真问题（不是我们造的缺口）

### 2.1 现有药物全部败于选择性，不是败于效价

| 药物 | 临床问题 |
|---|---|
| 苯溴马隆 | 肝毒性（线粒体抑制）；多国限用 |
| lesinurad | 肾脏不良反应；须与 XOI 联用 |
| 丙磺舒 | 广谱抑制 → DDI |
| **dotinurad** | 唯一以"选择性"为卖点获批 |

**SURI（selective urate reabsorption inhibitor）** 已是日本痛风指南（第 3 版 2022 补充版）收录的正式药理分类，定义为"强效 URAT1 抑制 + 对 ABCG2 与 OAT1/3 影响极小"。所以"按 SURI 画像挑分子"是领域公认目标，不是自创指标。

### 2.2 结构生物学家自己指出这个问题未解，且给了明确设计假设

- *Cell Discovery* 2025：「designing inhibitors for URAT1 is **challenging due to the high sequence homology among OAT family members**. Probenecid… also binds to OAT1. Aligning the URAT1 structure with probenecid-bound rOAT1 reveals a **high similarity in binding residues**, indicating a shared binding mode.」
- Dai & Lee *Cell Res.* 2024（结尾）：「**Engaging Arg477 can be explored and incorporated into the design of future compounds** targeting URAT1.」

第二条是一个**可检验的设计假设**，而且已发表结构里就有支持它的效价梯度：**verinurad 与 dotinurad 都与 Arg477 结合，二者正是最强的两个**；而低亲和力的 lesinurad 把转运体稳定在 apo 样构象，verinurad/dotinurad 稳定在收缩开放态。**没有人在计算上检验过这个假设。**

### 2.3 关键窗口：hOAT1 结构 2025 年才有

选择性计算过去做不了，因为只有 URAT1 一边有结构。现在两边都有：

| 蛋白 | 结构 | 分辨率 | 来源 |
|---|---|---|---|
| URAT1 | 9DK9 apo / 9DKA 苯溴马隆 / 9DKB lesinurad / 9DKC TD-3 | **2.68 / 3.00 / 2.74 / 2.55 Å** | Suo, Fedor & Lee, *Nat. Commun.* 2025, 16:5178 |
| URAT1 | 9IRY verinurad、9JE1 dotinurad、9J5X lingdolinurad(humanized rat) | 3.0–3.6 Å | Guo 2025；Cell Discov 2025；Fan & Lei *JACS Au* 2025, 5:1308 |
| **hOAT1** | apo / olmesartan+Cl / olmesartan+Br / **probenecid** | 3.33–3.88 Å | Jeon 等, *Structure* 2025, 33:1856 |
| **hOAT1** | apo / cidofovir 3.15 Å / glibenclamide 3.68 Å | 3.15–3.68 Å | Wu, Luo 等, *Sci. Adv.* 2025 |
| hOAT3 | **无实验结构** | — | 只能同源建模 |

这是真正的时间窗口：**hOAT1 的人源结构 2025 年才公开**，URAT1 与 OAT1 的结构性选择性分析在此之前不可能做。

---

## 3. 定量锚点：4 药 × 4 转运体

Taniguchi 等 *JPET* 2019, 371:162（同一实验体系，IC50 / µM）：

| 药物 | URAT1 | ABCG2 | OAT1 | OAT3 | **URAT1/OAT1 选择性比** |
|---|---|---|---|---|---|
| dotinurad | **0.0372** | 4.16 | 4.08 | 1.32 | **110×** |
| 苯溴马隆 | 0.190 | 0.289 | 3.14 | 0.967 | 16.5× |
| lesinurad | 30.0 | 26.4 | 6.99 | 1.07 | **0.23×**（反选择） |
| 丙磺舒 | 165 | 433 | 10.9 | 2.37 | **0.066×**（反选择） |

**选择性比跨越约 1600 倍，且方向在 dotinurad 与丙磺舒之间反转。** 这是一个有真实动态范围、可做回归与分类的目标量——比 469 个同骨架 IC50 有信息量得多，而且**不受标签病理影响**（这 4 个药的数据来自同一篇、同一体系）。

可扩展的补充数据：verinurad（URAT1 25–40 nM，Phe365/Ser35/Cys32 依赖）、lingdolinurad/ABP-671（URAT1EM ~70 nM，Phase 2b/3）、TD-3（9DKC）、AR882/pozdeutinurad（Phase 2）、CDER167 与 verinurad 类似物（URAT1/GLUT9 双抑）、Zhao 等 *J. Med. Chem.* 2020, 63:10829 系列。**整理成一份 URAT1-vs-抗靶选择性数据集本身就是论文的一个可交付资源。**

---

## 4. 主张与产品

**主张（不是方法学发明，是首次可行的分析 + 候选）**
> 利用 2025 年才同时具备的 URAT1 与 hOAT1 实验结构，建立 URAT1-over-OAT 选择性的结构判据，用已发表的 4 药 × 4 转运体选择性矩阵做回顾性校准，并据此提名具备 SURI 画像的候选分子；同时在计算上检验 Dai & Lee 提出的 Arg477 结合假设。

**产品**
1. 一份**候选分子表**（库筛 + 定向设计两路），每个含 URAT1 姿态、OAT1 对照姿态、选择性判据数值、MD 稳定性、ADMET；
2. 一份**URAT1/抗靶选择性数据集**（文献整理，可复用）；
3. **Arg477 假设的计算检验结果**（正或负都有意义，因为它是别人提出的公开假设，不是我们的退路）；
4. NLRP3 二级注释（沿用已跑通的结构门）。

---

## 5. 实施阶段

### S1 — 选择性数据集与结构对齐（低算力，先做）
- 整理 URAT1 / OAT1 / OAT3 / ABCG2 /（GLUT9）IC50 与选择性比，标注测定体系；
- URAT1（9DKA/9DKB/9DKC/9IRY/9JE1）与 hOAT1（probenecid / cidofovir / glibenclamide / apo）结合位点残基逐一对齐，列出**差异残基**；
- 输出：差异残基表 + 选择性数据集。
- 分辨率纪律：**URAT1 2.55–3.00 Å 优于 hOAT1 3.15–3.88 Å**，非对称质量必须声明；OAT1 侧侧链级结论置信度低于 URAT1 侧。

### S2 — Arg477 假设检验（别人提出的公开假设）
- 在已有共结构上量化各药与 Arg477 / Lys393 的结合，与效价、选择性比作图；
- 检验：Arg477 结合是否同时解释**效价**（verinurad/dotinurad 最强）与**选择性**（OAT1 侧对应位置是否不同）；
- 这一步只用晶体坐标 + 已发表 IC50，几乎不需要算力，且结论无论正负都可发表。

### S3 — 差分判据标定（核心验证）
- 同一配体分别对接 URAT1 与 hOAT1，取**差分**（不是各自的绝对分数）；
- 判据候选：差分 IFP、差异残基接触、口袋形状互补差；
- **在 S1 的选择性数据集上回顾性标定**：能否复现 dotinurad ≫ 苯溴马隆 > lesinurad > 丙磺舒 的选择性方向（含丙磺舒/lesinurad 的反选择）；
- 门：选择性方向排序正确，且与分子大小无强相关（Spearman |ρ| ≤ 0.4）。

### S4 — 候选获取（两路并行）
- **库筛路**：临床/再定位池（1588 → 8319）+ taosu 库，用 S3 判据取 URAT1 强 / OAT1 弱者；
- **设计路**：以 TD-3（9DKC，2.55 Å）、verinurad、lingdolinurad 化学型为母核，做**同系列**类似物设计，目标是增强 Arg477 结合并放大 URAT1/OAT1 差异。
  - 同系列（congeneric）是关键：相对自由能计算在**同系列小扰动**上才可靠，这正好避开了历史上跨化学型排序的失败模式。

### S5 — MD / 相对自由能确证
- 对照先行：苯溴马隆@9DKA、TD-3@9DKC、dotinurad@9JE1、lesinurad@9DKB、丙磺舒@hOAT1；
- 力场：Amber14SB / GAFF2 / POPC / TIP3P、PME、vdW 9 Å、2 fs、Parrinello–Rahman，生产 3 × 150 ns；
- 相对自由能（可选、同系列内）：膜蛋白 RBFE 已有验证记录（多个 GPCR 上 RMSE ≈ 0.80 kcal/mol、Spearman ρ ≈ 0.55，含完整脂双层与显式水），且在同源模型上表现与晶体结构相当；
- **先做回顾性**：在 §3 的已知药对上验证，再用于设计分子；
- 退路：若 RBFE 不收敛，只报 MD 接触占据与位移的**差分**，仍然给出选择性读出。

### S6 — NLRP3 二级注释
沿用已跑通的结构门（正例 10/10、背景 11/20、Fisher p = 0.0134；IFP Jaccard 参考 NP3-146 0.842），作为候选的**附加价值标注**，不作为主轴。

---

## 6. 风险与缓解（诚实清单）

| 风险 | 缓解 |
|---|---|
| hOAT3 无实验结构 | 由 hOAT1 同源建模，明确标注为模型、置信度低于 OAT1；OAT3 只作辅助 |
| ABCG2 是 ABC 家族、折叠不同 | 不做结构差分，仅作配体基/文献注释，或排除出主判据 |
| hOAT1 分辨率低于 URAT1 | 非对称质量声明；OAT1 侧只做骨架/口袋级判读 |
| RBFE 对带电配体 + 部分水化腔收敛困难 | 先回顾性验证；不收敛则退到 MD 差分读出 |
| 校准样本量（完整 4×4 只有 4 药） | 用 §3 的补充数据扩充；主张限于"方向正确 + 分类"，不在 n=4 上声称相关系数 |
| 设计分子需合成才能验证 | 明确标注为设计假设；库筛路同时提供可直接购买/已上市的候选 |
| 双靶身份变化 | 主轴改为 URAT1 选择性，NLRP3 降为二级注释——这是范围变更，已在此显式记录 |

## 7. 与历史资产的关系

| 资产 | 角色 |
|---|---|
| 冻结 P2 表、C1 全部产物 | 只读；作为"绝对幅值排序失败"的对照基线 |
| TrueDecoy / RandomDecoy | 仅用于 URAT1 侧召回检查，不再承担主判据 |
| 469 活性 + 80 弱活 | 建 URAT1 参考指纹；不再用于跨化学型排序 |
| NLRP3 结构门（p = 0.0134） | S6 二级注释，直接沿用 |
| C2 底物/反向阴离子排除 | 保留为候选安全过滤（吡嗪酸样 = 抗尿酸排泄风险） |
| C3 的 9DK 四联结构与突变体校准 | **保留并入 S1/S2**，只是不再作为论文主张 |
| URAT1 配体 ML（R² 0.53、命名药 2/4） | 已判失败，不复活 |

## 8. 禁止事项

- 不覆盖 `data/repurposing/p2/` 与 `data/campaigns/c1/`；
- 不把差分判据或对接分数称作 IC50 / 选择性比的预测值（只作方向与分类）；
- 不在 hOAT1 3.15–3.88 Å 模型上做侧链级断言；
- 不在 n = 4 的校准集上报相关系数；
- 不做跨化学型的绝对强弱排序（这是已失败的任务形态）；
- 相对自由能只在**同系列**内使用；
- 不声称双靶抑制剂；不把 substrate-like（吡嗪酸样）分子写进候选；
- 不声称"发明"了选择性计算——这是标准任务形态，新的是首次两边都有结构。

## 9. 参考锚点

- Suo, Fedor & Lee, *Nat. Commun.* 2025, **16:5178**（9DK9/9DKA/9DKB/9DKC；apo + 3 抑制剂，2.55–3.00 Å）
- Dai & Lee, *Cell Res.* 2024, **34:776–787**（三态、吡嗪酸三位点、gating hijack；**Arg477 设计假设**）
- Guo 等, *Nat. Commun.* 2025, **16:1512**（9IRW/9IRX/9IRY；F449/R447 跨态冲突）
- *Cell Discovery* 2025, **11**（四药 + 尿酸三态；**OAT 同源性导致选择性设计困难**；IC50 8 nM / 40 nM / ~200 nM / ~12 µM）
- Fan & Lei, *JACS Au* 2025, **5:1308–1319**（9J5X lingdolinurad；humanized rat URAT1）
- Jeon 等, *Structure* 2025, **33:1856–1866**（hOAT1 apo / olmesartan / **probenecid**，3.33–3.88 Å）
- Wu, Luo 等, *Sci. Adv.* 2025（hOAT1 + cidofovir 3.15 Å / glibenclamide 3.68 Å）
- Taniguchi 等, *JPET* 2019, **371:162**（4 药 × 4 转运体 IC50 矩阵；SURI 定义）
- Tan 等, *Sci. Rep.* 2017, **7:665**（定量突变体 fold-change；verinurad 依赖 Cys32/Ser35/Phe365/Ile481）
- Zhao 等, *J. Med. Chem.* 2020, **63:10829**；CDER167, *Acta Pharmacol. Sin.* 2021（同系列 SAR 与选择性数据）
- 膜蛋白 RBFE 验证：4 个 GPCR 上 RMSE ≈ 0.80 kcal/mol、Spearman ρ ≈ 0.55（FEP+，含脂双层）；同源模型上表现与晶体相当
