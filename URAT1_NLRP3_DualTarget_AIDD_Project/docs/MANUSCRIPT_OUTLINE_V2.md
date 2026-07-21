# 文章大纲 V2（2026-07-21）— 推翻旧结论后的重规划

> **状态**：取代 `MANUSCRIPT_OUTLINE_REVISED.md` 作为当前投稿规划主文档  
> **首投（非 OA / Hybrid）**：*Journal of Computer-Aided Molecular Design*（Springer Hybrid，可选不开 OA）  
> **备选**：*Chemical Biology & Drug Design*（Wiley；subscription + 可选 OA）  
> **三线**：《Journal of Molecular Graphics and Modelling》（Hybrid）  
> **明确不投（当前数据包）**：*Journal of Molecular Modeling*（MD/商业对接硬门槛）；Gold-OA-only 首投；JCIM / J. Med. Chem.（证据高度不够）

---

## 0. 推翻什么、保留什么

### 0.1 明确推翻 / 降级的旧结论

| 旧说法 | 新口径 |
|--------|--------|
| 主文 =「Glide XP 双靶漏斗发现 canagliflozin」 | 主文 =「**URAT1 对接协议如何选 + 不对称双节点漏斗如何用**」；canagliflozin 是审计后假说案例，不是发现叙事中心 |
| 8973 上 A vs D（近 RandomDecoy）富集 ≈ 已证明对接能力 | A vs D / RandomDecoy **不够**；必须以 **TrueDecoy** 为主判决，RandomDecoy 作对照 |
| 默认采用 Glide XP / 或默认 Vina 即可 | 必须先在 True/Random 上完成 **P0–P5（±Glide）协议筛选**，再锁定生产排序分 |
| redock top-1 RMSD&lt;2 Å「通过」或可放宽算通过 | **严格 top-1 未过关**（Vina/gnina≈4.5–5 Å）；采样可达近晶体姿 → 只证明盒子/采样，**不证明原生排序可靠** |
| 富集高 = 口袋/结合模式正确 | 富集 = **标签回收 / 排序协议考试**（9DKB 口袋假设下）；位点结论靠共晶 +（可选）RTMScore 选姿 / MD |
| 必须复刻 Hou 文全工具表才算按文章做 | 对齐 Gu/Hou **范式**（双诱饵、搜索∥打分、EF）；工具用开源 **Vina/gnina + RTMScore（±Glide）** |
| EGCG 可作主推 / 主文 MD lead | EGCG = **Pareto 盲筛 + PAINS 降级** 方法学负例；主推假说仍为更干净的 **canagliflozin**（机制边界写清） |
| 首投 JMM | **取消**；改 JCAMD（非 OA 路径） |
| 对接分数 ≈ 亲和力 / 双口袋直抑 | **禁止**；仅池内百分位与假说生成 |

### 0.2 仍然保留的硬资产

- 临床库 8319 → NLRP3 ML 缩库 → 双靶对接 → Pareto → 模块 A–F 审计提名  
- URAT1 不宜用回归 ML 做主排序（benchmark 回收不足）的不对称论证  
- EGCG 入前沿但被降级；canagliflozin 经清洁提名居前  
- 与 PLK1/NLRP3、湿法双靶路线的 Differentiation  
- 开源可复现脚本（Vina/gnina）+ True/Random 基准集  

---

## 1. 一句话 spine（新）

在痛风相关 **URAT1–NLRP3** 双节点上，先按 Gu/Hou 双诱饵框架为 URAT1（9DKB）选定开源可复现的对接排序协议，再将该协议嵌入 **NLRP3 ML 缩库 + 双靶对接 + Pareto + 成药性审计** 的不对称临床库漏斗，展示协议如何回收对照、降级 PAINS 型命中并提名 **canagliflozin** 类可检验假说——**不声称**已验证双口袋抑制剂或临床推荐。

---

## 2. 贡献分层（审稿人可读）

| 层级 | 贡献 | 证据 |
|------|------|------|
| **C1 方法（主）** | URAT1 对接协议筛选：TrueDecoy vs RandomDecoy；搜索与打分解耦 | P0–P5 表；AUC/EF；四药百分位；redock 双指标诚实报告 |
| **C2 系统** | 不对称双节点临床库漏斗（NLRP3-ML 缩库；URAT1 对接主导） | 8319→1588→对接合并→Pareto；8973/新基准对照 |
| **C3 审计** | Pareto ≠ 提名；EGCG 降级案例 | 模块 A–F；PAINS/ADMET；τ 提名表 |
| **C4 假说** | canagliflozin 双节点计算提名 + 构象讨论 | 对接百分位 +（P0）MD；机制边界 Discussion |

**Novelty 类型**：`new analysis/benchmark + new system/protocol`（非新打分函数、非新药发现）。

---

## 3. 非 OA 期刊定位

### 3.1 推荐排序

| 排名 | 期刊 | 模式 | 为何匹配 V2 |
|------|------|------|-------------|
| **#1** | **J. Comput.-Aided Mol. Des. (JCAMD)** | Hybrid → 可选非 OA | 协议/回顾验证友好；接受开源对接+有限 MD；奖励失败模式与负结果 |
| **#2** | **Chem. Biol. Drug Des. (CBDD)** | Subscription + 可选 OA | 疾病双节点与 canagliflozin 假说可读；需把 claim 压在假说生成 |
| **#3** | **J. Mol. Graph. Model.** | Hybrid | 偏建模图示；作备胎 |
| ❌ | J. Mol. Model. | Hybrid 但 Aims 硬 | ≥500 ns / 多构象 / 商业对接 discouraged → 当前包 desk-reject 风险高 |
| ❌ | JCIM / JMC / Nat Commun | — | 高度与湿实验预期超出数据包 |
| ❌ | Gold-OA-only 首投 | — | 违反非 OA 偏好 |

### 3.2 JCAMD 写作契约（投稿合同）

- **卖**：可复现协议、True vs Random 差异、排序失败（redock top-1）诚实写、Pareto≠Module F  
- **不卖**：首个双靶药、分数=Ki、口袋已由富集证明、未过门控却宣称 pose 可靠  
- **证据底线**：协议筛选主表 + redock 双指标 + 漏斗数字 +（强烈建议）50–100 ns 级代表 MD；开源参数进 SI  
- **APC**：走 subscription / 非 OA 路径（作者选择不开 Open Choice 即可）

### 3.3 标题方向（英）

*Protocol selection for URAT1 docking under TrueDecoy and RandomDecoy benchmarks and an asymmetric NLRP3-informed clinical-library funnel for gout dual-node hypothesis generation*

（可缩短；避免 “discovery of dual inhibitors / Glide XP lead compound” 旧标题气质。）

---

## 4. 数据与实验规划（按文章幕次）

### 4.1 Part A — URAT1 协议筛选（新主柱）

**输入**

- `data/benchmarks/urat1_true_decoy/true_decoy_benchmark.csv`  
- `data/benchmarks/urat1_true_decoy/random_decoy_benchmark.csv`  
- 受体：9DKB（同一搜索盒）

**候选协议（锁定）**

| ID | Pose | Score | 角色 |
|----|------|-------|------|
| P1 | Vina | Vina affinity | 物理基线 |
| P2 | gnina | CNNaffinity | AI 主候选 |
| P3 | gnina | gnina affinity (kcal) | 读出对照 |
| P4 | Vina ensemble | RTMScore | 搜索∥打分 |
| P5 | gnina ensemble | RTMScore | 同上 |
| P0 | gnina | CNNscore | 负对照 |
| P6–P8 | Glide ± RTMScore | （有许可） | 与 Hou 物理强项对照 |

**门控与选姿**

- Redock（lesinurad）：报告 **Top-1 RMSD** + **Best-in-ensemble RMSD** + **RTMScore 选姿 RMSD**  
- 严格门控：Top-1 ≤2 Å；未过关则协议不得宣传为 pose-accurate，仍可用于 enrichment 比较  
- 结构用姿 = RTMScore `struct_pose` 或晶体坐标；富集用姿 = 各协议 `rank_pose`

**主指标**：TrueDecoy EF@1%/5%、AUC；RandomDecoy 同指标作否决；四药百分位平局。

**产出**：选定 **生产协议 Π\***（可含 hierarchical：粗排引擎 + RTMScore 精排）。

### 4.2 Part B — 不对称漏斗（应用 Π\*）

1. 临床库 8319 → NLRP3 ML（P≥0.5）→ 1588  
2. Π\* @ 9DKB +（同引擎或既定）@ 7ALV  
3. Pareto（S_U, S_N）→ 短名单  
4. 模块 A–F：过滤器、适用域、稳健性 → **EGCG 降级**；**canagliflozin 提名**  
5. 可选：用 Π\* 重算/对照旧 8973 或 Glide 历史表（SI），证明协议迁移

### 4.3 Part C — 构象讨论（辅助，非发现证明）

- 五组 MD（若算力允许）：benz/dot/@9DKB；MCC950@7ALV；**canagliflozin**@9DKB+7ALV  
- 时长目标：JCAMD 可接受的 **50–100 ns** 级（不冲 JMM 500 ns）  
- 初始姿：优先晶体或 RTMScore struct_pose，**不用**已失败的原生 top-1 装可靠  

---

## 5. 正文结构蓝图

### Abstract（五句）

1. 痛风双节点；临床库缺尊重数据不对称且对接协议未经双诱饵筛选的可复现流程。  
2. 先在 URAT1 True/Random 上比较开源对接/重打分协议并诚实报告 redock 排序失败。  
3. 将选定协议嵌入 NLRP3-ML + 双靶对接 + Pareto + 审计漏斗。  
4. 展示对照药行为、EGCG 降级与 canagliflozin 提名。  
5. 输出可检验假说；需摄取/IL-1β 实验，不作临床推荐。

### Introduction（4 段）

1. 疾病双轴（URAT1 代谢 / NLRP3 炎症）与临床分轴治疗。  
2. 计算挑战：转运体结构新但数据噪声；对接富集 ≠ 位点证明；Random 诱饵易虚高。  
3. 相关工作：Gu/Hou VSDS-VD；湿法双靶；PLK1/NLRP3 不对称 VS（**必须 Differentiation**）。  
4. 本文贡献 C1–C4；claim 边界。

### Methods

```
2.1 数据集：临床库；URAT1 活性集；TrueDecoy/RandomDecoy 构建规则
2.2 受体/配体准备；9DKB/7ALV 搜索盒；redock 定义（双指标）
2.3 候选协议 P0–P5（±Glide）；RTMScore ensemble 重打分
2.4 富集指标与协议选优规则（预先锁定）
2.5 NLRP3 ML 缩库；双靶对接与百分位；Pareto
2.6 化学过滤器、适用域、y-scramble、提名规则（模块 A–F）
2.7 MD（若纳入主文）：体系、力场、时长、分析量
```

### Results

| 节 | 内容 | 图/表 |
|----|------|-------|
| R1 | True vs Random 协议比较；选定 Π\* | Fig 2；表 1 |
| R2 | Redock：top-1 失败 vs ensemble/RTMScore | Fig 3 或表 S |
| R3 | 不对称必要性（URAT1 ML 不足；NLRP3 ML 可用） | Fig 1 流程 + SI |
| R4 | 漏斗压缩与对照药 | Fig 4 |
| R5 | Pareto 6 与 EGCG 降级；Module F → canagliflozin | Fig 5；表 2–3 |
| R6 | MD 构象讨论（若完成） | Fig 6–7 |

### Discussion

1. 为何 TrueDecoy 改变协议选择（相对旧 A vs D）。  
2. Redock 排序失败对「pose 叙事」的限制；富集仍可指导排序。  
3. 与 Gu/Hou：范式对齐、工具裁剪、hierarchical 启示。  
4. EGCG vs canagliflozin：审计改变提名。  
5. 机制边界与实验证伪路径。  
6. Limitations：诱饵池规模、开源 vs Glide、单态口袋、无湿实验。

### Conclusions

协议筛选 + 不对称漏斗 → 可审计假说列表；非双靶药发现。

---

## 6. 图件规划（主文）

| 图 | 内容 |
|----|------|
| Fig 1 | 总流程：协议筛选 → 漏斗 → 审计提名 |
| Fig 2 | TrueDecoy vs RandomDecoy：各协议 EF/AUC |
| Fig 3 | lesinurad redock：RMSD 分布 / top-1 vs best / RTMScore |
| Fig 4 | 临床库漏斗数字与对照 |
| Fig 5 | Pareto 与 Module F 分离（EGCG 降级箭头） |
| Fig 6 | canagliflozin / 基准药 MD（若有） |

表：协议定义；选优结果；Pareto 6；提名 Top；redock 数值。

---

## 7. 投稿前 P0 / P1 清单

### P0（无则不投 JCAMD）

- [ ] True/Random 上 P0–P5 跑完并出主表  
- [ ] RTMScore ensemble 补齐；`rank_pose` / `struct_pose` 分离  
- [ ] lesinurad redock 双指标写入 Methods/Results（诚实）  
- [ ] 锁定 Π\* 并用于（或对照）漏斗 URAT1 轴  
- [ ] Claim 全文清除：双靶药发现 / 富集=位点证明 / top-1 已验证  

### P1（强烈建议）

- [ ] 五组 MD 或至少 canagliflozin + 1–2 基准  
- [ ] Differentiation vs PLK1/NLRP3 写入 Introduction  
- [ ] 开源参数、SMILES、分数 SI 打包  
- [ ] 中英稿数字与 Π\* 对齐  

### P2（加分）

- [ ] 有许可则加 Glide 对照行  
- [ ] EquiScore 可选  
- [ ] 扩大 decoy 池逼近 1:30–50 后敏感性分析  

---

## 8. 与旧大纲文件关系

| 文件 | 状态 |
|------|------|
| **`MANUSCRIPT_OUTLINE_V2.md`（本文件）** | **当前主规划** |
| `MANUSCRIPT_OUTLINE_REVISED.md` | 归档参考（Glide/JMM 气质；结论已部分推翻） |
| `MANUSCRIPT_OUTLINE_BENCHMARK.md` / `FAST_9DKB` | 仅作分轨素材，不再当主文 spine |
| `docs/paper_spine_ars_analysis/` | 期刊非 OA 论证仍有效；贡献句按 V2 的 C1 升权修订 |

---

## 9. 作者一句话备忘

**先选对接协议（TrueDecoy 主判），再用协议跑漏斗，用审计改提名；redock 没过 top-1 就别吹 pose；期刊走 JCAMD 非 OA，不走 JMM。**
