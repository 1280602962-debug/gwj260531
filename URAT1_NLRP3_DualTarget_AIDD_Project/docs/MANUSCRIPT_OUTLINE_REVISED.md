# 文章思路定稿（2026-07 修订）

> **首投**：*Journal of Molecular Modeling*（Subscription / 非 OA）  
> **定位**：痛风 URAT1–NLRP3 双节点下的 **临床库计算重定位漏斗**（假说生成，非新药发现）  
> **主推候选**：**canagliflozin**（模块 F）  
> **方法学案例**：**EGCG**（Pareto 盲筛 → PAINS/成药性降级）  
> **对接**：Glide XP @ **9DKB** + **7ALV**

---

## 0. 一句话 spine

在 ChEMBL 临床库上建立 **NLRP3 分类缩库 + 双靶 Glide XP + Pareto + 非对接审计** 的不对称漏斗；用 8973 独立轨证明 URAT1 应对接主导；漏斗能捞出 EGCG 类信号但会主动降级，并提名 **canagliflozin** 等更干净的双节点假说，辅以代表药 MD 作构象讨论——**不声称**双口袋直接抑制剂或临床推荐。

---

## 1. 已有数据 vs 缺口（按投稿必需度）

### 1.1 已完成（可直接写 Methods / Results）

| 模块 | 关键数字 / 产物 | 路径 |
|------|----------------|------|
| 临床库 | n=**8,319** | `data/repurposing/repurposing_manifest.csv` |
| NLRP3 ML | P≥0.5 → **1,588**；AUROC≈**0.89** | `data/repurposing/screening/` |
| 双靶对接合并 | n=**1,451** | `pareto_merged_scores.csv` |
| Pareto 前沿 | **6** 分子 | `pareto_shortlist.csv` |
| 8973 URAT1 回顾 | n=**8,973**；开发参考 AUC≈0.705 | `data/docking/8973_9DKB_with_manifest.csv` |
| 非对接 A–F | y-scramble、PAINS、ADMET、稳健性、提名 | `results/cheminformatics/`、`results/candidates/` |
| 提名首选 | **canagliflozin**（τ=90；clean+结构轴支撑） | `candidate_nomination_summary.json` |
| 主图 | Fig 2–4 PNG | `figures/generated/main/` |
| PPT 素材 | 结构/文献/示意图 | `figures/ppt_assets/` |

**Pareto 6**：SLV-334、LANPROSTON、LASALOCID、**EGCG**、FOSIGOTIFATOR、FOSRAVUCONAZOLE  
**模块 F Top（结构轴支撑）**：canagliflozin → caficrestat → nelutroctiv → SLV-334 → fulimetibant → PF-06952229

### 1.2 缺口（投稿前必须补）

| 优先级 | 缺口 | 用途 |
|--------|------|------|
| **P0** | 五组 MD（改：canagliflozin 替代 EGCG 作 lead） | Fig 5–6、Results R6 |
| **P0** | lesinurad@9DKB + 7ALV 共晶类似物 **redock RMSD** | Methods 可信度 |
| **P0** | Methods：Glide/AMBER 版本、网格、力场、时长 | 可复现性 |
| **P1** | 统一替换全文【待填】为 Glide XP 终版数字 | 一致性 |
| **P1** | Fig 1 流程图；Fig 7 非对接汇总 | 完整图集 |
| **P2** | canagliflozin 对接 2D 相互作用图（若尚未归档） | Discussion 结构段 |

### 1.3 明确不写 / 不混用

- 不写：发现首个双靶抑制剂、1+1>2 协同、临床推荐用药  
- 不混：8319 主筛、8973 仅 URAT1 回顾、benchmark/MD 平行表征  
- 不把：EGCG 写成主推；canagliflozin 写成 Pareto 前沿成员；URAT1 直接抑制剂（对 canagliflozin）

---

## 2. 核心叙事逻辑（五幕）

| 幕 | 问题 | 答案（已有数据） | 图 |
|----|------|------------------|-----|
| **1 为何不对称** | 为何不能双靶统一 ML？ | URAT1 回归 benchmark 2/4；NLRP3 分类 AUROC 0.89；两套库零重叠 | Fig 1 / SI |
| **2 缩库** | 临床库怎么压？ | 8319 → NLRP3 P≥0.5 → 1588；痛风药/colchicine 对照 | Fig 2 |
| **3 URAT1 辩护** | 为何 URAT1 对接主导？ | 8973 A vs D 富集；dotinurad 对接高/ML 低 | Fig 3 |
| **4 双轴筛选** | 谁在双轴上不被支配？ | 1451 → Pareto 6；lesinurad/colchicine 对照 | Fig 4 + 表1 |
| **5 审计与提名** | 6 人里跟谁？跟谁才干净？ | EGCG=盲筛案例但 PAINS 降级；模块 F → **canagliflozin** | Fig 7 + 表2 |
| **6 构象讨论** | pose 稳不稳？ | MD：benz/dot/canagliflozin@9DKB；MCC950/canagliflozin@7ALV | Fig 5–6 |

---

## 3. 两个候选的角色（全文统一口径）

| | **EGCG** | **Canagliflozin** |
|---|----------|-------------------|
| 如何进入故事 | Pareto 盲筛（6 人之一） | 模块 F 提名（非 Pareto 前沿） |
| 计算地位 | 方法学案例：漏斗能捞信号 | **主推计算假说** |
| 化学审计 | PAINS+Brenk，类药性差 → **降级** | 干净 + Lipinski/Veber |
| 文献 | 动物：URAT1 表达/MSU–NLRP3 | CANVAS post-hoc：SUA↓、痛风事件↓ |
| MD | **不进主文**（可选 SI 短轨迹） | **主文 lead**（9DKB + 7ALV） |
| 机制表述 | dual-node modulator（上游/表达） | pathway-adjacent；**非** lesinurad 式直接抑制 |
| 禁止说法 | 最终临床候选 | 已验证双口袋抑制剂 / 痛风适应症获批 |

---

## 4. 建议 MD 五组（与叙事对齐）

| # | 体系 | 角色 |
|---|------|------|
| 1 | benzbromarone @ 9DKB | URAT1 阳性基准 |
| 2 | dotinurad @ 9DKB | URAT1 阳性基准 |
| 3 | **canagliflozin @ 9DKB** | 主推假说（URAT1 结构轴） |
| 4 | MCC950 @ 7ALV | NLRP3 直接抑制基准 |
| 5 | **canagliflozin @ 7ALV** | 探索性 pose（与间接抗炎对照） |

可选 SI：`9DKB_EGCG` 短轨迹，仅说明盲筛代表的构象探索，不进主文结论。

---

## 5. 文章结构蓝图

### Title（英文）

*Clinical drug repurposing for gout-related URAT1 and NLRP3 nodes: an asymmetric NLRP3 machine-learning and Glide XP dual-target funnel with computational nomination of canagliflozin*

### Abstract（五句）

1. 痛风双轴；临床库缺可复现双节点漏斗。  
2. 8319 → NLRP3 ML → 双靶 XP → Pareto；8973 独立 URAT1 轨；模块 A–F。  
3. 1588 / 1451 / 6；EGCG 入前沿但 PAINS 降级；提名 canagliflozin。  
4. MD（完成后）：benchmark + canagliflozin pose 稳定性。  
5. 可检验假说；需 URAT1 摄取 + MSU–IL-1β 验证。

### Introduction（3 段）

1. 疾病：代谢轴 URAT1 vs 炎症轴 NLRP3；临床常分轴用药。  
2. 靶点与数据困境：转运体 vs 炎性小体；不对称可计算性。  
3. 空白 + 贡献：① 不对称漏斗 ② Pareto+审计（EGCG 降级）③ canagliflozin 提名 + MD。

### Methods

```
2.1 临床库 8319
2.2 NLRP3 分类（URAT1 回归仅回顾，不作主筛）
2.3 Glide XP @ 9DKB + 7ALV + redock
2.4 Pareto（S_U, S_N=max(ML,dock); min-su/sn=0）
2.5 8973 独立回顾轨
2.6 模块 A–F（PAINS/ADMET/适用域/稳健性/提名）
2.7 MD 五组（canagliflozin 为 lead）
```

### Results

| 节 | 内容 | 图/表 | 数据状态 |
|----|------|-------|----------|
| R1 | 数据不对称与漏斗设计 | Fig 1 | 待画 Fig1；逻辑已有 |
| R2 | NLRP3 ML：8319→1588 | Fig 2 | ✅ |
| R3 | 8973 URAT1 富集 | Fig 3 | ✅（终版数字核对） |
| R4 | 双靶 Pareto 6 + 对照药 | Fig 4 + 表1 | ✅ |
| R5 | 非对接审计：EGCG 降级；canagliflozin 提名 | Fig 7 + 表2 | ✅ 数字；Fig7 待拼 |
| R6 | MD：benz/dot/cana + MCC950/cana | Fig 5–6 | ❌ 待跑 |

### Discussion（固定 6 点）

1. 不对称双证据为何合理  
2. Pareto ≠ 临床最优（lesinurad、canagliflozin）  
3. EGCG：盲筛价值 + PAINS 降级（漏斗自我纠错）  
4. canagliflozin：通路邻近假说；对接≠直接抑制  
5. 对接用途边界：结构相容性排序 + 可检验 pose  
6. 局限：无湿实验、7ALV 非 MCC950 共晶、CANVAS 为 post-hoc、MD 参数待填

### 禁止出现在全文的表述

- “发现双靶抑制剂 / dual-pocket inhibitor discovery”  
- “协同 / synergistic / 1+1>2”  
- “真实世界已证实痛风疗效”（应写 CANVAS **事后分析**）  
- “EGCG 为最终推荐候选”  
- “canagliflozin 是 URAT1 直接抑制剂（lesinurad-like）”

---

## 6. 主图 / 表一览

| 编号 | 内容 | 状态 |
|------|------|------|
| Fig 1 | 双节点 + 漏斗（8319→1588→1451→Pareto→审计→MD） | 待画 |
| Fig 2 | NLRP3 ML | ✅ |
| Fig 3 | 8973 URAT1 | ✅ |
| Fig 4 | Pareto 散点（标 lesinurad、colchicine、EGCG） | ✅ |
| Fig 5 | URAT1 MD（benz、dot、**canagliflozin**） | 待 MD |
| Fig 6 | NLRP3 MD（MCC950、**canagliflozin**） | 待 MD |
| Fig 7 | 非对接：PAINS/ADMET/适用域/稳健性/提名 | 待拼图 |
| 表 1 | Pareto 6 + 文献裁决 | ✅ 可写 |
| 表 2 | 模块 F 干净候选（canagliflozin 居首） | ✅ |
| 表 3 | MD / MM-PBSA 汇总 | 待 MD |

---

## 7. 写作顺序（基于现有数据）

```
现在就能写
  1. Methods 2.1–2.6（对接/ML/Pareto/模块 F）
  2. Results R1–R5（含 EGCG 降级 + canagliflozin 提名）
  3. Discussion 1–4、6（局限可先写“MD 进行中”占位）
  4. Introduction + Abstract 骨架

等 MD + redock
  5. Methods 2.7 + Results R6 + Fig 5–6
  6. Abstract / Discussion 定稿
  7. 替换全部【待填】→ 投 JMM（Subscription）
```

---

## 8. 期刊

| 优先级 | 期刊 | 理由 |
|:---:|------|------|
| 1 | **Journal of Molecular Modeling** | 与现稿匹配；对接+MD 案例友好；非 OA |
| 2 | Chemical Biology & Drug Design | 疾病故事 + 重定位；可选非 OA |
| 3 | J. Comput.-Aided Mol. Des. | 方法学漏斗叙事 |

不上：JCIM / J Med Chem / Nat Commun（高度不够或需湿实验）。

---

## 9. 与旧稿差异（务必改口）

| 旧表述 | 新表述 |
|--------|--------|
| EGCG = 双节点代表 lead / MD 主对象 | EGCG = Pareto 盲筛案例，PAINS 降级 |
| MD：9DKB_EGCG + 7ALV_EGCG | MD：**canagliflozin** @ 9DKB + 7ALV |
| “真实世界证据”笼统说法 | CANVAS **事后分析** + 机制文献 |
| Pareto 6 = 最终跟进名单 | Pareto 6 = 计算前沿；跟进靠模块 F |

---

*本大纲取代早期“EGCG 作 MD lead”的默认叙事；与 `MANUSCRIPT_DRAFT_CN.md` 后续修订应对齐。*
