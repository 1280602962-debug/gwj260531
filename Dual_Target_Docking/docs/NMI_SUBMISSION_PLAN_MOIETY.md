# NMI 投稿规划：双药效团 / Passenger 污染范式

> **主线定调**  
> 单分子双靶配体含两个药效团；独立对接 + 整分子打分会把「乘客药效团」算进单口袋结合，系统误评。  
> 本文提出 **moiety-resolved dual-target docking evaluation**，配开放基准与非循环案例验证。  
> 朴素分数融合仅作 **基线对照**，不作封面创新。

总览见 [`PROJECT_MASTER_PLAN.md`](PROJECT_MASTER_PLAN.md)；诊断操作见 [`EGFR_HER2_DIAGNOSTIC_DEMO.md`](EGFR_HER2_DIAGNOSTIC_DEMO.md)。

---

## 1. 目标期刊与卖什么

| 项 | 内容 |
|----|------|
| 目标 | **Nature Machine Intelligence** |
| 备选 | JCIM / Chemical Science / Briefings in Bioinformatics（主结果不够「纠偏」时） |
| NMI 买的点 | **任务纠偏 + 可复现协议 + 开放基准**，不是新对接采样器 |
| 一句话 | Dual-target ligands are mis-scored by single-pocket whole-molecule docking because the second pharmacophore acts as a passenger; moiety-resolved scoring + per-target calibration restores dual-vs-single ranking |

### 1.1 与近年 NMI 药学相关文的对齐

| 对标模式 | 例子 | 本文对应 |
|----------|------|----------|
| 评测打假 / 资源 | VSDS-VD、CleanSplit | Dual-VSDS-Moiety：按药效团分解的双靶评测 |
| 姿态/分数可靠性 | PoseBusters 类精神 | passenger 污染诊断表 + 几何可行性门控 |
| 不宣称新引擎 | 多数 NMI 方法文复用已知零件 | GNINA/Vina + 自研 scoring protocol |

### 1.2 Cover letter 主张（可直接改写）

> We do not introduce a new docking sampler. We show that independent docking of dual-pharmacophore ligands systematically misattributes the passenger moiety to single-pocket scores, that naive cross-target fusion cannot fix this, and that moiety-resolved scoring with per-target calibration recovers dual-vs-single discrimination on a leakage-controlled public benchmark and an externally anchored NLRP3/JNK1 case study.

---

## 2. 论文主张（Claim ladder）

只承诺能用实验支撑的层级：

| 层级 | Claim | 必须交付的证据 |
|------|-------|----------------|
| C1 诊断 | 整分子对接分在双药效团分子上系统性抬升假双靶 / 错排 A-only | 诊断表：dual vs A-only vs B-only，whole-mol vs moiety |
| C2 机制 | passenger 体积/柔性/连接方式加重污染；fused vs linked 是严重度轴，非主因果 | 架构分层 + passenger 描述符相关 |
| C3 方法 | moiety-resolved + 分靶校准 + 几何门控 + shortfall 双靶决策优于整分子融合 | 主表：AUROC / EF / dual-vs-single |
| C4 资源 | Dual-VSDS-Moiety 公开 + 复现脚本 | 数据、划分、标签、协议 YAML |
| C5 案例 | NLRP3/JNK1 非循环：先外部锚点，再看实验室双靶分子 | 已知抑制剂/共晶 → 协议 QC → 私有分子解释 |

**不做的 claim：** 新对接引擎；协同预测；PROTAC 三元；细胞表型 = 双靶直接结合。

---

## 3. 文章结构（NMI 友好）

| 板块 | 篇幅感 | 内容 |
|------|--------|------|
| Abstract | 150–200 词 | 问题 → passenger → 方法一句 → 主结果数字 → 开放资源 |
| Intro | ~1.5 页 | 双靶药需求 → 现有 VS 把双靶当两次单靶 → 我们指出 passenger → 本文贡献 3–4 条 |
| Results 1 | 诊断 | 整分子打分失败模式；图：双药效团示意 + 分数污染示意 |
| Results 2 | 方法有效 | moiety vs whole-mol vs naive fusion；多靶点对 |
| Results 3 | 分层与消融 | fused/linked、校准有无、几何门控、shortfall vs 乘积 |
| Results 4 | 案例 | NLRP3/JNK1 外部锚点 + 实验室分子（holdout） |
| Discussion | 局限 | 药效团标注成本、NLRP3 表型≠结合、对接天花板 |
| Methods | 详 | 数据、对接协议、moiety 定义、校准、指标、泄漏控制 |
| Data/Code | 必 | GitHub + Zenodo；Dual-VSDS-Moiety |

图件建议（最少）：

1. **概念图**：一分子两药效团，对接 A 时 B 为 passenger  
2. **诊断主图**：whole-mol 抬高 A-only；moiety 纠正  
3. **基准主表/图**：三对公开靶点，指标对比  
4. **消融**  
5. **NLRP3/JNK1 案例**（结构锚点 + 分子映射，避免循环叙事）

---

## 4. 科学设计（实验骨架）

### 4.1 公开靶点对（已冻结）

见 [`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md) / `FROZEN_PUBLIC_PAIRS.yaml`：

1. **PIK3CA / mTOR** — 主规模  
2. **EGFR / HER2** — 两端共晶（TAK-285）  
3. **Mcl-1 / Bcl-xL** — PPI + 两端共晶（LC6）

标签：pChEMBL ≥ 6 = active；测得 < 6 = weak；未测 ≠ inactive。  
四类：dual / A-only / B-only / neither（硬负样本必须进主评测）。

### 4.2 方法流水线（投稿版最小完备）

```
分子 → 药效团/片段标注（人工+规则，可复现 schema）
     → 每靶独立对接（固定引擎与网格协议）
     → 整分子分（基线）vs moiety 分（主方法）
     → 分靶校准（rank / Z / decoy 分布）
     → 几何可行性（clash / linker 跨度 / PoseBusters 类检查）
     → shortfall 双靶决策（softmin / 短板）≠ 朴素平均
     → dual-vs-single 排序与 EF
```

### 4.3 必比基线

| 基线 | 作用 |
|------|------|
| 单靶对接分（A 或 B） | 说明不是「只会单靶」 |
| 朴素 mean / sum / min 融合 | 证明融合修不好 passenger |
| 仅校准的整分子融合 | 隔离「校准」与「moiety」贡献 |
| 随机 / 脚手架泄漏划分对照 | 证明不是泄漏 SOTA |

### 4.4 核心指标（主文只强调少数）

- **Dual-vs-single AUROC / AUPRC**（dual 对 A-only∪B-only）  
- **EF@1%、EF@5%**  
- **Calibration gap**：跨靶原始分 vs 校准后  
- **Passenger severity**：按 fused/linked 或 passenger 体积分层  
- 辅助：单靶回收率（证明 moiety 未毁掉真结合端）

### 4.5 NLRP3 / JNK1（非循环）

1. **外部锚点先**：JNK1 已知 ATP 抑制剂 + 共晶；NLRP3 以 MCC950/7PZC 等为锚（明确口袋假设）  
2. **协议 QC**：已知活性在正确端回收；passenger 诊断可复现  
3. **再解释实验室双靶分子**：moiety 是否分别复现两端特征  
4. **细胞数据角色**：时间盲测 holdout / 表型一致，**不写**「对接证明双靶结合」  
5. 若实验室读出仅为通路表型：Discussion 必须写清局限

---

## 5. 工作包与交付物（执行顺序）

| WP | 内容 | 交付物 | 过关标准 |
|----|------|--------|----------|
| WP0 | 主张冻结 + 旧文档对齐 | 本文档；MASTER 一句话更新 | 组内同意 claim ladder |
| WP1 | Moiety schema | 标注规范、示例分子、一致性检查脚本 | 三对靶点各 ≥N 标注分子 |
| WP2 | 公开数据导出 | ChEMBL 四类集、划分、scaffolds | 无 random-split 主结论 |
| WP3 | 对接协议固化 | YAML：引擎、盒子、exhaustiveness、种子 | 共晶自对接 RMSD 合格 |
| WP4 | 打分与校准 | whole-mol / moiety / fusion 实现 | 单元测试 + 小样本 sanity |
| WP5 | Dual-VSDS-Moiety | 基准包 + leaderboard 脚本 | 他人可一键复现主表 |
| WP6 | 主实验与消融 | 三对靶点全表 | C1+C3 数字达到预设门槛 |
| WP7 | NLRP3/JNK1 案例 | 锚点 QC + 私有 holdout 分析 | 叙事非循环；局限写清 |
| WP8 | 写作与投稿 | 主文、Methods、Code、Cover letter | 内审 + 预印本可选 |

**投稿前最低门槛（Go）：**

- [ ] 至少 **2–3 对**公开靶点上，moiety 在 dual-vs-single 上稳定优于 whole-mol 与 naive fusion  
- [ ] 一张 **不可辩驳的诊断图/表**（passenger 污染）  
- [ ] 泄漏控制划分 + 硬负样本  
- [ ] 代码与基准可公开  
- [ ] NLRP3/JNK1 **不**拖垮主张（锚点先行）

**No-Go → 改投 JCIM/Chem Sci：**

- 仅 1 对靶点有效；或 moiety 增益只来自标注泄漏；或细胞数据被迫写成结合证明

---

## 6. 写作与审稿预期

### 6.1 审稿人会打哪里

| 攻击点 | 预设答复 |
|--------|----------|
| 「不就是片段打分？」 | 贡献在双靶决策任务定义 + 污染诊断 + 校准/短板协议 + 基准，非整分子融合补丁 |
| 「药效团标注主观」 | schema 公开、双人一致率、规则优先、敏感性分析 |
| 「对接本身不可靠」 | 不宣称绝对亲和力；相对排序 + 共晶 QC + 几何门控 |
| 「NLRP3 口袋存疑」 | 明确假设与锚点；表型 ≠ 结合；案例降级为 illustration |
| 「和 VSDS-VD 比呢？」 | 单靶虚拟筛选评测 ≠ 双药效团乘客污染；任务不同 |

### 6.2 标题方向（暂定，可改）

- *Passenger pharmacophores confound dual-target docking scores*  
- *Moiety-resolved evaluation of dual-target ligands*  
- *Rethinking dual-target virtual screening under single-pocket docking*

### 6.3 作者叙事纪律

- Intro 前两段必须出现 **passenger / dual pharmacophore**，不要先写 fusion  
- Results 先诊断、后方法  
- 私有数据永远 holdout，不参与调参叙事  

---

## 7. 与仓库其他文档的关系

| 文档 | 关系 |
|------|------|
| [`PROJECT_MASTER_PLAN.md`](PROJECT_MASTER_PLAN.md) | 精简总览 |
| [`EGFR_HER2_DIAGNOSTIC_DEMO.md`](EGFR_HER2_DIAGNOSTIC_DEMO.md) | C1 诊断怎么跑 |
| [`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md) | 公开对已冻结 |
| [`NMI_REFERENCE_PAPER_PLAYBOOK.md`](NMI_REFERENCE_PAPER_PLAYBOOK.md) | 高分文流程对标 |
| [`DUAL_TARGET_SCORING_IMPLEMENTATION.md`](DUAL_TARGET_SCORING_IMPLEMENTATION.md) | 校准/短板组件 |
| 本文 | **现行 NMI 投稿执行规划** |

---

## 8. 近期两周执行清单（工程）

1. 冻结 moiety 标注 JSON schema（字段：target_A_moiety, target_B_moiety, linker, architecture）  
2. 从 EGFR/HER2 共晶双靶配体做 **最小可行诊断 demo**（whole-mol vs moiety）  
   → 操作细则：[`EGFR_HER2_DIAGNOSTIC_DEMO.md`](EGFR_HER2_DIAGNOSTIC_DEMO.md)  
3. 实现分靶 rank/Z 校准 + softmin 短板  
4. 导出一对靶点的 dual/A-only/B-only 评测表  
5. 若 demo 诊断表成立 → 扩到三对；若不成立 → 先修协议再谈写作  

---

*文档状态：现行投稿规划。主创新 = passenger 污染诊断 + moiety-resolved 双靶评测协议 + Dual-VSDS-Moiety，不是分数融合。*
