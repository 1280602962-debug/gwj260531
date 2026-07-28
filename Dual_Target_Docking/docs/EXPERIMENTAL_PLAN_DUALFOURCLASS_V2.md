# DualFourClass 实验总规划（v2）— 探索 / 开发 / 验证分离

> **本文件取代** [`EXPERIMENTAL_PLAN_DUALFOURCLASS_V1.md`](EXPERIMENTAL_PLAN_DUALFOURCLASS_V1.md) 作为权威实验规划。  
> **⚠ 2026-07-28 修订：** 红队审计发现本 v2 修的是统计纪律，但**测量效度**问题未修（主指标自相抵消、标签贴阈值、体积基线未被超过、prep 未统一）。  
> **Track B 之前必须先过 Stage M（测量审计）** — 见 [`PLAN_V2_REDTEAM_AND_REDESIGN.md`](PLAN_V2_REDTEAM_AND_REDESIGN.md)，其 §4 的 6 项重设计对本文 §2 / §5 / §7 具有优先效力。  
> 约束：可大批对接；**无湿实验**。  
> 期刊预期：JCIM / J. Cheminform.（不做 NMI 默认）。  
> 日期：2026-07-28  
> 触发：EGFR/HER2 panel120 **S1 No-Go** + 对 v1「小 panel 探索 → 冻结 → 同 pair 扩样确认」范式的方法论纠偏。

---

## 0. 一句话

交付 **DualFourClass-Bench + 可复现决策协议**：在多对公开靶上，用**事先列死的候选臂**与**跨 pair 选规则、留出验证**的流程，报告相对 `vina_mean` 的 Dual vs 硬负例判别能力——或诚实降级为诊断/基准论文。

**不是：** 新采样器、乘客/moiety 主线、在已坏臂上叠门控补丁、用「同 pair 更大 N」假装独立确认。

---

## 1. 为什么需要 v2（方法论失败，不只是 S1 阴性）

### 1.1 v1 实际做了什么

```
EGFR panel40 探索 → 见 rtm_min_z 好看 → 冻结 YAML → 同 pair 扩至 ~110 求显著
```

### 1.2 这为什么不科学

| 问题 | 具体表现 |
|------|----------|
| **探索与确认混同** | 规则在 panel40 上被选出；panel120 不是独立确认，而是同化学空间扩样 |
| **功效假设依赖效应量外推** | `POWER_ESTIMATION_V0` 假设 Δ≈0.17 可外推；新 70 上 Δ≈−0.19，假设被证伪 |
| **同 pair 扩样 ≠ 跨 pair 外推** | 旧 40 仍正向、新 70 显著负向 →「协议在该系列上成立」不成立 |
| **配体准备混杂** | 旧 LigPrep vs 新 RDKit+meeko；未先排除 confound 就判方法 |
| **阴性后仍易滑向补丁堆叠** | shortfall/clash/ML 在小面板已阴性或未预注册，不得再拧 |

### 1.3 已冻结的科学事实（仍有效）

- **C1（诊断）：** 朴素 `vina_mean` 在 EGFR/HER2 四类任务上接近随机 / 硬负例污染重。  
- **C2（机制线索）：** cognate 上 Top-1 姿态可错；RTM 可改姿态排序（panel40 证据）。  
- **C3（未成立）：** 不得再宣称 `rtm_min_z` 为通用主决策臂。  
- **C4（未成立）：** 跨 pair 协议增益未通过正确验证设计。  
- **分型资产：** T1/T2/T5 failure typology；决策消融阴性；两对靶 bootstrap CI。

这些资产进入 **Track A（诊断）** 或作为 **Track B 的探索池**，**不得**再当作「已冻结协议已确认」。

---

## 2. 科学原则（所有后续实验必须遵守）

### P1 — 三池分离

| 池 | 用途 | 禁止 |
|----|------|------|
| **Exploration** | 看现象、做分型、形成假设、列候选臂 | 用本池点估计定终身协议 |
| **Development** | 在预注册候选臂中做 **leave-one-pair-out（LOTO）** 选规则 | 看完 holdout 再改臂或改阈值 |
| **Confirmation** | 唯一选出的协议 **一次性** 在预留 pair / 预留化学块上报告 | 迭代调参后再报「最终」数字 |

### P2 — 候选臂预注册

在 Development 开工前，把候选决策臂写进 YAML（见 §5）。  
允许的事后分析：**仅**探索性、明确标注 *post hoc*，不得改主结论。

### P3 — 规则选择用跨 pair，不用同 pair 调优

- 选规则的目标函数：例如 LOTO 下各 fold 的 `ΔAUROC(arm − vina_mean)` 的 **中位数**，或「正向 fold 数」。  
- **禁止**在单个好看的 pair 上选臂再外推。  
- 同 pair 内可用 nested CV 估方差，但 **主选规则必须跨 pair**。

### P4 — 同 pair 扩样的合法用途

仅允许用于：

1. 估计该 pair 内效应量与 CI（描述性）；  
2. 检查 chemotype / prep 混杂；  
3. 功效复估（*prospective*，不得回溯改已报告结论）。

**不允许**作为「协议已确认」的证据。

### P5 — 统一配体准备与结构冻结

- Development / Confirmation 全库：**同一 prep 管线**（推荐全程 RDKit+meeko *或* 全程 LigPrep，二选一写死）。  
- 已有 LigPrep 旧姿态：要么重 prep 重对接，要么仅进 Exploration，并在 Limitations 标明。  
- 每对靶：结构冻结文件 + cognate QC（对特异 exhaustiveness）在面板对接前完成。

### P6 — 主标签与分层

- 主标签：**仅** Dual / A_only / B_only / neither（四类活性）。  
- 架构粗分层 `compact_ATP` / `clear_linker` / `unknown`：**只报告，不进主分数门控**。  
- 骨架配额：同 Murcko ≤3–5 / 面板；硬负例必须足额。  
- 归因 L1/L2/L3 与 ML 指纹消融：同 v1 D1/D2，仍强制。

### P7 — 统计报告规范

> **修订（Stage M / M1）：** 池化 `Dual vs A∪B` 已被证明会抵消自身信号（EGFR/HER2：D/A=0.69 与 D/B=0.31 平均成 0.52）。  
> **新主指标 = 方向分解**：`AUROC(dual vs A_only)` 与 `AUROC(dual vs B_only)` **分别报告**；需要单一数时取两者最小值或平均，且必须同时给出分项。  
> 平凡基线（重原子数 / MW / cLogP）为**必报对照**。详见 [`PLAN_V2_REDTEAM_AND_REDESIGN.md`](PLAN_V2_REDTEAM_AND_REDESIGN.md) §4.1、§4.3。

- 主指标（历史写法，已被上述修订取代）：AUROC Dual vs A∪B（或 Dual vs rest；**每对靶预先写死一种**）。  
- 并列：Top10 dual 数、Top10 hardneg（A_only+B_only）数。  
- 不确定性：配对 ligand-bootstrap 95% CI（B≥2000）。  
- 跨 pair：森林图 +（若 K≥4）随机效应汇总 *或* 仅报 LOTO 分布，不假装独立同分布强假设。  
- **禁止**裸点估计当显著。

---

## 3. 双轨产品（先选轨，再开工）

两条轨可串联：先完成 Track A 写作骨架，再决定是否投入 Track B 算力。

### Track A — 诊断 / 基准论文（默认保底）

**卖点：** 任务定义（四类硬负例）+ 朴素融合失败（C1）+ 姿态/重打分必要性线索（C2）+ **同系列扩样失败**（panel40≠panel120）+ failure typology T1/T2/T5 + 诚实 CI。

**Claim 天花板：** 诊断与可复现评测框架；**不**声称通用决策尺子已验证。

**还需补（低对接成本）：**

1. 论文主图：C1 + 旧40/新70 分裂 + T2/T5 对照。  
2. 可选：~20 新配体 LigPrep vs RDKit prep 对照（排除/量化 confound）。  
3. 打包：现有两对靶面板 + 协议草稿 YAML（标注 *exploratory*）+ 脚本。

**目标刊：** J. Cheminform. / Digital Discovery / JCIM short。

### Track B — 方法 / 多 pair 嵌套验证（冲 JCIM full）

**卖点：** 跨 pair 可迁移的 DualFourClass 决策协议；用 LOTO 选规则 + 一次确认，而不是单 pair 调参。

**Claim 天花板：** 「在预注册臂与 LOTO 流程下，选出的协议相对 `vina_mean` 在 held-out pair 上 …」——写到数据允许的强度为止。

**最低数据规模（方法轨）：**

| 项 | 最低要求 | 说明 |
|----|----------|------|
| 靶对 K | **≥4**（**修订：目标 6–8**） | 观测 per-pair ΔAUROC SE≈0.064–0.10；K=4 + 单 pair 确认几乎必然「点估计正、CI 含 0」 |
| 每对 N | **≥80–120**（四类配额） | 硬负例每类建议 ≥15–20 |
| Prep | 全库统一 | 旧 EGFR40 LigPrep 姿态不直接混入 Development |
| 引擎 | Vina 必做；**≥1** 第二引擎（GNINA 优先）作预注册臂 | 证明不单吃 Vina+RTM |
| 选规则 | LOTO 中位数 / 正向折数 | 预注册目标函数 |
| 确认 | **主估计量 = 跨 pair 汇总 Δ**；预留 pair 作复现检查 | 单 pair holdout 功效不足，不作唯一判据 |

**目标刊：** JCIM full（仍无湿实验；靠新任务+基准+正确验证设计对冲）。

### 选轨决策树

```
S1 已 No-Go
    ├─ 近期要成文、算力紧          → Track A（立即）
    ├─ 愿投多 pair 大批对接        → Track B（按下节阶段走）
    └─ 想先低成本排除 prep 混杂    → 先做 A 的 prep 对照，再决定 B
```

**明确禁止的第三条路：** 同一套 Vina→`rtm_min_z` 再扩 EGFR 赌显著；或在 RTM 臂上叠 shortfall/clash/ML 补丁当「新方法」。

---

## 4. Track B 阶段路线图

> **前置门（Stage M，2026-07-28 新增）：** M1 方向分解 / M2 margin 标签+阈值敏感性+噪声天花板 / M3 平凡基线组 / M4 统一 prep 重跑 / M5 清洗预注册清单。  
> **M1–M5 未全过，不得进入 B0。** Stage M 的 No-Go 条件：margin 标签 + 方向分解 + 体积对照下无任何臂在 ≥2 对靶上超过体积基线 → 直接收 Track A。

```
B0  预注册：候选臂 YAML + 选规则目标函数 + 靶对名单角色（dev/holdout）
  ↓
B1  靶对与面板建设（统一 prep；每对结构冻结 + cognate QC）
  ↓
B2  全库对接 + 预注册臂打分（不调参）
  ↓
B3  Development：LOTO 选唯一主臂（可并列报告次优作敏感性）
  ↓
B4  Confirmation：holdout 一次性评估 → Go / 降级
  ↓
B5  诊断层与归因（L1 全表；L2 配额分型；可选 ML 头仅 LOTO）
  ↓
B6  打包 DualFourClass-Bench + 写作投稿
```

任一关卡不达标：降级 Track A，不硬冲 claim。

---

## 5. 预注册：候选决策臂与选规则协议

> 细节清单可拆文件维护；**开工前必须冻结一版**。建议路径：`data/protocols/CANDIDATE_ARMS_V0.yaml`（尚未创建则在 B0 创建）。

### 5.1 基线（必报）

| ID | 定义 |
|----|------|
| `vina_mean` | 两端 Vina top1（或协议写死的 top-K 代表分）均值 |
| `vina_min` | 两端 Vina 代表分取 min（弱端短板） |

### 5.2 重打分族（预注册，不得看完再发明）

| ID | 定义 |
|----|------|
| `rtm_mean` | 两端 RTM best-of-K 均值 |
| `rtm_min` | 两端 RTM best-of-K 取 min |
| `rtm_min_z` | 分靶 z 校准后 min（*若进 Development，校准参数只能在 LOTO 训练折内估计*） |

### 5.3 第二引擎 / 共识族（Track B 强烈建议）

| ID | 定义 |
|----|------|
| `gnina_mean` / `gnina_min` | 同面板 GNINA（或等价） |
| `rank_consensus_*` | 预定义的 rank 聚合（如 mean rank / min rank）；**公式写死** |

### 5.4 明确排除出主臂竞赛的（可作 *post hoc* 诊断）

- 在 Exploration 小面板上拧出的 shortfall / clash 阈值门控  
- `warning_flags` / 架构标签灌进分数  
- 未 LOTO 的任意 ML 头  

### 5.5 LOTO 选规则（写死）

**输入：** Development 集中的 K_dev 个靶对（建议 K_dev≥3；另 1 对进 Confirmation）。  

**对每个候选臂 a：**

1. for fold i = 1..K_dev：在除 pair_i 外的对上**不调参**（或仅允许训练折内估计的校准参数），在 pair_i 上算 `Δ_i = AUROC_a − AUROC_vina_mean`。  
2. 汇总：`median(Δ)` 为主；次要：`#(Δ_i>0)`。  
3. **选出唯一主臂** `a*` = median(Δ) 最大者；若并列，取更简单者（奥卡姆：少校准 < 多校准 < 共识 < ML）。  

**Confirmation：** 仅 `a*` vs `vina_mean`（及预注册的 1–2 个敏感性对照）在 holdout pair 上评一次。

### 5.6 Confirmation Go / No-Go

| 结果 | 判定 |
|------|------|
| holdout `Δ` 的 CI 排除 0 且点估计正向 | **Go**：可写「LOTO 选出的协议在 held-out pair 上 …」 |
| 点估计正向但 CI 含 0 | **弱 Go**：报告效应与不确定性；降低外推措辞 |
| 点估计负向或 LOTO 中位数 ≤0 | **No-Go**：主文降级 Track A；方法章节报告阴性 |

---

## 6. 靶对角色建议（基于已有公开对审计）

依据 [`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md)：

| 角色 | 靶对 | 备注 |
|------|------|------|
| Exploration（已完成） | EGFR/HER2 panel40/120；PIK3CA/mTOR panel48 | 仅假设生成与分型；**不**再当确认 |
| Development 候选 | EGFR/HER2（**统一 prep 重建或重对接子集**）；PIK3CA/mTOR（扩或重建至配额）；AChE/BChE | 异质酶对加分 |
| Development / 异质加分 | Mcl-1/Bcl-xL | PPI 沟槽难，但对「非激酶」外推有价值 |
| Confirmation（预留） | **事先指定 1 对**（建议：未参与选臂的 AChE/BChE *或* Mcl-1/Bcl-xL） | 选定后禁止改换 |

**NLRP3/JNK1：** 若存在私有数据，可作最终外部 holdout；**不得**在选臂阶段使用。

面板建设硬规则：

1. 四类标签规则与 ChEMBL 审计一致（已测才进 A-only/B-only；未测≠阴）。  
2. 共价药 / PROTAC / 金属螯合特例：另表或排除，写入面板 README。  
3. 每对 N 与配额在对接前冻结 CSV。

---

## 7. Track A 与 Track B 的共用工作包（立即可做）

无论最终选哪条轨，下列工作都提高科学可信度：

| # | 工作包 | 目的 | 对接量 |
|---|--------|------|--------|
| W1 | 写作：S1 分裂图 + typology + C1 | Track A 主干 | 0 |
| W2 | Prep 对照：~20 新配体 LigPrep vs RDKit | 量化 panel120 confound | 小 |
| W3 | 预注册 YAML：候选臂 + LOTO 目标函数 | 防止事后发明 | 0 |
| W4 | 第三/四对靶选型与结构冻结 | 为 Track B 铺路 | 0→中 |
| W5 | 失败分型库扩样（配额，不调分） | 诊断深度 | 0–小 |

**优先级建议：** W1 + W3 立刻完成 → W2（若怀疑 prep）→ 再批准 Track B 的 W4 大批对接。

---

## 8. 与旧文档 / 旧关卡的关系

| 文档或关卡 | v2 地位 |
|------------|---------|
| `EXPERIMENTAL_PLAN_DUALFOURCLASS_V1.md` | **废止为权威**；保留作历史；S1 No-Go 记录仍有效 |
| v1 的 S1「扩 EGFR 求显著」 | **已关闭**；不得重启同协议赌局 |
| v1 的 S2/S3 | 由 Track B 的 B3–B5 **重新定义**（必须 LOTO，不得「≥3 对点估计正向即可宣称外推」） |
| `JCIM_GAP_TWO_SCENARIOS.md` | 方案甲 ≈ Track A；方案乙须按 v2 嵌套验证重写，**不可**再「先扩 EGFR 再加靶」 |
| `RESEARCH_DIRECTION_REFREEZE.md` | 方向（架构无关决策尺子）仍有效；claim 强度按本 v2 降级/重验 |
| `POWER_ESTIMATION_V0.md` | 仅作历史；新功效须在 **预注册臂 + 多 pair** 设定下重估 |
| `FAILURE_TYPOLOGY_V0.md` / decision ablation / bootstrap | Exploration 资产，继续引用 |

---

## 9. 每阶段强制检查清单

进入下一阶段前：

- [ ] 当前分析属于 Exploration / Development / Confirmation 哪一池，文中已标明  
- [ ] 主指标仅四类活性，未按架构过滤  
- [ ] 候选臂在打分前已写入预注册清单  
- [ ] Development 选臂使用 LOTO（或等价跨 pair），非单 pair 调优  
- [ ] Confirmation 只评一次，无阈值回拧  
- [ ] Prep 管线在 Development/Confirmation 内统一  
- [ ] 报告含 bootstrap CI；跨 pair 含森林图或 LOTO 分布  
- [ ] ML（若有）含 pair-held-out + 分数-only vs +指纹消融  
- [ ] 未把 warning flags / 架构标签灌进主分数  
- [ ] 未把 S1 No-Go 的 `rtm_min_z` 当作已验证主臂  

---

## 10. 资源与分工

| 工作 | 谁做 |
|------|------|
| 大批对接 / LigPrep / GNINA 等引擎 | 本地 agent |
| 预注册 YAML、面板规则、QC、LOTO 分析、CI、文档 | 云端 / 文档 agent |
| 架构粗标、分型个案 | 半自动 + 人工抽查 |
| 湿实验 | **不做** |

---

## 11. 当前状态与下一步（2026-07-28）

| 项 | 状态 |
|----|------|
| 方向（决策尺子，非 moiety） | 冻结 |
| v1 S1 EGFR120 | **No-Go** |
| `rtm_min_z` 通用主臂 | **撤销** |
| 权威规划 | **本文 v2** |
| Track 选择 | 待定（默认保底 Track A；Track B 需明确批准算力） |

**立即下一步（文档/分析，默认）：**

1. ~~冻结本文为权威；更新总览与红队引用。~~  
2. ~~创建预注册臂清单草稿 `CANDIDATE_ARMS_V0.yaml`（B0）。~~  
3. **执行 STEP0：** [`AGENT_COMMAND_STEP0_TRACK_A_STARTER.md`](AGENT_COMMAND_STEP0_TRACK_A_STARTER.md)（W1+W3）。  
4. 可选 W2 prep 对照设计一页 SOP（另开命令）。  
5. **在用户批准 Track B 之前，不启动新的同协议 EGFR 扩样对接。**
