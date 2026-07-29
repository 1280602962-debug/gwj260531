# JCIM 加厚规划 v1（评测文 · 尽可能多补实验）

> 日期：2026-07-29  
> 前提：K=4 对接 + CI 分析包已完成；claim = **Evaluation / Benchmark Article**，不卖通用决策臂。  
> 用户目标：提高投中把握；认为当前厚度不够撑 JCIM full paper → **在合理范围内尽量补实验**。  
> 权威约束：`CLAIM_CEILING.md`；禁止再扩 EGFR 赌显著；禁止 Track B 选臂胜利叙事。

---

## 0. 先把话说死（否则补实验会白做）

### 0.1 JCIM full paper 缺的不是「再算一轮同样的 AUROC」

对照近三年 JCIM / J. Cheminform. 对接评测文，一篇能站住的评测/基准文通常至少占两项：

| 厚度轴 | 你们现在 | JCIM 期望感 | 本规划怎么补 |
|--------|----------|-------------|--------------|
| A. 新任务 / 新数据集可复用 | 有雏形，K=4 小 | 要别人能下载、能提交、有协议 | Zenodo + 标准 split + 参考实现 |
| B. 多条件稳健性 | prep/θ/CI 初具 | 要结构、化学型、聚合、噪声、缺失 | Wave 1–2 分析 + 结构敏感对接 |
| C. 规模感（对/分子） | 4 对，~350 有效配体 | 评测文可小于 VS 文，但不能像两案例 | **扩 PM + 尽量加第 5 对 + 外部 holdout** |
| D. 「对接本身没坏」对照 | 弱 | 审稿常问 docking broken? | 单靶 enrichment / cognate QC 专节 |
| E. 方法新颖 | 无新打分 | 用任务+基准+失败学对冲 | 系统失败分型 + 基线门控协议 |

### 0.2 供给铁律（决定你能补多大）

公开 ChEMBL 上 **strict 双侧硬负≥50** 的非金属对极少（J0：约 3 对可用 + EGFR 作供给案例）。  
→ **不可能**做成 DiffDock 那种 15–43 靶规模的双靶四类基准。  
→ JCIM 路径只能是：**小而极严的 DualFourClass-Bench + 厚诊断 + 可复现发布**，不是「靶越多越好」。

### 0.3 成功标准（投稿 Go）

同时满足才建议冲 JCIM full（否则降 Molecular Informatics / JCAMD）：

1. **主文至少 3 个「厚」结果块**：供给审计；方向分解×基线×CI 森林；失败分型+化学型/结构敏感。  
2. **PM 扩面后** vina−heavy 的 Δ 要么 CI 不再跨 0，要么仍跨 0 但 N 与功效分析写死（两种都可投，措辞不同）。  
3. **至少 1 个外部 holdout**（文献双靶/临床双靶回顾，不进调参）。  
4. **Zenodo DOI + 一键复现 + 禁止泄漏的评测协议**。  
5. 全文无「通用 scorer / 显著通吃」表述。

---

## 1. 总路线：四波（按顺序，有闸门）

```text
Wave 0  写作骨架 + 缺口清单冻结          （1 次会议级，零算力）
Wave 1  零/低对接：把「薄」变「厚」        （云端为主）
Wave 2  中对接：扩面 + 结构敏感 + 单靶对照 （本地为主）
Wave 3  外推：第5对(若有) + 文献holdout   （有供给才开）
Wave 4  发布：Zenodo + 英文稿 + 预投稿检查
```

**闸门：**
- Wave 1 未完成（尤其 scaffold / 失败分型 / 聚合敏感）→ **不开** Wave 2 的扩面以外项目。  
- J0 补抓 22 靶后若 **0 个新 Y 对** → Wave 3 不做第 5 对，改加厚 holdout + SI。  
- 任何结果若诱使「复活 Track B 选臂」→ **拒绝改 claim**，只进诊断附录。

---

## 2. Wave 0 — 文章骨架冻结（先做，避免乱补）

**产出：** `docs/JCIM_MS_OUTLINE_V1.md`（可与本文合并执行时再写）

建议主文结构（对应实验）：

1. Intro：双靶 VS 报告习惯的缺口  
2. DualFourClass 任务 + 指标 + claim ceiling  
3. 供给审计（49→补全后）  
4. Bench 建设（K 对、协议、受体 QC）  
5. 主结果：森林 × 基线门控 × CI  
6. 稳健性：prep / θ / 聚合 / 噪声 / scaffold / 结构  
7. 失败学：端不对称、描述符捷径、同工酶、Top10 硬负  
8. 外部 holdout  
9. 单靶 sanity（对接没坏）  
10. Limitations + 开放基准使用说明  

**Wave 0 完成标志：** 每个 Results 小节都映射到下表 WP 编号。

---

## 3. 工作包（WP）明细

### Wave 1 — 零/低对接（必须做满；性价比最高）

| ID | 内容 | 回质疑 | 产出 | 谁做 |
|----|------|--------|------|------|
| **W1.1** | 补抓 J0 的 22 个缺失 ChEMBL 靶 → 重跑供给表 | cherry-pick / 供给故事不完整 | 更新 `jcim_j0j1_v0`；Y/T/− 终表 | 云端 |
| **W1.2** | 聚合敏感：pooled / mean / harmonic / summary_min 对照 | 「min 任意」 | 表+主文 1 panel | 云端 |
| **W1.3** | Cluster / Murcko bootstrap（系列去相关） | 「CI 虚高」 | 与配体 bootstrap 并列 | 云端 |
| **W1.4** | Scaffold 分层 AUROC + AChE 支架重叠 | TPSA 捷径=构建假象 | AChE 专节核心图 | 云端 |
| **W1.5** | 协变量控制：对 heavy/TPSA 分层或残差后再比对接 | 基线不公平 | SI 表 | 云端 |
| **W1.6** | 标签噪声天花板（扩到 K=4） | pChEMBL 噪声 | 复用 Stage M 脚本 | 云端 |
| **W1.7** | Dock-fail 缺失机制表（class×MW×torsdof） | MNAR | SI + Methods | 云端 |
| **W1.8** | **系统失败分型** T1–T5 扩到四对（每对≥2 案例） | 贡献薄 | 主图级 taxonomy | 云端+人工 |
| **W1.9** | 简单 ML 配体基线：ECFP4+RF / logistic（嵌套 CV，**禁止看 holdout**） | 「没跟现代 ligand-based 比」 | 与对接同表 | 云端 |
| **W1.10** | 评测协议成文：train/dev 禁用规则、LOTO 精神（即使本文不做选臂） | 可复用基准 | `PROTOCOL_BENCH_V1.md` | 云端 |

**对接次数：** 0。  
**完成标志：** 主文 Fig 2–4 与 SI S1–S6 的数据齐；AChE 捷径有支架级解释。

---

### Wave 2 — 中对接（加「规模感」与「结构稳健」；本地）

预算按 **Vina 端次数** 计（每配体每口袋 1 次）；RTM/GNINA 复用姿态另计但便宜。

| ID | 内容 | 回质疑 | 规模（约） | 优先级 |
|----|------|--------|------------|--------|
| **W2.1** | **PM strict 扩面** panel48→**~110–120**（四类定额；统一 RDKit） | N/CI/功效 | ~70 新配体 ×2 ≈ **140 Vina** + RTM + GNINA | **P0** |
| **W2.2** | AChE/BChE 从 100→**~120–140**（补 fail、补边界硬负） | 面板完整性 | ~20–40×2 ≈ **40–80 Vina** | P1 |
| **W2.3** | PIK3CA/PIK3CB 同理补到 **~120** | 同 | ~20×2 ≈ **40 Vina** | P1 |
| **W2.4** | **结构敏感（关键加分）**：PM 换 1 套替代受体（或换盒子）× **子集 40 配体×2** | 结构依赖 | **80 Vina** + RTM | **P0** |
| **W2.5** | AChE **或** BChE 二选一做换 PDB 子集 40 | 同上 | **80 Vina** | P1 |
| **W2.6** | Cognate / near-cognate **RMSD–分数–活性** 三角图（已有 QC 则补配体子集） | docking broken? | 小 | P0（偏分析） |
| **W2.7** | **单靶 enrichment sanity**：每受体用 ChEMBL actives vs property-matched decoys（或 LIT-PCBA/DUD-E 同源靶若可得）报 EF/AUROC | 「四类任务失败因对接废了」 | 每受体 ~50–100 act + decoy 对接 | **P0** |
| **W2.8** | （可选）PM 子集 GNINA **重采样/最小化** vs rescore | 引擎耦合 | 40×2 | P2 |

**Wave 2 合计粗算：**  
- 必做 P0：W2.1+W2.4+W2.7 ≈ **140 + 80 +（4 受体 × ~150）≈ 140+80+600 ≈ 820 Vina 量级**（单靶 decoy 是大头，可先做「每受体 50 act + 200 decoy」压到 ~1000 以内，或只做 PM 两端+EGFR 一端作代表以砍半）。  
- **务实砍法（推荐）：**  
  - **必做：** W2.1（PM 扩面）+ W2.4（PM 换结构子集）+ W2.7 **仅 PM 两受体**（证明主对上单靶对接可用）+ W2.6。  
  - **强烈建议：** W2.2 或 W2.5 二选一。  
  - **可进 SI：** W2.3、全受体 decoy、W2.8。  

**推荐必做对接预算：约 300–500 Vina**（不含全库 decoy）；若上全 4 受体 decoy 则 **+600–1000**。

**闸门 W2：** PM 扩面后更新 bootstrap；若 Δ CI 仍跨 0 → 正文写「功效仍不足」，**不要**再无脑扩 EGFR。

---

### Wave 3 — 外推层（有则极大加分）

| ID | 内容 | 条件 | 规模 |
|----|------|------|------|
| **W3.1** | 第 5 对靶（Tier-Y/T，非金属）建面板 + 对接 | W1.1 发现新 Y/T | ~100×2 ≈ **200 Vina** + RTM/GNINA |
| **W3.2** | **文献/临床双靶 holdout**：从 JMC/批准药目录抽 20–40 已知双靶 + 匹配单靶负例，**一次性**评测冻结分数规则 | 始终可做 | 视分子数；建议 **≤80 配体×现有受体** 或独立小面板 |
| **W3.3** | 跨对「协议可迁移」叙述：同一脚本出全部表 | 与 W1.10 绑定 | 0 对接 |

**不要：** 把 holdout 用于选阈值/选臂。

---

### Wave 4 — 发布与投稿（与实验并行收尾）

| ID | 内容 |
|----|------|
| **W4.1** | Zenodo：assembled 长表、面板 SMILES、受体 PDB/盒子、分数、脚本、CLAIM_CEILING |
| **W4.2** | GitHub release + `pip`/`README` 最小复现（至少重算 AUROC/CI） |
| **W4.3** | 英文稿（JCIM 版式）+ cover letter 明确 evaluation article |
| **W4.4** | 内部红队：用 `CRITIQUE` 清单逐条勾「正文何处已答」 |
| **W4.5** | 备降刊：同一稿可投 Mol. Inf. / JCAMD（若 Wave 2 缩水） |

---

## 4. 推荐执行顺序（白话第 1/2/3…步）

### 第 1 步（立刻，云端，~Wave 1 全集）
把「分析厚度」拉到 JCIM SI 级：W1.1–W1.10。  
**不停本地对接。**

### 第 2 步（本地，Wave 2 必做核）
1. PM → ~110–120 strict 扩面并对接三通道  
2. PM 换结构子集敏感  
3. PM 两受体单靶 enrichment sanity  
4. 更新 CI / 基线门控 / 主图  

### 第 3 步（本地或云端，Wave 2/3 加厚）
- AChE 补面板或换 PDB 子集（二选一优先）  
- 文献双靶 holdout  
- 若供给审计冒出第 5 对 → 再开 200 Vina  

### 第 4 步（发布）
Zenodo DOI → 英文稿 → 投稿。

---

## 5. 明确不做（避免工作量膨胀却不加分）

| 不做 | 原因 |
|------|------|
| 再扩 EGFR 赌 D/B 显著 | 供给上限 7 个 strict B_only；已 No-Go |
| 复活 Track B 主文选臂 | Stage M=Weak；与 claim ceiling 冲突 |
| HDAC 金属对进主基准 | 金属酶对接协议另一套故事 |
| 全 K=4 都做换受体全量 | 边际递减；PM+1 足够主文 |
| 湿实验 | 评测路线非必须；一开就变发现文资源模型 |
| 乘客/linker 设计线 | 另一篇文章 |
| 无新供给时硬凑第 5–8 对 θ=6 松标签 | 会毁「严格硬负」卖点 |

---

## 6. 工作量与「是否够撑 JCIM」对照

| 方案 | 内容 | 相对现况增量 | 估测把握 |
|------|------|--------------|----------|
| **现况 only** | K=4 + CI 包 | 0 | JCIM full **偏低**；Mol. Inf. **中** |
| **A：Wave1+W4** | 分析加厚+发布 | 零对接 | JCIM **仍紧**；短稿/中档刊 **中高** |
| **B：A + Wave2 核**（推荐） | +PM 扩面+结构敏感+单靶 sanity | ~300–500 Vina | JCIM full **中**（诚实评测定位） |
| **C：B + Wave2 全 + Wave3** | +补面板+第5对+holdout+多受体 decoy | ~800–1500 Vina | JCIM full **中高**（仍无湿实验则非稳赢） |

**结论：**  
你要「把握大一点 + 工作量够撑 JCIM」，应打 **方案 B 为下限、方案 C 为冲刺**；  
**只写现况不够**——这个判断成立；  
**但正确加厚是「扩主对 + 稳健性 + 单靶对照 + 发布 + holdout」，不是无限新靶对。**

---

## 7. 与现有文件关系

| 文件 | 关系 |
|------|------|
| `JCIM_NEXT_ROUTE_PLAIN.md` | 对接阶段已完成；**本文件接管「加厚阶段」** |
| `CLAIM_CEILING.md` | 仍有效；Wave 2/3 不得突破 |
| `BENCHMARK_ANALYSIS_V1.md` | Wave 1 在其上扩展，不推翻 |
| `JCIM_GAP_TWO_SCENARIOS.md` | 方案乙里「扩 EGFR」作废；本文件是修正后的方案乙 |

---

## 8. 下一步请求（执行时）

默认授权执行顺序：

1. **先跑 Wave 1 全套脚本包**（云端可做）  
2. 用户本地确认后开 **W2.1 PM 扩面**（需 ChEMBL 抽样名单 + 对接）  
3. 并行准备 **W3.2 holdout 分子清单**（可不对接先定名单）

若只选一个最贵但最值的对接：**W2.1 PM 扩面**；  
若只选一个最值的「防 docking broken」：**W2.7 缩略版（PM 两端）**；  
若只选一个最值的零对接：**W1.4+W1.8+W1.9**。

---

## 9. 审稿审计后的执行入口（优先用这个）

深度重算后发现主指标需改为**口袋匹配**，且 exhaustiveness 对照优先于盲目扩面。  
**复制给本地 agent 的总命令：**  
[`AGENT_COMMAND_JCIM_STRENGTHEN_SUPPLEMENT.md`](AGENT_COMMAND_JCIM_STRENGTHEN_SUPPLEMENT.md)  

审计正文：[`../data/jcim_bench_v0/analysis/REVIEWER_AUDIT_V1.md`](../data/jcim_bench_v0/analysis/REVIEWER_AUDIT_V1.md)
