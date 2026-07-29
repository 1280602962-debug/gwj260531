# 冲 JCIM：下一步路线（白话版）

> 依据：Stage M（Track B=Weak）+ `JCIM_ROUTE_ASSESSMENT_V1.md` + 已完成的 J0/J1（`data/jcim_j0j1_v0/`）  
> 日期：2026-07-29  
> **目标刊：** JCIM（评测/基准 Article，不要求方法赢）  
> **对接授权：** 本文批准后可开；范围只限下文「该对接什么」。

---

## 0. 论文要卖什么（先钉死）

**标题级主张（评测文，不是发现药）：**

> 双靶对接决策应建成「四类硬负」任务；常用池化分数会掩盖端不对称；公开数据里严格硬负供给极稀缺；在统一 prep + 多打分通道下报告方向分解与平凡基线。

| 贡献块 | 现在状态 |
|--------|----------|
| 任务定义 + 池化会抵消 | 已有（EGFR/PM） |
| 公开数据供给审计（49 对） | **已有**；待抓 22 靶可后补 |
| 4 对靶对接证据 + 3 打分通道 + 基线 | **缺** ← 下面步骤补这个 |
| prep 敏感性 | 已有（EGFR M4） |
| 公开基准打包 | 对接完成后做 |

**不卖：** 通用决策尺子已验证；再扩 EGFR 赌显著；乘客/湿实验。

---

## 1. 冻结的四对靶（K=4）

| 座位 | 靶对 | 角色 | 对接吗 |
|------|------|------|--------|
| 1 | PIK3CA / mTOR | 主开发对（strict 够厚） | **要**：先改成统一 RDKit；可选扩到 ~110 |
| 2 | EGFR / HER2 | **供给受限案例**（B_only_strict=7） | **不再新对接**；用现有统一 prep EH110 |
| 3 | AChE / BChE | 新主对（strict 够厚） | **要**：结构冻结 → 面板 → 对接 |
| 4 | PIK3CA / PIK3CB | 同工酶对照（叙事写「过近」） | **要**：结构冻结 → 面板 → 对接 |

不要进预算：HDAC 金属对、再扩 EGFR、NLRP3/JNK1、Mcl-1 等薄对（除非你另批）。

---

## 2. 路线：第一步 → 第二步 → …

### 第一步 — 冻结授权与协议（文档，几乎零算力）

- 把 K=4 从「草案」标成 **已批准对接**  
- 写死：全库配体准备 = **RDKit ETKDG + meeko**；主指标 = 方向分解；必报平凡基线；打分通道 = Vina + RTM +（稍后）GNINA  
- 产出：更新 `PAIR_ROLES_DRAFT_J1.yaml` → `PAIR_ROLES_APPROVED_JCIM.yaml`（或同文件改 `docking_authorized: true`）

**谁做：** 云端 agent 即可。  
**命令：** 见 `AGENT_COMMAND_JCIM_DOCKING_PHASE.md` 的「第一步」段。

---

### 第二步 — 统一旧数据 + 准备新受体（本地为主）

并行两件事：

1. **PIK3CA/mTOR panel48**：同一套 RDKit 重准备并对接两端（约 96 次 Vina）+ RTM；可选再按 strict 定额扩到 ~110。  
2. **新受体冻结**：AChE、BChE、PIK3CB（PIK3CA 继续用 4L23）— 选 PDB、盒子、cognate QC，**先 QC 过关再批量对接**。

EGFR：**跳过**（已齐）。

**谁做：** 本地对接 agent。  
**完成标志：** PM 有全 RDKit 分数表；三套新受体有冻结文件 + cognate QC 通过记录。

---

### 第三步 — 新两对靶建面板并对接（本地）

对 **AChE/BChE**、**PIK3CA/PIK3CB** 各建约 100–120 分子四类面板：

- **按 strict（6.5/5.5）定额抽样**，不要先按 θ=6 堆满再改标签  
- 协议对齐：E、seed、n_modes 与现有冻结协议一致（写进每对 `protocol.yaml`）  
- Vina → RTM best-of-K  

预算大约各 220 次 Vina（合计 ~440；加上第二步 ~96，总计常落在 **540–760** 若含 PM 扩面）。

---

### 第四步 — 第二打分通道（本地，尽量轻）

对全部已对接姿态跑 **GNINA CNN rescore**（复用姿态，不重新大采样）。  
目的：证明结论不绑死 Vina+RTM。

---

### 第五步 — 分析成文 + 打包（云端可做）

对 K=4 统一出：

- 方向分解 AUROC（D/A、D/B）+ 配对 bootstrap CI  
- 平凡基线（重原子 / MW / cLogP / TPSA / 可选 Morgan）  
- 跨对靶森林图  
- 供给审计图（已有 49 对）+ prep 敏感性（已有）  
- Zenodo/GitHub：**DualFourClass-Bench**  

然后按 JCIM **Evaluation / Benchmark** 叙事投稿。

---

## 3. 明确不要做的事

- 同协议再扩 EGFR「赌显著」  
- 未做受体 cognate QC 就开数百个 ligand  
- 把 PIK3CA/PIK3CB 写成「与 AChE 同等外推力」的主对（它是同工酶对照）  
- 把 Track B 方法竞赛混进这篇（Weak 未解除）  
- 等 ChEMBL 那 22 个待抓靶齐了才开工（可后补进供给表，不挡对接）

---

## 4. 退路（任何一步不够仍有文）

若 AChE 或 PIK3CB 结构/面板做崩：至少用 **PM（RDKit）+ EGFR（已有）+ 供给审计** 仍可投 Mol. Inf. / JCAMD；JCIM 可降为 short / 或补第三对后再投。

---

## 5. 状态（2026-07-29）

**已批准对接阶段。** `docking_authorized: true` → [`../data/protocols/PAIR_ROLES_APPROVED_JCIM.yaml`](../data/protocols/PAIR_ROLES_APPROVED_JCIM.yaml)。

| 步骤 | 状态 |
|------|------|
| 1 授权 + 目录 | **完成** |
| 2A PM48 RDKit Vina+RTM | **完成**（含 LigPrep Δ） |
| 2B 新受体 cognate QC | **完成**（AChE=4EY7、BChE=4BDS、PIK3CB=2WXF） |
| 3 新面板 Vina+RTM | **完成**（AChE 191/200；PIK3 199/200；超柔配体 skip） |
| 4 GNINA CNN rescore | **完成**（mode_01；AChE/PAB/PM/EGFR） |
| 5 基准包 | **主表 v0 已写入** `data/jcim_bench_v0/tables/directional_forest_v0.csv`；姿态大文件在本地 results |

执行命令：[`AGENT_COMMAND_JCIM_DOCKING_PHASE.md`](AGENT_COMMAND_JCIM_DOCKING_PHASE.md)。

**当前：** 对接打分阶段完成；PLK1/NLRP3 已另进程恢复。下一步可写作投稿或补平凡基线/bootstrap。
