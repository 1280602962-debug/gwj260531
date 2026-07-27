# 决策：PI-103 cognate QC No-Go 后的下一步

> 日期：2026-07-27  
> 依据：本地 Phase A/B 结果（`results/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/`）  
> 对照：EGFR/HER2 `SENSITIVITY_VERDICT.md` 选 E 规则

---

## 1. 事实

| 端 | E=8 mode1 / best9 | E=16 best9（同 seed） | 判读 |
|----|-------------------|----------------------|------|
| 4L23 | 0.624 / 0.624 | 0.624 | E=8 已够 |
| 4JT6 | 7.118 / **5.003** | **0.445**（mode3；mode1 仍~7.1） | **E=8 采样失败；E=16 采样成功** |

协议硬门槛：两端 `best_of_9 < 2 Å` → 当前 **No-Go**，禁止 48×2。

这与 EGFR 侧「不必升 E」的条件**正好相反**：EGFR 是 E=8 已有 best9&lt;2；这里是 **E=8 失败、E=16 成功**。

---

## 2. 三选一怎么选

| 选项 | 是否推荐 | 理由 |
|------|----------|------|
| 1 只修蛋白/盒子再硬撑 E=8 | **不作为主路径** | E=16 已找回 0.445 Å → 盒子/坐标系大体成立；死磕 E=8 可能反复空转 |
| **2 书面确认后本对靶用 E=16** | **推荐（主路径）** | 满足「E=8 采样失败且 E=16 成功」才升 E 的规则；**仅对本 pair 冻结，不回改 EGFR 的 E=8** |
| 3 只归档不对接 | 不推荐 | 会卡住 C4 外推；除非短期无算力 |

**推荐执行顺序 = 短核查（选项1的轻量版）→ 书面锁定 E=16（选项2）→ 重做 Phase B → Go 再全面板。**

---

## 3. 立刻要做的清单（按序）

### Step 0 — 轻量 sanity（半天内，不要大改）

只确认没有「明显准备错误」；**不要**为了硬过 E=8 重做一轮蛋白：

1. 4JT6 prepared 是否仍是 **ATP 位点激酶域**（勿混入 FRB/FKBP 逻辑）  
2. 盒子是否来自 **同一坐标系** 的 X6K（与 docking receptor 一致）  
3. 配体 PM48_01 是否确为 PI-103 / 与 X6K 同构型输入  
4. 受体是否丢掉关键辅因子/不该质子化的铰链（快速目视即可）

若发现硬错误 → 修正后 **仍先用 E=8 重跑 Phase B**。  
若无硬错误（预期如此）→ 进入 Step 1。

### Step 1 — 书面冻结本对靶 exhaustiveness=16

在 `protocol/protocol.yaml` / `SEED_POLICY` 旁增加：

```yaml
pair_id: DTPAIR-01_pik3ca_mtor
seed_fixed_global: 20260727
n_modes: 9
energy_range: 3
exhaustiveness_v0_1: 16   # pair-specific; EGFR/HER2 remains 8
reason: >
  PI-103@4JT6: E=8 best_of_9=5.003 Å fail; E=16 best_of_9=0.445 Å pass
  (seed 20260727). Raise E only because cognate sampling failed at 8.
```

**禁止**把 EGFR panel40 的 E=8 改成 16。

### Step 2 — Phase B′：E=16 确认（仍只 PM48_01）

| 作业 | 内容 |
|------|------|
| 主 seed | 4L23 + 4JT6，E=16，seed=`20260727` |
| 噪声抽查 | **仅 4JT6×PM48_01**，E=16，额外 seed `{7, 42}` |

通过标准（与 SOP 一致）：

- 两端（主 seed）`best_of_9 < 2.0 Å` → **Go**  
- 4JT6 三 seed 的 best9 均 &lt;2（或至少主 seed + 2/3 seed &lt;2）→ 认为 E=16 采样稳定  
- mode1 仍可失败（与 TAK-285@3POZ 同类）→ **仍算采样 QC 通过**，全面板必须保留 top-9 + 计划 RTM

写出：`analysis/cognate_redock_v0/COGNATE_QC_VERDICT_E16.md`

### Step 3 — Go 后全面板

- 96 作业：48 × {4L23,4JT6}  
- **一律 E=16，seed=20260727，n_modes=9**  
- 不要中途混用 E=8/16  
- 汇总后照旧做 Vina mean/min + RTM 消融

### Step 4 — 若 E=16 仍不稳

| 情况 | 动作 |
|------|------|
| 主 seed 过、噪声 seed 偶发 best9≥2 | 可试 E=32 **仅 PM48_01@4JT6**；仍不稳则换 mTOR 备用结构（如 4JT5）做对照 QC，不改面板分子 |
| E=16/32 主 seed 仍 best9≥2 | 回到深度蛋白/盒子排查；或暂停该对靶主表，改第三对 |

---

## 4. 给本地 agent 的一句话

> 保持全面板 No-Go。先做 4JT6 轻量 sanity；若无硬错误，**书面冻结本项目 exhaustiveness=16**（EGFR 仍为 8），用 seed `20260727`（外加 7/42 噪声）重做 PM48_01 双端 Phase B；两端 best9&lt;2 后再开 48×2（全 E=16）。

**详细逐步操作（可粘贴任务块）：** [`PIK3CA_MTOR_AGENT_OPS_E16.md`](PIK3CA_MTOR_AGENT_OPS_E16.md)

---

## 5. 科学表述（写进 Methods 的口径）

- EGFR/HER2：cognate 在 E=8 已满足 → `exhaustiveness=8`  
- PIK3CA/mTOR：mTOR(4JT6, 3.6 Å) 在 E=8 无法恢复 PI-103 近晶构象，E=16 可恢复 → **对该靶对**使用 `exhaustiveness=16`  
- 升 E 的理由是 **采样 QC**，不是刷分；Vina mode1 在 4JT6 上仍可能错排，故协议保留 ensemble + RTM
