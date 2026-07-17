# JNK2 对接前初筛 — ARS 重规划方案（v1）

**依据：** ARS-Codex deep-research（OpenAlex Q1–Q7 + 核心文献 A1–C4）  
**输入：** ~527,779 唯一胺（ChEMBL）  
**下游：** 松阈值 Schrödinger 共价对接 → AF3 mPAE 门控（Phase0 已校准）  
**贡献对齐：** Framework 必达；Chemistry 升级  

---

## 0. 规划原则（先接受再跑库）

1. **不把 Step0–7 原样称为“活性预测”** → 改称 **chemotype triage / enrichment**。  
2. **不把固定 Tc=0.25/0.55 当真理** → 先跑特征，再做阈值敏感性，再切 dock_ready。  
3. **两仓分榜到底**（Sim / Novel）；pan（JNK-IN-8）第三仓。  
4. **Glide 不主排序**；AF3 mPAE 才是采购/合成门控。  
5. 每个分子必须可解释：为何进/出 + 最近锚点 + core-Tc + 反应性桶。

---

## 1. 修订后的漏斗（MVS = 最小可行科学版）

```
L0  胺库体检 + 标准化
L1  胺预筛（控爆炸）
L2  丙烯酰胺枚举（主 warhead；变体第二轮）
L3  产物 QC + 反应性分桶（flag，慎硬删）
L4  多表征特征（core-Tc / full-Tc / ErG / Murcko / hinge）
L5  回顾性校准探针（小集）→ 定阈值与配额
L6  双轨+pan 分仓排序
L7  dock_ready 5k–15k 交付（分文件）
    → 对接粗筛 → AF3 → 人工/选择性负筛设计
```

相对原 Step0–7 的关键变化：

| 原方案 | 新方案 |
|--------|--------|
| 装弹头后只算 ECFP | **并行 core-Tc（去弹头）+ full-Tc** |
| 固定 Tc 窗直接分仓 | **先特征全库 → 校准/敏感性 → 再分仓** |
| hinge SMARTS 硬门 Novel | Novel 改为 **低 core-Tc + 新 Murcko + (ErG↑ 或 hinge)** |
| PAINS 硬删 | **flag**；过反应 SMARTS 进 reactivity_bucket |
| 混成“活性初筛”叙事 | 明确 enrichment；成功标准可证伪 |
| Glide 隐含核心 | 文档写死：**松阈值粗筛 only** |

---

## 2. 分层细节

### L0–L1 胺侧（基本保留）
- 去盐、最大片段、可酰化胺位数  
- sites==1 优先；sites==2 且胺 MW≤350；≥3 丢弃  
- 胺 MW 120–450  

### L2 枚举
- 主：标准丙烯酰胺  
- 第二轮（仅 Sim 高价值子集）：Me2N–CH2–CH=CH–C(=O)–（YL/JNK-IN 型）  
- 记录 attachment atom；产物 InChIKey 去重  

### L3 QC + 反应性
保留：单 Michael、MW 350–650、clogP −1～6、RB≤12、SA 记录  
新增 `reactivity_bucket`：
- `ok`：单丙烯酰胺，无额外强电泳  
- `watch`：多取代/延长共轭等提高反应性特征  
- `bad`：双亲电、氯乙酰胺、乙烯砜等 → **硬排除**  
PAINS：`pains_flag` 不直接删（进 watch 配额限制）

### L4 特征（核心升级）
对每个产物：
1. `tc_full_*` 对四锚点（ECFP4）  
2. `tc_core_*`：剥离丙烯酰胺/二甲胺丙烯酰胺后对锚点 **core** 的 ECFP4  
3. `max_tc_core`, `nearest_anchor_core`  
4. Murcko + `same_scaffold_as`  
5. `erg_max`（RDKit ErG 对四锚点）  
6. `hinge_hits`（SMARTS 计数，作辅）  
7. 可选：父分子 ChEMBL kinase/JNK 证据位  

### L5 校准（宣称理性前必做；可小集快速）
探针集（不必大）：
- 阳性：四锚点 + 若有同系列近邻  
- 难例：56d（正交）  
- 阴性对照：饱和酰胺（YL5084R/56a 类，不应靠共价几何进）  
- 无关共价丙烯酰胺（其他激酶）若干  

报告：
- leave-one-anchor-out 回收  
- core-Tc vs full-Tc 对排序影响  
- Tc 网格热图（仓规模稳定性）  
- 选定 **工作阈值** 写入 `thresholds.json`（可复现）

### L6 分仓规则（建议初值；以 L5 为准）

**Track-Sim（活性优先）**
- `0.22 ≤ max_tc_core ≤ 0.55`（初值；L5 调整）  
- `nearest_anchor_core ∈ {YL5084,YL2056,56d}` 优先  
- `max_tc_full > 0.70` → near_duplicate，不进对接主名单  
- 子标签：`sim_YL` / `sim_56d`  
- 排序：chembl_evidence ↓, max_tc_core ↓, SA ↑, reactivity ok 优先  

**Track-Novel（骨架探索）**
- `max_tc_core < 0.22` 且 Murcko 不同于四锚点  
- **且**（`erg_max` 高于分位数阈值 **或** `hinge_hits≥1`）  
- `reactivity_bucket != bad`；`watch` ≤20% 配额  
- 排序：erg_max ↓, hinge_hits ↓, SA ↑；偏好 core-Tc 落在 0.08–0.22（避免极端异类垃圾）  

**Track-Pan（对照）**
- nearest 为 JNK-IN-8 的单独文件；不作选择性主交付  

### L7 dock_ready 配额（默认总计 ~10k）
| 仓 | 数量 | 备注 |
|----|------|------|
| Sim (YL) | 4200 | 70% of sim |
| Sim (56d) | 1800 | 30%；后续 AF3 配额提高 |
| Novel | 3500 | ErG/hinge 过线 |
| Pan control | 500 | 对照 |
| **合计** | **~10k** | 可按算力改为 5k/15k |

不足则如实报告，**禁止用 discard 填数**。

---

## 3. 成功标准（可证伪）

对接前：
- [ ] L5 校准表与 `thresholds.json` 存在  
- [ ] core-Tc 与 full-Tc 分仓 Jaccard < 1（证明去弹头有影响）  
- [ ] 56d leave-one-out 不被系统性丢弃  
- [ ] Novel 随机 50 样例人工：明显非激酶配体比例有记录  

对接/AF3 后：
- [ ] Glide 仅作几何过滤；采购名单以 mPAE 为主  
- [ ] Sim 与 Novel 分别报告 mPAE 通过率  
- [ ] 设计 JNK1 counterscreen（湿实验或计算 Δ）  

湿实验（升级贡献）：
- [ ] 分层抽样（Sim高 / Novel高 / 边界 / 随机）→ labeling / kinact/KI / C116S  

若分层抽样命中率不优于“随机合格丙烯酰胺”基线 → **证伪本 triage**，回退改 ErG/3D 药效团权重。

---

## 4. 给本地 WSL Agent 的执行指令（可复制）

```text
使用已安装思路：ARS 重规划的 JNK2 对接前初筛 v1。
工作目录：~/jnk2_amine_screen/
输入：unique_amines.csv（~527779）
禁止打开整文件；chunksize=50000；每步写日志与计数。

锚点 SMILES（full）：
YL5084 / YL2056 / JNK-IN-8 / 56d（与项目 phase0_compounds_seed.csv 一致）
同时生成各锚点 warhead-stripped core SMILES 用于 tc_core。

按顺序执行并产出：
L0_report.md
L1_amines_prescreen.csv
L2_acrylamides_raw.csv          # 仅标准丙烯酰胺第一轮
L3_acrylamides_qc.csv           # 含 reactivity_bucket, pains_flag, sa
L4_features.csv                 # tc_full_*, tc_core_*, erg_*, murcko, hinge_hits
L5_calibration.md + thresholds.json
L6_track_sim.csv / L6_track_novel.csv / L6_track_pan.csv / L6_discard.csv
L7_dock_ready_{sim_yl,sim_56d,novel,pan}.csv + L7_summary.json

硬规则：
1) 先跑 L0–L4 全库特征，再 L5 定阈值，再 L6–L7（不要跳过校准）
2) Sim/Novel/Pan 分文件；禁止合并成单一总分榜
3) PAINS 不硬删；reactivity=bad 硬删
4) 第二 warhead 变体留到 L7 之后只对 Sim Top 做第二轮
5) 最终 summary 必须含漏斗计数、阈值敏感性摘要、两仓 ID 交集=0 证明

完成后打印验收清单（成功标准对接前四项）。
```

---

## 5. 明确不做什么（本阶段）

- 不训练深度活性模型当主筛（标签少、易泄漏）  
- 不把生成模型当作 hit  
- 不把疾病药效写进初筛目标  
- 不对 527k 全库跑 AF3  

---

## 6. Editorial one-liner

> 合理的初筛不是“ECFP 双轨脚本”，而是：**去弹头相似 + ErG 补新骨架 + 反应性分桶 + 校准阈值 + 分仓配额**，只为把 52 万胺压到可被 Glide粗筛 / AF3精排验证的几千分子，并用可证伪指标约束自己。
