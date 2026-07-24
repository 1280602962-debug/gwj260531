# EGFR/HER2：朴素双端融合失败诊断 Demo

> 回答三问：**为什么选这对** · **诊断表要定什么** · **怎么跑**  
> 对应决策尺子规划的 C1/C2（诊断 + 姿态/重打分必要性），不是全量基准。  
> **不再以 passenger/moiety 为必做前提**；moiety 仅可选附录。

---

## 1. 为什么先选 EGFR / HER2

不是因为它「最有药学故事」，而是因为它最适合 **用最少分子做姿态 QC + 证明朴素双端融合失败**。

| 条件 | EGFR/HER2 | 对比 |
|------|-----------|------|
| **双端同配体共晶** | TAK-285：EGFR [3POZ](https://www.rcsb.org/structure/3POZ) + HER2 [3RCD](https://www.rcsb.org/structure/3RCD) | PIK3CA/mTOR **无**双端同配体共晶 → 无法先做姿态金标准 |
| **口袋可对接** | 经典激酶 ATP，协议成熟 | Mcl-1/Bcl-xL 是 PPI，对接方差大，不宜做第一个诊断 |
| **四类标签够用** | dual 1182 / A-only 207 / B-only 46（ChEMBL 审计） | 有硬负样本，诊断表不只看 dual |
| **demo 成本** | 共晶 + 小面板（~20–50 分子）即可出诊断 | PIK3CA/mTOR 更适合 **规模主表**，不适合第一枪 |

**角色分工（不要搞反）：**

| 对 | 角色 |
|----|------|
| **EGFR/HER2** | **姿态 QC + 朴素融合失败 / 重打分纠偏诊断**（先跑通） |
| PIK3CA/mTOR | 主规模 dual-vs-single 外推 |
| Mcl-1/Bcl-xL | 异质口袋外推；linked 仅分层附录 |
| NLRP3/JNK1 | 外部锚点后的 holdout，不参与调参 |

**已知风险（写进 Methods）：** EGFR/HER2 同源，交叉活性高 → dual 比例偏高属预期；主证据是 **A-only/B-only 是否被朴素融合抬进前列**，以及协议能否压下它们。

---

## 2. 诊断表要先定死的内容

诊断表不是随便对接打分。下面 **全部要在跑数前冻结**；改一项就重跑并记版本。

### 2.1 分子面板（谁进表）=「验证对接列表」的组成

诊断阶段的对接列表 **不是** 海量虚拟筛选库，也 **不是** 实验室私有 NLRP3/JNK1 分子。  
它是一张 **带四类标签的小面板**，专门用来看 whole-mol 是否错排。

| 类 | 定义 | demo 建议 N | 来源 | 在诊断里的角色 |
|----|------|-------------|------|----------------|
| **Dual** | 两端 pChEMBL ≥ 6 | 10–20 | ChEMBL 配对 + 文献双 TKI；**必含 TAK-285** | 真双靶；看 moiety 是否仍保留两端信号 |
| **A-only** | EGFR≥6 且 HER2 **测得** <6 | 10–20 | ChEMBL 硬负 | **主攻击对象**：whole-mol 是否把弱端抬成「像双靶」 |
| **B-only** | HER2≥6 且 EGFR **测得** <6 | 尽量多（审计约 46） | ChEMBL 硬负 | 对称检查 |
| （可选）**Neither** | 两端都测且都 <6 | 5–10 | ChEMBL | 假阳性底噪 |
| （可选）**Decoy** | 未测 / 物化匹配 decoy | 后期 | DUD-E 类 | 全基准再用；demo 可省略 |

**组成公式（demo）：**

```text
对接列表 = {TAK-285}
         ∪ 抽样 Dual
         ∪ 抽样 A-only
         ∪ 可用 B-only
         ∪（可选 Neither）
```

每个分子要对 **EGFR(3POZ) 与 HER2(3RCD) 各对接一次**（整分子）→ 列表长度 N 时，对接任务数 = **2N**。

**冻结规则：**

- 活性类型优先 IC50/Ki/Kd/EC50；每分子每靶取 **max pChEMBL**  
- 阈值：pChEMBL ≥ 6 = active；测得 <6 = weak；**未测 ≠ inactive**（未测不进 A-only/B-only）  
- 记录 ChEMBL ID、SMILES、scaffold（Bemis–Murcko）防系列全是同一骨架  
- **不做：** 随机买库、只对接 dual、把私有细胞分子混进第一张诊断表

### 2.2 Moiety 标注（passenger 定义）

每个分子必须有可复现标注（JSON schema 试点）：

```json
{
  "mol_id": "TAK-285",
  "architecture": "merged",
  "smarts_or_atom_maps": {
    "moiety_EGFR": [/* atom indices or SMARTS */],
    "moiety_HER2": [/* 可与 EGFR 重叠 partial */],
    "linker_or_shared_core": [/* 共享核 / 连接原子 */]
  },
  "passenger_when_docked_to_EGFR": "moiety_HER2 (+ linker if outside pocket)",
  "passenger_when_docked_to_HER2": "moiety_EGFR (+ ...)",
  "annotator": "initials",
  "notes": "based on 3POZ/3RCD contacts"
}
```

**必须事先约定：**

| 决策 | 建议默认 |
|------|----------|
| merged 分子如何切 | 以 **共晶接触残基** 为主（TAK-285），文献 SAR 为辅；允许 shared core |
| moiety 分怎么算 | **方案 A（推荐 demo）：** 对接仍用整分子；分数只对 moiety 原子与口袋的相互作用重计 / 掩蔽 passenger 贡献  
| | **方案 B：** 切出 moiety 片段再对接（改变构象搜索，与「乘客污染」叙事不完全同构，作消融） |
| 共享核算哪边 | 计入 **当前对接靶的 moiety**；不双计进 passenger |
| 无清晰二分时 | 标 `ambiguous`，**主诊断表排除**，附录报告 |

Demo 主结论必须来自 **方案 A**（整分子采样 + moiety 重打分），否则审稿人会说你换了搜索问题。

### 2.3 对接协议（两端各自固定）

| 项 | 冻结值（demo 建议） |
|----|---------------------|
| 结构 | EGFR: **3POZ**（去配体、保辅因子按团队惯例）；HER2: **3RCD** |
| 引擎 | **普通 AutoDock Vina（或 smina）可以，且推荐作为 demo 主引擎**；全文只用一个，不要与 GNINA 混跑 |
| 盒子 | 以共晶配体几何中心 ± 扩展（如 20–25 Å）；两端记录中心与尺寸 |
| 构象 | 同一 SMILES → 同一 3D 生成种子；`exhaustiveness`（建议 ≥8）与 `cpu`/种子写死 |
| 自对接 QC | TAK-285 重对接 RMSD vs 3POZ / 3RCD；**RMSD 门槛先定**（如 ≤2.0 Å 重原子）不过关则改盒子/质子化，不改分子面板 |
| 输出 | 每分子每靶：top1 pose + whole-mol score（Vina 亲和力）；可选 top-k |

#### 普通 Vina 可不可以？

**可以。** 本诊断要证明的是「整分子分数被 passenger 污染」，不是「我们的采样器更准」。

| 选择 | 适用 |
|------|------|
| **AutoDock Vina / smina（推荐 demo）** | 足够；易复现；与 Zhou 2013 等双靶 VS 文献同族；先过 TAK-285 RMSD 即合格 |
| **GNINA**（Vina 搜索 + CNN 分） | 主文可作 **同一姿态上的重打分消融**，或第二引擎重复诊断；不要当成必须先有的条件 |
| DiffDock 等 AI 对接 | **仅消融**；不作第一张诊断表的唯一引擎 |

注意：Vina 分「越负越好」；做 rank/fuse 前统一符号或全程保持原生并在表头注明。  
Moiety 分在 Vina 姿态上用 **方案 A 重计/掩蔽**，不必为了 moiety 换引擎。

### 2.4 分数定义（表的列）

对每个分子、每个靶 \(T \in \{\mathrm{EGFR},\mathrm{HER2}\}\)：

| 列 | 含义 |
|----|------|
| `S_whole(T)` | 整分子对接分（引擎原生；注意符号方向统一为「越大越好」或全程保持原生并注明） |
| `S_moiety(T)` | 仅计「对 T 的药效团」与口袋的相互作用（方案 A） |
| `S_passenger(T)` | 乘客部分贡献（可定义为 whole − moiety，或显式乘客原子分） |
| `rank_whole(T)` / `rank_moiety(T)` | 在 **同一面板内** 的分位数或秩（跨靶勿直接比原始分） |
| `fuse_mean_whole` | \((\mathrm{rank}_A+\mathrm{rank}_B)/2\) |
| `fuse_min_whole` | \(\min(\mathrm{rank}_A,\mathrm{rank}_B)\) |
| `fuse_min_moiety` | \(\min(\mathrm{rank}^{\mathrm{moi}}_A,\mathrm{rank}^{\mathrm{moi}}_B)\) |

校准：demo 阶段用 **面板内 rank/percentile** 即可；全基准再上 decoy Z-score。

### 2.5 主诊断读数（表要回答的假设）

先验假设（写进笔记，跑前登记）：

1. **H1（污染）：** 在 A-only 分子上，`S_whole(HER2)` 或 `rank_whole(HER2)` 被不合理抬高（相对 moiety）。  
2. **H2（错排）：** `fuse_mean_whole` / `fuse_min_whole` 把 A-only/B-only 排进 dual 前列的比例，高于 `fuse_min_moiety`。  
3. **H3（真端保留）：** 在 A-only 上，`S_moiety(EGFR)` 仍应高于弱端 moiety；moiety 不应毁掉真结合端。

**主表（分子×方法）建议列：**

`mol_id | class | arch | S_w_EGFR | S_w_HER2 | S_m_EGFR | S_m_HER2 | fuse_mean_w | fuse_min_w | fuse_min_m | dual_rank_w | dual_rank_m`

**汇总表（一类一行）：**

| class | median fuse_mean_w | median fuse_min_m | EF@10% dual (whole) | EF@10% dual (moiety) | A-only 进入 top10% 比例 (whole vs moiety) |
|-------|--------------------|-------------------|---------------------|----------------------|-------------------------------------------|

**过关（demo Go）：**

- 自对接 RMSD 达标  
- A-only 的 **假双靶抬升** 在 whole 上可见，且 moiety 融合降低 A-only 进入 top 分位的比例  
- dual 的真端 moiety 分未系统性崩溃  

**不过关：** 先查标注/盒子/分数方向，再谈叙事。

---

## 3. 怎么跑（操作顺序）

### Phase 0 — 冻结（半天文档）

1. 复制本协议版本号（如 `egfr_her2_diag_v0.1`）  
2. 锁定 PDB、引擎、盒子、RMSD 门槛、pChEMBL 阈值、方案 A 打分公式  
3. 列出面板分子 ID 清单（含 TAK-285）

### Phase 1 — 结构与配体

1. 下载 3POZ / 3RCD；准备受体（去水策略、质子化一致）  
2. 抽出共晶配体 03Q/TAK-285，核对原子编号  
3. 面板 SMILES → 3D（固定种子）→ 质子化/互变异构规则写死  

### Phase 2 — 姿态 QC（只做 TAK-285）

1. 整分子重对接 EGFR(3POZ)、HER2(3RCD)  
2. 算 RMSD vs 晶体；不达标则只调协议，不改假设  
3. 按晶体接触完成 TAK-285 的 moiety 标注，作为 schema 金标准  

### Phase 3 — 面板对接

1. 每个分子 → EGFR、HER2 各对接一次（整分子）  
2. 存 pose + `S_whole`  
3. 按标注算 `S_moiety` / `S_passenger`（方案 A）  
4. 面板内 rank → 三种 fuse  

### Phase 4 — 出诊断表

1. 导出 CSV（分子表 + 汇总表）  
2. 画两张图：  
   - dual vs A-only 的 `fuse_min_whole` vs `fuse_min_moiety` 分布  
   - A-only：`S_whole(弱端)` vs `S_moiety(弱端)`  
3. 写 5 行结论：H1/H2/H3 成立与否 + 失败原因  

### Phase 5 — 决策

| 结果 | 下一步 |
|------|--------|
| Go | 扩面板到审计全量四类子集；再上 PIK3CA/mTOR |
| 弱信号 | 加 linked 对照分子；检查 passenger 体积相关 |
| No-Go | 公开 moiety 打分定义与失败案例；必要时改方案 B 作消融，不硬写 NMI |

### 最小命令流（实现后应对齐的接口）

仓库对接脚本尚未齐；目标接口形态：

```bash
# 1) 导出面板
python scripts/export_egfr_her2_panel.py --out data/diag_egfr_her2/panel.csv

# 2) 对接（YAML 协议）
python scripts/dock_panel.py --protocol configs/egfr_her2_vina.yaml --panel data/diag_egfr_her2/panel.csv

# 3) moiety 重打分
python scripts/score_moiety.py --poses results/... --annotations data/diag_egfr_her2/moiety.json

# 4) 诊断表
python scripts/make_diagnostic_table.py --scores results/... --out results/diag_egfr_her2_table.csv
```

当前仓库仅有配对审计脚本 `scripts/audit_public_target_pairs.py` 与 `FROZEN_PUBLIC_PAIRS.yaml`；**上表命令为待实现目标**，跑数前先完成 Phase 0–2 的人工冻结与 TAK-285 QC。

---

## 4. 一页清单（打印用）

- [ ] 为何这对：双端共晶 + 可对接 + 有 A-only  
- [ ] 面板：dual / A-only / B-only 名单冻结  
- [ ] 标注：方案 A + TAK-285 金标准  
- [ ] 协议：3POZ/3RCD、引擎、盒子、RMSD  
- [ ] 列：whole / moiety / passenger / fuse_*  
- [ ] 假设：H1 污染 · H2 错排 · H3 真端保留  
- [ ] 先 QC TAK-285，再跑面板，再出汇总表  

---

*与 [`NMI_SUBMISSION_PLAN_MOIETY.md`](NMI_SUBMISSION_PLAN_MOIETY.md) WP1–WP2 对齐。*
