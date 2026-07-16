# 无 Schrödinger 条件下：待补实验清单 + Gnina 结果不一致怎么办

**约束：** 无正版 Schrödinger；采购已定为 **690 + 2231**；湿实验仅 JNK1/2/3 IC50；主叙事 Option A。  
**原则：** Glide 漏斗作为**历史筛选记录**保留；投稿 Methods 的可复现对接/MD **以开源为准**（Vina / Gnina / OpenMM 或 GROMACS）。

---

## A. 完整待补实验清单（按优先级）

### A0. 写作/合规（非计算，但必须先定）

| # | 事项 | 状态 | 做法 |
|---|------|------|------|
| W1 | Methods 软件表述 | 待定稿 | 写明：短名单来自既往 Glide 管线；本文 pose/能量验证使用 Vina±Gnina；MD 使用 OpenMM/GROMACS |
| W2 | 不把 Glide 分数改写成 Vina/Gnina | 规则 | 禁止分数“换马甲”；可并列报告 |
| W3 | PaperSpine Intro/Methods/RQ-C | 有草稿 | 按开源路径改对接段 |

---

### A1. P0 — 到货前必须完成

| ID | 实验 | 开源工具 | 目的 | 状态 | 交付物 |
|----|------|----------|------|------|--------|
| **C1** | Chemotype novelty（690, 2231） | RDKit ECFP4/Murcko | 回答是否已知 JNK 近邻 | **已完成** | `results/chemotype_novelty/` |
| **C2** | 多 seed 对接共识 | **AutoDock Vina**（已跑）；建议再加 **Gnina CNN rescore** | pose 是否偶然 | Vina **已完成**；Gnina **待补** | `results/pose_consensus/` |
| **C3** | 无约束 MD replicas | **OpenMM 或 GROMACS** + GAFF2/OpenFF；≥2 seed × 3 亚型 × 20–50 ns | 替代单次/带约束 MD；检验 2231 grade C | **未做（最高优先级）** | `results/md_replicas/` |
| **C4** | IC50/SI 预注册分析 | 已有脚本 | 防 HARKing | **已锁 v2（690+2231）** | `results/assay_analysis/` |
| **C5** | 选择性预测器尸检表 | 已有归档数据 | RQ-C 主贡献 | **已完成** | `results/selectivity_autopsy/` |
| **C7** | PAINS/理化风险 | RDKit FilterCatalog | 排除assay 假象 | **已完成** | `results/purchase_risk/` |
| **C11** | 690 vs 2231 购买理由 | 表+叙事 | Discussion | **已完成** | `results/c11_2231_comparison/` |

**P0 真正还欠的核心：C3（MD）+ C2 的 Gnina 层（可选但强烈建议）。**

---

### A2. P1 — 强烈建议（到货前或到货后立刻）

| ID | 实验 | 开源工具 | 目的 |
|----|------|----------|------|
| **C2b** | 对 **690/2231/E1/CC-90001** 做 Gnina 局部精修 + CNN 打分（固定 Vina 盒子） | Gnina | 与 Vina 交叉验证 pose；不重筛百万库 |
| **C2c** | 共晶 redock：3ELJ / 3E7O / 3TTI | Vina+Gnina | 证明开源协议口袋可用（审稿高频） |
| **C8′** | 2231 在 JNK1 Ile106 vs JNK2/3 Leu 接触频率 | 对接+MD 分析 | 给 RQ-B 可证伪结构假说 |
| **C9** | 回顾富集：随机库 vs 短名单 score 分布 | RDKit + 已有分数 | 解释 n=2 的先验 |
| **C10** | SEA / SwissTargetPrediction | 网页/本地 | 脱靶假设（非 kinome） |
| **C5+** | 把 C5 主文表/图注按 JCIM 语气定稿 | 写作 | 与开源 Methods 对齐 |

---

### A3. P2 — 有算力再做 / 到货后

| ID | 实验 | 说明 |
|----|------|------|
| **C6** | 按 IC50 结果做 MM-GBSA/短 MD 能量一致性 | **事后解释**，不作选择性证明 |
| **C12** | 活性化合物更长无约束 MD（100 ns+） | 机制图用 |
| **C3-ext** | 2231 旧 200 ns（有配体约束）与新无约束轨迹对比说明 | 明确旧结果局限 |

---

### A4. 明确不做

- 用 Gnina/Vina **重跑整库**并改掉已订的 690/2231  
- 再堆 Δsel / Gly87 作为采购门  
- 无授权却在正文写 Glide/Prime/Desmond 细节装可复现  
- 用 MD hinge 不对称直接宣称“已证实 JNK1 选择性”

---

## B. 推荐开源技术栈（写进 Methods 的一套）

| 步骤 | 工具 | 备注 |
|------|------|------|
| 指纹/新颖性/PAINS | RDKit | 已用 |
| 粗对接 / 多 seed | AutoDock Vina 1.2.x | C2 已用 |
| 精修 + CNN 打分 | **Gnina** | 补 C2b；固定与 Vina 相同盒子 |
| 蛋白准备 | meeko / PDB2PQR | 已部分使用 |
| MD | OpenMM 或 GROMACS | C3 主战场 |
| 配体力场 | GAFF2 / OpenFF | 与 AmberTools `antechamber` 或 OpenFF toolkit |
| 能量（可选） | OpenMM 或 gmx_MMPBSA | 相对排序 |

历史 Glide：在 Methods 写一句  
> “Initial commercial-library triage used Glide XP in a prior institutional campaign; purchased IDs and funnel counts are archived. Independent open-source docking/MD were used for pose consensus and stability analysis reported here.”

---

## C. 若换成 Gnina，筛选结果不一样——怎么办？

### C.1 先分清“哪一种不一样”

| 情况 | 含义 | 正确处理 |
|------|------|----------|
| **① 仅对已购分子**：Gnina 分数/排序与 Vina 或历史 Glide 不一致 | 正常，评分函数不同 | **报告三者并列**；以 pose RMSD 共识为主，不以分数决胜 |
| **② Gnina 认为 2231 在某亚型 pose 很差**（如 JNK2 已 Vina 不稳） | 与 grade C 一致 | 强化“2231 高风险假说”；C3 MD 看是否站得住；**不改采购**（已订货） |
| **③ 若重筛整库，Top 列表与 Glide 短名单重叠很低** | 预期内（不同引擎） | **不要据此换分子**；可在 SI 做小规模重叠分析作局限讨论 |
| **④ Gnina CNN 分与实验 IC50 事后不符** | 常见 | 讨论评分局限；回归 Option A（方法失败/校准） |

### C.2 硬规则（避免审稿灾难）

1. **已订 690/2231 = 锁定实验对象**。开源结果再差，也是对这套分子做验证/证伪，不是重新选苗。  
2. **Gnina 只做“验证层”**，不做“新一轮发现层”（除非你愿意开新项目、新采购）。  
3. 正文区分：  
   - *Selection*（历史 Glide 漏斗）  
   - *Confirmation*（Vina/Gnina/MD）  
4. 若 Gnina 与 Vina 对同一分子 top pose 重原子 RMSD **> 2 Å** 且多 seed 不一致 → 该亚型 **不画“确定结合模式”主图**，改写“pose ambiguous”。  
5. 若将来整库 Gnina 重筛：最多作为 SI“引擎敏感性”，**结论不得改成“真正 hit 是另一批未购分子”**（除非补购并补测）。

### C.3 实用决策树

```
对 690/2231 跑 Gnina（同盒子）
├─ 与 Vina top pose RMSD ≤ 2 Å（多数 seed）
│   → 写“开源共识 pose”；可进主图
├─ RMSD 大 / 分数大乱 / 仅 JNK2 乱（2231 已现）
│   → 标明亚型依赖不稳；靠 C3 MD；IC50 仍测三亚型
└─ 想看“若当初用 Gnina 会买谁”
    → SI 小样本回顾即可；不改主线、不改已购集
```

### C.4 和你当前数据的衔接

- Vina C2：**690 三亚型共识均过**；**2231 在 JNK2 未过**。  
- 预期 Gnina 很可能同样挑出 2231 的 off-isoform 不稳 → 这是**支持风险披露**，不是否定购买。  
- 购买逻辑本来就是：690 = 稳的家族锚；2231 = 高风险高信息量的偏好假说。

---

## D. 建议执行顺序（无 Schrödinger，从现在到到货）

1. **装 Gnina**，对 690/2231/E1/CC-90001 + 三 PDB 做 C2b（固定 Vina 盒子；多 seed）  
2. **C2c redock** 共晶配体（Vina+Gnina RMSD）  
3. **C3 MD replicas**（OpenMM/GROMACS）——本月最大块  
4. 定稿 Methods“历史 Glide + 开源验证”双层表述  
5. 扩写 PaperSpine Intro/RQ-C  
6. 到货 → 填 C4 → 按预注册出 SI；Gnina/MD 只解释、不改终点  

---

## E. 一句话

**没有 Schrödinger 完全可以发**：把 Glide 当历史筛选，用 Vina/Gnina/OpenMM 把验证做扎实。  
**Gnina 和 Glide/Vina 不一致是常态**——用来做敏感性与风险披露，**不要用来推翻已订的 690/2231**，除非你准备新开一轮发现+新湿实验。
