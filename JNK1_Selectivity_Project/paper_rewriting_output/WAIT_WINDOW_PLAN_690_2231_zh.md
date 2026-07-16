# 到货前一个月：计算实验与写作规划（采购 = 690 + 2231）

**采购更新（相对此前 690+2157）：** 已订 **HIT103871685 (690)** + **HIT100544184 (2231)**；约 1 个月后到货。  
**湿实验约束不变：** 仅 JNK1/2/3 酶促 IC50；阳性对照 E1、CC-90001 自备。  
**叙事锚点：** PaperSpine **Option A**（家族富集管线 + 选择性预测器失败基准）仍为主贡献；2231 入购后 **RQ-B（JNK1 偏好假说）明显增强**，但仍不能升格为“选择性发现”主标题，除非 IC50 满足预注册 SI 规则。

**Skills 对齐：**
- **ARS experiment-agent / experiment-team**：实验设计、预注册终点、可复现、防 HARKing  
- **PaperSpine**：确认门已过（Option A）→ 本月可写 Intro/Methods/RQ-C Results；IC50 Results 等货  

---

## 0. 采购变更对科学问题的影响（先读）

| RQ | 690+2157（旧） | **690+2231（现）** |
|----|----------------|---------------------|
| **RQ-A** 家族活性富集 | 两枚 grade A，最稳 | 690 = grade A；**2231 = grade C / pass_md_overall 否** → 家族富集叙事变“混合风险” |
| **RQ-B** JNK1 偏好 | 弱（2157 Δsel 为负） | **明显更好**：2231 MD bias #1、hinge J1≫J2、Δsel_dock=+3.37 |
| **RQ-C** 预测器失败 | 已完成（与采购无关） | **不变**，仍是主可发表计算贡献 |
| 风险 | n=2 偏保守 | 2231 可能酶活阴性或 pose 不稳 → 必须预注册“可接受阴性/混合结局” |

**预注册成功标准（建议立刻改锁到 C4）：**
- **Primary (RQ-A):** {690, 2231} 中 ≥1 个在任亚型 IC50 ≤ 10 µM  
- **Secondary (RQ-B):** 仅当 SI_J2≥3 **且** SI_J3≥3 才称 JNK1 preference；**优先看 2231**  
- **Assay validity:** E1 方向大致保留；CC-90001 可测多亚型  
- **Allowed negatives:** 2231 无活性 / 无选择性 / MD 假阳性 → 全部可写进 Option A（方法校准 + 管线边界）

---

## 1. 现状盘点（昨天已完成 vs 因采购变更需重做）

| ID | 内容 | 状态 | 对 690+2231 的动作 |
|----|------|------|-------------------|
| C1 | Chemotype novelty | 做过 690/2157 | **补跑 2231**；主文表改为 690+2231 |
| C2 | Vina 多 seed | 做过 690/2157 | **补跑 2231×3 亚型×≥3 seed** |
| C3 | MD 多 seed 无约束 | 仅协议；2231 有 **带配体约束的 200 ns 单轨** | **本月最高优先级计算**：690+2231 无约束 replica |
| C4 | IC50/SI 预注册 | 锁的是 690/2157 | **改锁化合物列表为 690+2231**（到货前必须完成） |
| C5 | 选择性尸检 | 完成 | 保持；与采购解耦 |
| C7 | PAINS | 已含 2231 | 完成 |
| C11 | 2231 机会成本 | 按“未买”写的 | **改写为“已购 + 为何与 690 配对”** |
| PaperSpine | 蓝图+RQ-C 草稿 | 有 | 本月推进正文；更新购买描述 |

---

## 2. 四周执行日历（到货前）

### Week 1 — 锁定采购叙事 + 轻量补齐（可马上做）

| 任务 | Skills | 产出 | 优先级 |
|------|--------|------|--------|
| 更新 C4：NEW={690,2231}；写清 RQ-A/B 与允许阴性 | ARS | `c4_analysis_lock.json` 新版本 | **P0** |
| C1 加入 2231；汇总 vs E1/CC-90001/SP600125/ChEMBL | ARS C1 | `results/chemotype_novelty/` 更新 | **P0** |
| C2 Vina 多 seed 补 2231 | ARS C2 | `results/pose_consensus/` 更新 | **P0** |
| 重写 C11：从“未买理由”→“购 2231 的 RQ-B 假说 + grade C 风险” | ARS | `C11_*.md` | **P0** |
| PaperSpine：更新 `confirmed_*` 购买集；Intro 一句反映 2231 | PaperSpine | config / draft | P0 |
| 软件合规确认：Glide 是否有授权；无则 Methods 固定 Vina/OpenMM 路径 | 合规笔记 | Methods 决策备忘 | P0 |

### Week 2–3 — 重计算：MD replicas（核心缺口）

| 任务 | 说明 | 产出 |
|------|------|------|
| **C3a 690**：JNK1/2/3 × ≥2 seeds × 20–50 ns，**无配体生产约束** | 补齐 grade A 的可重复性 | `results/md_replicas/690_*` |
| **C3b 2231**：同上；**必须与旧 200 ns 带约束轨迹分开报告** | 旧轨迹仅作“受限探索”；新轨迹才可谈 pose 稳定性 | `results/md_replicas/2231_*` |
| 汇总 hinge occupancy / ligand RMSD mean±SD；重算 pass_md_overall 在 replica 上的翻转率 | 回答“grade C 是否偶然” | `C3_MD_REPLICA_REPORT.md` |
| （可选）短 MM-GBSA per isoform on replica windows | 相对排序，不作选择性证明 | 表 |

**若无 AMBER/Desmond 授权：** 用 OpenMM 或 GROMACS + GAFF2/OpenFF，协议写进 Methods；不要把旧 Schrödinger/带约束结果改名。

### Week 2–3 并行 — 结构解释与富集先验（P1）

| ID | 任务 | 为何现在做 |
|----|------|------------|
| **C8′** | 2231 在 Ile106(JNK1) vs Leu(JNK2/3) 口袋的交互对比（对接+MD 接触频率） | 给 RQ-B 一个可证伪的结构假说 |
| **C9** | 回顾富集：随机抽库 vs 短名单的 score/p_family 分布 | 到货后解释 n=2 hit/miss |
| **C10** | SEA / SwissTargetPrediction（690、2231） | Discussion 脱靶假设；标明非 kinome |
| **C6 预计算（可选）** | 对接 pose 上单点 MM-GBSA 三亚型 | 到货后只做“一致性检查”，不预注册为成功标准 |

### Week 4 — 写作与到货准备（PaperSpine）

| 任务 | 产出 |
|------|------|
| 完成 Intro + Methods + Results(RQ-C 选择性失败) 英文草稿 | 扩写 `draft_intro_methods_rqc_en.md` |
| Figures：漏斗、C5 尸检表、C1 新颖性、C2/C3 QC、690 vs 2231 对照 | 图注清单 |
| 预注册 SI 分析 notebook 试跑空表 + 假数据冒烟测试 | 防止到货后脚本挂掉 |
| 订板/assay 后勤清单：溶剂、浓度梯度、同板 E1/CC-90001、重复孔 | 实验备忘（非计算） |
| contribution_check：确认主贡献仍是 Option A，2231 只升级 **次级假说** | PaperSpine 关卡 |

---

## 3. 计算实验优先级总表（到货前）

### 必须完成（P0）

1. **C4 改锁**到 690+2231  
2. **C1/C2 覆盖 2231**  
3. **C3 无约束 MD replicas（690+2231×3 亚型×≥2 seed）** — 本月最大块  
4. **C5 保持主文**（已有）  
5. **C11 叙事改写**（已购）  
6. **PaperSpine 正文前半**（Intro/Methods/RQ-C）

### 强烈建议（P1）

7. C8′ 2231 异构口袋接触分析  
8. C9 回顾富集模拟  
9. C10 配体脱靶预测  
10. 软件 license 决策写入 Methods  

### 可到货后或并行（P2）

11. C6 按 IC50 结果做能量一致性（**事后**）  
12. C12 活性化合物更长无约束 MD  
13. 若 2231 出 preference：再考虑是否加买 2157 作 chemotype 对照（非必须）

### 明确不要做（浪费窗口）

- 再堆 Δsel / Gly87 库筛选（已证伪）  
- 无湿实验就宣称 kinome 干净  
- 用旧 **带配体约束** 的 2231 200 ns 当选择性铁证  

---

## 4. 与 Skills 的分工（本月怎么用）

| Skill | 本月用途 | 输出物 |
|-------|----------|--------|
| **ARS experiment-agent** | 采购变更后的设计审计、预注册、防 HARKing | 更新版 `ARS_EXPERIMENT_AUDIT.md`；C4 lock |
| **ARS code-runner 思路** | 跑 C1/C2/C3/C9 脚本、冒烟测试 C4 | `results/*` |
| **PaperSpine** | 在确认门之后写 Intro/Methods/RQ-C；IC50 节留空 | `section_blueprints` → 正式草稿 |
| **文献对照（已有 LITERATURE_COMPARISON）** | Discussion 对标 JCIM Kinase-Bench / JAK IFP | 引用不漂移到低分 OA |

---

## 5. 到货当周检查清单（预告）

- [ ] 化合物 ID / CAS / SMILES 复核（2231=HIT100544184；690=HIT103871685）  
- [ ] 填 `results/assay/ic50_raw.csv`  
- [ ] 跑 C4 → 出 SI 表  
- [ ] **禁止**看到数据后改 SI 阈值或把 MD hinge 改称为“证实选择性”  
- [ ] 若仅 690 有活性、2231 无：主文强调 RQ-C + 管线边界；2231 = MD 假阳性案例（很有价值）  
- [ ] 若 2231 满足 SI 规则：报告为 **secondary finding**，仍保留预测器失败主线  

---

## 6. 一句话执行令

**到货前把时间花在：改锁终点（C4）→ 补齐 2231 的 novelty/pose（C1/C2）→ 无约束 MD replicas（C3）→ 写完 Option A 前半正文；不要再做已被证伪的选择性过滤器。**
