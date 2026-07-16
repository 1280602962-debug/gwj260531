# 下一步行动清单（基于 AF3 vs 薛定谔 Phase 0）

**依据：** `jnk2_af3_vs_schrodinger/AF3_vs_薛定谔共价对接对比总结.md`  
**贡献约束：** Framework 必达 + Chemistry 升级（PaperSpine confirmed B）  
**日期：** 2026-07-16

---

## 预实验已回答的问题

| 问题 | 答案 |
|------|------|
| Phase 0 AF3 gate 过了吗？ | **过了**：mPAE AUC=1.0，EF@1%=36.8，0% decoy 压过 active |
| Glide 能否单独当主筛？ | **不能**：GScore AUC=0.957 但 Cohen's d≈0.09；Top15 中 13 个 decoy；56d anchor 假阳性高 |
| AF3 与 Glide 能否互换排序？ | **不能**：mPAE–GScore Spearman ρ≈0.27；Top15 Jaccard≈0.20 |
| 框架主指标选谁？ | **AF3 mPAE 主门控**；Glide 只做 pose/相互作用辅助 |

→ **校准阶段可以收束**，下一步转向 **固化 SOP + Layer1 化学库**，不要继续堆 decoy 做“更好看的 AUC”。

---

## 立即做（本周，Framework 必达）

### 1. 冻结筛选 SOP（写进项目，不再改主规则）

```
候选库
  → AF3 Cys116 共价共折叠，按 mPAE 排序
  → Gate：mPAE ≤ 阈值（建议先用“优于 decoy 中位数 / 或 ≤~1.05–1.20 Å”试点）
  → 取 Top N（如 100–500）
  → 薛定谔 Covalent Docking：只做几何/相互作用 QC（非唯一排名）
  → 剔除：无合理 Cys 几何、或仅 Glide 好看但 mPAE 差的 56d 样假阳性
  → 可采购/合成清单
```

**硬规则**
- 主排序 = **mPAE**，不用 ipTM/ranking_score 单独做 gate（有 ~18/143 decoy ipTM≥0.85 假阳性）
- **禁止**单独用 Glide 筛 56d 邻域
- 联合策略可用：`Glide Top ∩ mPAE 过线`

### 2. 归档可复现包（半页清单即可）

- [ ] 最终 AF3 输入表（修过 YL2056 bonded pairs）  
- [ ] decoy 生成参数（性质匹配 + Tc 阈值）  
- [ ] 本目录对比表 + 本地原始路径说明（已有）  
- [ ] 一页 `SCREENING_SOP.md`（上面漏斗）  
- [ ] Cys 编号表：UniProt 116 vs 激酶域 FASTA 113  

### 3. 把本对比挂到贡献证据

在报告/开题中写一句：

> Phase 0 在 8ELC/Cys116 上，AF3+mPAE 达到 COValid 风格富集（AUC=1，EF@1%=36.8），薛定谔共价对接 enrichment 较弱且与 mPAE 弱相关，故框架采用 AF3 门控 + Glide 辅助。

---

## 下一步主战场（1–4 周，Chemistry 升级开始）

### 4. 建 Layer1 邻域库（500–2000）

来源优先级：
1. YL5084 / YL2056 类似物（flag methyl、吡咯烷立体）  
2. 56d / 4WHZ ligand-first 氨基吡唑邻域（**必须走 AF3，勿只靠 Glide**）  
3. JNK-IN-8 骨架扩展（作 pan/对照，非选择性主线）  

过滤：丙烯酰胺 SMARTS、去过反应、Ro5/可采购优先（Enamine REAL/库存）。

### 5. 漏斗跑通一次小规模（证明框架可复现）

- 先 200–500 分子（WSL 可准备库；AF3 可分批上云/server）  
- 输出：mPAE 排序表 + Glide QC 表 + Top20 可采购清单  
- 成功标准：流程能从头跑通并产出表格，不要求一次出神药  

### 6. 湿实验并行规划（有 Top 候选后再启动）

- kinact/KI（JNK2 vs JNK1）  
- C116S  
- washout / 细胞 engagement  
- **不做** DSS 疗效作主 claim  

---

## 明确不要做

| 不要 | 原因 |
|------|------|
| 再扩大 Phase 0 decoy 到几千刷 AUC | 校准已充分，边际收益低 |
| 用 Glide GScore 单独决定采购 | 假阳性高，尤其 56d |
| 用 ipTM 单独 gate | 假阳性多于 mPAE |
| 百万库盲筛叙事 | 与确认贡献不符；算力与可复现性差 |
| 疾病模型当下一里程碑 | 贡献边界已禁止作主 claim |

---

## 建议的“本周完成定义”

1. 仓库内有 `SCREENING_SOP.md`（AF3 主 / Glide 粗筛分层）  
2. 已定稿完善漏斗：`LIBRARY_FUNNEL_COMPLETE.md`（ChEMBL 胺装弹头 + 双轨 + Glide→AF3）  
3. Layer1 种子列表 ≥200 SMILES（可先虚拟）或商业 acrylamide 试点集就绪  
4. 开题/报告写入 Phase 0 对比结论 + claim boundary  
5. 列出 Top 想采购/合成的 5–10 个化学型方向（不必已有活性）  

完成以上，即可进入 **第一次真实缩库筛选**（建议：商业 ~4k 快车道并行 + ChEMBL 枚举主轨道）。

## 流程完善结论（2026-07-16）

用户提出的「ChEMBL 类药 → 胺 → 丙烯酰胺 → 合成性/YL5084·新骨架过滤 → 薛定谔粗筛 → AF3 精排」**可以采用**，且比“只下载几千商业库”更适合新骨架目标。约束是：

- 2.9M 是**底物池**，不是对接规模；
- Glide 只做**松阈值缩库**，最终以 **mPAE** 定名单；
- 不建 BTK 库，建 **JNK2 Cys116 丙烯酰胺发现库**。
