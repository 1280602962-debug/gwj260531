# SCI 论文大纲 — URAT1 三态对接 Benchmark（修订版）

> **状态**：2026-06 主投稿大纲  
> **替代**：`MANUSCRIPT_OUTLINE.md`（TAPE-GATE 双靶算法版，已降级）  
> **战略说明**：`PAPER_PIVOT_BENCHMARK.md`  
> **执行计划**：`URAT1_THREE_STATE_BENCHMARK_PLAN.md`

---

## 文章类型与目标期刊

| 项目 | 内容 |
|------|------|
| 类型 | Original Research — **Computational benchmark / critical evaluation** |
| 首选 | *Journal of Chemical Information and Modeling* |
| 备选 | *Journal of Cheminformatics* |
| 字数 | 5500–7500 words + SI |

---

## 标题（推荐）

**Benchmarking URAT1 inhibitor docking across inward, occluded, and outward cryo-EM conformations: limitations of rigid Glide ensembles and a pose-transfer rescoring alternative**

备选：
- *When does three-state docking work for SLC transporters? A URAT1 case study*
- *Critical evaluation of cryo-EM conformational ensembles for URAT1 virtual screening*

---

## Abstract（~250 words 结构）

1. **Background**：URAT1 为 alternating-access 转运体；2024–2025 年多构象 cryo-EM 结构可用；抑制剂应对 inward/occluded 有偏好、outward 相对不利。
2. **Gap**：领域缺乏对「刚性三态 Glide 系综」在药物分子上可行性的系统检验；PDB 映射存在 9JDZ 等误用风险。
3. **Methods**：比较 Protocol A（单态 9DKB）、B（刚性三态）、C（inward dock + 叠合 pose 转移 + Prime/MM-GBSA）、D（IFD 敏感性）；面板含 4 个共晶抑制剂、文献 benchmark 与 property-matched decoys。
4. **Results**（据实填写）：9DKB 对接成功；**9B1K/9B1L 刚性对接对药物零 pose**；Protocol C 下四药 $S_\pi$ 方向性 X/4；A vs C enrichment …
5. **Conclusions**：刚性三态 Glide 不能直接推广至 URAT1 药物筛选；提供开源子集与可复现 rescoring 流程；讨论 transporter vs kinase 对接范式差异。

**禁止写法**：dual-target discovery、novel AI algorithm、Teacher distillation。

---

## 1. Introduction（~1200 words）

### 1.1 URAT1 与痛风治疗背景
- 高尿酸血症；URAT1 底物转运与抑制剂机制

### 1.2 Alternating access 与三态结构药理学
- Dai 2024；Suo 2025；Wu 2025 与 PDB 公开现状
- **明确贡献预告**：9DKB / 9B1K / 9B1L 正确映射

### 1.3 虚拟筛选如何误用转运体结构
- 激酶式单 PDB / 刚性系综的常见假设
- 9JDZ 被误标为 occluded/outward 的问题

### 1.4 本文问题与贡献（4 bullets）
1. 三态 PDB 映射与对接协议澄清  
2. 刚性 vs pose-transfer rescoring 系统比较  
3. Gate 四药 + decoy enrichment 定量  
4. 开源 benchmark 子集与 Schrödinger 工作流  

### 1.5 文章结构说明

---

## 2. Materials and Methods（~2500 words）

### 2.1 Structures and grid preparation
- 9DKB, 9B1K, 9B1L 来源与预处理（Protein Prep Wizard）
- Grid 中心：lesinurad vs urate 位点
- 备用结构 9JDZ 仅作 SI 讨论

### 2.2 Compound panels
- Gate 四药 SMILES（`teacher_gate_qc_panel_b_direction.csv`）
- `literature_benchmarks.csv` Tier1a
- Decoy 生成：从 `distill_subset_d.csv` 子采样 + 理化性质匹配规则

### 2.3 Docking protocols A–D
- **A**：9DKB SP→XP  
- **B**：三态刚性 SP→XP，记录 pose viability  
- **C**：叠合、pose 复制、Prime、MM-GBSA、clash penalty  
- **D**：9B1K IFD（4 药）

### 2.4 Scoring functions
- GlideScore（A/B 成功 case）  
- MM-GBSA $\Delta G$（C）  
- Boltzmann $\pi_s$；$S_\pi = \pi_{in}+\pi_{occ}-\pi_{out}$

### 2.5 Evaluation metrics
- Redock RMSD；pose viability；$S_\pi$ sign；enrichment@k；ROC-AUC  
- 统计：效应量 + bootstrap CI（jcim.5c01609）

### 2.6 Software and reproducibility
- Schrödinger 版本；GitHub 仓库；`utils_three_state_scoring.py`

---

## 3. Results（~2000 words）

### 3.1 Structure alignment and pocket consistency
- **Fig 1**：三态示意图 + 叠合（标注 B1K/B1L global RMSD）
- 讨论 pocket 视觉重叠 vs 全局构象差

### 3.2 Protocol B：刚性三态 Glide 可行性
- **Table 1**：四药 × 三态 pose 数、GlideScore、失败 log 摘要
- **关键句**：药物在 9B1K/9B1L **0 poses**；urate 对照成功

### 3.3 Protocol A vs C：已知抑制剂
- **Fig 2**：四药 $\Delta G$ 三态柱状 + $S_\pi$  
- Gate 2：X/4 满足 $S_\pi>0$

### 3.4 Enrichment against decoys
- **Fig 3**：enrichment 曲线（A vs C）  
- **Table 2**：四药排名、AUC

### 3.5 IFD sensitivity (Protocol D)
- **Fig S2**：4 药 occluded IFD 是否挽救 pose

### 3.6 Gate 3：活性集 vs decoy 分离
- median($\pi_{in}+\pi_{oc}$) 比较 + 效应量

---

## 4. Discussion（~1200 words）

### 4.1 为何刚性三态 Glide 对 URAT1 药物失败
- 口袋几何、clash、诱导契合需求  
- 与 P-gp / 其他 transporter 文献对比（IOMemP 等多构象思路）

### 4.2 Pose-transfer rescoring 的合理性与局限
- 文献支持：ensemble rescoring 常见；**URAT1 专用协议无先例**  
- Sindt 2025：rescoring 在 hit list 上常失效 — 你的 decoy 结果如何呼应

### 4.3 对虚拟筛选实践的启示
- 何时单态 9DKB 足够  
- 何时必须柔性/MD/实验

### 4.4 PDB 策展与社区 benchmark 需求

### 4.5 局限性
- 无湿实验；叠合 RMSD 大；decoy 规模；MM-GBSA 近似  
- **明确不写**：NLRP3 双靶、Teacher ML

### 4.6 Future work
- 全库 rescoring；诱导契合 MD；细胞 IC50 验证

---

## 5. Conclusions（~200 words）

- 刚性三态 Glide 不是 URAT1 药物筛选的即插即用方案  
- Protocol C 作为可检验的工程替代  
- 开源 benchmark 供社区复现与扩展  

---

## 图表清单

| ID | 内容 |
|----|------|
| **Fig 1** | URAT1 alternating access + 三 PDB 叠合 |
| **Fig 2** | Protocol B pose viability 热图 / 条形图 |
| **Fig 3** | Protocol C：四药 $\Delta G$ 与 $S_\pi$ |
| **Fig 4** | Enrichment@k：A vs C |
| **Table 1** | 协议定义与计算成本对比 |
| **Table 2** | 四药 + benchmark 化合物完整分数表 |
| **Fig S1** | lesinurad redock（9DKB） |
| **Fig S2** | IFD 敏感性（4 药） |
| **Table S1** | Decoy 子集 SMILES 与性质分布 |
| **Table S2** | Glide 失败 log 摘录 |

---

## 投稿前 Checklist

- [ ] Gate 1：lesinurad RMSD ≤ 2.0 Å
- [ ] Gate 2：Protocol C 四药 $S_\pi>0$（目标 4/4）
- [ ] Protocol B 失败率已表格化（含 urate 阳性对照）
- [ ] A vs C enrichment 至少 50 decoys
- [ ] GitHub 子集与 `utils_three_state_scoring.py` 可运行
- [ ] 摘要与 title **无** dual-target / novel algorithm 措辞
- [ ] 引用 jcim.5c00730、jcim.5c01609
- [ ] 旧 TAPE-GATE 图表 **未** 混入主文

---

## 明确不纳入主文的内容

| 内容 | 处理 |
|------|------|
| NLRP3 ML / assay-conditioned | 删或独立第二篇 |
| OAT 迁移 (+0.004 ρ) | 不写或脚注 |
| Path B 生成式 | Future work |
| Teacher M-CPDL 8973 | Gate 全过后 SI 预告 |
| PLK1-style baseline | 非 URAT1 benchmark 必要 |
