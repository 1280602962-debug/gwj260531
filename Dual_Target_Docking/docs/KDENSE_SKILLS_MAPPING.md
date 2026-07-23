# K-Dense Scientific Agent Skills × 本课题可用映射

> 源仓库：[K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills)（~148 skills）  
> 课题：双靶兼容性打分 + Dual-VSDS（冲击 NMI / 保底 JCIM）  
> 日期：2026-07-23

**一句话：** 这套 skills 强在 **造数据、跑化学信息学、对接消融、统计与写稿辅助**；**不能替代**你们已定的「任务重定义 + 打假表」主叙事，也 **没有** GNINA/PoseBusters/RTMScore/Dual-VSDS 专用 skill。

与已有审计分工：

| 工具包 | 擅长 |
|--------|------|
| ResearchStudio | idea 质量 / scoop / 可证伪 |
| academic-research-skills (ARS) | RQ / 方法蓝图 / 诚信闸门 / 审稿路由 |
| **本包 (K-Dense)** | **落地执行：库查询、分子处理、对接批跑、统计、图表、写作** |

---

## 1. 对你当前阶段（Phase 1 数据 + Phase 2 对接）最有用

| Skill | 对本课题的具体帮助 | 优先级 |
|-------|-------------------|--------|
| **`database-lookup`** | 可复现查询 ChEMBL / PubChem / PDB / UniProt；带 provenance，适合建 D2 配对活性骨干 | **P0** |
| **`bioservices`** | 同一脚本跨 ChEMBL↔UniProt↔PDB ID 映射 | **P0** |
| **`rdkit` / `datamol`** | SMILES 标准化、InChIKey、scaffold、相似度、去重、描述符；防泄漏划分的指纹/骨架 | **P0** |
| **`medchem`** | PAINS / Lipinski 等过滤 RandomDecoy 或库筛选；**不当主创新** | P1 |
| **`diffdock`** | 仅作 **消融对照采样器**（你们主引擎仍是 GNINA/Vina）；skill 自己也写明不做亲和力 | P1（消融） |
| **`modal`** | 大规模对接 / DiffDock 批跑需要 GPU 时云端并行 | P1（有预算时） |
| **`pytdc`** | scaffold split、DTI/ADMET 基准习惯；可借鉴划分 API，**勿把单靶 DTI 当双靶主任务** | P1 |
| **`statistical-power`** | 补 ARS 指出的缺口：dual 正例最小 n、效应量、power 曲线（打假表前必做） | **P0** |
| **`statistical-analysis`** | 打假表后的置信区间、配对检验、多重比较 | P1（有结果后） |

---

## 2. 有结果之后（图表 / 写作 / 审稿）

| Skill | 帮助 | 注意 |
|-------|------|------|
| **`scientific-schematics`** | 方法总图、校准→softmin 流程图 | 投稿图需人工把关 |
| **`infographics` / `pptx` / `latex-posters`** | 组会 / 壁报 | 非论文正文 |
| **`markdown-mermaid-writing`** | Dual-VSDS 数据流 mermaid | 已有文档可继续用 |
| **`scientific-writing`** | IMRaD 段落化 | **有 CSV 再写**；禁幻觉数字 |
| **`citation-management` / `paper-lookup` / `bgpt-paper-search` / `exa-search` / `literature-review`** | Related Work 补全、引用核验 | 与 ARS citation-check 重叠，可二选一为主 |
| **`peer-review` / `scholar-evaluation` / `scientific-critical-thinking`** | 投前自审 | 与 ARS reviewer / ResearchStudio 重叠；可作「第三方清单」 |
| **`hypothesis-generation`** | 把打假观察写成可检验假设 | 你们 RQ 已较清晰，收益递减 |
| **`what-if-oracle`** | 「打假失败怎么办 / 只做一对靶点怎么办」情景树 | 战略决策用，不产出实验 |

---

## 3. 慎用 / 易把课题带偏

| Skill | 风险 |
|-------|------|
| **`diffdock` 当主引擎** | 与 VSDS-VD 教训及你们主协议冲突；只作对照 |
| **`deepchem` / `torchdrug` 训通用亲和/DTI SOTA** | 易滑回「再做一个 scorer」——已锁定 **不当主创新** |
| **`deepchem` ADMET 当 PK 主结论** | 你们明确不做 PK 预测模型；最多附录描述符对照 |
| **`primekg` 网络药理学主线** | 变成多靶网络故事，偏离对接兼容性 |
| **`molecular-dynamics`** | 成本高、非 Dual-VSDS 最小充分集 |
| **`literature-review` 强制 AI 示意图** | 可能与你们已有调研文档重复劳动 |

---

## 4. 明确帮不上 / 本包没有的

- **没有** GNINA / AutoDock Vina / RTMScore / PoseBusters / EquiScore 专用 skill  
- **没有**「双靶配对标签 / dual-vs-single / softmin 融合」现成流水线  
- **不能**自动证明 NMI 可中；执行层工具 ≠ 科学主张成立  
- **不能**替代湿实验亲和力（你们也已接受计算主线）

对接主协议仍需你们自己写：YAML + GNINA subprocess + PoseBusters + 校准融合脚本。

---

## 5. 建议调用顺序（贴合现有 Phase）

```text
1. statistical-power     → 定 dual 最小 n / MDE（写进方案）
2. database-lookup
   + bioservices         → 抽 ChEMBL 配对活性 → D2-public
3. rdkit / datamol       → 标准化、InChIKey、scaffold split
4. medchem（可选）       → decoy / 库过滤
5. 自管 GNINA 主对接
   + diffdock（可选消融）+ modal（可选加速）
6. statistical-analysis  → 打假表 + 消融推断
7. scientific-schematics → 方法图
8. scientific-writing
   + peer-review         → 有表之后再写/自审
```

---

## 6. 和「现在该不该装整包」

| 情况 | 建议 |
|------|------|
| 立刻建 D2 + 划分 | 装 **`database-lookup` `bioservices` `rdkit`/`datamol` `statistical-power`** 即可 |
| 要云端大批对接 | 再加 **`modal`**；DiffDock 仅消融加 **`diffdock`** |
| 已有 ResearchStudio + ARS 审计 | **不必**再为「想清楚课题」装 critical-thinking / hypothesis；重复度高 |
| 冲写作 | 等 K1 打假表落地后再开 writing / peer-review |

**结论：** K-Dense 包是你们从「方案」跨到「可复现实验」时最对口的执行技能库；对「课题是否成立 / 是否被 scoop」帮助有限——那两件事已由 ResearchStudio 与 ARS 覆盖。
