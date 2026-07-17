# ARS × PaperSpine 对照整理（对接前初筛）

**日期：** 2026-07-17  
**目的：** 回答「两套 skill 结论是否一样、何处互参」

---

## 1. 一句话对照

| | **PaperSpine** | **ARS-Codex deep-research** |
|--|----------------|------------------------------|
| **角色** | 立项/论文贡献骨架 + 索赔边界 | 文献系统检索 + 方法可证伪 + 实验计划 |
| **对 Step0–7** | 工程 SOP，**非方法创新** | **REVISE**，不能当理性“活性初筛” |
| **主缺口** | 位点校准 triage（H）+ JNK1 选择性（I） | core-Tc / ErG / 反应性 / **先校准再阈值** |
| **该听谁的操作** | 贡献怎么写、别吹什么 | **本地 WSL 怎么跑**（L0–L7） |
| **该听谁的边界** | 疾病不作主 claim；Framework 必达 | Devil’s advocate：无校准=不理性 |

**结论：核心科学判断高度一致；产品形态不同，值得合并成一套「贡献叙事 + 执行 SOP」。**

---

## 2. 相同点（已交叉验证，可信）

1. **ECFP + 固定 Tc + Murcko 双轨 alone ≠ 创新 / ≠ 活性预测**  
2. **胺→丙烯酰胺→QC→缩库** 是领域常见模块，不是你们的卖点  
3. **Sim / Novel 必须分仓分榜**；JNK-IN-8 作 pan 对照  
4. **AF3 mPAE = 最终门控；Glide = 松阈值粗筛**（绑定 Phase0）  
5. **要升级就补：** 去弹头相似、scaffold-hop 表征（ErG/药效团）、反应性、阈值敏感性/UQ、湿实验闭环  
6. **文献锚点一致：** Lu/YL5084·8ELC、Wydra/56d、JNK-IN-8、COValid/AF3、反应性/warhead 精度  

---

## 3. 不同点（互补，不是打架）

| 维度 | PaperSpine 更强 | ARS 更强 |
|------|-----------------|----------|
| **问题框定** | 贡献类型 B/C、claim boundary、疾病叙事禁令 | FINER RQ、inclusion/exclusion、OpenAlex 可复现检索 |
| **文献广度** | 精选高权重地图 + 项目内 Phase0 | 7 查询系统检索 + 去重语料 JSON |
| **批判强度** | “勿作主创新 / gap H·I” | 正式 **REVISE** + 证伪条件清单 |
| **可执行性** | 原则与升级包（P0–P3） | **L0–L7 配额、thresholds.json、WSL 粘贴指令** |
| **选择性** | 明确列为缺口 **I**（JNK2 vs JNK1） | 强调需负筛，但执行细节少于 PaperSpine 的贡献对齐 |
| **论文写作** | `confirmed_contribution` 直接可写 Intro/边界 | 偏方法附录与实验设计，不组织全文 |

---

## 4. 值得互相借鉴的具体条目

### 从 PaperSpine → 补进 ARS 执行计划
1. **贡献分层：** L7 名单 = Framework 验证材料；kinact/KI = Chemistry 升级——写进 L7 summary 的“用途标签”。  
2. **缺口 I：** 在 L7 后或 AF3 过线名单上加 **JNK1 负筛设计**（即使首轮只是计划表）。  
3. **索赔话术：** 对外/对本地 agent 禁止说“活性预测”，统一 enrichment/triage（与 claim boundary 一致）。  
4. **双源库：** 商业丙烯酰胺快车道 vs ChEMBL 胺枚举——ARS 计划可并行一条小对照，避免只绑枚举叙事。

### 从 ARS → 补进 PaperSpine 框架叙事
1. **L5 校准门：** 把“thresholds.json + leave-one-out + core vs full Tc”写成 Framework 的 **Layer0 扩展**（对接前校准，不只 AF3 回顾校准）。  
2. **可证伪成功标准：** 分层抽样不优于随机合格丙烯酰胺 → 证伪 triage——可写进 methods 的 negative control 逻辑。  
3. **ErG 进 Novel 仓：** 把 gap H 从口号落成 ARS 的分仓规则。  
4. **reactivity_bucket：** 支撑“共价专属”而不只是 kinase-likeness。

---

## 5. 合并后的推荐工作流（两套 skill 合一）

```text
PaperSpine 定：做什么贡献、不说什么话
        ↓
ARS 定：L0–L7 怎么跑、何为校准通过
        ↓
本地 WSL：先 L0–L5 → 你确认 thresholds
        ↓
L6–L7 → Glide 松筛 → AF3 mPAE
        ↓
PaperSpine claim：Framework 报告 +（可选）Chemistry 湿实验
```

**统一优先级（双方一致）：**  
P0 执行 ARS L0–L5（core-Tc + 校准）→ P0.5 ErG Novel → P1 反应性分桶 → P2 JNK1 负筛（PaperSpine I）→ P3 湿实验分层验证。

---

## 6. 文件索引

| 来源 | 路径 |
|------|------|
| PaperSpine 贡献 | `confirmed_contribution.md` |
| PaperSpine SOTA/缺口 | `sota_gap_map.md`（含初筛 addendum） |
| PaperSpine 初筛深研 | `prefilter_pipeline_deep_research.md` |
| ARS 摘要 | `ars_deep_research/EXECUTIVE_BRIEF.md` |
| ARS 重规划 SOP | `ars_deep_research/phase4_plan/prefilter_replan_v1.md` |
