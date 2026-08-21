# 定位与“框架”用语边界（Intro / Abstract / Cover letter 用）

> 目的：避免把常规评测流程包装成“新算法 Framework”，看起来像 AI 生成论文。  
> 与 `CLAIM_CEILING.md`、`JCIM_PREWRITING_CHECKLIST_V1.md` 配套。

---

## 1. 用户批评是否合理？

**合理，且应采纳。** 原因：

1. **步骤过于工整**（Data audit → Hard-negative → Pocket-matched AUROC → Confounder → Robustness）每一步都是领域已有概念的组合，不是可发表的新算法步骤。写成 “Step 1…Step 5 Framework” 会让审稿人问：这是真方法，还是重新包装？
2. **宏大缩写名**（如 D-DRAF）在没有独立软件产品、标准 pipeline API、大规模社区采用时，容易显得包装；JCIM 评测文通常不轻易给评价流程起缩写。
3. **真正需要的是概念框架，不是算法声明。** 贡献是“系统评价怎么做”，不是“发明了能赢的新对接分”。

---

## 2. 允许 vs 禁止

### 允许（推荐）

| 英文 | 中文意图 |
|------|----------|
| systematic benchmarking framework | 系统评测/基准框架（评价流程） |
| evaluation protocol / benchmarking protocol | 评测协议 |
| DualFourClass-Bench | **公开基准资源名**（数据集 + 标签 + 分数 + 脚本） |
| We established a systematic benchmarking framework to evaluate the reliability of docking-based dual-target recognition. | 建立评价体系 |
| Evaluating the reliability and limitations of docking-based dual-target recognition | 全文定位句 |

### 禁止 / 强烈不推荐

| 写法 | 为什么不行 |
|------|------------|
| Dual-target Docking Reliability Assessment Framework (D-DRAF) | 宏大缩写；像 AI 起名；无独立软件/标准产品支撑 |
| We developed a novel dual-target docking framework named … | “novel + framework + named” = 算法声明；审稿人立刻问 novelty |
| Step 1…Step 5 作为主文方法图的正式“Framework 步骤” | 工整包装；把常规分析写成发明流程 |
| new algorithm / new scoring function / decision framework validated | 越出 claim ceiling |

### 灰色区（可用，但要克制）

- Methods 里用普通编号小节（2.1 Data curation, 2.2 …）**可以**——这是论文结构，不是“Framework Steps”。
- TOC / 示意图若展示流程，用 **evaluation workflow** 或 **analysis workflow**，不要标题成 D-DRAF / Reliability Assessment Framework。
- “ATP recognition framework” 等结构生物学用语保留（指口袋识别模式，不是本文方法名）。

---

## 3. 推荐贡献句（可直接进 Abstract / Intro / Cover letter）

**Preferred:**

> We established a systematic benchmarking framework and the DualFourClass-Bench resource to evaluate the reliability and limitations of docking-based dual-target recognition under pocket-matched directional metrics, hard-negative selective ligands, and physicochemical confounder controls.

**Avoid:**

> We developed a novel Dual-target Docking Reliability Assessment Framework (D-DRAF) consisting of five steps…

**Cover letter 一句版：**

> This Article presents a multi-pair evaluation benchmark (DualFourClass-Bench), not a new docking score; the contribution is a systematic protocol for assessing when docking retains residual dual-target signal and when apparent signal is explained by ligand properties.

---

## 4. 与现有命名的关系

| 名字 | 角色 | 是否保留 |
|------|------|----------|
| DualFourClass-Bench | 基准资源（panels + scores + scripts） | **保留**；这是数据集/基准名，不是算法框架缩写 |
| DualFourClass | 四类任务定义的简称 | 保留 |
| D-DRAF 或同类新缩写 | — | **不要引入** |
| pocket-matched directional AUROC | 主指标名 | 保留（描述性，非产品名） |

---

## 5. Intro 结构建议（防 AI 感）

不要写成五步 Framework 清单。推荐叙事弧：

1. 双靶对接在实践中被使用，但缺少严格四类硬负评测。
2. 公开数据供给限制了可平衡的双靶基准规模。
3. 池化分数可能掩盖方向失败；需要口袋匹配指标。
4. 表观双靶信号常被配体属性/化学型解释。
5. 生成式双靶方法（DualDiff / FuseDiff）用两端 Vina 相对参考配体的 Dual High Affinity 评测成功，与“硬负分臂”不是同一问题；本基准可作其下游诚实评测（零新实验写作点）。
6. 因此本文建立 **systematic benchmarking protocol + DualFourClass-Bench**，在 K=4 上报告可靠边界与局限。

流程细节放 Methods；不要在 Intro 用 “our Framework Step 1–5” 推销。
