# JCIM Article 写作前注意事项清单 + 逐项核对

> DualFourClass-Bench · 目标类型：**Articles**（全文研究论文）  
> 核对日期：2026-07-30  
> 依据：JCIM Author Guidelines（2026-04-15）、2015 editorial ([acs.jcim.5b00180](https://pubs.acs.org/doi/10.1021/acs.jcim.5b00180)）、2021 Data/Software Availability ([acs.jcim.0c01389](https://pubs.acs.org/doi/10.1021/acs.jcim.0c01389)）、本仓库 `CLAIM_CEILING.md` / `JCIM_P0_COMPLETION_GUIDE.md`  
> 状态：✅ 已就绪 · 🟡 材料有、稿未成 · ❌ 缺口 · ➖ 不适用

---

## 总览

| 区块 | ✅ | 🟡 | ❌ |
|------|---:|---:|---:|
| A. 类型与主张 | 5 | 1 | 0 |
| B. 2015 编辑门槛 | 4 | 2 | 1 |
| C. 2021 数据/软件 | 3 | 2 | 2 |
| D. Author Guidelines 格式 | 1 | 3 | 5 |
| E. 章节与证据链 | 3 | 4 | 4 |
| F. 图/表定稿 | 1 | 2 | 4 |
| G. Claim ceiling 安全 | 4 | 1 | 0 |
| **合计（去重前粗计）** | **21** | **15** | **16** |

**写作前阻塞项（必须先清或边写边盯）：**  
1. ❌ Zenodo DOI + Data and Software Availability 定稿句  
2. ❌ 完整英文 IMRaD（现有 Results + Introduction + Methods 英稿；Discussion 仍以中文工作稿为主）  
3. ❌ Cover letter  
4. ❌ Fig1 schematic / Fig2 supply / Fig4–5 定稿；现有 forest/gate 图 DPI 偏低（~180）  
5. 🟡 Methods 中明确 **Data curation** 小节（2.1 已写；英文稿已有）

非阻塞但建议开写前钉死：主指标措辞、K=4 边界、禁止“通吃 scorer”。

---

## A. 稿件类型与主张边界

| ID | 事项 | 官方/内部依据 | 状态 | 核对记录 |
|----|------|---------------|------|----------|
| A1 | 投稿类型定为 **Articles** | Author Guidelines · Manuscript Types | ✅ | 评测/基准全文研究；非 Letter/App Note/Review/Viewpoint |
| A2 | Cover letter 写明 evaluation/benchmark，非新打分函数、非无湿实验单靶筛药 | 2015 RAER + Author Guidelines Cover Letter | ❌ | 尚无 cover letter 文件 |
| A3 | 主张钉死：口袋匹配方向 AUROC + 平凡基线；K=4 冻结；EGFR=供给案例 | `CLAIM_CEILING.md` | ✅ | ceiling 文件齐全；Results 英稿遵守 |
| A4 | 禁止：通用决策臂 / RTM 通吃 / LigPrep 混主表 / Track B 选臂胜利 | `CLAIM_CEILING.md` | ✅ | Results 无越界措辞 |
| A5 | 期刊范围：不是“单靶对接无实验验证”应用稿 | Scope + 2015 | ✅ | 叙事为 multi-pair metric evaluation |
| A6 | 题名草案 ≤12 词；摘要 3–4 句 | Author Guidelines Title/Abstract | 🟡 | 未写 Title/Abstract；开写时按此限 |

**建议题名方向（草稿，未定）：**  
`DualFourClass-Bench for Dual-Target Docking Evaluation`（约 7 词）  
备选（更贴定位句）：`Evaluating Docking Reliability for Dual-Target Ligand Recognition`（约 8 词）

**全文定位句（Introduction/Abstract 用，避免绝对化）：** 不写成 "Docking can/cannot identify dual-target ligands"；改写为 *Evaluating the reliability and limitations of docking-based dual-target recognition*——强调这是对现有对接分数的可靠性评测，而非提出新方法或做全有全无判决，与 Results 开篇的定位句（`RESULTS_DRAFT_ZH_JCIM_V1.md` 顶部）一致。问题定义句（Intro 末段已写入；Abstract / Cover letter 围绕同一思想改写，勿三处逐字重复）：*The central question is not whether a ligand can obtain favorable docking scores at two targets, but whether docking can distinguish experimentally characterized dual-active ligands from target-selective hard negatives across both target directions.* Introduction **不写死 K = 4**。

**“框架”用语（防 AI 包装感；详见 [`POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`](POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md)）：**
- ✅ 可以说：*systematic benchmarking framework* / *evaluation protocol* + 资源名 **DualFourClass-Bench**
- ✅ 推荐句：*We established a systematic benchmarking framework to evaluate the reliability of docking-based dual-target recognition.*
- ❌ 不要：D-DRAF / “Dual-target Docking Reliability Assessment Framework” / *We developed a novel … framework named …*
- ❌ 不要：Intro/Abstract/TOC 把常规分析写成 Framework Step 1–5（Methods 普通小节编号可以）

---

## B. 2015 Editorial（审稿门槛，非排版）

| ID | 事项 | 状态 | 核对记录 |
|----|------|------|----------|
| B1 | 有 **data curation** 专节（标签规则、过滤、面板建造） | 🟡 | J0 审计 + panel 建造材料齐；**Methods 未成文** |
| B2 | 至少一套公开数据集 | 🟡 | 分数/面板在 git；**Zenodo 未发布** |
| B3 | 方法可理解、可复现（不必全开源，须写清） | ✅ | 协议 YAML、脚本、ENV_PIN 在仓 |
| B4 | 多任务 / 多靶对 benchmark，非单终点套路 | ✅ | K=4 + 供给审计 49 对 |
| B5 | 与平凡基线 / 既有通道对照，不只报对接分 | ✅ | heavy/MW/cLogP/TPSA + Vina/RTM/GNINA |
| B6 | 专有数据若用须论证；否则公开 | ✅ | 主线用 ChEMBL 公共数据 |
| B7 | Cover letter 回应 RAER 点（范围、推进、非套路） | ❌ | 同 A2 |

**2015 Table 1 对本稿的读法：** 写成 *methodological best practices + multi-pair benchmark + public data*，不要写成 *new scorer that beats everything*。

---

## C. 2021 Data and Software Availability（自 2021-01-01）

| ID | 事项 | 状态 | 核对记录 |
|----|------|------|----------|
| C1 | 文末独立 **Data and Software Availability** 节 | ❌ | 英文全稿未组；Results 无该节 |
| C2 | 数据机器可读（结构、标签、分数） | 🟡 | `assembled_all_pairs_long.csv`、panel CSV、主表齐；待 Zenodo 打包 |
| C3 | 验证用筛选条件可精确复现 | ✅ | J0 规则、θ、strict 配额有文档/表 |
| C4 | 新指标/流程：伪代码或工作流 + 参数机器可读 | ✅ | `PRIMARY_METRIC_V2.md` + build 脚本 |
| C5 | 第三方软件：名称、版本、参数 | 🟡 | `ENV_PIN.md` 有；Methods 须抄成正文句 |
| C6 | Zenodo（或同等）DOI 写入正文与 README | ❌ | README：`Zenodo DOI：（发布后填这里）` |
| C7 | 姿态：可后补；分数包可先发 | 🟡 | `POSE_UPLOAD_CHECKLIST.md` 已列；姿态在本地未进 git |

**复现命令（已就绪，待写入 Zenodo README / Methods）：**

```bash
cd Dual_Target_Docking
python3 data/jcim_bench_v0/scripts/build_pocket_matched_diagnostics_v1.py
python3 data/jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py
python3 data/jcim_bench_v0/scripts/plot_forest_ci_v1.py
```

---

## D. Author Guidelines 格式（Articles）

| ID | 事项 | 状态 | 核对记录 |
|----|------|------|----------|
| D1 | Fast Format：标准章节齐全、无批注高亮 | ❌ | 仅有 Results 草稿，无完整稿 |
| D2 | 图/表嵌在相关正文旁（首投） | ❌ | 待组稿 |
| D3 | Title ≤12 词；Abstract 3–4 句 | ❌ | 未写 |
| D4 | TOC / Abstract graphic | ❌ | 未做 |
| D5 | References 完整且含标题 | ❌ | 未建参考文献库 |
| D6 | SI 单独文件 + 文末非完整句描述（含扩展名） | 🟡 | 材料可进 SI；描述句未写 |
| D7 | 图形分辨率：彩色 ≥300 dpi；单栏 ≤3.33 in；双栏 4.17–7 in | ❌ | 现有 PNG 元数据 ~**180 dpi**（1890×1476 / 1530×936），定稿前须重导 ≥300 dpi |
| D8 | 图注自洽；表有短标题 | 🟡 | Table 2 已在 Results；图注未定稿 |
| D9 | ACS Research Data Policy Level 1（鼓励公开 + DAS） | ✅ 路径明确 | 与 C 节绑定；执行靠 Zenodo |

**Articles 无硬字数上限。** 字数折算（单栏图≈300 词等）主要用于 Letter / App Note / Perspectives / Reviews，不是 Articles 硬约束。

---

## E. 章节与证据链（开写顺序）

按 P0 建议：**Results 数字 → Methods → Intro → Discussion → Abstract**。

| ID | 章节 | 状态 | 材料指针 |
|----|------|------|----------|
| E1 | Abstract（3–4 句） | ❌ | 待写；须低于 Results 主张强度 |
| E2 | Introduction | 🟡 中英稿已按五段论证重构 | [`INTRODUCTION_DRAFT_ZH_JCIM_V1.md`](INTRODUCTION_DRAFT_ZH_JCIM_V1.md) + [`INTRODUCTION_SECTION_JCIM_EN_V1.md`](INTRODUCTION_SECTION_JCIM_EN_V1.md)；引用核验 [`INTRODUCTION_REFS_JCIM_V1.md`](INTRODUCTION_REFS_JCIM_V1.md)。K=4 不在 Intro 写死；DualDiff/FuseDiff 为互补评测而非竞争对象 |
| E3 | Methods 2.1–2.13 协议化重构 | 🟡 中英稿 | [`METHODS_DRAFT_ZH_JCIM_V1.md`](METHODS_DRAFT_ZH_JCIM_V1.md) + [`METHODS_SECTION_JCIM_EN_V1.md`](METHODS_SECTION_JCIM_EN_V1.md)。结果数字已移出 Methods；four-state + pairwise primary；wrong-pocket = falsification；holdout ≠ external validation；换晶 = structure sensitivity |
| E6 | Results 3.1–3.6 | ✅ 草稿 | 已按证据链压缩为 6 节：[`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) + [`RESULTS_SECTION_JCIM_EN_V1.md`](RESULTS_SECTION_JCIM_EN_V1.md)。旧 3.8 面板重抽样已移出 Results |
| E7 | Discussion / Limitations | ✅ 中英稿 | [`DISCUSSION_DRAFT_ZH_JCIM_V1.md`](DISCUSSION_DRAFT_ZH_JCIM_V1.md) + [`DISCUSSION_SECTION_JCIM_EN_V1.md`](DISCUSSION_SECTION_JCIM_EN_V1.md)；引用 [`DISCUSSION_REFS_JCIM_V1.md`](DISCUSSION_REFS_JCIM_V1.md)；正文局限仅 5 条 |
| E8 | Conclusions | ❌ | 评测主张收束，勿发明新 claim |
| E9 | Data and Software Availability | ❌ | 同 C1/C6 |
| E10 | Keywords ~5–8（Articles 惯例；Perspectives 要求 8–10） | ❌ | 待定 |

**Results 草稿已覆盖：** 3.1 供给约束 → 3.2 有限且对靶的 docking 判别 → 3.3 ligand/chemotype 混淆 → 3.4 配体层持续 vs 受体层崩溃 → 3.5 错口袋 holdout 反转 → 3.6 探索性结构线索。θ 网格、E8/E16、PM110 通道、GNINA mode-1 vs best-of-9、单靶 enrichment 留 SI。

---

## F. 图 / 表定稿

| ID | 显示项 | 状态 | 备注 |
|----|--------|------|------|
| F1 | Fig 1 任务 + 口袋匹配示意 | ❌ | 需新画 |
| F2 | Fig 2 供给审计（49 对） | ❌ | 有 J0 表；无投稿级图 |
| F3 | Fig 3 森林图 pocket-matched ± CI + 基线 | 🟡 | `forest_summary_min_ci_v1.png` 有；须确认标注为 pocket-matched 并升 DPI |
| F4 | Fig 4 混淆（错口袋 / LE / 匹配） | ❌ | 有 CSV；无定稿图 |
| F5 | Fig 5 配体层 vs 受体层 | 🟡 | 脚本已改面板标题；bitmap 需重跑 `plot_jcim_article_figures_v1.py` |
| F6 | Table 1 靶对清单 / N / 协议 | 🟡 | inventory 材料有 |
| F7 | Table 2 主结果（口袋匹配） | ✅ 文内已有 | 与 `PRIMARY_METRIC_V2` 一致 |
| F8 | TOC graphic | ❌ | 同 D4 |

建议主文 **3–5 张图**（与 Perspectives 鼓励量同级即可）；其余进 SI。

---

## G. Claim ceiling 安全核对（写作全程）

| ID | 检查项 | 状态 |
|----|--------|------|
| G1 | 主指标 = pocket-matched directional AUROC | ✅ Results 已用 |
| G2 | PM 仅作 exploratory positive control；Δ vs baseline CI 含 0 须写明 | ✅ Results 3.2–3.3 |
| G3 | 不写 RTM/GNINA 通吃或通用决策臂 | ✅ |
| G4 | EGFR 不作新对接扩面结论 | ✅ |
| G5 | LigPrep 不进投稿正文/SI（无正式权限；主协议仅 RDKit/meeko） | ✅ 已从 Methods/Results/SI 删除 |
| G6 | Abstract/Intro 主张 ≤ Results/Discussion 证据 | 🟡 | 后写章节时再核 |

---

## H. 开写前最小动作清单（按优先级）

### 立刻可开写（不阻塞 Results→Methods）

1. 以英文 Results 为锚，起草 **Methods 2.1–2.3**（含 **Data curation** 标题）  
2. 从 `ENV_PIN.md` 抽出软件版本句进 Methods  
3. 起草 **Limitations** 子弹（讨论用）

### 与写作并行（投稿前必须完成）

4. 打包并发布 **Zenodo** 分数包 → 填 DOI  
5. 重导 Fig3（及后续图）至 **≥300 dpi**；补 Fig1/2/4/5  
6. 写 **Title / Abstract / Cover letter / Data and Software Availability**  
7. 组完整英文稿 + SI 描述句  
8. 全文搜禁语：universal scorer、decision arm validated、显著通吃四对、D-DRAF、novel framework named、Framework Step 1–5

### 可不挡首投、但加分

9. top1 姿态进 Zenodo v1.1  
10. Cover art（修回时再交亦可）

---

## I. 一页 Cover letter 必写点（模板要点）

- Manuscript type: **Article**  
- Contribution: systematic benchmarking protocol + DualFourClass-Bench resource for dual-target docking reliability evaluation  
- Explicitly **not** a new scoring function / named method framework (no D-DRAF-style acronym)  
- Not a wet-lab validated single-target VS campaign  
- Primary metric: pocket-matched directional AUROC with trivial baselines and confound controls  
- Data/code: Zenodo DOI + GitHub （填空）  
- Suggest: fits JCIM molecular modeling / cheminformatics evaluation audience  

---

## 核对结论

| 问题 | 答案 |
|------|------|
| 现在能不能开始写正文？ | **能。** 以 Results 英稿为核写 Methods → Intro → Discussion。 |
| 现在能不能投稿？ | **不能。** 缺完整 IMRaD、Zenodo DOI、DAS 节、cover letter、定稿图/TOC。 |
| 还要不要新对接？ | **不必**（P0 明确）。 |
| 最大风险 | Cover letter / Intro 写成“新方法赢了”；或无 Data curation / 无公开数据被 2015+2021 标准卡住。 |

**下一步建议顺序：** Methods（含 data curation）→ Discussion/Limitations → Abstract → Zenodo → 主图定稿 → Cover letter → 组稿投稿。

---

## 参考链接

- [JCIM Author Guidelines](https://researcher-resources.acs.org/publish/author_guidelines?coden=jcisd8)  
- [2015 Letter from the Editors](https://pubs.acs.org/doi/10.1021/acs.jcim.5b00180)  
- [2021 Method and Data Sharing](https://pubs.acs.org/doi/10.1021/acs.jcim.0c01389)  
- 本仓：[`CLAIM_CEILING.md`](../data/jcim_bench_v0/CLAIM_CEILING.md) · [`JCIM_P0_COMPLETION_GUIDE.md`](JCIM_P0_COMPLETION_GUIDE.md) · [`RESULTS_SECTION_JCIM_EN_V1.md`](RESULTS_SECTION_JCIM_EN_V1.md)
