# JCIM 投稿前全项目审计与行动计划（2026-08-27）

> 本文件取代此前分散的“下一步计划”作为当前投稿决策入口。历史计划用于追溯，不再代表待办优先级。

## 一、结论先行

当前课题已经是一项可复现、诚实且有明确问题意识的 **四靶对 formulation / failure-mode audit**，但还不是可以直接投稿 JCIM 的完整 Research Article。最主要的障碍不是“没有湿实验”这一点本身，而是：

1. 证据宽度仍是 K = 4，且两个靶对共享 PIK3CA、三个靶对涉及激酶 ATP 位点；
2. 最有辨识度的 formulation gap 主要来自 EGFR/HER2，另外三对没有形成一致规律；
3. 四个主 `summary_min` 的 95% CI 都包含 0.5；
4. 预冻结的时间外和 BindingDB-native 路线均未形成可对接的双靶外部验证集；
5. MCL1/Bcl-xL LC6 QC 的 RMSD 实现不是拓扑约束的 symmetry-corrected RMSD，且没有完成预声明的 physical-validity / PLIF / 多种子门槛，不能称为正式 pose-gold gate；
6. 英文主稿约 11,062 个 word-like tokens，摘要 331 词、15 句，且 SI 仍是中文工作稿并扩张到 53 张表，投稿叙事明显过载；
7. 数据仍指向移动中的 GitHub 分支，尚无冻结 release 和 Zenodo DOI。

如果现在投稿，我的编辑判断是：**作为“一般双靶 docking benchmark”不成立；作为严格限界的范式批判/最佳实践 Article，处于 borderline，仍有较高 desk-reject 或大修后拒稿风险。** 不建议靠继续堆 docking engine、短 MD 或更多事后阈值来“做大”。最佳策略是先修复 QC 硬伤、压缩叙事、冻结可复现制品，再决定是否把 MCL1 完整保留在 SI。

## 二、已读取和核验的项目状态

- 当前分支：`cursor/pik3ca-mtor-structure-freeze-0b1a`，审计基线提交 `3b36df0`。
- 主稿：`docs/MANUSCRIPT_JCIM_EN.md`。
- 主评价：EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB、PIK3CA/mTOR 四对；MCL1/Bcl-xL 不进入 Table 2。
- BindingDB-native 202608：预冻结规则后没有靶对通过 primary gate，也没有 EGFR thin replication；停止对接是合规的阴性 stop-rule 结果，不是未完成的外部验证。
- MCL1/Bcl-xL：24/24/24/24 冻结面板，93/96 配体双端完成；结果只能作为 applicability stress-test。
- 主结论仍由 EGFR/HER2 驱动：Vina Dual-vs-neither 0.756，而 directional `summary_min` 0.430；独立 GNINA pose generation 保留该描述性差距。其他靶对没有同样模式。
- ECFP4 加 docking 的 scaffold-grouped CV AUROC 绝对变化不超过 0.020；受体替换可提高或降低结果，支持“receptor realization matters”，不支持单向 collapse。
- 现有 claim ceiling 对外部验证、K = 4、置信区间、complete-case selection 和 assay heterogeneity 的限制总体诚实。

## 三、哪些是合理的、必须做的、错误的、画蛇添足的

### A. 已经合理，应保留

| 项目 | 裁决 | 原因 |
|---|---|---|
| dual / A-only / B-only / neither 四状态定义 | 保留为核心 | 这是文章真正的新意与审稿价值，不是新打分函数 |
| 两条 pocket-matched directional AUROC + 两臂完整报告 | 保留为主终点 | 直接对应“双端同时成立”的证据要求；`summary_min` 只能作描述性摘要 |
| ECFP4、四个描述符、caliper、wrong-pocket、文献阻断和受体替换 | 保留但压缩 | 它们共同排除最明显的化学、来源和受体实现捷径 |
| 失败配体和 rank-extreme bounds | 保留 | 避免只在 protocol-processable subset 上无说明地报乐观结果 |
| BindingDB-native 预冻结 stop rule | 保留为供给审计 | 阴性结果有价值；必须继续明确“未对接、非 external validation” |
| 全部 CI 包含 0.5 的如实陈述 | 必须保留 | 这是当前证据强度的客观边界 |

### B. 投稿前必须完成（P0）

1. **修复 MCL1 LC6 QC。** 用保持分子图同构/对称性的 RMSD（例如 RDKit `GetBestRMS` 或 graph-isomorphism RMSD）重算；检查晶体 LC6 与输出 PDBQT 的化学图一致性。按原预声明补 physical-validity、关键相互作用/PLIF 恢复和至少第二随机种子。若无法完整执行，则删除“pose-gold pass/fail”的正式表述，只保留“preliminary coordinate screen was inadequate; pair is nonconfirmatory”。
2. **把 SI 改成英文并做减法。** 53 张表不适合直接送审。建议保留约 15–22 个直接支撑主张的表/图；运行日志、候选受体尝试、完整流水表、历史 REST 供给核对放入 repository archive，不占 SI 主叙事。
3. **重写摘要。** 官方要求 3–4 句且简洁。目标 180–230 词、4 句：问题；设计；最关键的 EGFR 与跨对照结果；严格结论。BindingDB 失败供给、AND 工作点、17 对普查和 MCL1 不应同时挤进摘要。
4. **压缩主稿并建立一条主线。** 建议正文从约 11k 降到约 7–8k（不含参考文献），Results 由 9 节合并成 5–6 节。主线只回答：负类定义是否改变结论；这种变化能否被化学/来源/受体解释；结论边界在哪里。
5. **冻结可复现版本。** 清洁环境运行全部 manuscript-facing 回归检查，生成 checksum；创建只读 Git tag/GitHub Release；把代码、表、必要受体/pose 或明确的可再生下载清单存入 Zenodo 并获得 DOI。正文 Data Availability 必须引用永久版本，而不是 moving branch。
6. **补齐投稿文件。** 英文 SI PDF、TOC graphic、cover letter、Supporting Information 内容说明、Data Availability、Author Contributions、Funding、Conflict of Interest，以及 AI 使用披露（ACS 要求说明文本/图像/文献组织中的 AI 使用）。
7. **做最终一致性审计。** 主文、SI、MASTER_RESULTS_TABLE、claim ceiling、图注中的 n、CI、PDB、阈值、版本必须逐项一致；禁止将 full-map ECFP、BindingDB supply freeze 或 MCL1 stress-test写成额外 docking validation。
8. **使 checksum 跨平台稳定。** 当前 manifest 直接哈希工作树字节，在 Windows `core.autocrlf` 下会让多数已提交 CSV/Markdown 与原 manifest 不一致。发布前应在脚本中对文本统一 LF 后再计算 SHA-256（或对 Git blob 计算），随后在同一干净环境重建整份 manifest；不能只手改个别 hash 宣称全库通过。

### C. 建议完成但不是决定性门槛（P1）

| 项目 | 建议 |
|---|---|
| 161 个尚未解析的 ChEMBL development documents | 为审计完整性补齐最好；但剩余 n 是上界且所有靶对已失败，补齐只会减少样本，不能把失败变成通过，因此不是继续 docking 的前提 |
| 主 AUROC 的 document/scaffold cluster uncertainty | 若当前主表仍以 ligand bootstrap 为主，建议并列提供 cluster bootstrap/leave-document-out 可估计结果；不可估计的单元格明确标记，不做插补 |
| 文献更新 | 增加 data leakage / assay-derived benchmarking / property-unmatched decoy 的直接讨论；不要把文献列表扩成泛泛综述 |
| 报告清单 | 增加一个面向读者的最小复现实验表：输入快照、预处理、受体/盒、随机种子、软件版本、失败规则、主指标、冻结文件 hash |

### D. 当前错误或必须禁止的做法

1. 把按元素全局 Hungarian 坐标匹配称为 symmetry-corrected RMSD 或 topology-checked RMSD；它允许同元素原子跨拓扑任意配对，可能低估 RMSD。
2. 把 BindingDB-native 剩余候选称为外部验证集，或在看到计数后放松 θ、Tanimoto、来源数和类别 n 门槛。
3. 把 MCL1/Bcl-xL 称为第五个主靶对、异质折叠扩展、正式 screening evidence 或 pose-gold passed。
4. 用 `summary_min` 替代两条方向 AUROC，或把它解释为生物双靶活性/新 scoring function。
5. 把同一 ChEMBL 完整病例图、unused pool、PubChem/BindingDB 镜像记录称为真正数据库外或前瞻验证。
6. 用没有匹配实验任务的短 MD/MM-GBSA 作为“生物验证”。它最多是姿态稳定性敏感性，不能验证真实活性。

### E. 画蛇添足，原则上不要做

- 再加多个 docking/rescoring engine 做无终点的 bake-off；已有 Vina、RTMScore、GNINA 足够说明 engine dependence 没有消除核心不确定性。
- 为了让 MCL1 过门槛，在看到 3WIZ 结果后更换盒子、受体或参数；这会形成 post-hoc receptor shopping。
- 继续增加相似激酶对；这会增加 n，却不能解决靶域覆盖和真正独立验证。
- 将 17 对 θ = 6.0 标签普查全部对接。供给存在不等于每对都有合格受体、来源独立性和可比的 hard negatives。
- 对每个诊断都设立正文小节和 SI 表。项目现在的问题已经从“控制不够”转为“信息架构过载”。
- 把失败的 BindingDB 和 MCL1 结果同时放进摘要。这些适合 Results/Limitations 或 SI，用于证明研究边界，不是摘要主发现。

## 四、推荐的执行顺序与停止规则

### Gate 1：科学完整性（先做）

- 重算 MCL1 topology-aware RMSD，并补足预声明 QC；否则将其降级为 repository-only exploratory stress-test。
- 检查主 AUROC 的 cluster uncertainty 与所有失败敏感性是否可由已提交表重建。
- **停止规则：** 若 MCL1 仍失败，不再换受体“救结果”；若 document-blocked 单元不可估计，报告不可估计，不改变分组。

### Gate 2：稿件收敛

- 先写 4 句摘要和一段 120–180 词的 central claim，再据此删正文。
- 主文只保留 4–6 个主要 display items；SI 英文化并合并重复表。
- MCL1、BindingDB、17-pair census 各自最多承担一个作用：边界/供给/可扩展性，不能被写成三项新增验证。

### Gate 3：可复现发布

- 在干净环境重跑零 docking 分析和全部验证脚本；核对生成文件与 committed checksums。
- tag → GitHub Release → Zenodo DOI → 回填稿件固定链接。
- **停止规则：** moving branch、中文 SI、回归测试未通过或 DOI 未生成时，不提交 JCIM。

### Gate 4：投稿决策

- 若 P0 全部完成且文章被压缩为“well-grounded provocative evaluation/best-practice paper”，可以尝试 JCIM。
- 若不愿进一步修复 MCL1/英文 SI/永久归档，或仍想强调“benchmark performance”，应转投 Journal of Computer-Aided Molecular Design 或 Molecular Informatics，而不是靠扩大措辞冲 JCIM。

## 五、按 JCIM 要求核实稿件水平

JCIM 最新官方指南（更新于 2026-07-03）说明：期刊发表 chemical informatics/molecular modeling 新方法及“有实验验证的应用”；同时，编辑政策明确把“充分扎实、批判/否定既有范式并提出解决思路的论文”和 best-practice 方法学列为更可能成功的类型。Article 摘要应为简洁的 3–4 句；SI 必须随稿单独提交；JCIM 执行 ACS Research Data Policy Level 2，投稿时必须有 Data Availability Statement。ACS 也要求披露 AI 工具使用。参见 [JCIM Author Guidelines](https://researcher-resources.acs.org/publish/author_guidelines?coden=jcisd8)。

这意味着：**没有新湿实验并非形式上绝对禁止，但本稿不能走“应用发现”路线，只能走严格的 computational evaluation / paradigm critique 路线。** 现有 ChEMBL 实验标签是 retrospective experimental grounding，不等于 prospective experimental validation。

### 审稿维度评分（当前，不是接受概率）

| 维度 | 当前水平 | 判断 |
|---|---:|---|
| 问题重要性 | 8/10 | 双靶筛选中 selective hard negatives 的缺失确实重要 |
| 概念新意 | 6.5/10 | formulation audit 有价值，但四状态与 hard-negative 思想不是全新算法 |
| 内部严谨性 | 7/10 | 控制很多且保留阴性结果；MCL1 QC 硬伤必须修 |
| 外部有效性 | 3/10 | K = 4、无可评估 time/database-external docking |
| 统计说服力 | 5/10 | 所有主 `summary_min` CI 包含 0.5，主要强叙事由一对驱动 |
| 可复现性 | 8/10（冻结前） | 脚本/表/claim ceiling 很强；缺永久 DOI 和干净环境最终复现 |
| 写作与投稿完成度 | 4/10 | 过长、摘要不合要求、SI 非英文且过度膨胀 |
| JCIM 当前就绪度 | **约 55/100** | 不建议现在提交；完成 P0 后才进入可尝试区间 |

### 我作为审稿人会提出的核心问题

1. 为什么四靶对、且主要现象只在 EGFR/HER2 明显，足以支持超出个案的评价建议？
2. 当全部主置信区间包含 0.5 时，作者的 confirmatory claim 究竟是什么？哪些分析是预设，哪些是看到结果后增加？
3. complete-case 仅覆盖有任一端可用测量结构的 14.5%–34.0%，文献/系列集中和 assay heterogeneity 会怎样改变四状态标签？
4. 为什么不能提供至少两个真正 time/database-external、文献和结构独立的 docked pairs？BindingDB gate 失败后，文章应如何收缩结论？
5. ligand-only ECFP4 已解释大量类别结构，docking 的独立结构信息在哪里？增量 ≤0.020 是否意味着主要发现是数据集构造而非 docking？
6. 受体选择、质子化/互变异构、单构象准备和单随机种子对结论有多大影响？
7. cognate redocking 是 search coverage 还是 top-ranked pose recovery？为什么部分受体 top-1 失败？
8. MCL1 的 RMSD 原子映射是否保持化学图和对称性？physical plausibility 与关键相互作用是否通过？
9. `summary_min` 作为两个噪声估计的最小值存在下偏，为什么不始终把两臂作为主结果？
10. 53 张 SI 表中哪些是验证核心主张所必需，哪些只是项目日志？读者能否从固定 DOI 一键重建主表？

## 六、投稿前最终判定标准

只有同时满足以下条件，才建议点击 JCIM submission：

- MCL1 QC 修复或从正式证据链中移除；
- 英文 SI 完成并显著合并；
- 摘要压缩为 3–4 句，全文主线收敛；
- 所有主结论均能在 MASTER_RESULTS_TABLE 和固定 checksum 中定位；
- GitHub Release + Zenodo DOI 完成；
- cover letter 明确承认 K = 4、全部主 CI 包含 0.5、BindingDB 未通过 supply gate，并把贡献定位为 evaluation standard / failure-mode audit；
- 不再新增事后“救结果”的受体、阈值或引擎。

达到这些条件后，文章会从“数据很多但未收敛的内部审计”提升为“可被严肃评审的、范围受限的计算评价论文”。它仍不是 JCIM 稳收稿，但会有清晰、诚实且可辩护的投稿理由。
