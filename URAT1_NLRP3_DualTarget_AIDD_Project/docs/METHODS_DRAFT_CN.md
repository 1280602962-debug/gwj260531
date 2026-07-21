# Methods 中文大纲与正文草稿（V2 · 先验协议框架）

> **这是当前中文 Methods 写作底稿**（取代旧稿 `MANUSCRIPT_DRAFT_CN.md` 中的 Glide XP Methods）。  
> 接在引言 [`INTRO_DRAFT_CN.md`](INTRO_DRAFT_CN.md) 之后。  
> **叙事合同**：本节只陈述**事先规定的设计、规则与操作**；协议胜出结果、富集终值、Pareto 名单属 Results，不在此预告。  
> 公式与脚本：`build_urat1_true_decoy.py`、`utils_ml.py`、`merge_docking_pareto.py`、`run_vina_batch.py` / `run_gnina_batch.py`、模块 A–F。  
> Π\* 与 MD 时长待实验回填。

---

## 叙事合同（作者自用，可不进正文）

| 应写 | 不写 |
|------|------|
| 为何两靶证据权重不同（数据不对称） | “我们本来打算用某某，后来发现不行才改” |
| 协议候选集与选优规则**事先锁定** | 按时间线复盘换引擎 |
| 阴性对照与否决规则（如 P0、RandomDecoy 否决） | 把开发失败写成主创新 |
| 富集=排序考试；pose 另报 | 富集=位点已证明 |

---

## 一、Methods 大纲（建议二级标题）

| 节号 | 标题 | 本课题涉及什么 |
|------|------|----------------|
| 2.0 | 总体设计（可并入 2.1 首段） | 三阶段：协议筛选 → 不对称漏斗 → 审计提名 |
| 2.1 | 数据来源与预处理 | ChEMBL URAT1/NLRP3；临床库 8319；蒸馏子集 D |
| 2.2 | TrueDecoy / RandomDecoy 基准 | 性质匹配、近邻过滤、round-robin |
| 2.3 | 受体与配体准备、自对接报告 | 9DKB / 7ALV；三类 RMSD |
| 2.4 | 候选对接与重打分协议 | 预先定义的 P0–P5 |
| 2.5 | 协议选优规则（锁定 Π\*） | EF@f、AUC、四药百分位、否决条件 |
| 2.6 | NLRP3 机器学习缩库 | 骨架 GroupKFold；进入对接池阈值 |
| 2.7 | 双靶百分位与 Pareto | \(S_U\)、\(S_N\)、非支配前沿 |
| 2.8 | 成药性审计与提名 | PAINS/Brenk、适用域、y-scramble、提名规则 |
| 2.9 | 分子动力学（若纳入主文） | 体系、力场、时长、分析项 |
| 2.10 | 软件与可复现性 | 版本、种子、配置路径 |

---

## 二、本课题方法一览

1. **化学信息学预处理**：RDKit 标准化/去重；Murcko；物化描述符。  
2. **双诱饵基准**：TrueDecoy（性质匹配）+ RandomDecoy（等量随机对照）。  
3. **协议筛选（阶段 A）**：预先定义 P0–P5；按锁定规则输出 Π\*。  
4. **自对接报告**：lesinurad@9DKB；Top-1 / Best-in-ensemble / RTMScore 姿 RMSD。  
5. **NLRP3 分类 ML（阶段 B 缩库）**：assay-conditioned；骨架 CV。  
6. **URAT1 回归 ML**：按设计仅作基准回收对照；回收不足则**不进入**主库排序。  
7. **不对称漏斗**：P≥0.5 → Π\* 双靶对接 → 百分位 → Pareto。  
8. **审计提名（阶段 C）**：结构警报、适用域、稳健性；Pareto ≠ 最终提名。  
9. **MD（可选）**：代表复合物定性讨论。

---

## 三、关键公式

### 3.1 活性标度

\[
\mathrm{pActivity} = 9 - \log_{10}(\mathrm{IC}_{50,\mathrm{nM}})
\]

（若 IC50 以 µM 计：\(6-\log_{10}(\mathrm{IC}_{50,\mu\mathrm{M}})\)。活性集默认 \(\mathrm{pActivity}\ge 6\)。）

### 3.2 TrueDecoy 性质距离与近邻过滤

\[
d(a,x)=\sqrt{\frac{1}{|\mathcal{D}|}\sum_{d\in\mathcal{D}}\left(\frac{a_d-x_d}{w_d}\right)^2},\quad
\mathrm{TC}_{\max}(x)=\max_{a\in\mathcal{A}}\ \mathrm{Tanimoto}\big(fp(x),fp(a)\big)\le 0.5
\]

\(\mathcal{D}=\{\mathrm{MW},\log P,\mathrm{TPSA},\mathrm{HBD},\mathrm{HBA},\mathrm{NRot}\}\)，窗口半宽 \(w_d\) 见脚本；不足时可启用 \(1.5w_d\) 放宽轮。

### 3.3 排序分方向统一

对接能越低越好时：\(s=-E_{\mathrm{dock}}\)。CNNaffinity 等“越高越好”指标可直接作 \(s\)。

### 3.4 富集因子

\[
\mathrm{EF}@f=\frac{\big(\tfrac{1}{n_f}\sum_{i\in\mathrm{Top}f} y_i\big)}{\big(\tfrac{1}{N}\sum_{i=1}^{N} y_i\big)}
\]

主文预设报告 \(f=0.01,0.05\)。

### 3.5 百分位与双轴

\[
S=100\times\mathrm{rankpct}(-E),\quad
S_U=S_{\mathrm{URAT1\text{-}dock}},\quad
S_N=\max\big(S_{\mathrm{NLRP3\text{-}ML}},\,S_{\mathrm{NLRP3\text{-}dock}}\big)
\]

（取 max 的设计理由：两证据近正交时避免单轴埋没；相关分析属 Results。）

### 3.6 Pareto 非支配

不存在 \(j\) 使 \(S_U^{(j)}\ge S_U^{(i)}\) 且 \(S_N^{(j)}\ge S_N^{(i)}\)，并至少一维严格更优。

### 3.7 分类 / 回归指标

AUROC、AUPRC；回归侧 RMSE、\(R^2\)、Spearman \(\rho\)（URAT1 仅对照）。

### 3.8 骨架分组交叉验证

Murcko 为组，GroupKFold（\(K=5\)）。

### 3.9 适用域（简化）

\[
\mathrm{AD}=\mathbf{1}\big[\mathrm{TC}_{\max}\ge 0.35\big]
\]

（阈值以脚本为准，敏感性分析进 SI。）

### 3.10 自对接 RMSD

\[
\mathrm{RMSD}=\sqrt{\frac{1}{n}\sum_{k=1}^{n}\|\mathbf{r}_k^{\mathrm{dock}}-\mathbf{r}_k^{\mathrm{xtal}}\|^2}
\]

预先规定报告三类：Top-1、Best-in-ensemble、RTMScore 选姿。严格几何门控：Top-1 ≤ 2 Å；未达标时，富集结论仅解释为**排序协议表现**，结构讨论改用 ensemble 近晶体姿或 RTMScore 选姿。

### 3.11 MM-PBSA（若写主文）

\[
\Delta G_{\mathrm{bind}}\approx \langle E_{\mathrm{MM}}\rangle+\langle G_{\mathrm{PB/GB}}\rangle+\langle G_{\mathrm{SA}}\rangle-T\langle S\rangle
\]

仅作定性/半定量，不换算实验 \(K_i\)。

### 3.12 提名分（示意）

\[
\tau = \alpha\,S_U + \beta\,S_N + \gamma\,\mathbf{1}_{\mathrm{clean}}
\]

权重与过滤器以 `candidate_nomination` 最终设定为准。

---

## 四、Methods 正文草稿

### 2.0 总体设计

本研究将计算流程预先分为三个衔接阶段，并在分析前锁定各阶段规则。**阶段 A（协议筛选）**：在 URAT1（PDB 9DKB）上构建 TrueDecoy 与 RandomDecoy 两套基准，对预先定义的开源对接/重打分协议（P0–P5）计算富集与对照回收，按 §2.5 规则选定生产用排序协议 Π\*。**阶段 B（不对称漏斗）**：对临床库以 NLRP3 分类概率缩库；对进入对接池的分子，在 9DKB 与 7ALV 上按 Π\* 生成双靶结构分，再以百分位定义 \(S_U\)、\(S_N\) 并求 Pareto 非支配前沿。**阶段 C（审计提名）**：对合并结果施加结构警报、类药性、适用域与稳健性检查；Pareto 命中不等于最终提名。URAT1 回归模型仅用于基准回收对照：若回收不满足预设门槛，则不进入阶段 B 的主排序。全文对接分数仅用于池内相对比较。

### 2.1 数据来源与预处理

URAT1 活性数据取自 ChEMBL，经关系符清洗、单位统一与分子级聚合后得到 \(\mathrm{pActivity}\)；默认活性集取 \(\mathrm{pActivity}\ge 6\)。NLRP3 采用多 assay 条件下的分类标签构建训练表。临床重定位库由 ChEMBL 临床阶段化合物整理得到（n=8319）。诱饵池来自蒸馏子集 D（未标记多样性负样本）。所有分子经 RDKit 解析、消毒与规范 SMILES；骨架采用 Murcko 定义。

### 2.2 TrueDecoy 与 RandomDecoy 构建

以 URAT1 活性集为锚，从子集 D 中按物化窗口做性质匹配，并要求 \(\mathrm{TC}_{\max}\le 0.5\)。匹配采用 round-robin，避免少数活性垄断诱饵；不足时用 1.5 倍窗口补齐。得到 TrueDecoy 集后，从同一池中抽取等量随机分子构成 RandomDecoy。两套基准共用同一批活性分子，仅诱饵不同，以便在同一活性标签下对照“性质匹配”与“随机”两种负样本设定。

### 2.3 结构准备与自对接报告

URAT1 采用 inward-open 结构 9DKB；NLRP3 采用 NACHT 结构 7ALV（共晶配体为 MCC950 类类似物 NP3-146）。受体去水、加氢并导出对接格式；配体经同一套 3D 嵌入与准备流程。搜索盒以共晶/参考配体几何为中心（边长约 22 Å，以配置文件为准）。自对接以 lesinurad–9DKB 为主，按 §3.10 报告三类 RMSD，用于区分采样可达性与原生 top-1 排序可靠性。

### 2.4 候选对接与重打分协议

在分析前预先定义下列协议（同一受体、同一搜索盒、同一配体准备）：P1 AutoDock Vina affinity；P2 gnina CNNaffinity；P3 gnina affinity（kcal）；P4 Vina 姿 + RTMScore；P5 gnina 姿 + RTMScore；P0 gnina CNNscore（**预先设定的负对照读出**，不作为生产排序候选）。每个分子保留多构象（num_modes≥9）。富集与临床库排序使用各协议定义的 rank_pose 分数；需要讨论结合几何时，另行报告 RTMScore 最高姿（struct_pose）或晶体坐标，二者角色在设计中分离。主文生产路径限定为上述开源协议集；商业对接软件不纳入默认流程，亦不与开源分数做绝对值混比。

### 2.5 协议选优规则

在 TrueDecoy 与 RandomDecoy 上分别计算 ROC-AUC 与 EF@1%/5%。选优规则在对接完成前锁定如下：（i）以 TrueDecoy 的 EF@1% 为主判据，并列时比较 EF@5% 与 AUC；（ii）若某协议在 RandomDecoy 上相对 TrueDecoy **明显变差**，予以否决；（iii）仍并列时，比较 lesinurad、benzbromarone、verinurad、dotinurad 四药的回收百分位。胜出协议记为 Π\*，固定用于阶段 B 的 URAT1 轴；NLRP3 轴采用同一开源引擎与准备流程，以保持漏斗内分数尺度一致。

### 2.6 NLRP3 机器学习缩库与 URAT1 对照模型

构建 assay-conditioned 分类模型，以 Murcko 骨架 GroupKFold（5 折）评估 AUROC、AUPRC。对临床库输出 \(P(\mathrm{active})\)，预设阈值 0.5 定义对接池。同步训练 URAT1 回归模型并报告骨架交叉验证与尿酸药基准回收；按总体设计，该模型**不参与**临床库主排序，其作用限于论证“URAT1 轴应对接主导”的数据前提。

### 2.7 双靶百分位与 Pareto

对接池分子在 9DKB 与 7ALV 上按 Π\* 对接。按 §3.5–3.6 计算 \(S_U\)、\(S_N\)，求 Pareto 非支配前沿。已知对照药（如 lesinurad、colchicine）同步进入同一百分位池，用于检验漏斗行为是否符合机制预期。

### 2.8 成药性审计与提名

对合并结果标注 PAINS/Brenk 等结构警报，计算类药性与适用域，并进行 y-scramble / 稳健性检查。预先规定：**Pareto 前沿不等于最终提名**；命中警报或成药性差的分子降级；仅在通过过滤器的清洁子集上，按提名规则给出跟进假说。

### 2.9 分子动力学

对代表体系（URAT1 基准药、NLRP3 工具药及审计后提名候选）运行 MD，分析配体/口袋 RMSD、关键接触，必要时计算 MM-PBSA 类评分。初始坐标优先晶体姿或 struct_pose。时长与力场在定稿时据实填写；MD 用于构象讨论，不作为活性证明。

### 2.10 软件与可复现性

RDKit、scikit-learn、AutoDock Vina、gnina、RTMScore 及 MD 软件版本写入正文或 SI；随机种子与配置见 `config/docking_open_source.yaml` 与相应脚本。

---

## 五、写作注意

- Methods = **规则与操作**；Results = **Π\* 是谁、EF 多少、名单如何**。  
- 讨论中的阴性结果（如 top-1 RMSD 未过关、某读出作为负对照失败）写成**框架内的预期检验**，不要写成路线变更史。  
- 明确：富集评估 ≠ 结合位点已证明。
