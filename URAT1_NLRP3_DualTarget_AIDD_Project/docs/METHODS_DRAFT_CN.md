# Methods 中文大纲与正文草稿（V2 流程）

> **这是当前中文 Methods 写作底稿**（取代旧稿 `MANUSCRIPT_DRAFT_CN.md` 中的 Glide XP Methods）。  
> 接在引言 [`INTRO_DRAFT_CN.md`](INTRO_DRAFT_CN.md) 之后。公式与步骤以仓库脚本为准：`build_urat1_true_decoy.py`、`utils_ml.py`、`merge_docking_pareto.py`、`run_vina_batch.py` / `run_gnina_batch.py`、模块 A–F。  
> 对接协议 Π\* 终值待服务器跑完回填；MD 时长按最终实验填写。

---

## 一、Methods 大纲（建议二级标题）

| 节号 | 标题 | 本课题涉及什么 |
|------|------|----------------|
| 2.1 | 数据来源与预处理 | ChEMBL URAT1/NLRP3；临床库 8319；蒸馏子集 D；标准化 SMILES |
| 2.2 | URAT1 TrueDecoy / RandomDecoy 基准构建 | 性质匹配、近邻过滤、round-robin 分配 |
| 2.3 | 受体与配体准备、自对接门控 | 9DKB / 7ALV；搜索盒；Top-1 / Best-in-ensemble / RTMScore 姿 RMSD |
| 2.4 | 对接与重打分候选协议 | P0–P5（Vina、gnina、RTMScore）；主文不开源 Glide |
| 2.5 | 富集指标与协议选优规则 | AUC、EF@f、四药百分位；预先锁定选优 |
| 2.6 | NLRP3 机器学习缩库 | 骨架 GroupKFold；AUROC/AUPRC；P≥0.5 |
| 2.7 | 双靶对接百分位与 Pareto 整合 | S_U、S_N、非支配前沿 |
| 2.8 | 成药性审计与提名（模块 A–F） | PAINS/Brenk、ADMET、适用域、y-scramble、提名分 |
| 2.9 | 分子动力学（若纳入主文） | 体系、力场、时长、RMSD/接触/MM-PBSA |
| 2.10 | 软件与可复现性 | 版本、种子、脚本路径 |

---

## 二、本课题方法一览（写给作者自己核对）

1. **化学信息学预处理**：RDKit 标准化/去重；Murcko 骨架；物化描述符（MW、logP、TPSA、HBD、HBA、可旋转键等）。  
2. **TrueDecoy 构建**：DUD-E 风格窗口匹配 + Morgan TC 过滤 + round-robin；RandomDecoy 等量随机对照。  
3. **对接协议筛选**：AutoDock Vina；gnina（CNNaffinity / CNNscore 对照）；RTMScore 对 ensemble 重打分。  
4. **自对接评估**：lesinurad@9DKB；报告三种 RMSD。  
5. **富集评估**：ROC-AUC、EF@1%/5%、临床四药百分位。  
6. **NLRP3 分类 ML**：assay-conditioned 集成分类；骨架分组交叉验证。  
7. **URAT1 回归 ML（仅作阴性对照/不对称论证）**：骨架 CV + benchmark 回收（2/4）→ 不用于主排序。  
8. **临床库漏斗**：P(active)≥0.5 缩库 → 选定协议 Π\* 双靶对接 → 百分位 → Pareto。  
9. **下游审计**：PAINS/Brenk/NIH；类药性；适用域（Max Tc）；y-scramble；Pareto 稳健性；提名规则。  
10. **MD（可选主文）**：代表复合物稳定性与定性自由能。

---

## 三、关键公式（投稿 Methods / SI 可用）

### 3.1 活性标度

由 IC50（nM）得：

\[
\mathrm{pActivity} = 9 - \log_{10}(\mathrm{IC}_{50,\mathrm{nM}})
\]

若单位为 µM：

\[
\mathrm{pActivity} = 6 - \log_{10}(\mathrm{IC}_{50,\mu\mathrm{M}})
\]

URAT1 TrueDecoy 活性集默认取 \(\mathrm{pActivity}\ge 6\)。

### 3.2 TrueDecoy 性质距离（归一化 L2）

对描述符集合 \(\mathcal{D}=\{\mathrm{MW},\log P,\mathrm{TPSA},\mathrm{HBD},\mathrm{HBA},\mathrm{NRot}\}\)，窗口半宽 \(w_d\)（如 MW±40，logP±1，TPSA±25，HBD±1，HBA±2，NRot±2）。活性分子 \(a\) 与候选诱饵 \(x\) 的性质距离：

\[
d(a,x)=\sqrt{\frac{1}{|\mathcal{D}|}\sum_{d\in\mathcal{D}}\left(\frac{a_d-x_d}{w_d}\right)^2}
\]

仅当 \(\forall d:\ |a_d-x_d|\le w_d\)（或放宽轮 \(1.5w_d\)）才允许匹配。近邻过滤：Morgan 指纹（半径 2）与任一活性分子的最大 Tanimoto 相似度

\[
\mathrm{TC}_{\max}(x)=\max_{a\in\mathcal{A}}\ \mathrm{Tanimoto}\big(fp(x),fp(a)\big)\le 0.5
\]

### 3.3 对接分到“越高越好”的排序分

Vina / gnina affinity 等为越低越好（kcal/mol）。富集计算前统一为：

\[
s = -E_{\mathrm{dock}}
\]

gnina CNNaffinity（pK 越高越好）可直接用作 \(s\)，或按协议定义取负亲和力以统一方向。

### 3.4 富集因子与 ROC-AUC

设标签 \(y\in\{0,1\}\)（1=活性），分数 \(s\) 越高越像活性。取前比例 \(f\) 的分子（\(n_f=\lfloor N f\rfloor\)）：

\[
\mathrm{EF}@f=\frac{\big(\tfrac{1}{n_f}\sum_{i\in\mathrm{Top}f} y_i\big)}{\big(\tfrac{1}{N}\sum_{i=1}^{N} y_i\big)}
\]

主文报告 \(f=0.01,0.05\)（即 EF@1%、EF@5%）。ROC-AUC 按标准定义，以 \(s\) 对 \(y\) 计算。

### 3.5 百分位排名（池内相对排序）

对对接能 \(E\)（越低越好），在池 \(\mathcal{P}\) 内：

\[
S = 100\times \mathrm{rankpct}(-E)
\]

即能量越低，百分位越高。对 NLRP3 ML 概率 \(P\)（越高越好）：

\[
S_{\mathrm{ML}}=100\times \mathrm{rankpct}(P)
\]

### 3.6 双轴得分与 Pareto

URAT1 轴（对接主导）：

\[
S_U = S_{\mathrm{URAT1\text{-}dock}}
\]

NLRP3 轴（默认 both）：

\[
S_N=\max\big(S_{\mathrm{NLRP3\text{-}ML}},\ S_{\mathrm{NLRP3\text{-}dock}}\big)
\]

（因 ML 与 7ALV 对接近正交，Spearman 约 −0.04，取 max 避免单轴埋没。）

Pareto 非支配：分子 \(i\) 在前沿，当且仅当不存在 \(j\) 使

\[
S_U^{(j)}\ge S_U^{(i)}\ \wedge\ S_N^{(j)}\ge S_N^{(i)}
\]

且至少一维严格更大。短名单可再加阈值 \(S_U\ge \tau_U,\ S_N\ge \tau_N\)（默认 \(\tau_U=\tau_N=0\)）。

### 3.7 分类与回归评价指标

\[
\mathrm{AUROC},\quad \mathrm{AUPRC}
\]

回归侧（URAT1，仅论证用不主筛）：

\[
\mathrm{RMSE}=\sqrt{\tfrac{1}{n}\sum(y_i-\hat y_i)^2},\quad
R^2=1-\frac{\sum(y_i-\hat y_i)^2}{\sum(y_i-\bar y)^2},\quad
\rho=\mathrm{Spearman}(y,\hat y)
\]

### 3.8 骨架分组交叉验证

Murcko 骨架为组，GroupKFold（\(K=5\)）划分，避免同一骨架泄漏到测试折。

### 3.9 适用域（简化）

相对训练集最大 Morgan Tanimoto：

\[
\mathrm{AD}=\mathbf{1}\big[\mathrm{TC}_{\max}\ge 0.35\big]
\]

（阈值以项目脚本为准，可在 SI 敏感性分析。）

### 3.10 自对接 RMSD

相对晶体配体重原子：

\[
\mathrm{RMSD}=\sqrt{\frac{1}{n}\sum_{k=1}^{n}\|\mathbf{r}_k^{\mathrm{dock}}-\mathbf{r}_k^{\mathrm{xtal}}\|^2}
\]

报告：Top-1 RMSD、Best-in-ensemble RMSD、RTMScore 所选姿 RMSD。严格门控：Top-1 ≤ 2 Å。

### 3.11 MD / MM-PBSA（若写主文）

结合自由能近似：

\[
\Delta G_{\mathrm{bind}}\approx \langle E_{\mathrm{MM}}\rangle+\langle G_{\mathrm{PB/GB}}\rangle+\langle G_{\mathrm{SA}}\rangle-T\langle S\rangle
\]

具体项按所用软件输出填写；主文只作定性/半定量，不宣称实验 \(K_i\)。

### 3.12 提名分（模块 F，示意）

在通过过滤器的集合上，综合结构轴百分位与清洁标签，例如：

\[
\tau = \alpha\,S_U + \beta\,S_N + \gamma\,\mathbf{1}_{\mathrm{clean}}
\]

权重与阈值以 `candidate_nomination` 脚本最终设定为准（开发参考 canagliflozin 居前）。

---

## 四、Methods 正文草稿（中文，可直接改）

### 2.1 数据来源与预处理

URAT1 活性数据取自 ChEMBL，经关系符清洗、单位统一与分子级聚合后得到 \(\mathrm{pActivity}\)；默认活性集取 \(\mathrm{pActivity}\ge 6\)。NLRP3 采用多 assay 条件下的分类标签构建训练表。临床重定位库由 ChEMBL 临床阶段化合物整理得到（n=8319）。诱饵池来自蒸馏子集 D（未标记多样性负样本）。所有分子经 RDKit 解析、消毒与规范 SMILES；骨架采用 Murcko 定义。

### 2.2 TrueDecoy 与 RandomDecoy 构建

以 URAT1 活性集为锚，从子集 D 中按物化窗口做性质匹配，并要求 \(\mathrm{TC}_{\max}\le 0.5\)。匹配采用 round-robin，避免少数活性垄断诱饵；不足时用 1.5 倍窗口补齐。得到 TrueDecoy 集后，从同一池中抽取等量随机分子构成 RandomDecoy。两套基准共用同一批活性分子，仅诱饵不同。

### 2.3 结构准备与自对接

URAT1 采用 inward-open 结构 9DKB；NLRP3 采用 NACHT 结构 7ALV（共晶配体为 MCC950 类类似物 NP3-146）。受体去水、加氢并导出对接格式；配体经同一套 3D 嵌入与准备流程。搜索盒以共晶/参考配体几何为中心（约 22 Å 立方盒，以配置文件为准）。自对接以 lesinurad–9DKB 为主，报告 Top-1、Best-in-ensemble 与 RTMScore 选姿三类 RMSD。

### 2.4 对接与重打分协议

候选协议包括：P1 AutoDock Vina affinity；P2 gnina CNNaffinity；P3 gnina affinity（kcal）；P4 Vina 姿 + RTMScore；P5 gnina 姿 + RTMScore；P0 gnina CNNscore（负对照）。每个分子保留多构象（num_modes≥9）。富集排序使用各协议原生 top-1（rank_pose）；结构分析优先使用 RTMScore 最高姿（struct_pose）或晶体坐标。**主文不以 Schrödinger Glide XP 为默认对接引擎**；历史 Glide 表仅作开发参考或 SI，不与开源分数混比、不写入 Methods 主流程。

### 2.5 协议选优

在 TrueDecoy 与 RandomDecoy 上分别计算 ROC-AUC 与 EF@1%/5%。预先锁定规则：以 TrueDecoy 的 EF@1%（并列 EF@5%、AUC）为主判；RandomDecoy 明显变差则否决；平局比较 lesinurad、benzbromarone、verinurad、dotinurad 的回收百分位。胜出协议记为 Π\*，用于后续临床库 URAT1 轴（NLRP3 轴沿用同一开源引擎或既定准备流程）。

### 2.6 NLRP3 机器学习缩库

构建 assay-conditioned 分类模型，以 Murcko 骨架 GroupKFold（5 折）评估 AUROC、AUPRC。对临床库输出 \(P(\mathrm{active})\)，默认阈值 0.5 进入对接池。URAT1 回归模型仅用于展示基准回收不足，不参与主库排序。

### 2.7 双靶百分位与 Pareto

对接池分子在 9DKB 与 7ALV 上按 Π\* 对接。按 §3.5–3.6 计算 \(S_U\)、\(S_N\)，求 Pareto 非支配前沿。对照药（如 lesinurad、colchicine）同步记录百分位行为。

### 2.8 成药性审计与提名

对合并结果标注 PAINS/Brenk 等结构警报，计算类药性与适用域，并进行 y-scramble / 稳健性检查。Pareto 前沿不等于最终提名：命中警报或成药性差的分子降级；在清洁子集上按提名规则给出跟进假说。

### 2.9 分子动力学

对代表体系（URAT1 基准药、NLRP3 工具药及提名候选）在相应复合物上运行 MD，分析配体/口袋 RMSD、关键接触，必要时计算 MM-PBSA 类评分。初始坐标优先晶体姿或 RTMScore struct_pose。时长与力场在定稿时据实填写。

### 2.10 软件与可复现性

RDKit、scikit-learn、AutoDock Vina、gnina、RTMScore 及 MD 软件版本写入正文或 SI；随机种子与配置见 `config/docking_open_source.yaml` 与相应脚本。分数仅用于池内相对比较，不换算为实验亲和力。

---

## 五、写作注意

- Methods 写**做了什么**，不要提前写 EF 终值、Pareto 6、canagliflozin 等 Results。  
- 公式放正文关键几条即可，其余进 SI。  
- 明确：富集评估 ≠ 结合位点已证明；Top-1 RMSD 未过关时不得写 pose 已验证。
