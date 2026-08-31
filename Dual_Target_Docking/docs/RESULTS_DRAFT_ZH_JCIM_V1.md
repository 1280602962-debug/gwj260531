# Results（中文工作稿 · JCIM Articles）

## 3. 结果

### 3.1 实验数据供给限制了严格双靶基准的构建

为确定公开生物活性数据是否能够支持严格的双靶点识别评测，我们首先对 49 对有 ChEMBL 缓存的候选靶标进行供给审计（Figure 2）。一端达到活性阈值、对端明确低活性的配体定义为方向性选择性硬负样本。

在严格标签规则下（dual：两端 pChEMBL ≥ 6.5；选择性类：活性端 ≥ 6.5 且对端 ≤ 5.5），能够同时提供足量 A-only 与 B-only 硬负样本的靶对十分有限。两端严格硬负均不少于 50 的厚面板条件仅有 4 对满足。排除金属依赖 HDAC1/HDAC6 后，PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 构成三个规模相对充足的靶对；EGFR/HER2 仅有 7 个严格 B-only 配体，因此被保留为供给受限案例（Table 1）。BindingDB 与 PubChem 的零对接计数核对支持同一供给稀缺结论（Table S12）。

严格 6.5/5.5 规则用于量化供给并记录面板构建，而 θ = 6.0 定义全部主 AUROC 的实验状态标签（Methods 2.1）。对同一 49 对在主规则 θ = 6.0 下重计后，有 17 对 dual/A-only/B-only 均 n ≥ 10（Table S44）。对接评价仍是原来的四对。完整病例覆盖为 14.5%–34.0%（Table S37）。186 个优先分子的元数据审核为 179 include / 7 uncertain / 0 exclude，未改变冻结类别（Table S42）。PIK3CA/mTOR Dual versus B-only 按文献阻断后无法稳定估计（Table S40）。在 2026-08-26 API 重拉快照内部，采用 max 与 median 聚合时 EGFR/HER2 的最弱臂 AUROC 分别为 0.417 和 0.424；该比较独立于冻结主分析中 0.430 的 Table 2 估计（Table S29）。

### 3.2 基准设定改变了表观双靶判别

在冻结的四对靶标上，采用统一 θ = 6.0 标签规则和口袋匹配方向 AUROC 对 Vina 对接分数进行评价（Figure 1B；Methods 2.4）。EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 的最弱臂 AUROC 分别为 0.430、0.606、0.500 和 0.692（Table 2；Figure 4A）。四条配体 bootstrap 95% 区间均包含 0.5。算术、几何与调和平均下四对排序不变（Table S26）。

同一套冻结分数再按 Dual versus neither 计分（Table 3；Figure 3）。EGFR/HER2 上 Dual versus neither 的 AUROC 为 0.756 [0.562, 0.920]（n_neg = 12），而方向性最弱臂 AUROC 仍为 0.430 [0.282, 0.578]；Dual versus all non-duals 降至 0.551。在 110 个 EGFR/HER2 配体的混合库中按 `vina_mean` 取 Top-10，含 1 个 dual 与 9 个实验选择性配体（硬负比例 0.90；Table S25）。固定口袋 A 分数后，neither 与 B-only 负类相差 0.378 [0.205, 0.547]（Table S34）。AChE/BChE 与 PIK3CA/PIK3CB 的 Dual-versus-neither 增量很小，区间与方向性臂重叠。PIK3CA/mTOR Dual versus neither 因 neither n = 4 而效能不足。

独立 GNINA 1.3.2 姿态生成仍保留该配方差距：EGFR/HER2 Dual versus neither 0.783 [0.610, 0.922]，方向性最弱臂 AUROC 0.220 [0.109, 0.343]（Table S32）。PIK3CA/mTOR 最弱臂 AUROC 为 0.633。五个预先规定的 Vina 种子上，方向性最弱臂估计保持了类似的靶对特异格局（Table S54）。

**Table 2.** 冻结 K = 4 评价集上的口袋匹配方向 AUROC（Vina，统一 θ = 6.0），并列出四个预先指定描述符的 `summary_min`。表中类别样本量为 n_scored（dual / A-only / B-only）。最高描述符是最佳单一描述符参考。

| 靶对 | n_scored (dual / A-only / B-only) | dual 对 A_only（口袋 B） | dual 对 B_only（口袋 A） | summary_min [95% CI] | heavy | MW | cLogP | TPSA |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|------:|---:|------:|-----:|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.282, 0.578] | 0.369 | 0.416 | 0.482 | 0.427 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.437, 0.730] | 0.582 | 0.579 | 0.467 | 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.350, 0.650] | 0.622 | 0.620 | 0.595 | 0.418 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.470, 0.813] | 0.463 | 0.448 | 0.310 | 0.260 |

**Table 3.** 同一套 Vina 分数在 Dual-versus-neither 与方向性设定下的 AUROC（统一 θ = 6.0）。Dual-versus-neither 使用实验 inactive（`vina_mean`）。PIK3CA/mTOR Dual versus neither 效能不足（n_neg = 4）。

| 靶对 | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|------|--------------------------------:|------------------------------:|----------:|----------------------:|
| EGFR/HER2 | 0.430 [0.282, 0.578] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.437, 0.730] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/PIK3CB | 0.500 [0.350, 0.650] | 0.559 [0.373, 0.746] | 16 | 0.556 [0.437, 0.672] |
| PIK3CA/mTOR | 0.692 [0.470, 0.813] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |

当前样本量更容易分辨较大的方向性效应（Table S31）。CI 未能排除 0.5 并不能建立与随机等价。

### 3.3 配体化学是竞争性解释

口袋匹配对接与四种预先定义的物化性质及支架分组 ECFP4 比较（Figure 4B–C；Tables 2, S19–S20, S24）。相对最佳单一描述符，最弱臂差值在四对上的区间均包含 0。AChE/BChE 上 TPSA 单独即可超过对应方向的 Vina；加入重原子数和 TPSA 后 dual-versus-B-only AUROC 从 0.606 增至 0.807，而对接分数的优势比接近 1。支架分组 ECFP4 在多个方向上高于对接，例如 EGFR/HER2 dual-versus-B-only 约 0.89 对 0.43。把口袋匹配对接分数加到 ECFP4 后，AUROC 最大绝对变化为 0.020（Table S24）。物化 caliper 匹配后，有足够样本的 Dual-versus-B-only 臂仍接近随机（Table S45）。

按文献阻断后，EGFR/HER2 弱臂仍为 0.430（document-cluster 95% CI [0.321, 0.617]；Table S39）。

### 3.4 受体实现改变表观判别的幅度与方向

一端受体冻结、只替换另一端时，表观判别向相反方向变化（Figure 5；Table S30）。在 PIK3CA/mTOR 上，将 PIK3CA 4L23 替换为 4JPS 或 5DXT、mTOR 4JT6 保持不变，最弱臂 AUROC 由 0.692 降至 0.486 [0.259, 0.692] 和 0.505 [0.292, 0.696]。同一套 PIK3CA 晶体用于 PIK3CA/PIK3CB、2WXF 保持冻结时，最弱臂 AUROC 由 0.500 升至 0.691 和 0.685。受体替换因此同时改变了表观判别的幅度与方向。

对接失败集中于大或柔性配体（Table S27）。AChE/BChE 上 rank-extreme lower bounds 与完整病例方向一致；PIK3CA/PIK3CB 上一个失败 A-only 配体改用可用口袋分数后最弱臂 AUROC 仍为 0.500。

### 3.5 BindingDB 外部切片

BindingDB 202608 原生归档按对接前冻结的合约重建后，文献、结构与 ECFP4 < 0.70 过滤没有一对达到预先冻结的主外部门槛；该切片未对接，也不作为外部验证（Tables S48–S49；Figure S8）。[16] 预先冻结的 2018 文献年份分割同样未通过样本量门槛。

### 3.6 双靶筛选的实际后果

EGFR/HER2 上以 Dual 中位 `vina_worst` 做 AND-like dual filter 时，保留 14/28 个 Dual，同时留下 33 个选择性配体（precision 0.298；硬负比例 0.702；Table S46）。四对完整 ChEMBL 图上的配体层 ECFP4 仍比 Dual versus 选择性更容易分开 Dual versus neither（EGFR/HER2：0.921 对 Dual versus B-only 0.864；Table S47）。这些诊断描述双口袋过滤在实验标签化学上的行为。
