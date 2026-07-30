# Results（中文工作稿 · 与英文 JCIM 稿对齐）

> 对应英文定稿：[`RESULTS_SECTION_JCIM_EN_V1.md`](RESULTS_SECTION_JCIM_EN_V1.md)  
> 已按 GitHub nature paper skills / deslop 去 AI 化（见英稿 header）。  
> 投稿以**英文**为准；本稿供中文讨论与校对数字。

---

## 3. Results

### 3.1 Strict hard-negative supply on public target pairs

双靶对接评测需要四类配体：双靶活性、A-only 硬负、B-only 硬负、以及 inactive/neither。我们在 49 对有 ChEMBL 缓存的靶对上，按严格标签规则审计（dual：两端 pChEMBL ≥ 6.5；A_only / B_only：活性端 ≥ 6.5 且对端 ≤ 5.5）。两端严格硬负均 ≥ 50 的只有 4 对。排除金属依赖的 HDAC1/HDAC6 后，适合常规对接厚面板的剩 3 对：PIK3CA/mTOR、AChE/BChE、PIK3CA/PIK3CB（Table 1）。作为对照，文献上常见的 EGFR/HER2 在同一规则下仅有 7 个严格 B_only，达不到厚四类面板门槛。冻结的 K=4 评测集由该审计结果确定（构建规则见 Methods）。

### 3.2 Pooled scores versus pocket-matched directional AUROC

双靶分数必须同时压住两条单靶硬负臂。若两臂共用同一池化分数（如 `vina_mean`），强臂会掩盖弱臂。EGFR/HER2 上池化汇总 AUROC 接近 0.50，而较弱方向臂（dual vs B_only）落到约 0.28（Figure 1）。因此主指标采用口袋匹配方向 AUROC：dual vs A_only 用口袋 B 分数，dual vs B_only 用口袋 A 分数，并以两臂 AUROC 的最小值 `summary_min` 汇总。池化分数仅作对照。相对池化，口袋匹配抬高了四对的点估计，但排序不变：仍只有 PIK3CA/mTOR 明显高于随机，其余三对 `summary_min` ≤ 0.61（Table 2）。

### 3.3 Directional discrimination on the frozen K=4 set

各面板在统一 RDKit/meeko 配体协议下，用 AutoDock Vina、RTMScore best-of-K 与 GNINA CNN（mode_01 rescore）打分。主报告为 Vina 口袋匹配 `summary_min`；RTM 与 GNINA 作通道对照。Bootstrap 95% CI 使用 B = 2000、seed = 20260729（Table 2；Figure 2）。

**Table 2.** Pocket-matched directional AUROC (Vina), K=4.

| Target pair | n (D / A / B) | D vs A (pocket B) | D vs B (pocket A) | summary_min [95% CI] | Wrong-pocket min | LE-PM min | Best trivial baseline |
|---|---:|---:|---:|---|---:|---:|---|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.281, 0.576] | 0.260 | 0.311 | cLogP 0.482 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.442, 0.737] | 0.444 | 0.413 | TPSA 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.340, 0.648] | 0.349 | 0.332 | heavy atoms 0.622 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.457, 0.813] | 0.602 | 0.657 | heavy atoms 0.463 |

PIK3CA/mTOR 是唯一 Vina `summary_min` 点估计同时高于 0.5 与 heavy-atom 基线（0.463；Δ ≈ +0.23）的一对。该 CI 下界仍接近 0.5。EGFR/HER2 与 PIK3CA/PIK3CB 的 `summary_min` 为 0.430 与 0.500，均低于各自最优平凡基线。AChE/BChE 对接（0.606）低于 TPSA（0.733）。RTM 与 GNINA 未扭转这一格局：在 PIK3CA/mTOR 上弱于 Vina，在阴性对上同样未过基线门。

排名读出与弱臂一致。EGFR/HER2 按 Vina 池化分取 Top-10 时，硬负占 9/10（bootstrap 均值 ≈ 8.9；CI ≈ 7–10）。

### 3.4 Baseline gate against trivial descriptors

对每对用对接 `summary_min` 减去最优平凡描述符（heavy atoms / MW / cLogP / TPSA）的 `summary_min` 得到 Δ。在池化分数的 bootstrap 门控下，EGFR/HER2 与 PIK3CA/PIK3CB 的 Δ CI 整体落在 0 以下。AChE/BChE 点估计未过门，Δ 区间刚跨过 0。PIK3CA/mTOR 点估计高于 heavy atoms，但 Δ 的 95% CI 仍包含 0（Figure 3）。升格为口袋匹配后，同类差距仍在：EGFR 落后于 cLogP，AChE 落后于 TPSA，PIK3CB 落后于 heavy atoms。这些门控把可写主张限制在评测诊断，而不是通用双靶决策分。

### 3.5 Wrong-pocket, ligand-efficiency, polarity, and matched-subset controls

若对接分主要跟随配体属性，错口袋对照应偏离 0.5。错口袋 `summary_min` 为 0.260（EGFR/HER2）、0.444（AChE/BChE）、0.349（PIK3CA/PIK3CB）、0.602（PIK3CA/mTOR）。四对口袋匹配相对错口袋的 gap 均 > 0.09（Table 2）。配体效率归一后，仅 PIK3CA/mTOR 仍高于 heavy-atom 基线（0.657 vs 0.463）；其余三对在该归一下不再支持方向信号。

AChE/BChE 上，dual 平均 TPSA ≈ 75，硬负 ≈ 51。TPSA 单独区分 dual 与硬负的 AUROC ≈ 0.769，高于同对比下的 Vina（≈ 0.56）。将 heavy-atom 数与 TPSA 作为协变量后，口袋匹配 dual vs B 的 AUROC 从 0.606 升到 0.807（Δ ≈ +0.20）。PIK3CA/mTOR 上相应升幅较小（约 +0.07 至 +0.11）。

在效价匹配（|ΔpChEMBL| ≤ 0.5）或尺寸匹配（|Δheavy| ≤ 2）子集上，EGFR/HER2 与 PIK3CA/PIK3CB 的 dual vs B 仍弱或近随机（约 0.45–0.52）。PIK3CA/mTOR 在匹配子集上方向仍在，但每臂 n 常 < 15、区间宽，分层结果放 Supporting Information。

### 3.6 Scaffold-grouped ligand fingerprint baseline

配体侧基线采用 ECFP4 逻辑回归，并用 Murcko 支架 GroupKFold，使同一支架不跨训练/测试折。支架折 AUROC 多数在 0.78–0.91，高于对应口袋匹配对接臂。EGFR/HER2 dual vs B_only：指纹模型 0.85，对接 0.43。随机按配体分折平均只比支架折高约 0.01（最大约 0.10），与大量近 singleton 支架、泄漏控制有限一致。正文报告支架折数字；随机折仅作泄漏诊断。指纹基线表明标签与化学型相关，本身不能证明口袋物理特异性。

### 3.7 Exhaustiveness, single-target enrichment, and PM110 expansion

将 PIK3CA/mTOR 从 E = 16 重对接到 E = 8（同配体、盒子、seed）后，Vina `summary_min` 从 0.692 变为 0.660（E16 相对 E8 的 Δ = +0.032）。该差距不足以解释该对相对其他对的优势，故两种设置并列报告。

单靶 enrichment 使用同靶 ChEMBL 弱效已测分子（pChEMBL ≤ 5.5）作 property-matched decoy，而非随机无关物。4L23（PIK3CA）与 4JT6（mTOR）的 enrichment AUROC 分别为 0.603 与 0.629。

PM110 保留 PM48 全部 48 个配体并按配额扩样，不是独立复制。PM110 上 Vina 口袋匹配 `summary_min` = 0.648 [0.51, 0.76]，相对 PM48 的 0.692（Δ ≈ −0.04），区间更窄、方向不变。同面板 RTM 为 0.576，GNINA 为 0.522。PM48 配体准备对照中，池化 Vina `summary_min` 在 RDKit 协议下约 0.671，旧 LigPrep 姿态约 0.597；主表仅用 RDKit。
