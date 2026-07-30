# Results（中文工作稿 · 与英文 JCIM 稿对齐）

> 对应英文定稿：[`RESULTS_SECTION_JCIM_EN_V1.md`](RESULTS_SECTION_JCIM_EN_V1.md)  
> 与 [`METHODS_DRAFT_ZH_JCIM_V1.md`](METHODS_DRAFT_ZH_JCIM_V1.md) 的 2.1–2.7 一一对应（见每节末尾标注）。  
> 全部数字均可追溯至 `data/jcim_bench_v0/` 与 `data/jcim_strengthen_t0t1_v0/` 下的表与脚本，未做的分析不写入本节。  
> 投稿以**英文**为准；本稿供中文讨论与校对数字。

---

## 3. Results

### 3.1 公开靶对上严格硬负的供给

双靶对接评测需要四类配体：dual、A_only 硬负、B_only 硬负与 neither。我们在 49 对有 ChEMBL 缓存的靶对上按严格标签规则审计（dual：两端 pChEMBL ≥ 6.5；A_only / B_only：活性端 ≥ 6.5 且对端 ≤ 5.5）。两端严格硬负均 ≥ 50 的只有 4 对；排除金属依赖的 HDAC1/HDAC6 后，剩余 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 三对适合建成规模均衡的严格四类面板（Table 1）。作为对照，文献中常见的 EGFR/HER2 在同一规则下仅有 7 个严格 B_only，达不到该门槛。本文的四对评价集（下文记为 K = 4）由该审计结果确定，构建细节见 Methods 2.1–2.3。

### 3.2 池化分数与口袋匹配方向 AUROC 的差异

双靶分数需要同时压住两条单靶硬负臂。若两臂共用同一池化分数（如两端 Vina 分数的均值），强臂可能掩盖弱臂：在 EGFR/HER2 上，池化 AUROC 接近 0.50，而较弱一臂（dual 对 B_only）单独计算时降到约 0.28。这是本文采用口袋匹配方向 AUROC 作主指标的依据：dual 对 A_only 用口袋 B 的分数，dual 对 B_only 用口袋 A 的分数，summary_min 取两臂较小值（定义见 Methods 2.6）。相对池化，口袋匹配普遍抬高了四对的点估计，但排序未变：仍只有 PIK3CA/mTOR 明显高于随机，其余三对 summary_min ≤ 0.61（Table 2）。

### 3.3 冻结 K = 4 评价集上的方向判别

各面板在同一 RDKit/meeko 配体协议下，用 AutoDock Vina、RTMScore best-of-9 与 GNINA CNN（mode_01 重打分）打分。正文以 Vina 口袋匹配 summary_min 为主报告；RTMScore 与 GNINA 作通道对照。Bootstrap 95% 区间使用配体层重采样，B = 2000，种子 20260729（Table 2）。

**Table 2.** 冻结 K = 4 评价集上的口袋匹配方向 AUROC（Vina）。

| 靶对 | n (dual / A_only / B_only) | dual 对 A_only（口袋 B） | dual 对 B_only（口袋 A） | summary_min [95% CI] | 错口袋 min | 配体效率归一 min | 最优描述符基线 |
|---|---:|---:|---:|---|---:|---:|---|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.281, 0.576] | 0.260 | 0.311 | cLogP 0.482 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.442, 0.737] | 0.444 | 0.413 | TPSA 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.340, 0.648] | 0.349 | 0.332 | 重原子数 0.622 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.457, 0.813] | 0.602 | 0.657 | 重原子数 0.463 |

PIK3CA/mTOR 是唯一 summary_min 点估计同时高于 0.5 与重原子数基线（0.463；Δ ≈ +0.23）的靶对，但其置信区间下界仍接近 0.5。EGFR/HER2 与 PIK3CA/PIK3CB 的 summary_min 分别为 0.430 与 0.500，均低于各自最优描述符基线；AChE/BChE 的对接结果（0.606）低于 TPSA（0.733）。RTMScore 与 GNINA 未改变这一格局：二者在 PIK3CA/mTOR 上弱于 Vina，在其余三对上同样未过描述符基线。

排名读出与弱臂结果一致：EGFR/HER2 按池化 Vina 分数取 Top-10 时，9 个为硬负配体（配体层 bootstrap 均值 ≈ 8.9；95% CI 为 7–10）。

（对应 Methods 2.3、2.5、2.6）

### 3.4 相对描述符基线的门控

对每对靶标，用对接 summary_min 减去表现最好的描述符基线（重原子数、分子量、cLogP 或 TPSA）的对应值得到 Δ。EGFR/HER2 与 PIK3CA/PIK3CB 的 Δ 置信区间整体落在 0 以下；AChE/BChE 点估计未过门，区间刚跨过 0；PIK3CA/mTOR 点估计高于重原子数基线，但 Δ 的 95% 区间仍包含 0。这一差距在口袋匹配前后保持一致：EGFR/HER2 落后于 cLogP，AChE/BChE 落后于 TPSA，PIK3CA/PIK3CB 落后于重原子数。上述结果将本文可支持的主张限定在评测诊断层面，而非通用的双靶决策分数。

（对应 Methods 2.6）

### 3.5 错口袋、配体效率、协变量与匹配子集对照

若对接分数主要跟随配体本身的物理化学属性而非结合口袋，错口袋对照的 summary_min 应偏离正确口袋结果不远、甚至同样偏离 0.5。四对的错口袋 summary_min 分别为 0.260（EGFR/HER2）、0.444（AChE/BChE）、0.349（PIK3CA/PIK3CB）与 0.602（PIK3CA/mTOR）；口袋匹配相对错口袋的差距在四对上均超过 0.09（Table 2）。配体效率归一（分数除以重原子数）后，仅 PIK3CA/mTOR 仍高于重原子数基线（0.657 对 0.463），其余三对在该归一下不再支持方向信号。

在 AChE/BChE 上，dual 配体的平均 TPSA 约为 75，硬负配体（A_only 与 B_only 合并）约为 51；TPSA 单独区分 dual 与硬负的 AUROC 约为 0.769，高于同一对比下的 Vina 分数（约 0.56）。将重原子数与 TPSA 作为协变量纳入逻辑回归后，口袋匹配 dual 对 B_only 的判别 AUROC 从 0.606 升至 0.807（Δ ≈ +0.20）；PIK3CA/mTOR 上相应升幅较小，约 +0.07 至 +0.11。

在效价匹配（|ΔpChEMBL| ≤ 0.5）或尺寸匹配（|Δheavy atoms| ≤ 2）子集上，EGFR/HER2 与 PIK3CA/PIK3CB 的 dual 对 B_only 仍偏弱或接近随机（约 0.45–0.52）；PIK3CA/mTOR 在匹配子集上方向仍保持，但各臂样本量常低于 15、区间较宽，完整分层结果见 Supporting Information Table S5。

（对应 Methods 2.7 第 1–3 条）

### 3.6 支架分组的二维结构基线

配体侧基线为 ECFP4 指纹加逻辑回归，交叉验证按 Murcko 支架分组（GroupKFold），使同一支架不跨训练/测试折。支架折 AUROC 多在 0.78–0.91，普遍高于对应的口袋匹配对接臂：EGFR/HER2 上 dual 对 B_only，指纹模型为 0.85，对接为 0.43。改用随机配体分折后，AUROC 平均仅比支架折高约 0.01（最大约 0.10），与面板中支架数目有限、多数接近单例的情况一致。正文以支架折数字为准，随机折仅作泄漏诊断。该指纹基线表明标签与二维化学型相关，但不能单独证明对接分数具有口袋特异性。

（对应 Methods 2.7 第 5 条）

### 3.7 Exhaustiveness 对照、单靶富集与 PM110 扩面

将 PIK3CA/mTOR 从 exhaustiveness = 16 改为 8（同配体、同盒子、同种子）重新对接后，Vina summary_min 从 0.692 降至 0.660（E16 相对 E8 的 Δ ≈ +0.03）。该差距远小于 PIK3CA/mTOR 与其余三对之间的差距，不足以解释该对相对更优的结果，故两种设置并列报告（对应 Methods 2.4）。

单靶富集分析以同靶已测定的弱效分子（pChEMBL ≤ 5.5）作性质匹配 decoy，而非随机无关分子。4L23（PIK3CA）与 4JT6（mTOR）的富集 AUROC 分别为 0.603 与 0.629，对应 EF1% 为 2.04 与 2.00，EF5% 为 1.22 与 3.20。两端富集能力均高于随机但幅度有限，与“对接本身没有完全失效，但不构成强单靶虚拟筛选”的判断一致（对应 Methods 2.7 单靶参照段）。

PM110 保留 PM48 的全部 48 个配体并按配额扩样，不是独立重复实验。PM110 上 Vina 口袋匹配 summary_min 为 0.648 [0.51, 0.76]，相对 PM48 的 0.692（Δ ≈ −0.04），区间更窄，方向未变；同一面板上 RTMScore 为 0.576，GNINA 为 0.522（对应 Methods 2.3 扩面段）。
