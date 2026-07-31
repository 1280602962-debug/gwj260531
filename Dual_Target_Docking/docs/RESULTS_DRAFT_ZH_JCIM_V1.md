# Results（中文工作稿 · JCIM 式结构修订）

> 结构按 JCIM 评测文习惯重排（供给发现 → 标签稳健 → 对接主结果 → 混淆主导 → 稳健性与案例依赖成功 → 结构线索）。  
> 全部数字可追溯至 `data/jcim_bench_v0/` 与 `data/jcim_strengthen_t0t1_v0/`；未做的全库 PLIF / 口袋叠合分析不写入。  
> 投稿以英文为准；本稿供中文审改。错口袋、配体效率、描述符明细见 Supporting Information Table S5–S6。

---

## 3. Results

### 3.1 双靶识别基准的构建：公开数据对硬负配体供给的限制

双靶对接评测需要四类配体：dual、仅 A 端强的选择性配体、仅 B 端强的选择性配体，以及两端均弱的 neither。我们将后两类实验定义的选择性配体作为硬负选择性配体（hard-negative selective ligands），用于检验对接分数能否同时压住两条单靶臂。

在 49 对有 ChEMBL 缓存的靶对上，按严格标签规则（dual：两端 pChEMBL ≥ 6.5；选择性类：活性端 ≥ 6.5 且对端 ≤ 5.5）做供给审计。尽管候选靶对数量不少，**可平衡构建的双靶基准仍受到实验表征硬负配体稀缺的严重约束**：两端严格硬负均 ≥ 50 的只有 4 对。排除金属依赖、不适合作为常规对接主对象的 HDAC1/HDAC6 后，剩余 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 三对适合建成规模较均衡的严格四类面板（Table 1）。文献中常见的 EGFR/HER2 在同一规则下仅有 7 个严格 B 端选择性配体，达不到该门槛，因而作为供给受限案例纳入，而不是严格厚面板。K = 4 评价集由该审计结果确定，而非事后挑选“对接好看”的靶对；构建细节见 Methods 2.1–2.3。

（对应 Methods 2.1–2.3）

### 3.2 统一标签规则下的跨对稳健性

面板建造时，AChE/BChE 与 PIK3CA/PIK3CB 使用严格规则，EGFR/HER2 与 PIK3CA/mTOR（PM48）因严格选择性配体不足而使用 θ = 6.0（Table 1）。为检验结论是否依赖建造阈值，我们在既有面板配体与既有 Vina 分数上，按 θ ∈ {5.5, 6.0, 6.5} 与严格 6.5/5.5 规则统一重标四类，并重算口袋匹配 summary_min（Supporting Information Table S4）。该统一重标作为跨对主稳健分析；建造规则下的点估计仅作 construction readout。

在严格规则下，AChE/BChE 与 PIK3CA/PIK3CB 的 summary_min 仍为 0.606 与 0.500；PIK3CA/mTOR 为 0.639（相对建造规则 θ = 6.0 的 0.692 略降）；EGFR/HER2 为 0.324（相对 0.430 下降）。EGFR/HER2 与 PIK3CA/mTOR 在严格规则下分别仅有 7 与 4 个 B 端选择性配体，标记为 underpowered。四对排序趋势保持一致：PIK3CA/mTOR 最高，其余三对不超过 0.61。因此，后文对接比较的排序结论不随统一重标而翻转。

（对应 Methods 2.1）

### 3.3 对接对真双靶配体的判别能力有限且高度依赖靶对

双靶分数需要同时压住两条单靶选择性臂。若两臂共用同一池化分数（如两端 Vina 分数的均值），较强一臂可能掩盖较弱一臂的失败。本文因此以口袋匹配方向 AUROC 为主指标：dual 对 A_only 用口袋 B 的分数，dual 对 B_only 用口袋 A 的分数，summary_min 取两臂较小值（Methods 2.6）。分数定义为 \(S=-E_{\mathrm{Vina}}\)（越大越好），dual 为正类。

需要分开两件事：其一，**方向特异的判别失败**——例如在 EGFR/HER2 上，dual 对 B_only 的口袋匹配 AUROC 仅为 0.430，而若单独审视池化协议下较弱一臂可读到约 0.28，说明 B 端方向本身接近或低于随机，并非“池化运算把 0.50 压成 0.28”；其二，**池化会掩盖上述弱臂**——同一对上池化 summary 可接近 0.50，给人以“尚可”的假象。相对池化，口袋匹配普遍抬高了四对的点估计，但排序未变（Supporting Information Table S6）。

各面板在同一 RDKit/meeko 协议下用 AutoDock Vina 打分；RTMScore（九姿态取最优）与 GNINA（仅对 Vina mode 1 rescore）作通道对照，二者姿态覆盖不对称，不作公平赛马解释。Bootstrap 95% 区间为配体层重采样（B = 2000，种子 20260729）。

**Table 2.** 冻结 K = 4 评价集上的口袋匹配方向 AUROC（Vina）。错口袋、配体效率归一与描述符基线见 Supporting Information Table S6。

| 靶对 | n (dual / A_only / B_only) | dual 对 A_only（口袋 B） | dual 对 B_only（口袋 A） | summary_min [95% CI] |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.281, 0.576] |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.442, 0.737] |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.340, 0.648] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.457, 0.813] |

对接信号总体偏弱，且高度依赖靶对。PIK3CA/mTOR 是唯一 summary_min 点估计同时高于 0.5 与最优平凡描述符基线（重原子数 0.463）的靶对，但其 95% CI 下界仍接近 0.5，不足以支持强判别主张。EGFR/HER2、PIK3CA/PIK3CB 的 summary_min 分别为 0.430 与 0.500，均不超过各自最优描述符；AChE/BChE 为 0.606，仍低于 TPSA 基线（0.733）。RTMScore 与 GNINA 未改变这一格局。排名读出与弱臂一致：EGFR/HER2 按池化 Vina 取 Top-10 时，9 个为硬负选择性配体（bootstrap 均值 ≈ 8.9；95% CI 为 7–10）。

（对应 Methods 2.5–2.6）

### 3.4 物理化学与结构混淆在多数靶对上主导表观信号

核心问题不是简单的“对接不好”，而是：**许多表观双靶信号可由配体属性与二维化学型解释。**

相对最优平凡描述符（重原子数、分子量、cLogP 或 TPSA），对接 summary_min 的差值 Δ 在 EGFR/HER2 与 PIK3CA/PIK3CB 上置信区间整体落在 0 以下；AChE/BChE 点估计未过门；PIK3CA/mTOR 点估计高于重原子数基线，但 Δ 的 95% CI 仍包含 0（Supporting Information Table S6）。错口袋对照的 summary_min 分别为 0.260、0.444、0.349 与 0.602；口袋匹配相对错口袋的差距均超过 0.09，提示分子层混淆仍显著。配体效率归一后，仅 PIK3CA/mTOR 仍高于重原子数基线（0.657 对 0.463）。

在 AChE/BChE 上，dual 与硬负选择性配体的平均 TPSA 分别约为 75 与 51；TPSA 单独区分 dual 与硬负的 AUROC 约为 0.769，高于同一对比下的 Vina（约 0.56）。将重原子数与 TPSA 纳入逻辑回归后，dual 对 B_only 的判别 AUROC 从 0.606 升至 0.807。该升幅表明：**表观对接贡献在很大程度上依赖于物理化学协变量**，而不是证明对接分数已提供强的口袋特异信息；对接分数的优势比接近 1（OR ≈ 1.18），不宜解读为“保留独立方向信息”。相比之下，PIK3CA/mTOR 在控制尺寸与极性后 AUROC 升幅较小（约 +0.07 至 +0.11），OR 约为 2.19 与 3.08，仅提示存在有限的残余口袋相关信号（residual pocket-specific signal），且须与仍跨 0 的 Δ 区间一并阅读。

二维结构基线进一步强化混淆叙事。ECFP4 加逻辑回归、按 Murcko 支架 GroupKFold 时，支架折 AUROC 多在 0.78–0.91，普遍高于对应口袋匹配对接臂；EGFR/HER2 上 dual 对 B_only，指纹为 0.85，对接为 0.43。随机配体分折仅平均高出约 0.01，与面板支架接近单例、泄漏控制有限一致。标签与化学型相关，不能单独证明对接分数具有口袋物理特异性。

效价匹配或尺寸匹配子集上，EGFR/HER2 与 PIK3CA/PIK3CB 的 dual 对 B_only 仍偏弱或接近随机（约 0.45–0.52）；PIK3CA/mTOR 的排序趋势保持一致，但各臂 n 常低于 15、区间较宽（Table S5）。

（对应 Methods 2.6–2.7）

### 3.5 稳健性核对与案例依赖的成功

将 PIK3CA/mTOR 的 exhaustiveness 从 16 改为 8（同配体、同盒子、同种子）后，Vina summary_min 从 0.692 降至 0.660（Δ ≈ +0.03）。该差距远小于靶对间差异，不足以单独解释 PIK3CA/mTOR 相对更优。

单靶富集以同靶已测定弱效分子（pChEMBL ≤ 5.5）作性质匹配 decoy。4L23 与 4JT6 的富集 AUROC 分别为 0.603 与 0.629，EF1% 为 2.04 与 2.00，EF5% 为 1.22 与 3.20。对接保留有限富集能力（limited enrichment capability），但不构成强单靶虚拟筛选引擎。

PM110 保留 PM48 全部 48 个配体并按配额扩样，用作**稳定性核对（stability check）**，不是独立验证实验，也不是用扩面“挽救”点估计。PM48 本身 dual / A_only / B_only 仅 18 / 14 / 12，功效有限。PM110 上 Vina summary_min 为 0.648 [0.51, 0.76]，相对 PM48 的 0.692（Δ ≈ −0.04），区间更窄，排序趋势保持一致；同面板 RTMScore 为 0.576，GNINA 为 0.522。

综合 §3.2–3.5：**仅 PIK3CA/mTOR 表现出可重复但幅度有限的口袋相关判别**；其余三对的表观信号在很大程度上可由配体属性或二维化学型解释。

（对应 Methods 2.3–2.4、2.7）

### 3.6 PIK3CA/mTOR 可重复信号的结构线索（个案级）

为避免将唯一较优靶对写成黑箱，我们基于已导出的姿态级诊断（failure typology；非全面板 PLIF 工程）归纳结构线索。PIK3CA（4L23）与 mTOR（4JT6）均为 ATP 竞争型激酶口袋；共晶配体 PI-103 / X6K 在协议检查中可回收近晶姿态（4JT6 需 elevating exhaustiveness 与 best-of-9，见 Table S3）。在代表性真 dual（如 PI-103）上，两端均可观察到 hinge 接触与对共晶配体的高占用；而若干 PIK3CA 选择性硬负（如 amino-triazine / morpholine–ATP 化学型）在**弱端 mTOR 上同样给出几何干净、hinge 阳性的姿态**，使两端分数同时偏高——这与“化学型同源假双靶”（T2）失败型一致，也说明即便在相对最好的靶对上，对接仍可能把 ATP 位点交叉化学型误判为 dual。

相反，部分经典 dual（如 Torin1、Omipalisib）在 Vina 上很强，但 RTM 优选姿态可偏离 PIK3CA hinge/共晶位，表现为重打分误伤（T5）。因此，PIK3CA/mTOR 的 summary_min 优势应解读为：**在共享 ATP 识别框架下、对部分化学型可出现有限方向信号**，而不是已验证的通用双靶决策规则；全面板残基级守恒分析与相互作用指纹定量比较超出当前冻结分析包，留待后续工作。

（对应已有 failure typology 与 cognate QC；非新对接）
