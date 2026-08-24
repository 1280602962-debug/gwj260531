# Results（中文工作稿 · JCIM 式结构修订）

> 结构按 JCIM 评测文习惯重排（供给发现 → 标签稳健 → 对接主结果 → 混淆主导 → 稳健性与案例依赖成功 → 跨对结构决定因素 → 个案结构线索）。  
> 全部数字可追溯至 `data/jcim_bench_v0/`、`data/jcim_strengthen_t0t1_v0/`、`data/jcim_bench_v0/analysis/structural_context_v1/`、`data/jcim_holdout_v0/`、`data/jcim_structure_robust_v0/analysis/pocket_mechanism_v1/` 与 `data/jcim_supply_crossdb_v0/`；未做的全面板残基级 PLIF 定量比较不写入。  
> 错口袋、配体效率、描述符明细见 Supporting Information Table S5–S6；holdout 错口袋匹配诊断见 Table S13；GNINA 全 9 姿态公平重打见 Table S14–S15。
> **本文定位（不用绝对化标题、不包装成新算法）：** 不是 "Docking can/cannot identify dual-target ligands"，也不是 "we developed a novel framework named D-DRAF"；而是 *Evaluating the reliability and limitations of docking-based dual-target recognition*——建立 systematic benchmarking framework / DualFourClass-Bench 评价体系，评价现有对接分数的可靠边界。详见 [`POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`](POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md)。

---

## 3. Results

### 3.1 双靶识别基准的构建：公开数据对硬负配体供给的限制

双靶对接评测需要四类配体：dual、仅 A 端强的选择性配体、仅 B 端强的选择性配体，以及两端均弱的 neither。我们将后两类实验定义的选择性配体作为硬负选择性配体（hard-negative selective ligands），用于检验对接分数能否同时压住两条单靶臂。

在 49 对有 ChEMBL 缓存的靶对上，按严格标签规则（dual：两端 pChEMBL ≥ 6.5；选择性类：活性端 ≥ 6.5 且对端 ≤ 5.5）做供给审计。尽管候选靶对数量不少，**可平衡构建的双靶基准仍受到实验表征硬负配体稀缺的严重约束**：两端严格硬负均 ≥ 50 的只有 4 对。排除金属依赖、不适合作为常规对接主对象的 HDAC1/HDAC6 后，剩余 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 三对适合建成规模较均衡的严格四类面板（Table 1）。文献中常见的 EGFR/HER2 在同一规则下仅有 7 个严格 B 端选择性配体，达不到该门槛，因而作为供给受限案例纳入，而不是严格厚面板。K = 4 评价集由该审计结果确定，而非事后挑选“对接好看”的靶对；构建细节见 Methods 2.1–2.3。

对该四对的 BindingDB / PubChem 计数核对（Supporting Information Table S12；零对接）表明：**厚面板门槛在等式测定下不随数据源翻转。** 与 pChEMBL 更接近的 `equal_only` 规则下，三对已冻结厚面板的 min HN 分别为 BindingDB 76 / 92 / 58、PubChem 86 / 97 / 61（ChEMBL 缓存为 80 / 78 / 56），仍全部 ≥ 50。EGFR/HER2 的 min HN 从 ChEMBL 的 7 升至 BindingDB 31、PubChem 30，达到薄面板（≥ 20）但仍达不到 ≥ 50 厚面板门槛。若把 `>` 截尾（如 IC50 > 10 µM）也当作弱端点估计（`as_is`），EGFR/HER2 的 BindingDB min HN 会升至 85（85/92）从而过 ≥ 50；但 92 个 as-is B_only 中有 49 个在 EGFR 端只有 `>` 记录。这改变的是“两端都有等式定量测定”的定义，不能写成 ChEMBL 漏掉了约 80 个干净的 HER2 选择性配体。PubChem 与 BindingDB 计数接近，符合沉积重叠，不作两次独立普查。K = 4 对接面板仍按 ChEMBL pChEMBL 冻结；本核对未重建面板、未做新对接。

（对应 Methods 2.1–2.3）

### 3.2 主分析采用单一统一标签规则（θ = 6.0），阈值敏感性作为支持性分析

为避免“不同靶对用不同阈值”的质疑，Table 2 的四对结果统一在**同一 θ = 6.0 规则**下报告（两端 ≥ 6.0 为 dual，一端 ≥ 6.0 且对端 < 6.0 为对应单靶类）。对 EGFR/HER2 与 PIK3CA/mTOR，该规则与建造时直接采用的规则相同；对 AChE/BChE 与 PIK3CA/PIK3CB，建造时因供给审计用的是更严格的 6.5/5.5 规则，但在本数据上 θ = 6.0 给出**完全相同**的配体分类与 AUROC（Supporting Information Table S4），即这两对的标签在所测阈值范围内不敏感——这是先核实、后报告的经验事实，不是选择性展示。

作为支持性稳健性分析（而非与 Table 2 竞争的第二套主标准），我们进一步在 θ ∈ {5.5, 6.5} 与严格 6.5/5.5 规则下重标四类并重算 summary_min（Table S4）。EGFR/HER2 与 PIK3CA/mTOR 对阈值更敏感：严格规则下二者分别降至 0.324（仅 7 个 B 端选择性配体，underpowered）与 0.639（仅 4 个，underpowered），相对 θ = 6.0 的 0.430 与 0.692 均下降但排序未变。四对排序趋势在整张阈值网格内保持一致：PIK3CA/mTOR 最高，其余三对不超过 0.61。

（对应 Methods 2.1）

### 3.3 对接对真双靶配体的判别能力有限且高度依赖靶对

双靶分数需要同时压住两条单靶选择性臂。若两臂共用同一池化分数（如两端 Vina 分数的均值），较强一臂可能掩盖较弱一臂的失败。本文因此以口袋匹配方向 AUROC 为主指标：dual 对 A_only 用口袋 B 的分数，dual 对 B_only 用口袋 A 的分数，summary_min 取两臂较小值（Methods 2.6）。分数定义为 \(S=-E_{\mathrm{Vina}}\)（越大越好），dual 为正类。

需要分开两件事：其一，**方向特异的判别失败**——例如在 EGFR/HER2 上，dual 对 B_only 的口袋匹配 AUROC 仅为 0.430，而若单独审视池化协议下较弱一臂可读到约 0.28，说明 B 端方向本身接近或低于随机，并非“池化运算把 0.50 压成 0.28”；其二，**池化会掩盖上述弱臂**——同一对上池化 summary 可接近 0.50，给人以“尚可”的假象。相对池化，口袋匹配普遍抬高了四对的点估计，但排序未变（Supporting Information Table S6）。

各面板在同一 RDKit/meeko 协议下用 AutoDock Vina 打分；RTMScore（九姿态取最优）与 GNINA 作通道对照。GNINA 最初仅对 Vina mode 1 rescore，与 RTM 姿态覆盖不对称；2026-08-24 已补做全 9 姿态公平重打（Supporting Information Table S14–S15），二者姿态覆盖对齐后 GNINA 结论未变，仍不作“三引擎已验证同一决策臂”的公平赛马主张。Bootstrap 95% 区间为配体层重采样（B = 2000，种子 20260729）。

**Table 2.** 冻结 K = 4 评价集上的口袋匹配方向 AUROC（Vina）。错口袋、配体效率归一与描述符基线见 Supporting Information Table S6。

| 靶对 | n (dual / A_only / B_only) | dual 对 A_only（口袋 B） | dual 对 B_only（口袋 A） | summary_min [95% CI] |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.284, 0.576] |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.440, 0.740] |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.347, 0.648] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.464, 0.802] |

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

PM110（历史命名；实际面板 n = 115，其中 dual / A_only / B_only 各 30）保留 PM48 全部 48 个配体并按配额扩样，用作**稳定性核对（stability check）**，不是独立验证实验，也不是用扩面“挽救”点估计。PM48 本身 dual / A_only / B_only 仅 18 / 14 / 12，功效有限。PM110 上 Vina summary_min 为 0.648 [0.51, 0.76]，相对 PM48 的 0.692（Δ ≈ −0.04），区间更窄，排序趋势保持一致；同面板 RTMScore 为 0.576，GNINA（全 9 姿态公平重打后，2026-08-24；取代此前仅重打 mode 1 的 0.522）为 **0.613** [0.46, 0.74]，PM48 同口径为 **0.655** [0.43, 0.81]（此前 mode-1 为 0.579）；GNINA 换用与 RTM 同等的 best-of-9 覆盖后数值上升，但仍不高于同面板 Vina 口袋匹配，排序结论不变（Supporting Information Table S15；`GNINA_POCKET_MATCHED_BEST9_VERDICT_V1.md`）。

综合 §3.2–3.5：**仅 PIK3CA/mTOR 的点估计在主面板与 PM110 稳定性核对上同时高于 0.5 与重原子数基线，幅度有限**；其余三对的表观信号在很大程度上可由配体属性或二维化学型解释。该方向随后在 unused-pool holdout 上保持（§3.9），但换晶体后显示受体依赖（§3.10），不得读成结构不变的可重复优势。

（对应 Methods 2.3–2.4、2.7）

### 3.6 对接判别力的跨对结构决定因素（探索性）

四对之间判别力差异悬殊，仅凭配体层混淆分析（§3.4）尚不能回答：这与两靶标口袋本身的结构相似度是否相关？我们从冻结受体 PDB 中提取每个靶标最长蛋白链的一级序列（Biopython `PDBParser`），用 BLOSUM62 全局两两比对（`PairwiseAligner`，gap open = −11、extend = −1）计算靶对内两靶标的**全链序列一致性**，作为整体结构相似度的粗粒度代理；脚本与原始输出见 `data/jcim_bench_v0/analysis/structural_context_v1/`（Supporting Information Table S7）。

**Table 3.** 靶对内全链序列一致性与判别力（summary_min）对照。

| 靶对 | 全链序列一致性（对齐长度为分母，%） | 全链序列一致性（较短链为分母，%） | summary_min（θ = 6.0 / 严格） |
|------|-----------------------------------:|---------------------------------:|-------------------------------|
| PIK3CA/mTOR | 18.1 | 21.0 | 0.692 / 0.639 |
| PIK3CA/PIK3CB | 40.5 | 43.3 | 0.500 / 0.500 |
| AChE/BChE | 51.9 | 53.1 | 0.606 / 0.606 |
| EGFR/HER2 | 71.4 | 76.6 | 0.430 / 0.324 |

在这四对上，判别力**最高**的 PIK3CA/mTOR 恰是全链序列一致性**最低**的一对；判别力**最低**的 EGFR/HER2 恰是序列一致性**最高**的一对（ErbB 家族激酶域高度同源，已知支持多种泛 ErbB/双靶抑制剂化学型，如拉帕替尼、阿法替尼）。这是 n = 4 的描述性对照，不是统计相关：不得写成“整体越相似就越难区分”的一般规则。

需要审慎解读的边界：这是全链一级序列一致性的粗粒度代理，不是经结构叠合限定的 ATP 口袋残基级 RMSD 或相互作用指纹（PLIF）相似度；后者需要经验证的结构叠合工具（如 TM-align、PyMOL align/super）并对齐口袋残基编号，本轮环境未配置、也未做残基级验证，不在此虚构口袋 RMSD 数值，留待后续工作。PIK3CA 与 mTOR 全链一致性最低，但二者同属 PIKK 相关家族，ATP 竞争位点存在已知的局部结构同源性——这正是文献中已有真实 PI3K/mTOR 双靶抑制剂化学型（如 PI-103、omipalisib）的结构基础；表中的低全链一致性不应解读为"两口袋不相似"，而应理解为"整体蛋白架构分化、但局部 ATP 位点保留可及重叠"，与 §3.7 的姿态级线索一致。四对样本量为 n = 4，本节为描述性对照，不做正式相关性检验或统计显著性主张。

（对应 Methods 2.6；Supporting Information Table S7）

### 3.7 PIK3CA/mTOR 有限方向信号的结构线索（个案级）

为避免将唯一较优靶对写成黑箱，我们基于已导出的姿态级诊断（failure typology；非全面板 PLIF 工程）归纳结构线索。PIK3CA（4L23）与 mTOR（4JT6）均为 ATP 竞争型激酶口袋；共晶配体 PI-103 / X6K 在协议检查中可回收近晶姿态（4JT6 需 elevating exhaustiveness 与 best-of-9，见 Table S3）。在代表性真 dual（如 PI-103）上，两端均可观察到 hinge 接触与对共晶配体的高占用；而若干 PIK3CA 选择性硬负（如 amino-triazine / morpholine–ATP 化学型）在**弱端 mTOR 上同样给出几何干净、hinge 阳性的姿态**，使两端分数同时偏高——这与“化学型同源假双靶”（T2）失败型一致，也说明即便在相对最好的靶对上，对接仍可能把 ATP 位点交叉化学型误判为 dual。

相反，部分经典 dual（如 Torin1、Omipalisib）在 Vina 上很强，但 RTM 优选姿态可偏离 PIK3CA hinge/共晶位，表现为重打分误伤（T5）。因此，PIK3CA/mTOR 的 summary_min 优势应解读为：**在共享 ATP 识别框架下、对部分化学型可出现有限方向信号**，而不是已验证的通用双靶决策规则；全面板残基级守恒分析与相互作用指纹定量比较超出当前冻结分析包，留待后续工作。

（对应已有 failure typology 与 cognate QC；非新对接）

### 3.8 尚未解决的稳健性缺口：面板构成层重抽样

现有 bootstrap（Methods 2.6；B = 2000）重采样的是**固定面板内**的配体，回答"这批已对接分子的不确定度有多大"，不回答"如果换一批同配额的分子，结论是否还成立"——即面板构成本身的抽样分布。要回答后者，需从各靶对严格供给池中按同一配额重复抽取大量独立 panel（如 1000 次），并对每次抽到的、尚未对接的池内分子跑 Vina/RTM/GNINA，再汇总 summary_min 的分布。这需要对当前冻结分数包之外的分子重新对接，超出本轮分析范围；`data/jcim_holdout_v0/` 已完成 unused-pool holdout 对接与评价（Results 3.9）；那是固定配额的一次外推检验，仍不等于对供给池做大量独立 panel 重抽的抽样分布。读者不应将§3.2–3.6 与 3.9 的结果等同于"面板构成已完全验证稳健"。

（对应 Methods 2.6 局限说明；`JCIM_SUPPLEMENTARY_EXPERIMENTS_PLAN_V2.md`）

### 3.9 面板外冻结验证集（unused-pool holdout）

在未参与面板构建与协议调优的 ChEMBL 严格池中，三对靶标各抽取 60 个配体（dual / A_only / B_only = 20 / 20 / 20；种子 20260731），按冻结协议对接后计算口袋匹配 summary_min（Supporting Information Table S8；`HOLDOUT_VERDICT.md`）。

| 靶对 | holdout summary_min [95% CI] | 主面板 summary_min | Δ(holdout−主面板) | 相对最强平凡基线 |
|------|-----------------------------:|-------------------:|------------------:|------------------|
| PIK3CA/mTOR | 0.765 [0.603, 0.891] | 0.692 | +0.073 | 点估计高于 heavy 0.555（Δ = +0.21） |
| AChE/BChE | 0.618 [0.422, 0.759] | 0.606 | +0.012 | 点估计略高于 cLogP 0.575（Δ = +0.043）；CI 仍含 0.5 |
| PIK3CA/PIK3CB | 0.425 [0.241, 0.618] | 0.500 | −0.075 | 点估计低于 heavy 0.691（Δ = −0.266） |

PIK3CA/mTOR 在 holdout 上点估计与主表同向：summary_min = 0.765，bootstrap 下界高于 0.5，并在点估计上继续高于最强平凡基线（heavy 0.555，Δ = +0.21）。AChE/BChE 点估计与主表接近，但 CI 仍跨越 0.5。PIK3CA/PIK3CB 仍未显示可用方向信号。含硼配体 HOAP_028 因 AutoDock 原子类型不支持而两端失败，已从 AUROC 装配中剔除（59/60 配体进入分析）。该 holdout 共享同一 ChEMBL 抓取批次，不能读成跨机构独立验证；其作用是检验协议在“建面未见过”的同规则配体上是否同向。

三对 holdout 上均出现 **Vina 错口袋对照不低于口袋匹配**（PM：0.788 对 0.765；AChE/BChE：0.643 对 0.618；PIK3CB：0.520 对 0.425）。为检验这是否只是 unused-pool 抽样在效价/尺寸上不同于主面板，我们按与 Table S5 相同的规则做最近邻匹配后重算错口袋对照（Supporting Information Table S13；零新对接）。效价匹配后三对仍是错口袋 ≥ 口袋匹配（AChE/BChE：0.642 对 0.593，n_min = 18；PIK3CA/PIK3CB：0.562 对 0.363，n_min = 11；PIK3CA/mTOR：0.734 对 0.715，n_min = 12）；尺寸匹配同样不翻转。holdout 相对主面板确有抽样偏移——最明显的是 PIK3CA/mTOR：holdout dual / A_only 的 pA 均值比主面板低约 1.1–1.3，B_only 的 pB 低约 1.8——但匹配后悖论仍在，故**不能**把错口袋 ≥ 匹配写成“holdout 抽得更强/更不平衡”的伪象。为检验这是否仅为打分函数伪象，我们在已冻结的 mode-1 姿态上计算不依赖打分的几何量 `contact_count`（配体重原子中与受体重原子距离 ≤4.0 Å 的原子数；粗粒度埋藏/接触计数，非 PLIF），并用它单独重复同一口袋比较（Supporting Information Table S11；`WRONG_POCKET_MECHANISM_VERDICT_V1.md`）。几何量在 **B 臂**上清晰高于随机（dual 对 B_only：0.698–0.714），与 dual 相对 B_only 的尺寸差更大一致（holdout 重原子均值：AChE/BChE 35.1 对 29.5；PM 33.5 对 31.0；PIK3CB 34.5 对 28.3）。**A 臂**上 dual 与 A_only 尺寸差很小（AChE/BChE 35.1 对 34.0；PM 33.5 对 32.3），contact_count AUROC 也更接近随机（0.552–0.622）。因此姿态层存在**真实的、不依赖打分函数的尺寸/埋藏混淆，主要出现在 B 臂**。但它**不能按幅度解释** Vina 的 holdout 错口袋，尤其是 PM：Vina 错口袋 summary_min = 0.788，而两臂 contact_count 的较小值仅为 0.552；PIK3CB 上该较小值（0.622）甚至高于 Vina 错口袋（0.520）。主面板在同一定义下口袋匹配仍高于错口袋（Supporting Information Table S6），与 holdout 的反差仍为开放问题，不得把接触计数或效价匹配写成“已解释 Vina 错口袋”。

### 3.10 替代晶体结构的稳健性（cognate QC + 面板重对接）

对替代晶体按 Methods 2.4 协议做共晶重对接：PIK3CA **4JPS**（1LT）best_of_9 = 0.607 Å、**5DXT**（5H5）= 0.624 Å；mTOR **4JSX**（Torin2/17G）= 0.515 Å，**三者均通过** &lt; 2 Å 门槛（`STRUCTURE_ROBUSTNESS_QC_V1.md`）。嵌合体 3T8M 已排除。

在冻结 PM48 配体上将口袋 A 换为 4JPS 或 5DXT（B 端仍用原 4JT6 分数）后，口袋匹配 summary_min 由主面板的 0.692 分别降至 **0.486** [0.259, 0.692] 与 **0.505** [0.292, 0.696]；D/A 臂（依赖未更换的 4JT6）保持 0.714，下降集中在依赖新 PIK3CA 口袋的 D/B 臂。将口袋 B 换为 4JSX（A 端仍用原 4L23）后，summary_min 为 **0.639** [0.418, 0.776]（Δ ≈ −0.05）。按方案预声明判据，PIK3CA 端记为**受体依赖**：4L23 上的有利信号不能自动外推到其他已过 QC 的 PIK3CA 晶体。mTOR 端换晶后点估计仍高于 0.5，但 95% CI 包含 0.5，只表明优势减弱、不能读成该端已稳健。4JPS / 5DXT 的 CI 同样包含 0.5。完整记录见 `STRUCTURE_ROBUSTNESS_VERDICT_V1.md`。

### 3.11 受体依赖的探索性结构对照：晶体间 Cα 差异

为对照 3.10 中 PIK3CA 端换晶即崩、mTOR 端换晶仅小幅下降这一不对称，我们在已有晶体坐标上做刚体叠合（探索性；PIK3CA 替代结构 n = 2，mTOR n = 1；`POCKET_MECHANISM_VERDICT_V1.md`）：以各自共晶配体重原子 ≤5 Å 界定口袋残基（PIK3CA 20 个、mTOR 18 个，均以残基编号+残基名精确匹配，匹配位点零错配），用全部匹配 Cα 做单次刚体叠合，分别报告**全域 Cα RMSD**（整链）与**同一变换下的口袋局域 Cα RMSD**。

| 参考（主面板） | 替代 | 匹配 Cα | 全域 Cα RMSD | 口袋局域 Cα RMSD | 共晶配体质心距离 |
|------|------|------:|------:|------:|------:|
| 4L23（PIK3CA） | 4JPS | 982 | **1.486 Å** | 0.867 Å | 2.566 Å |
| 4L23（PIK3CA） | 5DXT | 862 | **1.441 Å** | 0.343 Å | 2.072 Å |
| 4JT6（mTOR） | 4JSX | 1054 | **0.454 Å** | 0.467 Å | 2.196 Å |

结果与换晶后的方向**相一致，但不能定量解释因果**：**这批** PIK3CA 沉积结构彼此的整链 Cα RMSD（1.44–1.49 Å）大于**这批** mTOR 沉积结构彼此的差异（0.45 Å）。比较有边界：5DXT 仅匹配 862 个 Cα，少于 4JPS 的 982 个，全域 RMSD 不是同一覆盖度下的等长比较；替代结构各仅 1–2 个，不能外推为该靶标的一般晶体变异规律。局域口袋 Cα **并不**随 AUROC 崩溃幅度变化：5DXT 口袋局域 RMSD 仅 0.343 Å（小于全域 1.441 Å），ATP 位点骨架保守，但 summary_min 仍降至 0.505，与口袋局域偏移更大的 4JPS（0.867 Å，0.486）几乎相同。因此 **Cα 层面的口袋保守性不足以保证判别力可迁移**；侧链、质子化或搜索空间等未测因素仍可能参与——本轮无 PLIF/rotamer 证据，不作机制已解决主张。共晶配体质心距离 2.1–2.6 Å 只说明对接的仍是同一大类 ATP 竞争位点，不能证明结合模式相同。

（对应 `data/jcim_structure_robust_v0/analysis/pocket_mechanism_v1/`；零新对接，仅用已冻结晶体坐标）
