# Results（中文工作稿 · JCIM 式证据链）

> 结构按 *performance → confounding → robustness boundary → failure mode* 重排为 6 节，不再按实验清单写 3.1–3.11。  
> Results 只报告数字、预定义比较与直接可见的模式；解释、机制假说与使用边界放 Discussion。  
> 全部数字可追溯至 `data/jcim_bench_v0/`、`data/jcim_strengthen_t0t1_v0/`、`data/jcim_holdout_v0/`、`data/jcim_structure_robust_v0/` 与 `data/jcim_supply_crossdb_v0/`。未做的全面板残基级 PLIF 不写入。  
> **主图：** Fig 1 任务定义；Fig 2 供给；Fig 3 主森林图；Fig 4 弱臂/描述符；Fig 5 配体层 vs 受体层；Fig 6 错口袋反转；Fig 7 化学型/协变量检验。详见 [`FIGURE_AND_TOC_PLAN_JCIM_V1.md`](FIGURE_AND_TOC_PLAN_JCIM_V1.md)。  
> **定位：** DualFourClass-Bench 是 *constrained but experimentally grounded* 的四状态评测资源，不是“完整双靶对接排行榜”。禁止 D-DRAF；禁止 “docking can/cannot identify dual-target ligands”。见 [`POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`](POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md)。

---

## 3. Results

### 3.1 Experimental Data Supply Limits Strict Dual-Target Benchmark Construction

为确定公开生物活性数据是否能够支持严格的双靶点识别评测，我们首先对 49 对有 ChEMBL 缓存的候选靶标进行供给审计（Figure 2）。双靶对接评测需要四种实验标签状态：dual、A-only、B-only 与 neither（Figure 1A）。基准是四状态数据集；预先指定的主终点是两条方向 pairwise 判别（dual 对 A-only、dual 对 B-only）。我们将一端达到活性阈值、对端明确低活性的配体定义为方向性选择性硬负样本，用于检验对接分数能否同时压住两条单靶臂。

在严格标签规则下（dual：两端 pChEMBL ≥ 6.5；选择性类：活性端 ≥ 6.5 且对端 ≤ 5.5），尽管候选靶对数量较多，能够同时提供足量 A-only 与 B-only 硬负样本的靶对十分有限。两端严格硬负均不少于 50 的厚面板条件仅有 4 对满足。排除不适合作为常规小分子对接评测对象的金属依赖 HDAC1/HDAC6 后，PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB 构成三个规模相对充足的靶对；EGFR/HER2 仅有 7 个严格 B-only 配体，因此被保留为供给受限案例，而不是与前三对等价的厚面板（Table 1）。

这一供给限制并非 ChEMBL 单一数据库特有的计数现象。对最终四对靶标进行 BindingDB 与 PubChem 的零对接计数核对后（Supporting Information Table S12），与 pChEMBL 更接近的 `equal_only` 规则下，前三对的 min HN 分别为 BindingDB 76 / 92 / 58、PubChem 86 / 97 / 61（ChEMBL 缓存为 80 / 78 / 56），仍全部 ≥ 50。EGFR/HER2 虽可在其他数据库中达到约 30 个 B 端硬负样本，仍不足以满足 ≥ 50 的厚面板标准。将不等式活性记录作为点估计（`as_is`）会显著增加 EGFR/HER2 的表观供给（BindingDB min HN 升至 85），但 92 个 as-is B_only 中有 49 个在 EGFR 端只有 `>` 记录；这一处理改变了“两端具有等式定量测定”的标签定义，因此未用于冻结 benchmark。PubChem 与 BindingDB 计数接近，符合沉积重叠，不作两次独立普查。

因此，最终 benchmark 的规模并非根据对接表现事后筛选，而主要由公开实验数据中方向性选择性硬负样本的可获得性所约束。DualFourClass-Bench 是一套受供给约束、但由实验标签锚定的评价集（a constrained but experimentally grounded benchmark），不是覆盖全部双靶任务的完整抽样。构建细节见 Methods 2.1–2.3。

### 3.2 Docking Provides Limited and Target-Pair-Dependent Dual-Target Discrimination

在冻结的四对靶标上，采用统一 θ = 6.0 标签规则和口袋匹配方向 AUROC 对 Vina docking scores 进行评价（Figure 1B；Methods 2.8）。分数定义为 \(S=-E_{\mathrm{Vina}}\)（越大越好），dual 为正类；summary_min 取两臂较小值，使较强一臂不能掩盖较弱一臂。AChE/BChE 与 PIK3CA/PIK3CB 在建造时使用更严格的 6.5/5.5 规则，但在本数据上 θ = 6.0 给出完全相同的配体分类与 AUROC（Table S4）；EGFR/HER2 与 PIK3CA/mTOR 对阈值更敏感，严格规则下 B_only 过少并标记 underpowered，故严格规则只作支持性敏感性分析，不作第二套主标准。整张阈值网格内排序趋势保持一致（Figure S1A）。

EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 的 summary_min 分别为 0.430、0.606、0.500 和 0.692（Table 2；Figure 3）。两条方向性比较进一步显示，不同靶对的主要限制来自不同的弱臂：EGFR/HER2 的 dual-versus-B-only AUROC 为 0.430，而 PIK3CA/PIK3CB 的对应值为 0.500；相比之下，PIK3CA/mTOR 两个方向分别达到 0.714 和 0.692（Figure 4A）。相对池化协议，口袋匹配普遍抬高了四对的点估计，但排序未变（Table S6）。EGFR/HER2 按池化 Vina 取 Top-10 时，9 个为硬负选择性配体（bootstrap 均值 ≈ 8.9；95% CI 为 7–10）。

同一套冻结分数再按 Zhou 式 Dual-versus-neither 以及 Dual versus all non-duals 计分（Table 3；Figure S4）。EGFR/HER2 上 Dual versus neither（`vina_mean`）AUROC 为 0.756 [0.562, 0.920]（n_neg = 12），而方向性 summary_min 仍为 0.430 [0.284, 0.576]；Dual versus all non-duals 降至 0.551 [0.443, 0.666]。AChE/BChE 与 PIK3CA/PIK3CB 的 Dual-versus-neither 增量很小（0.649 与 0.559），区间与方向性臂重叠。PIK3CA/mTOR 的 Dual versus neither 因 neither n = 4 而 underpowered，不作解释；该对 Dual versus all non-duals 为 0.674，接近 summary_min 0.692。因此 benchmark formulation 在 EGFR/HER2 上改变了解释，但没有在四对上给出统一反转。

**Table 2.** 冻结 K = 4 评价集上的口袋匹配方向 AUROC（Vina，统一 θ = 6.0）。错口袋、配体效率归一与描述符基线见 Table S6。

| 靶对 | n (dual / A_only / B_only) | dual 对 A_only（口袋 B） | dual 对 B_only（口袋 A） | summary_min [95% CI] |
|------|---------------------------:|-------------------------:|-------------------------:|----------------------|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | 0.430 [0.284, 0.576] |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | 0.606 [0.440, 0.740] |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | 0.500 [0.347, 0.648] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | 0.692 [0.464, 0.802] |

**Table 3.** 同一套 Vina 分数在传统 formulation 与方向性 formulation 下的 AUROC（统一 θ = 6.0）。Dual-versus-neither 使用实验 inactive；Dual versus all non-duals 把 A-only、B-only 与 neither 都计为负类。方向性 summary_min 的 CI 来自 Table 2。PIK3CA/mTOR Dual versus neither underpowered（n_neg = 4）。

| 靶对 | directional summary_min [95% CI] | Dual vs neither (`vina_mean`) | n_neither | Dual vs all non-duals |
|------|--------------------------------:|------------------------------:|----------:|----------------------:|
| EGFR/HER2 | 0.430 [0.284, 0.576] | 0.756 [0.562, 0.920] | 12 | 0.551 [0.443, 0.666] |
| AChE/BChE | 0.606 [0.440, 0.740] | 0.649 [0.484, 0.812] | 15 | 0.579 [0.442, 0.716] |
| PIK3CA/PIK3CB | 0.500 [0.347, 0.648] | 0.559 [0.373, 0.746] | 16 | 0.556 [0.437, 0.672] |
| PIK3CA/mTOR | 0.692 [0.464, 0.802] | 0.514 [0.222, 0.806] | 4 | 0.674 [0.515, 0.817] |

因此，docking discrimination 并未表现出一致的跨靶对能力，而呈现明显的 target-pair dependence。PIK3CA/mTOR 是唯一一个 summary_min 点估计同时高于 0.5 和其最优简单物化性质基线（重原子数 0.463）的靶对，但其 95% bootstrap confidence interval 仍为 0.464–0.802，未能排除随机水平。AChE/BChE 的 summary_min 为 0.606，但低于 TPSA 基线（0.733）；EGFR/HER2 和 PIK3CA/PIK3CB 分别为 0.430 和 0.500，也未显示超过相应简单基线的明确优势。四个预先指定的物化描述符全部计算；“最强描述符”仅指该面板上 AUROC 最高者，不是预设的独立检验基准。

采用 RTMScore 或 GNINA 作为替代 scoring channel 未改变这一总体排序。GNINA 在统一 best-of-9 pose coverage 后，EGFR/HER2、AChE/BChE 与 PIK3CA/mTOR 的口袋匹配 summary_min 仍不超过同面板 Vina；PIK3CA/PIK3CB 上 GNINA best-of-9 为 0.533、Vina 为 0.500，二者均近随机且区间重叠（Table S14–S15；Figure S1B）。GNINA 仍只是单一 CNN 通道对照，不作三引擎公平赛马主张。

### 3.3 Ligand Properties and Chemotype Explain Much of the Apparent Signal

为判断 docking discrimination 是否超越简单的 ligand-level signal，我们首先将口袋匹配对接与四种预先定义的物化性质进行比较（Figure 4B）。相对于每个靶对中表现最佳的简单 descriptor，docking summary_min 的 paired difference 在 EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 和 PIK3CA/mTOR 中分别为 −0.052、−0.128、−0.122 和 +0.229；四个 95% confidence intervals 均包含 0（Table S19；Figure S3C）。由此可见，即使 PIK3CA/mTOR 的点估计表现出最大的正向差异，现有样本仍不足以将其与简单 ligand-property baseline 明确区分。该对照使用口袋匹配 summary_min，不是 pooled `vina_mean` 门控（EGFR/HER2 的 `vina_mean` 为 0.2824，≠ Table 2 的 0.4297）。

AChE/BChE 提供了一个较为直接的混淆案例。dual 配体平均 TPSA 约为 75，而选择性硬负配体约为 51（Figure 4C）；TPSA 单独获得约 0.769 的 AUROC，高于相同比较下的 Vina（约 0.56）。进一步加入 heavy-atom count 和 TPSA 后，dual-versus-B-only AUROC 从 0.606 增至 0.807，而 docking score 的 OR 仅约为 1.18（Figure 7C）。该结果表明，该方向上的 docking discrimination 很大程度上依赖与配体物化性质相关的信号，而不能直接解释为独立的 pocket-specific information。

PIK3CA/mTOR 的情况有所不同。加入 heavy-atom count 和 TPSA 后，AUROC 的变化约为 +0.07 至 +0.11，docking score 的 OR 约为 2.19 和 3.08，提示该靶对可能存在一定 residual pocket-related signal；然而，与 descriptor 的 paired difference 置信区间仍包含 0，因此这一残余信号不能被视为已确证的独立优势。配体效率归一后，仅 PIK3CA/mTOR 仍高于重原子数基线（0.657 对 0.463）。

二维化学结构 baseline 进一步说明了这一问题（Figure 7A）。ECFP4 + logistic regression 在 Bemis–Murcko scaffold GroupKFold 下多个方向获得约 0.78–0.91 的 fold AUROC，明显高于部分对应 docking contrasts，例如 EGFR/HER2 dual-versus-B-only 中 ECFP4 AUROC 为 0.85，而 docking AUROC 仅为 0.43。同一设定下随机 `StratifiedKFold` 相对支架折的平均差为 +0.011（八个方向对比；Table S20；Figure S3D），泄漏很小。该结果表明 dual/selective 标签与 chemotype 存在系统性关联，因此单独观察 docking score 的 AUROC 并不足以证明其识别来源于 pocket-specific physical interactions。

在同一支架折上把口袋匹配对接分数加到 ECFP4 后，CV AUROC 的变化至多约 0.01，若干方向为负（Table S24）。logistic docking AUROC 不是 Table 2 的 rank AUROC，且常常更低。ECFP4 Tanimoto ≥ 0.7 的 chemotype-matched A-only/B-only 子集为空。T ≥ 0.3 时，未匹配时最强的一臂（PIK3CA/PIK3CB dual versus A-only，0.691）降至 0.503（n_neg = 11），而远缘硬负（T < 0.3）升至 0.819（Table S23）。T ≥ 0.4/0.5 的格子常为 n_neg ≤ 7，不作为第二套主结果解释。

效价匹配或尺寸匹配子集上，EGFR/HER2 与 PIK3CA/PIK3CB 的 dual 对 B_only 仍偏弱或接近随机（约 0.45–0.52）；PIK3CA/mTOR 的排序趋势保持一致，但各臂 n 常低于 15、区间较宽（Table S5；Figure 7D）。全部四个平凡描述符见图 7B。

### 3.4 The PIK3CA/mTOR Signal Persists across Ligand Panels but Not Receptor Realizations

为判断 PIK3CA/mTOR 的较高 summary_min 是否仅由特定 panel 构成或 docking 搜索参数造成，我们进行了 ligand-panel 和 protocol-level sensitivity analyses（Figure 5A）。将 exhaustiveness 从 16 降至 8 后，summary_min 从 0.692 降至 0.660，变化约 0.03，明显小于不同 target pairs 之间的性能差异（Figure S1D）。

在包含 PM48 全部配体并扩展至实际 n = 115 的 PM110 面板中（分析用 dual / A_only / B_only 各 30），Vina summary_min 为 0.648 [0.51, 0.76]，相比 PM48 的 0.692 下降约 0.04，但排序趋势保持一致（Figure S1C）。该结果支持 PIK3CA/mTOR 的方向性信号并非完全由 PM48 的特定成员驱动，但 PM110 与 PM48 并非独立验证集，因此该结果应解释为 stability check。同面板 RTMScore 为 0.576；GNINA best-of-9 为 0.613 [0.46, 0.74]，PM48 同口径为 0.655 [0.43, 0.81]，仍不高于同面板 Vina。

更重要的是，在未参与主面板构建和协议调优的 unused-pool holdout 中（每对 20 / 20 / 20，种子 20260731；EGFR/HER2 不具备同等配额，记为 not eligible），PIK3CA/mTOR 的 summary_min 进一步达到 0.765 [0.603, 0.891]，高于主面板的 0.692；AChE/BChE 为 0.618 [0.422, 0.759]，与主面板接近但 confidence interval 跨越 0.5；PIK3CA/PIK3CB 则下降至 0.425 [0.241, 0.618]（Table S8 / Table S16）。含硼配体 HOAP_028 因 AutoDock 原子类型不支持而两端失败，已从 AUROC 装配中剔除（59/60 配体进入分析）。该 holdout 共享同一 ChEMBL 抓取批次，不能读成跨机构独立验证；其作用是支持所观察信号在未参与建面配体池中的持续性（persistence of the observed signal in an unused ligand pool）。

因而，PIK3CA/mTOR 的方向性 signal 在同一 ChEMBL 体系的未参与建面配体中仍然可观察到，而 PIK3CA/PIK3CB 的 signal 则未能保持。这进一步说明 docking performance 主要由 target-pair context 决定，而不是一个可在不同靶对之间稳定迁移的属性。

尽管 PIK3CA/mTOR 在 ligand-panel sensitivity analysis 中保持了方向性信号，我们进一步测试这一结果是否依赖于特定 receptor realization（Figure 5B）。三个替代晶体结构均通过 cognate redocking QC，best-of-9 RMSD 分别为 0.607 Å（4JPS）、0.624 Å（5DXT）和 0.515 Å（4JSX）；嵌合体 3T8M 已排除（Table S9）。

当 PIK3CA 4L23 替换为 4JPS 或 5DXT，而 mTOR 4JT6 保持不变时，PM48 的 summary_min 分别由 0.692 降至 0.486 [0.259, 0.692] 和 0.505 [0.292, 0.696]。变化主要发生在依赖替代 PIK3CA 结构的 D/B direction，而依赖原始 mTOR 结构的 D/A direction 保持 0.714。相比之下，将 mTOR 4JT6 替换为 4JSX 后 summary_min 为 0.639 [0.418, 0.776]。按方案预声明判据，PIK3CA 端记为受体依赖：4L23 上的有利信号不能自动外推到其他已过 QC 的 PIK3CA 晶体。mTOR 端换晶后点估计仍高于 0.5，但 95% CI 包含 0.5，只表明优势减弱，不能读成该端已稳健。

因此，PIK3CA/mTOR 的较高 summary_min 并不是 receptor-independent 的性质。尤其是 PIK3CA 端，在两个通过 cognate QC 的替代晶体上均接近随机水平。Cα structural comparison 进一步显示，5DXT 与 4L23 的口袋局域 Cα RMSD 仅为 0.343 Å，但 summary_min 仍降至 0.505，说明简单的 backbone similarity 并不足以解释判别性能的保持（Table S10）。这批 PIK3CA 沉积结构彼此的整链 Cα RMSD（1.44–1.49 Å）大于这批 mTOR 沉积结构彼此的差异（0.45 Å），与换晶后的不对称方向一致，但不能作为定量因果解释：5DXT 仅匹配 862 个 Cα，少于 4JPS 的 982 个；替代结构各仅 1–2 个。共晶配体质心距离 2.1–2.6 Å 只说明对接的仍是同一大类 ATP 竞争位点。姿态生成 QC 通过，并不等于 screening discrimination 可迁移。

### 3.5 Wrong-Pocket Controls Reveal an Unresolved Failure Mode

在主面板中，pocket-matched summary_min 均高于 wrong-pocket control，四对的 matched-minus-wrong differences 分别为 0.170、0.161、0.151 和 0.090；其中 EGFR/HER2 和 AChE/BChE 的差异置信区间排除 0，PIK3CA/PIK3CB 与 PIK3CA/mTOR 的区间包含 0（Table S6；Table S17；Figure 6A；Figure S3A）。错口袋对照的 summary_min 分别为 0.260、0.444、0.349 与 0.602。

然而，这一关系在 unused-pool holdout 中发生反转（Figure 6B）。PIK3CA/mTOR、AChE/BChE 和 PIK3CA/PIK3CB 的 wrong-pocket summary_min 分别为 0.788、0.643 和 0.520，而 matched-pocket 分别为 0.765、0.618 和 0.425。相应的 matched-minus-wrong point differences 均为负（−0.023 / −0.025 / −0.095），但其 95% confidence intervals 均包含 0（Table S17；Figure S3B）。holdout 反转是点估计格局，不是区间已排除零的统计结论。

为检验这一反转是否可以由 holdout 中的配体效价或分子大小差异解释，进一步进行 potency- and size-matched comparisons（Figure 6C；Table S13）。wrong-pocket ≥ matched-pocket 的关系仍未翻转（效价匹配后：AChE/BChE 0.642 对 0.593，n_min = 18；PIK3CA/PIK3CB 0.562 对 0.363，n_min = 11；PIK3CA/mTOR 0.734 对 0.715，n_min = 12）。holdout 相对主面板确有抽样偏移——最明显的是 PIK3CA/mTOR：holdout dual / A_only 的 pA 均值比主面板低约 1.1–1.3，B_only 的 pB 低约 1.8——但匹配后悖论仍在。

scoring-independent contact_count 在 B direction 获得 0.698–0.714 的 AUROC，提示 ligand size/burial 对该方向确实存在影响，但其幅度不足以解释全部 Vina wrong-pocket signal（Figure 6D；Table S11）。例如 PIK3CA/mTOR 中 Vina wrong-pocket summary_min 为 0.788，而 contact_count 的较弱一臂仅为 0.552。A 臂上 dual 与 A_only 尺寸差很小，contact_count AUROC 更接近随机（0.552–0.622）。

因此，holdout 中 wrong-pocket reversal 应被视为当前 benchmark 暴露出的 unresolved failure mode，而不是可以由单一尺寸或效价因素解释的现象。

### 3.6 Structural Context Provides Only Exploratory Clues

作为探索性分析，我们进一步比较了靶对内全链序列一致性与 summary_min（Table S7）。四对靶标中，PIK3CA/mTOR 具有最低的全链序列一致性（对齐长度分母 18.1%）但最高的 summary_min，而 EGFR/HER2 恰好相反（71.4%；ErbB 家族激酶域高度同源）。该现象与简单的“靶点越相似越难区分”假设并不一致。然而，该分析仅包含四个靶对，且全链序列一致性并不是 binding-pocket similarity 的直接度量，因此该观察仅作为结构背景线索，而不作为相关性证据。PIK3CA 与 mTOR 同属 PIKK 相关家族，ATP 竞争位点存在已知的局部结构同源性；表中的低全链一致性不应解读为“两口袋不相似”。

在已导出的姿态级诊断中，PIK3CA/mTOR 上可观察到两类代表性 failure typology（非全面板 PLIF）。T2：选择性硬负在两个口袋都形成几何干净、hinge 阳性的 ATP-like pose（例如 amino-triazine / morpholine–ATP 化学型在弱端 mTOR 上仍占用高、hinge 阳性），使两端分数同时偏高。T5：部分经典 dual（如 Torin1、Omipalisib）在 Vina 中较强，但 alternative rescoring 的优选姿态可偏离 PIK3CA hinge/共晶位。共晶配体 PI-103 / X6K 在协议检查中可回收近晶姿态（Table S3）。这些是观察到的姿态模式，不是已验证的残基级机制。
