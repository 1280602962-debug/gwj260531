# Results（中文稿 · JCIM Evaluation Article）

> DualFourClass-Bench 正文 Results 整合稿（v1）  
> 结构：短实验小节（Cieplinski / Vu 式）+ 完整主表；措辞按 claim ceiling。  
> 数字来源：`jcim_j0j1_v0`、`jcim_bench_v0/tables/pocket_matched_*`、`jcim_strengthen_t0t1_v0`。  
> 写作约束：话题式小标题；数字进句；阴性平铺；每节最多一句克制解读；禁“通吃 scorer / 显著优于基线通吃四对”。

---

## 3. Results

### 3.1 Strict hard-negative supply

在 49 对可审计公共靶对上，按 dual = 两端 pChEMBL ≥ 6.5、A_only / B_only = 一端 ≥ 6.5 且对端 ≤ 5.5 统计严格硬负。两端硬负均 ≥ 50 的只有 4 对；去掉金属酶 HDAC1/HDAC6 后，适合常规对接主面板的约 3 对：PIK3CA/mTOR、AChE/BChE、PIK3CA/PIK3CB（Table 1；审计表见 J0）。EGFR/HER2 的严格 B_only 仅 7 个，达不到厚面板门槛，因此保留为供给受限案例，不作新对接扩面。后续 K=4 冻结集即由此供给结果选出，而不是按文献热度铺开。

### 3.2 Pooled score versus pocket-matched directional AUROC

双靶决策需要同时压住两条单靶硬负臂。若两臂共用同一池化分数（如 `vina_mean`），弱臂会被强臂抬高：EGFR/HER2 上池化汇总 AUROC 约 0.50，而方向分解后弱臂（D vs B_only）可落到 0.28 附近（Figure 1；`asymmetry_pooled_vs_directional_v1`）。主指标因此定为口袋匹配方向 AUROC：D vs A_only 用口袋 B 分数，D vs B_only 用口袋 A 分数，并以两臂 AUROC 的最小值 `summary_min` 汇总。池化分数只作对照。升格后四对数值普遍上移，但排序不变：仍只有 PIK3CA/mTOR 明显高于随机，其余三对 ≤ 0.61（Table 2）。

### 3.3 Directional discrimination on the frozen K=4 set

AutoDock Vina、RTMScore best-of-K 与 GNINA CNN（mode_01 rescore）在同一 RDKit/meeko 配体协议下重打分。主报告以 Vina 口袋匹配为准；RTM/GNINA 作通道对照。Bootstrap 95% CI（B = 2000，seed = 20260729）见表 2 与 Figure 2。

**Table 2.** Pocket-matched directional AUROC (Vina), K=4.

| Target pair | n (D / A / B) | D vs A (pocket B) | D vs B (pocket A) | summary_min [95% CI] | Wrong-pocket min | LE-PM min | Best trivial baseline |
|---|---:|---:|---:|---|---:|---:|---|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 | 0.430 | **0.430** [0.281, 0.576] | 0.260 | 0.311 | cLogP 0.482 |
| AChE/BChE | 27 / 25 / 28 | 0.650 | 0.606 | **0.606** [0.442, 0.737] | 0.444 | 0.413 | TPSA 0.733 |
| PIK3CA/PIK3CB | 28 / 27 / 28 | 0.691 | 0.500 | **0.500** [0.340, 0.648] | 0.349 | 0.332 | heavy 0.622 |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 | 0.692 | **0.692** [0.457, 0.813] | 0.602 | 0.657 | heavy 0.463 |

PIK3CA/mTOR 是唯一 `summary_min` 点估计明显高于 0.5、且高于 heavy-atom 基线（0.463；Δ ≈ +0.23）的一对；其 CI 下界仍贴近 0.5，不足以写成“显著压过平凡基线”。EGFR/HER2 与 PIK3CA/PIK3CB 的 `summary_min` 分别落在 0.43 与 0.50，点估计低于各自最优平凡基线。AChE/BChE 的对接点估计（0.606）低于 TPSA（0.733）。RTM 与 GNINA 未把上述格局翻过来：在 PIK3CA/mTOR 上二者也弱于 Vina；在阴性对上同样通不过基线门。

排名式读出更直接。EGFR 面板按 Vina 池化分取 Top-10 时，硬负占 9/10（bootstrap 均值约 8.9，CI 约 7–10）。对接分数把许多单靶硬负排到前面，与方向性弱臂一致。

### 3.4 Baseline gate

对每对取最优平凡描述符（heavy / MW / cLogP / TPSA）与对接 `summary_min` 的差 Δ。在池化主表的 bootstrap 门控中，EGFR 与 PIK3CA/PIK3CB 的 Δ CI 整体落在 0 以下（对接显著劣于基线）；AChE 点估计 FAIL，Δ CI 刚跨 0；PIK3CA/mTOR 点估计胜过 heavy，但 Δ 的 95% CI 仍跨 0（Figure 3）。口袋匹配升格后，EGFR 相对 cLogP、AChE 相对 TPSA、PIK3CB 相对 heavy 的点估计差距并未消失。基线门的作用不是“证明对接无用”，而是把可发表主张钉在评测诊断，而不是通用决策臂。

### 3.5 Confounds: wrong pocket, size, polarity, and matched subsets

若分数只反映分子属性而非口袋匹配，错口袋对照应偏离 0.5。四对错口袋 `summary_min` 为 0.260（EGFR）、0.444（AChE）、0.349（PIK3CB）、0.602（PIK3CA/mTOR）；口袋匹配相对错口袋的 specificity gap 均 > 0.09（Table 2；`pocket_specificity_gap_v1`）。配体效率归一后，仅 PIK3CA/mTOR 仍高于 heavy 基线（LE-PM 0.657 vs 0.463）；其余三对 LE 后不再支撑方向信号。

AChE/BChE 上化学型捷径清楚：dual 平均 TPSA ≈ 75，硬负 ≈ 51；TPSA 区分 dual 与硬负的 AUROC ≈ 0.769，高于同设定下的 Vina（约 0.56）。把 heavy 与 TPSA 作为协变量进入 logistic 后，AChE 的 D vs B 口袋匹配 AUROC 从 0.606 升到 0.807（Δ ≈ +0.20），说明该臂上对接分与尺寸/极性高度共线。PIK3CA/mTOR 的协变量调整幅度较小（约 +0.07–0.11）。

效价匹配（|ΔpChEMBL| ≤ 0.5）或尺寸匹配（|Δheavy| ≤ 2）子集上，EGFR 与 PIK3CA/PIK3CB 的 D vs B 仍弱或近随机（约 0.45–0.52），不是简单用“dual 更大/更强效”就能解释完的。PIK3CA/mTOR 在匹配子集上方向仍在，但每臂 n 往往 < 15，CI 宽，只作 SI。

### 3.6 Two-dimensional ligand baseline under scaffold CV

用 ECFP4 + 逻辑回归作配体侧基线，交叉验证按 Murcko 支架分组（GroupKFold），避免同系物跨折。支架折 AUROC 多数落在 0.78–0.91，普遍高于对应口袋匹配对接臂；例如 EGFR D vs B：ML 0.85，dock 0.43。随机折相对支架折平均只高约 0.01（最大约 0.10），因为许多支架接近 singleton，分组折对泄漏的压制有限——这一点写进 Limitations，正文不以随机折数字作主张。2D 基线说明标签与化学型相关；它不证明对接“被碾压”，只限制把方向 AUROC 解读成口袋物理特异性。

### 3.7 Robustness: exhaustiveness, single-target enrichment, and PM110

PIK3CA/mTOR 主面板在 E = 16 下得到 `summary_min` = 0.692；同一配体/盒子/seed 重对接到 E = 8 为 0.660（Δ = +0.032）。exhaustiveness 差解释不了该对相对其他对的优势，正文并列报告 E = 8。

单靶 sanity 用同靶 ChEMBL 弱效已测分子（pChEMBL ≤ 5.5）作 property-matched decoy，而非随机无关物。4L23（PIK3CA）与 4JT6（mTOR）的 enrichment AUROC 分别为 0.603 与 0.629：有弱辨别，不是强虚拟筛选引擎。

PM110 是 PM48 的超集扩样（原 48 全部保留 + 严格配额新增），不是独立复制实验。Vina 口袋匹配 `summary_min` = 0.648 [0.51, 0.76]（相对 PM48 的 0.692，Δ ≈ −0.04），CI 收窄、方向不变。同面板 RTM 0.576、GNINA 0.522，均未增强。配体准备对照上，PM48 的 RDKit 主协议 Vina `summary_min`（池化）约 0.671，旧 LigPrep 姿态约 0.597；主表只用 RDKit，LigPrep 不进主表混报。

---

## 写作备忘（不进正文）

- 小节标题保持话题式（supply / metric / forest / baseline / confound / ML / robustness），勿改成结论句。  
- Discussion 再收：评测主张、claim ceiling、K 小与 singleton 支架、无湿实验。  
- 英译时保留 Table 2 数字与“CI spans / does not exclude 0”的克制语气；勿加成 “robustly outperforms”。  
- 环境中无独立 nature-paper skill；本稿按 Nature/JCIM 常见约束人工去 AI 化：短句与长句交错、少路标词、阴性不辩解、解读压到一句。
