# Methods（中文工作稿 · JCIM Articles）

> DualFourClass-Bench · 与 [`RESULTS_DRAFT_ZH_JCIM_V1.md`](RESULTS_DRAFT_ZH_JCIM_V1.md) 配套  
> 写法对照：Vu et al. JCIM 2025（Dataset → Docking → Scoring → Evaluation）式短小节；参数进句；过去时报告操作。  
> 审稿人视角已边写边改：标签规则不统一、EGFR 角色、PM110 超集、enrichment 参数以代码为准等已写入正文，不藏进脚注。

---

## 2. Methods

### 2.1 任务定义

本文将双靶对接评测建成四类配体判别任务。对每一对靶标 A/B，依据两端实测活性将配体分为：dual（两端均强）、A_only（仅 A 强）、B_only（仅 B 强）和 neither（两端均弱）。评价目标不是比较某一口袋上谁的对接分更高，而是判断对接分数能否把 dual 配体同时与两类单靶选择性配体（A_only 与 B_only）区分开。

为此，我们分别计算 dual 对 A_only、dual 对 B_only 的 AUROC，并以二者的最小值作为该靶对的汇总指标（记为 summary_min）。采用最小值是为了避免只报告较好一侧：若仅有一端可区分而另一端接近随机，则不足以支持“双靶判别”的主张。

### 2.2 数据策展与靶对选择

活性数据来自 ChEMBL。pChEMBL 为 ChEMBL 提供的统一活性标度，定义为相关摩尔活性值（如 IC50、Ki、Kd、EC50 等）的 −log10。同一配体–靶标若有多条记录，本文取可用的最大 pChEMBL（最强活性）作为该对的代表值；任一端缺少可用 pChEMBL 的配体不进入四类主分析。该聚合方式的局限（未按 assay 中位数或置信度过滤）在 Limitations 中说明。ChEMBL 结构条目常含盐形式；对接前将结构去盐并保留最大有机片段，以得到单一可对接分子。

四类标签的定义如下。**严格规则**：dual 为两端 pChEMBL ≥ 6.5；A_only 为 A ≥ 6.5 且 B ≤ 5.5；B_only 为对称定义；neither 为两端 ≤ 5.5；介于 5.5 与 6.5 之间的灰区配体不进入严格面板。**单阈值规则**：选定阈值 θ 后，两端均 ≥ θ 为 dual，一端 ≥ θ 且对端 < θ 为对应单靶选择性类，两端均 < θ 为 neither。

为判断公开数据能否支撑上述严格四类面板，我们在 49 对可审计靶对上统计严格规则下 A_only 与 B_only 的可用配体数。两端均 ≥ 50 的仅 4 对；去掉金属依赖、不适合作为常规对接主对象的 HDAC1/HDAC6 后，剩余 PIK3CA/mTOR、AChE/BChE 与 PIK3CA/PIK3CB。EGFR/HER2 在同一规则下 B_only 仅 7 个，不足以按严格规则建成规模均衡的四类面板。据此将评价集定为上述三对，并另纳入 EGFR/HER2 作为供给受限的案例对（沿用既有统一 RDKit 准备面板）。

各对建成面板时实际采用的标签规则不完全相同：AChE/BChE 与 PIK3CA/PIK3CB 按严格规则配额抽样；EGFR/HER2 与 PIK3CA/mTOR 主面板（PM48）因严格规则下单端选择性配体过少，改用 θ = 6.0 建成。为检查标签定义对结论的影响，我们在 θ ∈ {5.5, 6.0, 6.5} 与严格规则下统一重算口袋匹配 summary_min，作为敏感性分析（详见 Supporting Information）。

### 2.3 评价集面板

各对按类别配额抽样，并限制同类内相近骨架的重复（目标上限约每类 5；AChE/PIK3CB 建造脚本当时以 `chembl_id` 前缀作多样性代理，后续分析改用真 Murcko 支架清点）。面板规模与分析用有效 n 见表 1。对接失败的配体–受体组合从该受体分数中剔除；方向 AUROC 仅使用 dual / A_only / B_only 三类，neither 保留在面板中供描述，但不进入主方向对比。

**Table 1.** K=4 评价集（主报告规模）。

| 靶对 | 角色 | 标签规则 | 受体 (A/B) | 面板规模 | 分析用 n (D/A/B) | Vina exhaustiveness |
|------|------|----------|------------|----------|------------------|---------------------|
| PIK3CA/mTOR (PM48) | 主开发对 | θ = 6.0 | 4L23 / 4JT6 | 48 | 18 / 14 / 12 | 16 |
| AChE/BChE | 主开发对 | strict 6.5/5.5 | 4EY7 / 4BDS | 100（对接后有效约 95） | 27 / 25 / 28 | 8 |
| PIK3CA/PIK3CB | 同工酶对照 | strict 6.5/5.5 | 4L23 / 2WXF | 100（有效约 99） | 28 / 27 / 28 | 8 |
| EGFR/HER2 | 供给受限案例 | θ = 6.0 | 3POZ / 3RCD | 110 | 28 / 38 / 32 | 8（既有面板，无新对接） |

PIK3CA/mTOR 另建扩面面板（下文称 PM110）：保留 PM48 全部 48 个配体，并按严格规则配额追加新分子（目标 dual/A_only/B_only/neither = 30/30/30/25），实际面板 n = 115。PM110 是 PM48 的超集扩样，不是独立复制实验。

### 2.4 受体准备与对接盒

受体选自已发表共晶结构，优先含可对接的小分子配体、分辨率可接受、且口袋定义清晰者。PIK3CA 与 mTOR 分别冻结为 4L23 与 4JT6（共晶配体 PI-103/X6K）；AChE 与 BChE 为 4EY7 与 4BDS；PIK3CB 为 2WXF；EGFR 与 HER2 为 3POZ 与 3RCD。对接盒由共晶配体轴对齐包围盒外扩 5 Å 定义，边长下限 20 Å。受体处理为对接用 PDBQT；盒子坐标写入各面板 `boxes/` 与 `protocol.yaml`。

共晶配体重对接用于协议 sanity。PIK3CA/mTOR 上，mTOR（4JT6）在 exhaustiveness = 8 时 PI-103 重对接偏离较大，升至 16 后回到亚埃级，故该对主面板采用 E = 16；其余对主报告用 E = 8，并在 PIK3CA/mTOR 上另报 E = 8 对照。

### 2.5 配体准备

主协议冻结为 RDKit ETKDGv3 构象生成 + meeko 转 PDBQT。Schrodinger LigPrep 仅在 PM48 上作准备敏感性对照，不与 RDKit 姿态混入主结果表。构象生成、面板抽样与对接均使用固定随机种子；完整种子与参数取值见 Supporting Information Table S1。

### 2.6 对接与重打分

姿态采样使用 AutoDock Vina 1.2.7，默认输出 9 个模式，`energy_range = 3`；exhaustiveness 按表 1，PIK3CA/mTOR 主面板、其扩面面板与单靶 enrichment 用 16，其余主面板用 8。

同一组 Vina 姿态上运行两条重打分通道。RTMScore（`rtmscore_model1`）对 9 个模式取最优分（best-of-9）。GNINA v1.3.2 在 CPU 模式下，将 Vina 的 `mode_01` 经 Open Babel 转为 SDF 后，执行 `--cnn_scoring rescore --minimize`。主报告以 Vina 口袋匹配分为准；RTM 与 GNINA 作通道对照，不另选“优胜臂”改写主张。

### 2.7 评价指标与平凡基线

对每个配体，两端口袋各有对接分数。Vina 亲和力取负号使分数越高越好；RTM 与 GNINA CNN 分数保持越高越好。**口袋匹配方向 AUROC**定义为：dual 对 A_only 使用口袋 B 的分数；dual 对 B_only 使用口袋 A 的分数。两臂 AUROC 的最小值记为 `summary_min`。池化分数（两端平均）仅作对照，不作主指标。

平凡基线用同一套方向 AUROC 流程，但以配体描述符代替对接分：重原子数、分子量、cLogP、TPSA。基线门控报告对接 `summary_min` 相对最优平凡基线的差 Δ 及其 bootstrap 区间。

不确定度用配体层 bootstrap（B = 2000 次重采样），报告 `summary_min` 的 95% 百分位区间。支架层重采样另作对照，主文仍以配体 bootstrap 为准。

### 2.8 混淆与稳健性对照

为检验分数是否主要反映分子属性而非口袋匹配，设置下列对照。（i）错口袋：把方向对比用的口袋对调。（ii）配体效率：分数除以重原子数后再算口袋匹配 AUROC。（iii）匹配子集：在 \|ΔpChEMBL\| ≤ 0.5 或 \|Δheavy atoms\| ≤ 2 的近邻匹配子集上重算单臂 AUROC。（iv）协变量：以对接分为主效应，加入 heavy atoms 与 TPSA 的逻辑回归，比较仅分数与分数+协变量的判别 AUROC。（v）配体二维基线：ECFP4（Morgan 半径 2，2048 bit）逻辑回归；主用按 Murcko 支架的 GroupKFold，随机折仅作泄漏诊断。

单靶 enrichment 在 4L23 与 4JT6 上分别进行：活性分子 pChEMBL ≥ 6.5，decoy 为同靶已测弱效分子（pChEMBL ≤ 5.5），并按分子量 ±50、logP ±1.5、TPSA ±25 作属性匹配；目标规模约 50 活性 + 150 decoy，Vina E = 16。该设定严于随机无关 decoy，用于回答“单靶虚拟筛选信号是否存在”，不替代四类方向评测。

标签阈值、分数聚合方式（口袋匹配 / 错口袋 / 最差口袋 / 池化 / LE）与配体准备（RDKit vs LigPrep，仅 PM48）的敏感性结果写入 Supporting Information；主文只报告影响主张边界的摘要数字。

### 2.9 软件与复现

分析在 Python 3 环境完成。关键版本：RDKit 2026.3.1，meeko 0.7.1，AutoDock Vina 1.2.7，GNINA 1.3.2，RTMScore 使用公开权重 `rtmscore_model1.pth`。面板、分数表、分析脚本、协议文件与完整随机种子/参数表（Table S1）随仓库提供；公开数据包与 DOI 见 Data and Software Availability（Zenodo 发布后填入）。最小复现命令为重算口袋匹配主表、加强分析包与森林图脚本（见仓库 README）。

---

## 审稿人视角自审（不进正文）

| 审稿人可能问 | 本稿怎么处理 |
|--------------|--------------|
| 为什么 EGFR 也在 K=4，却又说供给不够？ | 角色在 **2.2** 定义为供给受限案例对；Results 3.1 只报审计数字 |
| 标签规则为何不统一？ | **2.2 先审计选对，再说明各对建面板时实际用的规则**，并统一重算作敏感性分析 |
| NLRP3/JNK1 为何出现？ | **已删**：正文未做该对实验，Methods 不再提私有 holdout |
| PM110 是不是第二次独立验证？ | **2.3 写明超集扩样**，n=115 |
| 为何 PM 用 E=16、别的用 E=8？ | **2.4 用 cognate 重对接说明**；Results 再报 E8 对照 Δ |
| 有没有 data curation？ | **2.2 整节**对应 2015 editorial |
| enrichment decoy 是不是 DUD-E 式随机物？ | **2.8 写明同靶弱效已测 + 属性匹配窗口**（以代码 ±1.5/±25 为准） |
| 主指标是不是又在挑对对接有利的聚合？ | **2.7 先定义口袋匹配**；池化降为对照 |
| GNINA/RTM 是不是事后选优胜？ | **2.6 写明通道对照，不改主张** |
| 多样性过滤用了 chembl_id 前缀？ | **2.3 诚实写入**；后续用真 Murcko 清点 |
| 随机种子数值有什么意义、为什么反复出现？ | **已改**：具体种子数值（如 20260727/20260729）不再逐节重复出现在正文，只在 2.5 提一次“固定种子，完整取值见 SI Table S1”；对结果解读有意义的参数（n_modes=9、energy_range=3、E=8/16、bootstrap B=2000）仍留在正文，仿 Vu et al. 2025 的处理方式 |

---

## 与 Results 3.1 的分工（已同步改 Results）

- **Methods 2.2–2.3**：谁进 K=4、EGFR 为何是案例、标签规则。  
- **Results 3.1**：只报告 49 对审计的供给数字与“厚面板稀缺”这一发现；EGFR 的 7 个 B_only 作为审计反例出现，不宣布实验决策。
