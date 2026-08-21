# Discussion（中文工作稿 · JCIM Articles）

> 与 [`DISCUSSION_LIMITATIONS_DRAFT_ZH_JCIM_V1.md`](DISCUSSION_LIMITATIONS_DRAFT_ZH_JCIM_V1.md) 配套：本文件写解释与外推边界；局限清单仍用 Limitations 稿。  
> 不把开放问题写成已解决。

---

本文的核心不是“对接能不能识别双靶配体”的全称判断，而是：在公开可构建的四类硬负面板上，对接分数何时还保留方向信号、何时被分子属性带走。K = 4 是供给审计冻结的评价集，不是对全部双靶任务的抽样。ChEMBL 作为对接主数据源是合理的：BindingDB/PubChem 计数核对未翻转厚面板门槛（Table S12），但 EGFR/HER2 在等式测定下仍达不到 ≥ 50，不能改写成厚面板。

口袋匹配指标把弱臂放到明处。池化可以让 EGFR/HER2 看起来接近随机，而分臂后 B 端失败是清楚的。这与 DualDiff / FuseDiff 一类生成式评测形成对照：那些工作用两端 Vina 相对参考配体的 Dual High Affinity 定义成功，并不检验实验硬负。本基准可以承接那类分子做下游诚实评测——要求报告口袋匹配 summary_min 相对 A_only / B_only，而不是只报两端 Dock 分都优于参考。我们未对 DualDiff 或 FuseDiff 的生成物重新对接，故这是用途陈述，不是方法对比表。

Holdout 上三对均出现错口袋对照不低于口袋匹配，与主面板（口袋匹配更高）方向相反。效价/尺寸匹配诊断（Table S13）表明：unused-pool 相对主面板确有抽样偏移（PIK3CA/mTOR holdout 更弱），但匹配后悖论仍在，故不能把反差归因于“holdout 抽偏了”。接触计数证明 B 臂存在真实的尺寸/埋藏混淆，但不能按幅度复现 Vina 错口袋。该反差保持开放，审稿时不应被包装成机制已阐明。

PIK3CA/mTOR 是唯一在点估计上同时高于随机与最强描述符基线的靶对，但其 Δ 的置信区间仍含 0，且换晶体后受体依赖。文章能主张的是：对接增量高度对靶、对结构；评测协议必须同时报告硬负、分臂与混淆对照。GNINA 仍仅为 mode-1 rescore，与 RTM 的 best-of-9 不对称（本环境无法补 9 姿态，见 `GNINA_NINE_POSE_SKIP_V1.md`），三引擎公平赛马不得写入。

局限见 Limitations 稿。
