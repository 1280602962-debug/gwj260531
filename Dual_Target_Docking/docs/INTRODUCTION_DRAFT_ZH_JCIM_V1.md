# Introduction（中文工作稿 · JCIM Articles）

> 结构按 [`POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md`](POSITIONING_AND_FRAMEWORK_LANGUAGE_V1.md)：叙事缺口，不包装成新算法。投稿以英文为准。  
> 可引用：Vu et al. 2025（大规模对接评测）；Hall 等虚拟筛选规模；DualDiff / FuseDiff（生成式双靶，评测口径与本基准互补）。

---

双靶（或多靶）小分子在耐药与通路代偿场景中被反复提出，对接也常被用来给候选分子两端同时打分。但“两端分数都不错”并不等于“能把真正的双靶配体与单端选择性硬负区分开”。若两条口袋的分数被平均或只与参考配体的对接分比较，较弱一端的失败会被掩盖；公开活性库里又很少有化合物在两个靶上都被定量测到、并拉开严格选择性间隔。因此，双靶对接需要的不是再一个打分函数，而是一套带硬负、带混淆对照、声明可靠边界的评测协议。本文建立 systematic benchmarking framework，并释放 DualFourClass-Bench 资源，用来评价对接分数在该任务上的可靠性与局限——不是提出名为新算法的对接框架。

公开数据使可平衡的四类面板很少。我们在 ChEMBL 上按严格规则（两端 ≥ 6.5 为 dual；活性端 ≥ 6.5 且对端 ≤ 5.5 为选择性）审计候选靶对：两端硬负均 ≥ 50 的非金属对只剩三对可作厚面板，EGFR/HER2 作为供给受限案例进入 K = 4。随后的 BindingDB/PubChem 计数核对表明，该厚面板门槛在等式测定下不随数据源翻转（Supporting Information Table S12）。主实验因此合理冻结在 ChEMBL 面板上，而不是声称已经穷尽一切公开库。

评测指标必须与任务同构。本文以口袋匹配方向 AUROC 为主：dual 对 A_only 用口袋 B 的分数，dual 对 B_only 用口袋 A 的分数，summary_min 取较弱一臂；池化、错口袋与平凡描述符并列报告。表观双靶信号经常被分子量、极性或化学型解释，必须把这些对照写进主结果，而不能只报对接 AUROC。

近两年的双靶**生成**工作把同一缺口推到了另一条线上。Zhou 等（NeurIPS 2024）的 DualDiff 与 Wu 等（2026）的 FuseDiff 用 AutoDock Vina 在两个口袋上报告 P-1 / P-2 Vina Dock、Max Vina Dock（两端中较弱的对接分）以及 Dual High Affinity（生成分子在两端都优于参考配体的对接分）。这是相对参考配体的对接双成功，**不是**相对实验标记的 A_only / B_only 硬负。Dual High Affinity 要求两端都过线，并不等于本文所批评的均值池化；但“两端 Vina 都优于参考”仍然把成功定义在对接分上，而本文显示硬负选择性配体也常常在两端给出看似合理的对接分（holdout 上错口袋对照不低于口袋匹配）。因此 DualFourClass-Bench 的四类硬负与口袋匹配指标，可被生成式方法用来做更诚实的下游评测：检查被标为 Dual High Affinity 的分子，在本基准上是更像 dual，还是更像硬负。本文**没有**重跑 DualDiff/FuseDiff 的生成分子；这是基准用途上的互补，不是对那些方法的实证打分。

下文报告 K = 4 冻结集上的供给、主指标、混淆对照、holdout 与结构稳健性。主张强度以实验结果为上限。
