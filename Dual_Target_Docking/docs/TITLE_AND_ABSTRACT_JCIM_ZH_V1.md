# 题名与摘要（JCIM Articles 中文工作稿）

## 题名

**对接双靶识别中的基准设定与化学混淆**

## 摘要

两端有利的对接分数是否构成双靶识别证据，尚未在实验定义的单靶选择性配体上得到充分检验。我们构建 DualFourClass-Bench，这是一套经策展的四对、四状态面板，含两条方向性主任务：dual-active 对 A-only 选择性配体在口袋 B 上打分，以及 dual-active 对 B-only 选择性配体在口袋 A 上打分；两臂中较弱者作为靶对汇总（`summary_min`）。在同一套冻结 AutoDock Vina 分数上，EGFR/HER2 的 Dual-versus-neither 对照给出明显更强的表观结果（AUROC 0.756；方向性 `summary_min` 0.430），而混合库 Top-10 中有 9 个实验选择性配体。AChE/BChE 与 PIK3CA/PIK3CB 只显示小且区间重叠的 formulation 增量；PIK3CA/mTOR 的 Dual-versus-neither 对照因 neither n = 4 而效能不足。在支架分组交叉验证下，把 docking 加入 ECFP4 后 AUROC 的绝对变化至多约 0.02。将最大 pChEMBL 换成重复测定的中位数只产生很小的靶对层变化。相比之下，替代受体实现使 PIK3CA/mTOR 从 0.692 变为 0.486/0.505，却使 PIK3CA/PIK3CB 从 0.500 变为 0.691/0.685。因此，表观双靶判别同时依赖基准设定、靶对、配体化学和受体实现；在本支架感知评价中，对接相对配体层化学基线只提供有限增量信息。这些结果支持采用实验定义的选择性硬负样本与混淆感知对照，但不构成普遍高估定律，也不证明 docking 缺乏口袋特异信息。

## 关键词

双靶对接；基准设定；选择性硬负样本；化学混淆；受体实现；虚拟筛选
