# 题名与摘要（JCIM Articles 中文工作稿）

## 题名

**基于对接的双靶识别：四靶对的评价设定审计**

## 摘要

评价集如何构建，决定了两个口袋中的有利分数能够为双靶活性提供多强的证据。我们将四个冻结靶对的 ChEMBL 配体划分为 dual、A-only、B-only 与 neither，并在两个口袋匹配方向评价 AutoDock Vina：dual 对 A-only 使用口袋 B，dual 对 B-only 使用口袋 A；较小的方向 AUROC（`summary_min`）作为保守的描述性摘要。EGFR/HER2 的 Dual versus neither AUROC 为 0.756，而 `summary_min` 为 0.430；独立 GNINA 姿态生成保留了这一差距（0.783 对 0.220），其余靶对未出现同样的设定效应。在 5 个冻结 Vina 种子中，四个靶对的 Dual-versus-neither 减 `summary_min` 差值符号均未改变；EGFR/HER2、AChE/BChE、PIK3CA/PIK3CB 与 PIK3CA/mTOR 的 `summary_min` 中位数分别为 0.373、0.599、0.478 与 0.704。将 docking 加入支架分组 ECFP4 模型后，交叉验证 AUROC 的最大绝对变化为 0.020；替代 PIK3CA 受体在两个相关靶对上引起方向相反的位移；四个主种子 `summary_min` 的 95% CI 均包含 0.5。预先冻结的 BindingDB 原生门槛没有产生可进入外部对接的靶对。因此，双口袋对接的解释需要选择性硬负、配体层对照以及受体和随机种子敏感性；当前四靶对面板不足以建立靶标通用性能。

## 关键词

双靶对接；基准设定；选择性硬负样本；化学混淆；受体实现；虚拟筛选
