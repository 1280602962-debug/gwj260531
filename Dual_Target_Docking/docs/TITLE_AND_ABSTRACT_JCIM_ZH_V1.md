# 题名与摘要（JCIM Articles 中文工作稿）

## 题名

**基于对接的双靶识别：三靶对的评价设定审计**

## 摘要

对接常被用来解释双靶识别，但负类通常是非结合配体或未匹配 decoy，而不是实验测定的单靶选择性配体。本文将 DualFourClass 冻结为该选择的四状态、三靶对配方审计；第四个候选靶对 PIK3CA/PIK3CB 在事后受体身份审计中被撤回——其对接所用的"PIK3CB"受体（PDB 2WXF）实为小鼠 PIK3CD，并非人源 PIK3CB，因此仅作为已记录的受体身份错误在 Supporting Information 中报告。在 AutoDock Vina 下，实验定义的双靶配体相对 A-only 与 B-only 选择性配体排序。最弱臂 AUROC 分别为 EGFR/HER2 0.430、AChE/BChE 0.606 与 PIK3CA/mTOR 0.692；三条 95% 置信区间均包含 0.5。同一套 EGFR/HER2 分数在双靶对 neither 时 AUROC 为 0.756，定向评价则降至 0.430。独立 GNINA 姿态生成再现了该配方差距（0.783 对 0.220）。支架分组的纯配体模型已捕获大部分表观排序，替换备选 PIK3CA 晶体结构则使 PIK3CA/mTOR 对比降至接近随机水平。双口袋分数过滤因此需要选择性感知的负类；两个口袋同时给出有利分数，本身并不是双靶识别的证据。

## 关键词

双靶对接；选择性；硬负例；AutoDock Vina；GNINA；虚拟筛选
