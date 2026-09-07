# 题名与摘要（JCIM Articles 中文工作稿）

## 题名

**基于对接的双靶识别：多靶对的评价设定审计**

## 摘要

对接常被用来解释双靶识别，但负类通常是非结合配体或未匹配 decoy，而不是实验测定的单靶选择性配体。本文将 DualFourClass 建成该选择的四状态配方审计。最初的主评价集是 2026-07-23 的四靶对冻结。事后 ChEMBL 37 宇宙普查之后，五个通过预先声明门槛的普通非共价靶对被加入同一分析栈——这是收集补全，不是预先注册的八靶对冻结。原先对接的第四个候选靶对 PIK3CA/PIK3CB 在受体身份审计中被撤回——其"PIK3CB"受体（PDB 2WXF）实为小鼠 PIK3CD，并非人源 PIK3CB，因此仅作为已记录的受体身份错误在 Supporting Information 中报告。主表因此为八行：PIK3CA/mTOR、AChE/BChE、供给受限的 EGFR/HER2，以及 F2/F10、JAK1/TYK2、JAK1/JAK2、PPARG/PPARA、PPARA/PPARD。在 AutoDock Vina 下，实验定义的双靶配体相对 A-only 与 B-only 选择性配体排序。最弱臂 AUROC 分别为 EGFR/HER2 0.430、AChE/BChE 0.606、PIK3CA/mTOR 0.692、F2/F10 0.345、JAK1/TYK2 0.365、JAK1/JAK2 0.588、PPARG/PPARA 0.649 与 PPARA/PPARD 0.446。仅 PPARG/PPARA 的 95% 置信区间完全高于 0.5，且该区间对面板成员与打分公式敏感。同一套 EGFR/HER2 分数在双靶对 neither 时 AUROC 为 0.756，定向评价则降至 0.430；JAK1/TYK2 再现同一配方差距（0.770 对 0.365）。独立 GNINA 姿态生成再现了这两处差距。支架分组的纯配体模型已捕获大部分表观排序，替换备选 PIK3CA 晶体结构则使 PIK3CA/mTOR 对比降至接近随机水平。双口袋分数过滤因此需要选择性感知的负类；两个口袋同时给出有利分数，本身并不是双靶识别的证据。

## 关键词

双靶对接；选择性；硬负例；AutoDock Vina；GNINA；虚拟筛选
