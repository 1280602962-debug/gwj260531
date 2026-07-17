# 漏斗一句话 + 计数

对 **527,779** 唯一胺跑通 ARS×PaperSpine **L0–L7** 对接前预筛选（chemotype triage，非活性预测）。

| 阶段 | 结果 |
|------|------|
| L0 | 527,779 胺入库体检 |
| L1 | **345,626**（胺 MW 120–450；sites=1 优先；sites=2 且 MW≤350） |
| L2–L4 | 标准丙烯酰胺枚举 + QC + 特征 → **239,990** 产物（ok 80,174 / watch 159,816；bad 硬删） |
| L5 | core↔full Tc Spearman ρ=**0.66**；Sim 窗 Jaccard(core vs full)=**0.27**（去弹头有影响）；ErG Novel 门=p75=**0.73** |
| L6 | Sim **8,063** / Novel **103,187** / Pan **9,204** / discard 余 |
| L7 dock_ready | **10,000**：sim_yl 4200 + sim_56d 1800 + novel 3500 + pan 500（仓间 ID 交集=0） |

主交付目录：`L7/L7_dock_ready_*.csv`  
阈值：`config/thresholds.json`  
校准：`L5_calibration.md`

| **L7b（收紧 watch + Novel keep=yes）** | **8,543**：sim_yl 4198 + sim_56d 1800 + novel 2245 + pan 300 |

主交付：`L7b/L7b_dock_ready_*.csv`；交接：`handoff_glide_af3/`

