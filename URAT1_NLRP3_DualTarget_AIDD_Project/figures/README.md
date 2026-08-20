# 论文图表（Arial 8 pt）

由 `scripts/plot_available_figures.py` 生成。生产协议 **gnina P2**。  
8973 回顾图与历史商业对接 Pareto 图已删除。P2 双靶百分位图为 `generated/main/fig03_p2_dual_percentiles.*`（n=1,580）。

## 当前保留

| 路径 | 内容 |
|------|------|
| `generated/main/` | Fig 2：NLRP3 ML 缩库；Fig 3：P2 双靶百分位 |
| `generated/nlrp3/` | NLRP3 单图 |
| `generated/si/` | 数据不对称、临床库分期、OOF（若有） |
| `ppt_assets/schematics/` | 痛风双轴、gnina P2 流程、与 Eurycoma 对照 |
| `ppt_assets/structures/` | 生产口袋 9DKB、7ALV（及 8ETR 文献对照） |

## 重新生成

```bash
python3 scripts/plot_available_figures.py
python3 scripts/build_ppt_assets.py
```
