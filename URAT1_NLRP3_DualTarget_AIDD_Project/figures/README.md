# 论文图表（Arial 8 pt · SciencePlots）

由 `scripts/plot_available_figures.py` 生成，遵循 [SciencePlots](https://github.com/garrettj403/SciencePlots) 出版规范。

## 样式规则（v2）

| 规则 | 说明 |
|------|------|
| 字体 | **Arial 8 pt** 全文统一 |
| 靶点标题 | 图顶居中一行：`Target: URAT1 …` 或 `Target: NLRP3 …`（不进绘图区） |
| 坐标轴 | 横纵轴均有完整英文标签 + `labelpad` |
| 统计量 | 放入 **legend**，不压在曲线/柱上 |
| 柱顶数值 | `ax.bar_label()`，自动留 headroom |
| 面板标记 | `a`–`d` 在轴外左上角，不与数据重叠 |

## 靶点配色

| 靶点 | 主色 | 辅色 |
|------|------|------|
| URAT1 | `#0072B2` | `#56B4E9` |
| NLRP3 | `#D55E00` | — |

## 目录

| 子目录 | 内容 |
|--------|------|
| `generated/main/` | Fig 2（NLRP3）、Fig 3（URAT1）组合图 |
| `generated/nlrp3/` | NLRP3 单图 |
| `generated/urat1/` | URAT1 单图 |
| `generated/si/` | 补充图 |

## 重新生成

```bash
pip install SciencePlots
python3 scripts/plot_available_figures.py
```
