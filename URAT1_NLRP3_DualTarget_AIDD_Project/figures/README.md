# 论文图表（Arial 8 pt）

由 `scripts/plot_available_figures.py` 从现有数据生成。

## 靶点配色（全文统一）

| 靶点 | 颜色 | 用途 |
|------|------|------|
| **URAT1** | `#0072B2` 蓝 | 8973 回顾、OOF parity |
| **NLRP3** | `#D55E00` 橙 | 8319 ML 筛库、ROC/PR |

每张图左上角有 **靶点横幅**（`URAT1 (SLC22A12)` 或 `NLRP3`）。

## 目录

| 子目录 | 内容 |
|--------|------|
| `generated/main/` | 主文组合图 Fig 2（NLRP3）、Fig 3（URAT1） |
| `generated/nlrp3/` | NLRP3 单图 |
| `generated/urat1/` | URAT1 单图 |
| `generated/si/` | 补充图（ROC/PR、parity、数据不对称） |

## 重新生成

```bash
python3 scripts/plot_available_figures.py
```

输出清单：`generated/figure_manifest.json`

## 待对接/MD 后补充

- Fig 4 Pareto（双靶对接）
- Fig 5–6 结构图与 MD（需 Maestro/GROMACS）
