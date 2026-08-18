# URAT1 蒸馏 / 回顾集（`data/distill/`）

**用途：** 仅 URAT1 回顾（可选 SI）。**不是**临床库主筛选，也不是 TrueDecoy/RandomDecoy 协议筛选池。

协议筛选诱饵见 `data/benchmarks/urat1_true_decoy/`。临床库见 `data/repurposing/`。

| 子集 | 文件 | 说明 |
|------|------|------|
| A | `distill_subset_a.csv` | URAT1 训练活性 |
| B | `distill_subset_b.csv` | Murcko 骨架代表 |
| C | `distill_subset_c.csv` | SLC22 邻域（URAT1 + OAT1/OAT3） |
| D | `distill_subset_d.csv` | 多样性负样本（**不是**现行协议筛选诱饵） |
| E | `distill_subset_e.csv` | Benchmark 参考 |

合并表：`distill_manifest.csv`。

```bash
python3 scripts/00b_build_distill_set.py
```

统计见 `distill_set_summary.json`。本目录 **不做** Teacher/MASFL 三态 Glide 对接。
