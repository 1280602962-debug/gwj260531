# MASFL v3.1 蒸馏集 (`data/distill/`)

## 子集定义

| 子集 | 文件 | 说明 |
|------|------|------|
| A | `distill_subset_a.csv` | URAT1 训练活性（pIC50） |
| B | `distill_subset_b.csv` | Murcko 骨架代表（每骨架 1 条） |
| C | `distill_subset_c.csv` | SLC22 邻域 FPS（URAT1 + OAT1/OAT3） |
| D | `distill_subset_d.csv` | 百万库随机无活性多样性负样本（**MASFL/蒸馏用**；**不是**现行 URAT1 TrueDecoy/RandomDecoy 协议筛选池——后者见 `data/benchmarks/urat1_true_decoy/`） |
| E | `distill_subset_e.csv` | Benchmark 参考 + 类似物边界集 |

合并去重表：`distill_manifest.csv`（按 A > E > B > C > D 优先级保留主标签）

## 构建命令

```bash
# 子集 D 需先抽样（或从外部拷贝 distill_subset_d.csv）
python3 scripts/sample_distill_subset_d.py --library /path/to/1M.csv --n 8000 --diversity

# 生成 A/B/C/E 并与 D 合并
python3 scripts/00b_build_distill_set.py
```

## 统计

见 `distill_set_summary.json`。ChEMBL SLC22 导出规模有限，子集 C 可能低于文档目标 2000–5000。

## 下一步（MASFL Stage 2）

对 `distill_manifest.csv` 全量运行 Teacher M-CPDL 三态 Glide 对接（不对百万库跑）。
