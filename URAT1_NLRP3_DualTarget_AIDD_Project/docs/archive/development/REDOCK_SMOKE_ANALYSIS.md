# lesinurad × 9DKB 重对接烟雾：结论与下一步

数据：`data/redock_smoke/redock_results_lesinurad_9DKB.csv`（2026-07-24）  
设置：同一搜索盒；`num_modes=9`；gnina `--cnn_scoring rescore --no_gpu`；RTMScore model1，口袋 10 Å。

## 1. 正式门控（exhaustiveness = 32）一览

| 协议 | 读出 | Top-1 RMSD (Å) | Best RMSD (Å) | 采样 ≤2 Å | Top-1 ≤2 Å |
|------|------|----------------|---------------|-----------|------------|
| **P1** | Vina affinity | **0.86** | **0.86** | ✅ | ✅ |
| **P2** | CNNaffinity | 4.16 | 0.99 | ✅ | ❌ |
| **P3** | gnina affinity | 5.69 | 0.99 | ✅ | ❌ |
| **P0** | CNNscore（负对照） | 4.16 | 0.99 | ✅ | ❌ |
| **P4** | Vina → RTMScore | **0.86** | **0.86** | ✅ | ✅ |
| **P5** | gnina → RTMScore | **0.99** | **0.99** | ✅ | ✅ |

## 2. 可下的结论

1. **流水线有效：搜索盒与配体/受体准备可用。** 全部协议的 Best-in-ensemble RMSD ≤ 1.0 Å，近晶体姿稳定出现在 9 构象集合中。
2. **可宣称 pose-accurate 的读出（exh=32）：P1、P4、P5。** Top-1 均 ≤ 1 Å。
3. **CNN 读出（P0、P2）不能作构象证据。** 集合里有近原生姿，但 Top-1 落在 ~4 Å；符合“采样成功、排序失败”。P0 作为负对照行为符合预期。
4. **gnina affinity（P3）不稳定。** exh=8 时 Top-1≈0.81 Å 通过；exh=32 时 Top-1≈5.69 Å 失败。不宜单独作为结构用姿读出；仍可进 TrueDecoy 富集比较。
5. **RTMScore 在 exh=32 下能从 Vina/gnina 集合中捞回近原生姿（P4/P5）。** exh=8 时 P4/P5 未过 Top-1 门控（2.6 / 4.6 Å），正式实验应固定 **exhaustiveness=32**。
6. **烟雾测试不等于协议锁定。** 本表只回答“能不能对接 / 谁能排对姿”；Π\* 仍须在 TrueDecoy（主）+ RandomDecoy（否决）上按 EF@1% 选定。

## 3. 对 Methods / 写作的即时用法

- 自对接三类 RMSD 数字可直接引用本 CSV（exh=32 行）。
- 正文表述建议：P1/P4/P5 可通过 Top-1≤2 Å 门控；P0/P2/P3 限制 pose 主张，但仍参与富集协议筛选。
- 生产搜索深度写死：**exhaustiveness=32，num_modes=9**（与本烟雾正式门控一致）。

## 4. 接下来做什么（按优先级）

| 优先级 | 任务 | 目的 |
|--------|------|------|
| **1** | 在 TrueDecoy + RandomDecoy 上跑 **P0–P5（exh=32）** | 按预定规则锁定 Π\* |
| **2** | 用锁定的 Π\* 对接临床池（1588）@ 9DKB + 7ALV | 进入不对称漏斗 / Pareto |
| **3** | 可选：benzbromarone、dotinurad 同盒交叉对接（SI） | 扩展对照，非门控 |
| **4** | MCC950 @ 7ALV 类似物对照（已完成，P2） | 药理学阳性对照，**不是**自对接；共晶仍是 NP3-146 |

MCC950 P2 结果：CNNaffinity **7.018**，gnina affinity −10.20 kcal/mol，CNNscore 0.9013（`data/redock_smoke/redock_results_mcc950_7alv.csv`）。相对 RM5 的 RMSD **不作**自对接门控。

**不必再重复 lesinurad@9DKB 烟雾**；门控数据已齐。下一步算力应投向基准库协议筛选。
