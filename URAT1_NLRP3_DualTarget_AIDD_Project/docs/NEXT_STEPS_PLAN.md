# 接下来做什么（面向 Mol Div 投稿的执行计划）

协议筛选已完成，Π* = **P2（gnina CNNaffinity）**（见 `PROTOCOL_SELECTION_RESULT.md`）。
剩余工作是把"协议选择"补成一篇完整的、假说级的双靶重定位计算文。

## 阶段总览

| 阶段 | 产出 | 脚本/工具 | 状态 |
|------|------|-----------|------|
| 0 | 归档协议筛选结果表 | 本地 EF 计算输出 → 提交 | 你已算，待提交 |
| 1 | 1588 双靶对接（P2，9DKB+7ALV） | `run_funnel_p2.sh` | 待跑 |
| 2 | 合并 + Pareto 短名单 | `merge_docking_pareto.py`（脚本已内含调用） | 随阶段 1 |
| 3 | 短名单 2–4 个 lead：MD + 相互作用 (+MM-GBSA) | GROMACS/AMBER（本地） | 待跑 |
| 4 | ADMET / PAINS / 化学空间审计 | `10_admet_druglikeness.py`, `11_chemical_space_novelty.py` | 待跑 |
| 5 | Results/Methods 写作 + 严格降调 claim | docs | 待写 |

## 阶段 1–2：双靶对接 + Pareto（先做，算力主项）

- 输入：`data/repurposing/screening/docking_pool_p05.csv`（~1588）
- 引擎：gnina，`score_mode=cnnaff`，`cnn_scoring=rescore`，exhaustiveness=32
- 只对 1588 各对接一次（9DKB、7ALV），不重复
- 产出：`results/repurposing/pareto_shortlist.csv`、`pareto_merged_scores.csv`

CPU-only 粗估：1588 × 2 靶，gnina exh=32 单分子约 20–60s（视核数），建议多核/夜间跑。

## 阶段 3：代表 lead 的 MD（Mol Div 观感关键）

- 从 `pareto_shortlist.csv` 选 2–4 个（优先：双轴百分位高、非 PAINS、非已知对照药）
- 每个 lead：9DKB 和/或 7ALV 复合物，100 ns 级 MD；报 RMSD/RMSF、关键氢键/盐桥
- 起始构象用晶体或 RTMScore 选姿，不用失败的 top-1
- 可选 MM-GBSA（仅同协议内相对比较，勿当绝对亲和力）
- 对照：lesinurad@9DKB、MCC950/类似物@7ALV 各跑一条作基线

## 阶段 4：审计（已有脚本）

```
python3 scripts/10_admet_druglikeness.py      # 类药性 / ADMET
python3 scripts/11_chemical_space_novelty.py  # 与已知配体的新颖性 / Tc
```
短名单标注 PAINS/Brenk、Lipinski/Veber、与已知 URAT1/NLRP3 配体的最近邻相似度。

## 阶段 5：写作要点（对齐 Mol Div 无实验虚筛文）

- 创新句放前面：**先在双诱饵上选协议**、**不对称双证据漏斗**、**Pareto ≠ 提名（审计降级）**
- Results：协议表（据实，含 P5 Random 失败）→ 漏斗计数 → 短名单 → 2–4 lead 的 MD/相互作用主图 → ADMET
- 主张固定：`computational dual-node repurposing hypotheses, pending experimental validation`
- 明确未来实验：URAT1 尿酸摄取抑制、NLRP3 IL-1β 释放
- Limitations：中等富集、单构象对接、诱饵为库分子而非实验无活、MM-GBSA 仅相对

## 不要做

- 不重复 lesinurad 烟雾（已完成）
- 不对 True/Random 两份 CSV 各对接一遍（协议筛选已完成）
- 不写 "identified dual inhibitors"
- 不把 P5 当生产协议
