# PRIMARY_METRIC_V2 — 口袋匹配方向 AUROC 升格

> 日期：2026-07-29  
> 依据：`data/jcim_bench_v0/tables/pocket_matched_directional_v1.csv`  
> 复现：`python3 data/jcim_bench_v0/scripts/build_pocket_matched_diagnostics_v1.py`

## 定义变更

| 对比 | 旧主指标 (v1) | **新主指标 (v2)** |
|------|---------------|-------------------|
| dual vs A_only | 池化 `vina_mean` | **口袋 B 分数**（该对比标签差仅在 B 端） |
| dual vs B_only | 池化 `vina_mean` | **口袋 A 分数** |
| 汇总 | `summary_min = min(AUROC_DA, AUROC_DB)` | 同左，但各臂用匹配口袋 |
| 对照 | — | 池化 / 错口袋 / worst / LE / heavy 并列报告 |

**Claim ceiling 同步更新**：主文 primary = pocket-matched directional AUROC；池化降为 Methods 对照。

## K=4 主结果表（Vina · pocket-matched ± bootstrap CI, B=2000, seed=20260729）

| 靶对 | n(D/A/B) | AUROC D vs A | AUROC D vs B | **summary_min** | **95% CI** | 错口袋 min | LE-PM min | heavy min |
|------|----------|--------------|--------------|-----------------|------------|------------|-----------|-----------|
| EGFR/HER2 | 28/38/32 | 0.666 | 0.430 | **0.430** | [0.281, 0.576] | 0.260 | 0.311 | 0.369 |
| AChE/BChE | 27/25/28 | 0.650 | 0.606 | **0.606** | [0.442, 0.737] | 0.444 | 0.413 | 0.582 |
| PIK3CA/PIK3CB | 28/27/28 | 0.691 | 0.500 | **0.500** | [0.340, 0.648] | 0.349 | 0.332 | 0.622 |
| PIK3CA/mTOR | 18/14/12 | 0.714 | 0.692 | **0.692** | [0.457, 0.813] | 0.602 | 0.657 | 0.463 |

## 与旧池化主表对比（结论变化）

| 靶对 | 池化 summary_min | 口袋匹配 summary_min | Δ | 解读 |
|------|------------------|----------------------|---|------|
| EGFR/HER2 | 0.311 | 0.430 | +0.12 | 池化**低估** B 端信号；仍 <0.5 基线门未过 |
| AChE/BChE | 0.530 | 0.606 | +0.08 | 略升；TPSA 捷径仍可见（错口袋 0.44） |
| PIK3CA/PIK3CB | 0.412 | 0.500 | +0.09 | 升至随机；B 端仍弱 |
| PIK3CA/mTOR | 0.671 | 0.692 | +0.02 | **唯一稳健对**；E=16；错口袋仍偏高 (0.60) |

## 特异性缺口（pocket − wrong-pocket）

见 `pocket_specificity_gap_v1.csv`：四对 specificity_gap 均 >0.09，说明错口袋对照偏离 0.5 → 分子层混淆（尺寸/效价/化学型）仍显著。仅 PM 在 LE 归一后仍 > heavy 基线 (`survives_le_normalisation=True`)。

## 尺寸分层（underpowered 警告）

`pocket_matched_size_strata_v1.csv`：多数 stratum 每臂 n<8，主文仅作 SI；扩面 PM110 旨在收窄 CI。

## 主文建议措辞

- Primary endpoint: pocket-matched directional AUROC with bootstrap CI.
- Pooled vina_mean retained as **legacy control** (Table Sx).
- Do not claim cross-pair generalization; PM is exploratory positive control only.
