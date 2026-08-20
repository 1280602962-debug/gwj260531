# 补充分析（不重对接临床池、不重锁 Π*）

复现：`python3 scripts/archive_p2_export.py`（漏斗归档与协议 bootstrap）；`python3 scripts/si_supplement_analyses.py --skip-dock`（assay 重叠）。  
生产对接池仍为 `data/repurposing/screening/docking_pool_p05.csv`（n=1588）；百分位表为 `data/repurposing/p2/pareto_merged_scores.csv`（n=1580）。

---

## 0. MCC950 要不要在 7ALV 上对接？

**要，作为药理学类似物对照；不是自对接。**

7ALV 共晶配体是 NP3-146（残基 RM5），不是 MCC950（Dekker 等，2021）。MCC950 是 NLRP3 工具抑制剂阳性对照，不在本临床库 8319 中（无 CHEMBL230208 / CHEMBL3183703 / 该规范 SMILES），因此生产 1588 池本来就不会给它排名。把它单独按生产 P2 接到 7ALV，只回答“口袋与配体准备能否给阳性对照一个可用姿”，**不重新选择协议**，也**不能**用相对 RM5 的 RMSD 当自对接门控。真正的结构自对接仍应是 NP3-146/RM5 晶体坐标重对接；本次未做。

ChEMBL 文献常用 ID `CHEMBL230208` 不在 NLRP3 清洗集；同一规范 SMILES 以 `CHEMBL3183703` 出现在训练集（15 条 IL-1 相关记录，pActivity 5.54–8.30）。因此 MCC950 的 **ML 阳性不是独立外推**；7ALV 对接才是与训练标签分开的结构对照。

本次 P2（gnina CNNaffinity，exh=32，`num_modes=1`，`--cnn_scoring rescore --no_gpu`）结果：

| 读出 | 值 | 方向 |
|------|----|------|
| CNNaffinity | **7.018** | 高优；生产 `dock_score` 存为 −7.018 |
| gnina affinity | −10.20 kcal/mol | 低优 |
| CNNscore | 0.9013 | 高优 |

文件：`data/si/mcc950_7alv/`，汇总 `data/redock_smoke/redock_results_mcc950_7alv.csv`。

---

## 1. Assay top-1 / top-3 / top-5 缩库重叠（不重对接）

对已归档的 8319 条临床库分数，用同一 `nlrp3_model.joblib` 在出现频率最高的 1/3/5 个测定条件下重打分，取 max，阈值仍为 \(q_N\ge0.5\)。

| n assays | n (\(q_N\ge0.5\)) | 与生产 1588 的 Jaccard | 仅生产有 | 仅该集合有 |
|----------|-------------------|------------------------|----------|------------|
| 1 | 1587 | 0.9994 | 1 | 0 |
| 3 | 1587 | 0.9994 | 1 | 0 |
| 5 | 1587 | 0.9994 | 1 | 0 |

top-1、top-3、top-5 **彼此 Jaccard = 1.0**（同一 1587 个分子）。与冻结生产池的唯一差异是 **CT-1578（CHEMBL2035185）**：归档 \(q_N\ge0.5\)，当前 XGBoost 反序列化后略低于 0.5。这是软件版本漂移，不是测定个数改变缩库。

**结论：** 1588 集合没有大变。**不替换生产池，不重对接。** 测定上下文取 5 个的生产设定予以保留；SI 表明缩库几乎由频率最高的单一测定（CHEMBL5549264）决定。

---

## 2. 协议表 EF / AUC 区间（ranking bootstrap）

P0–P5 分子级分数现已归档于 `data/benchmarks/protocol_selection/mol_protocol_scores.csv`（9,839 行）。下表在 True / Random 基准上对有分数的分子做有放回重采样（1,000 次），每次重算 EF 与 AUC，取 2.5–97.5 百分位。**不重新选择 Π\***。前 1% 宽度为 \(\lfloor 0.01N\rfloor\)，故 hits 分母多为 51 而非原文 52。

EF@1%（与锁定点估计一致至舍入）：

| 协议 | True hits@1% | True EF@1% (95% CI) | Random hits@1% | Random EF@1% (95% CI) | True AUC (95% CI) |
|------|--------------|---------------------|----------------|-----------------------|-------------------|
| P5 | 13/51 | 2.80 (1.51–4.36) | 0/51 | 0.00 (0.00–0.00) | 0.590 (0.563–0.617) |
| **P2** | **12/51** | **2.59 (1.31–4.07)** | **1/51** | **0.22 (0.00–1.04)** | **0.580 (0.548–0.609)** |
| P0 | 9/51 | 1.94 (0.89–3.33) | 9/51 | 1.94 (0.92–3.23) | 0.647 (0.619–0.671) |
| P4 | 3/50 | 0.65 (0.00–1.33) | 0/44 | 0.00 (0.00–0.00) | 0.625 (0.603–0.648) |
| P1 | 2/51 | 0.43 (0.00–1.21) | 4/50 | 0.86 (0.00–1.55) | 0.531 (0.508–0.553) |
| P3 | 2/51 | 0.43 (0.00–1.16) | 3/51 | 0.65 (0.00–1.30) | 0.503 (0.480–0.527) |

P2 的 True 超几何 p≈0.0016；P5 的 Random EF@1% 仍为 0。全文见 `data/si/protocol_enrichment_ci/protocol_ef_ci.csv`。

## 3. P2 完整案例 1,588 → 1,580

百分位只在双靶都有有效 gnina P2 分数的 1,580 个分子上计算。相对 NLRP3 缩库 1,588：

| 步骤 | n | 说明 |
|------|---|------|
| \(q_N\ge0.5\) 池 | 1588 | `docking_pool_p05.csv` |
| 配体 PDBQT | 1583 | 5 个脂质/核苷酸前药未进入 manifest |
| 单靶 docked | 1582 | `REP_05842`（tauroselcholic acid）空姿态 |
| 双靶完整案例 | **1580** | SMILES 内连接后再落 3 条（fostriecin 盐对、plocabulin） |

缺失分子表：`data/si/complete_case_drop/`。缺失集分子量中位约 721 Da，池中位约 480 Da。

## 4. 姿态质控（非 MD 数值）

7 个优选候选在 9DKB / 7ALV 生产构象上均 `both_in_pocket=True`（质心位移 ≤ 6 Å 或关键残基接触 ≥ 3；冲突截断 2.2 Å）。表：`data/si/pose_qc/pose_qc_dual.csv`。体系搭建记录含 GSK-3008348、Vecabrutinib、Zelenirstat；**不报告轨迹数值**。跟进假说仍为前两者。结果正文见 [`RESULTS_DRAFT_CN.md`](RESULTS_DRAFT_CN.md)。
