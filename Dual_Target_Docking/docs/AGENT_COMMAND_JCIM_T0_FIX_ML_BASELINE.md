# Agent 命令 — 修 ML 基线泄漏 + CI 展示 bug（零对接）

> 依据：[`../data/jcim_strengthen_t0t1_v0/analysis/CLOUD_QC_REVIEW_V1.md`](../data/jcim_strengthen_t0t1_v0/analysis/CLOUD_QC_REVIEW_V1.md)  
> 工作量：几十行代码、零对接、几分钟运行

```text
【任务】
修复 data/jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py 中两个问题：

1. ligand_ml_baseline()：
   当前用 StratifiedKFold(shuffle=True) 按配体随机分折，Morgan/ECFP4 指纹在
   支架高度聚集的小面板上会跨折泄漏同系物信息，AUROC 被夸大（观察到 0.89–0.92）。
   改为：
   - 用已算好的 Murcko scaffold（scaffold_inventory_v1.csv 或重算）做分组
   - sklearn GroupKFold（或等价的按 scaffold 手动分折），保证同一支架不跨折
   - 折数按最小类别样本数与支架簇数取 min，且写清 n_splits
   - 输出两个版本对照存档：
     ligand_ml_baseline_random_cv_v1.csv   （旧版，保留作对照，加列 note=potential_leakage）
     ligand_ml_baseline_scaffold_cv_v1.csv （新主用版本）
   - 在 analysis/ 下写 ML_BASELINE_LEAKAGE_CHECK.md：对比两版本 AUROC 差值，
     若差值大（如 >0.15）需在文中明确写「随机折高估，支架分组折为准」。

2. matched_subset_rows() 中 2 类子集（potency_matched_D_vs_A/B、size_matched_D_vs_A/B）
   目前复用 boot_pm_ci()（该函数假设三类都在，算 min(D/A, D/B)），
   导致缺失的第三类 AUROC 恒为 nan，min 也是 nan，bootstrap 全部被过滤，CI 输出为空。
   改为：
   - 新写 boot_single_contrast_ci(recs, score_key, pos_cls="dual", neg_cls)，
     直接对该单一对比做 AUROC 配体 bootstrap（B=2000），不经过 3 类 min。
   - 用它替换 potency_matched_*/size_matched_* 四行的 CI 计算；
     full_panel_pocket_matched 行继续用原 boot_pm_ci()（三类都在，不受影响）。
   - 重新生成 matched_subset_directional_v1.csv，确认不再有空 CI（除非 n<8 应该
     标 underpowered 而不是留空）。

3. 顺手：unified_threshold_sensitivity_v2.csv 增加一列 underpowered
   （min(n_dual,n_A_only,n_B_only) < 8 时标 1，如 PM θ=5.5 的 B_only=5）。

【约束】
- 不做任何新对接；只改统计/评测代码并重跑。
- 不改 Vina/RTM/GNINA 分数表。
- CLAIM_CEILING 不变。

【完成标准】
- [ ] 两版本 ML 基线表存在，leakage check md 写清差值与结论
- [ ] matched_subset 表里 D_vs_B 系列不再有异常空 CI（n≥8 时必须有数）
- [ ] threshold 表加 underpowered 列
- [ ] git add/commit/push；更新 PR #23

结束用中文总结：随机折 vs 分组折 AUROC 差多少、哪些结论需要改措辞。
```
