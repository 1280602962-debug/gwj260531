# A 组裁决 — 零对接加厚（T0.1–T0.10）

> 执行包：`data/jcim_strengthen_t0t1_v0/`  
> 脚本：`scripts/build_t0_strengthen_v1.py` + 既有 `build_pocket_matched_diagnostics_v1.py`

## 完成了什么

| ID | 产出 | 状态 |
|----|------|------|
| T0.1–0.2 | 口袋匹配主表升格 → `PRIMARY_METRIC_V2.md` | ✅ |
| T0.3 | `tables/matched_subset_directional_v1.csv` | ✅ |
| T0.4 | `tables/covariate_adjusted_v1.csv` | ✅ |
| T0.5 | `tables/aggregation_sensitivity_v1.csv` | ✅ |
| T0.6 | `tables/unified_threshold_sensitivity_v2.csv` | ✅ |
| T0.7 | `tables/chembl_aggregation_sensitivity_v1.csv` + SKIP 说明 | ✅（median/confidence 缺字段见 T0_SKIPS.md） |
| T0.8–0.9 | `scaffold_inventory_v1.csv`, `scaffold_bootstrap_ci_v1.csv` | ✅ |
| T0.8 加分 | `ligand_ml_baseline_scaffold_cv_v1.csv`（主）+ `…_random_cv_v1.csv` + `ML_BASELINE_LEAKAGE_CHECK.md` | ✅ 已修泄漏评测 |
| T0.10 | ENV_PIN.md, POSE_UPLOAD_CHECKLIST.md, MANIFEST 修复, CLAIM_CEILING 更新 | ✅ |

## 主结论是否因口袋匹配而改变？

**是，排序与数值均变，但 claim ceiling 不变。**

1. **指标修正**：池化 `vina_mean` 在 EGFR/PIK3CB 上系统性低估或扭曲 B 端对比；升格口袋匹配后主表数字与 `REVIEWER_AUDIT_V1.md` 一致。
2. **仍无通用决策臂**：四对中仅 **PIK3CA/mTOR** pocket-matched summary_min=0.69 [0.46, 0.81] 明显高于 0.5 且 LE 归一后仍优于 heavy；其余三对 summary_min ≤0.61。ECFP4 ML 基线须用 **支架 GroupKFold**（见 `ML_BASELINE_LEAKAGE_CHECK.md`）；随机折不可作正文数字。
3. **混淆未清除**：效价/尺寸匹配后（单对比 CI 已修复）EGFR/PIK3CB 的 D vs B 口袋匹配仍偏弱/近随机；错口袋对照远离 0.5 → 信号多为分子属性而非口袋特异性。
4. **协变量**：logistic 调整 heavy+TPSA 后，AChE D vs B 的 AUROC 从 0.61→0.81（Δ=+0.20），提示对接分数与尺寸/TPSA 共线；PM 调整幅度较小（+0.07–0.11）。
5. **Murcko bootstrap**：scaffold 重采样 CI 与配体 bootstrap 同量级；PM scaffold CI [0.46, 0.81] 仍宽。
6. **ML 基线**：支架折相对随机折平均 Δ≈+0.01、最大 Δ≈+0.10；但支架折 AUROC 仍常 >0.75，高于多数对接臂 → 正文可写「2D 化学型本身可部分区分标签」，**不得**写「随机折 0.89–0.92 全面碾压对接」。

## 对投稿的含义

- 正文必须写 **pocket-matched** 为主指标；池化/错口袋/LE 进 SI 作诊断。
- 不得写「对接通吃四类」；最多写「PM 对在严格面板下方向信号可检测，但特异性有限」。
- Wave 2 对接（exhaustiveness / enrichment / PM110）为 B 组任务，用于回应 E=16 与单靶 sanity 质疑。

## 建议下一步

B 组已完成（见 `B_GROUP_VERDICT.md`）。正文引用 ML 基线时只用 scaffold CV 表。
