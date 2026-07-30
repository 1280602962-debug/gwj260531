# ML baseline leakage check — random CV vs scaffold GroupKFold

> 同一 ECFP4 + LogisticRegression；仅交叉验证分折方式不同。
> 随机折：`StratifiedKFold(shuffle=True)`（易同系物泄漏）。
> 支架折：`GroupKFold` 按 Murcko，同一支架不跨折（主用）。

| pair | contrast | AUROC random | AUROC scaffold | Δ (rand−scaf) | dock PM |
|------|----------|--------------|----------------|---------------|---------|
| AChE/BChE | D_vs_A | 0.9096 | 0.9096 | +0.0000 | 0.6504 |
| AChE/BChE | D_vs_B | 0.8241 | 0.8426 | -0.0185 | 0.6058 |
| EGFR/HER2 | D_vs_A | 0.7961 | 0.8327 | -0.0366 | 0.6664 |
| EGFR/HER2 | D_vs_B | 0.8884 | 0.8527 | +0.0357 | 0.4297 |
| PIK3CA/PIK3CB | D_vs_A | 0.8042 | 0.7751 | +0.0291 | 0.6905 |
| PIK3CA/PIK3CB | D_vs_B | 0.8890 | 0.7857 | +0.1033 | 0.5000 |
| PIK3CA/mTOR | D_vs_A | 0.7262 | 0.7817 | -0.0555 | 0.7143 |
| PIK3CA/mTOR | D_vs_B | 0.9213 | 0.8889 | +0.0324 | 0.6921 |

## 结论

随机折相对支架折：平均 Δ=+0.011，最大 Δ=+0.103（PIK3CA/PIK3CB D_vs_B）。

**未达到「Δ>0.15 全面虚高」阈值**，但：
1. 仍以 `ligand_ml_baseline_scaffold_cv_v1.csv` 为主用版本；随机折仅作泄漏对照。
2. 支架折 AUROC 仍普遍 **0.78–0.91**，高于多数口袋匹配对接臂（尤其 EGFR D_vs_B：ML 0.85 vs dock 0.43）→ 可解释为「标签与 2D 化学型相关」，**不能**写成「对接被随机折 0.89–0.92 碾压」。
3. 面板支架簇小、多数支架接近 singleton 时，GroupKFold 对泄漏的压制有限；该限制应在 Limitation 中写明。
