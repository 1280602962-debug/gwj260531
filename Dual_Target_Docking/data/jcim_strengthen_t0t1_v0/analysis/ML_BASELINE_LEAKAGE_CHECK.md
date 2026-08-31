# ML baseline leakage check — random CV vs scaffold GroupKFold

> 同一 ECFP4 + LogisticRegression；仅交叉验证分折方式不同。
> 随机折：`StratifiedKFold(shuffle=True)`（易同系物泄漏）。
> 支架折：`GroupKFold` 按 Murcko，同一支架不跨折（主用）。

| pair | contrast | AUROC random | AUROC scaffold | Δ (rand−scaf) | dock PM |
|------|----------|--------------|----------------|---------------|---------|
| AChE/BChE | D_vs_A | 0.9096 | 0.8948 | +0.0148 | 0.6504 |
| AChE/BChE | D_vs_B | 0.8241 | 0.8214 | +0.0027 | 0.6058 |
| EGFR/HER2 | D_vs_A | 0.7961 | 0.7453 | +0.0508 | 0.6664 |
| EGFR/HER2 | D_vs_B | 0.8884 | 0.8895 | -0.0011 | 0.4297 |
| PIK3CA/PIK3CB | D_vs_A | 0.8042 | 0.7817 | +0.0225 | 0.6905 |
| PIK3CA/PIK3CB | D_vs_B | 0.8890 | 0.7691 | +0.1199 | 0.5000 |
| PIK3CA/mTOR | D_vs_A | 0.7262 | 0.7619 | -0.0357 | 0.7143 |
| PIK3CA/mTOR | D_vs_B | 0.9213 | 0.8889 | +0.0324 | 0.6921 |

## 结论

随机折相对支架折：平均 Δ=+0.026，最大 Δ=+0.120。
两版本差距不大；仍优先报告支架分组折。
若支架折 AUROC 接近对接或更低，应写：「表观易分是支架泄漏假象」。
