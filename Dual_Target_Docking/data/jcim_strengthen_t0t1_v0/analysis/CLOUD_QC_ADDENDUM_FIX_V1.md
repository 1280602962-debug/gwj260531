# Cloud QC addendum — ML baseline / matched-subset fix (`a03676a`)

> 复核对象：本地 agent 对 `AGENT_COMMAND_JCIM_T0_FIX_ML_BASELINE.md` 的修复推送  
> 结论：**修复请求已正确完成；另外意外修掉了匹配子集 D_vs_B 的错口袋点估计 bug。**

## 1. 请求项核对

| 要求 | 状态 | 证据 |
|------|------|------|
| 支架 GroupKFold 主用 ML 基线 | ✅ | `ligand_ml_baseline_scaffold_cv_v1.csv`；`ligand_ml_baseline_v1.csv` 已指向 scaffold 版 |
| 随机折保留作对照 | ✅ | `ligand_ml_baseline_random_cv_v1.csv`，`note=potential_leakage` |
| 泄漏对照 md | ✅ | `ML_BASELINE_LEAKAGE_CHECK.md` |
| 2 类子集单对比 bootstrap CI | ✅ | `boot_single_contrast_ci()`；`matched_subset_directional_v1.csv` 的 D_vs_B 行 CI 已填满 |
| θ 表 underpowered 列 | ✅ | PM θ=5.5、PM strict、EGFR strict 标 1 |

独立抽查：EGFR D_vs_B 随机折 AUROC **0.8884** 与表完全一致；GroupKFold 折间支架重叠 **0**。

## 2. 科学读数（比“修好了”更重要）

### 2.1 随机折 vs 支架折：泄漏没有大到 0.15，但原因要写清

| 统计 | 值 |
|------|----|
| 平均 Δ(rand−scaf) | **+0.011** |
| 最大 Δ | **+0.103**（PIK3CA/PIK3CB D_vs_B） |
| EGFR D_vs_B singleton 支架占比 | **~88%** |

**解读：** 多数支架是 singleton 时，GroupKFold 对泄漏的压制本来就有限——所以 AUROC 掉不下来，**不等于**“没有泄漏风险、随机折也可放心用”。正文仍应以 scaffold 版为准，并把 singleton 限制写进 Limitations。

### 2.2 支架折后 ML 仍普遍高于对接 → 这是真发现

即使 scaffold CV，ECFP4+LR 仍常在 **0.78–0.91**，且：
- EGFR D_vs_B：ML **0.85** vs dock **0.43**
- AChE：ML **0.84–0.91** vs dock **0.61–0.65**
- PM D_vs_A：ML **0.78** ≈ dock **0.71**（最接近）

**可写：** 四类标签与 2D 化学型相关，配体基线是必要对照。  
**不可写：** “随机折 0.89–0.92 证明对接全面失败”。

### 2.3 匹配子集点估计变化 = 旧 bug 被修掉（加分）

旧表 D_vs_B 点估计偏低（EGFR potency **0.316** / size **0.341**），是因为旧代码对两类子集误走 `directional_pm(..., vina_A, vina_B)`，把 **D_vs_B 算成了错口袋 vina_B**。  
新表改用单对比 + **vina_A**（EGFR potency **0.469** / size **0.519**）。独立复算：错口袋 vina_B ≈0.34–0.35，与旧表接近 → **确认旧值是错口袋伪影**。

修正后的匹配读数更合理：
- EGFR / PIK3CB：匹配后 D_vs_B ≈ **0.45–0.52**（近随机，CI 跨 0.5）→ 弱端识别仍不成立
- PM：匹配后 D_vs_A/B 仍约 **0.71–0.78**，CI 下界多仍 >0.4 → 方向信号在匹配后仍可见（但 N 小）

## 3. 对投稿状态的判决

| 问题 | 答案 |
|------|------|
| 上轮 CLOUD_QC 要求的返工是否完成？ | **是** |
| 还要不要再开对接？ | **不要** |
| A/B 组材料能否进正文/SI？ | **可以**，ML 只用 scaffold 表 |
| 还缺什么？ | 写作（英文稿）+ Zenodo 姿态；可选 C 组 holdout |

**一句话：** 这次补充把评测代码层面的洞补上了；科学故事更干净——PM 仍是唯一在匹配/LE 后站得住的对，其余对的表观信号更多是化学型/混淆，而不是口袋特异对接。
