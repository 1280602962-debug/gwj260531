# URAT1 ML 模型选型评估：为什么选 XGBoost？其他模型更好吗？

> 基于 `scripts/08_urat1_model_comparison.py` 在 **相同 5-fold 骨架 CV** 与 **相同 4 个 benchmark 药物** 上的对比实验（2026-06）。

---

## 一、为什么最初选 XGBoost + Morgan 指纹？

| 原因 | 说明 |
|------|------|
| **小数据常规 baseline** | 822 化合物、218 骨架，表格特征 + 树模型是 QSAR 领域最常用、最可复现的起点 |
| **训练快、无 GPU** | 相对 Chemprop/GNN 秒级完成，适合方法学流水线骨架 |
| **可解释 + conformal UQ** | 比深度模型更容易接 split conformal 不确定性区间 |
| **并非「最优选择」** | 是 **MVP 默认项**，不是文献论证后的 URAT1 专用最优模型 |

**结论**：选 XGBoost 是工程上的合理起点，**不是**因为实验证明它最适合 URAT1 转运体场景。

---

## 二、对比实验结果（同一数据、同一 CV、同一 benchmark）

### 2.1 骨架 OOF 交叉验证（越高越好）

| 模型 | RMSE↓ | R²↑ | Spearman↑ | EF@5%(p≥7)↑ | Benchmark(4药) |
|------|-------|-----|-----------|-------------|----------------|
| **XGBoost + Morgan + RDKit（当前）** | 0.649 | **0.530** | **0.732** | 3.10 | **2/4** |
| XGBoost + Morgan only | 0.661 | 0.512 | 0.715 | 3.00 | 2/4 |
| Random Forest | 0.658 | 0.516 | 0.732 | **3.39** | 2/4 |
| kNN Tanimoto (k=5) | 0.657 | 0.517 | 0.716 | 2.62 | 2/4 |
| SVR + PCA50（PLK1 式） | 0.721 | 0.419 | 0.658 | 3.00 | 2/4 |
| XGBoost assay-conditioned | 0.677 | 0.487 | 0.703 | 1.94 | 2/4 |
| **Chemprop D-MPNN（DL）** | **0.743** | **0.384** | **0.640** | 2.71 | **2/4** |
| Kernel Ridge RBF | 3.632 | -13.7 | 0.278 | 1.45 | 4/4 ⚠️ |

### 2.2 Benchmark 明细（全数据训练后预测）

**所有正常模型**（XGBoost / RF / SVR / kNN / Chemprop）失败模式一致：

| 药物 | XGBoost | RandomForest | Chemprop | 在训练集？ |
|------|---------|--------------|----------|-----------|
| lesinurad | 5.67 ❌ | ~6.6 ✅/❌ | 5.48 ❌ | **否**（assay 冲突被清洗剔除） |
| benzbromarone | 6.65 ✅ | ✅ | 6.02 ✅ | **否** |
| verinurad | 6.98 ✅ | ✅ | 6.57 ✅ | **是** |
| dotinurad | 5.10 ❌ | ❌ | 5.46 ❌ | **否** |

### 2.3 Kernel Ridge「4/4 满分」是假象

Kernel Ridge 在全训练集上预测 lesinurad=7.77、dotinurad=6.97（全过），但 **OOF R² = -13.7**，说明严重过拟合。

**教训**：只看 benchmark 回收、不看骨架 CV，会被过拟合模型欺骗。

---

## 三、核心结论：换模型解决不了 URAT1 的根本问题

```
问题不在「XGBoost vs Chemprop」，而在：

1. 数据覆盖：lesinurad / dotinurad 等上市药物骨架不在训练集
2. 标签质量：ChEMBL 中 lesinurad 多 assay 中位 pAct≈5.1，与文献 ~7.0 矛盾
3. 任务本质：URAT1 是转运体，构象依赖性强，指纹回归难以捕捉机制
4. 骨架外推：218 个 Murcko 骨架上要外推到全新药物系列，ML 本身就不稳
```

**实证**：7 种 ML/DL 模型 + assay-conditioned 变体，**没有任何一个**在严格 OOF 评估下同时做到：
- 显著优于 XGBoost（Chemprop 更差）
- benchmark 4/4 回收（正常模型全部 2/4）

Random Forest 在 EF@5%(p≥7) 上略好（3.39 vs 3.10），但 **benchmark 仍为 2/4**，不具备替换价值。

---

## 四、各模型路径的客观评价

| 方向 | 能否显著改善 URAT1？ | 说明 |
|------|---------------------|------|
| 换 XGBoost → Random Forest | ⚠️ 微小 | CV 略好，benchmark 不变 |
| 换 → SVR+PCA（PLK1 式） | ❌ | CV 更差，无优势 |
| 换 → Chemprop D-MPNN | ❌ | OOF R² 0.38 < 0.53，benchmark 仍 2/4 |
| 换 → kNN 相似性 | ❌ | 略差于 XGBoost |
| Assay-conditioned 回归 | ❌ | 66 个 assay 异质性低于 NLRP3，收益有限 |
| MiniMol / 预训练 GNN 微调 | ❓ 未测 | 可能略升 0.02–0.05 R²，**不太可能**单独解决 benchmark |
| **增加专利 SAR 数据** | 可能 | 补全 lesinurad 等骨架的训练覆盖 |
| **三态 \(S_{\mathrm{trap}}\) 系综对接** | **未实现** | 禁止写成当前方法 |
| **放弃 URAT1 纯 ML 筛** | **当前策略** | NLRP3 分类缩库；URAT1 由 **9DKB P2 对接** 排序 |

---

## 五、当前 URAT1 证据角色（已锁定）

```
URAT1：结构排序（Π* = P2，gnina CNNaffinity @ 9DKB inward-open）
       ML 回归只作对照 / SI，不主排临床库
NLRP3：assay-conditioned 分类缩库（P≥0.5）+ 7ALV 对接百分位门控
禁止：用 URAT1 ML 排名决定候选去留；把未做的三态对接写成生产证据
```

本对比表说明 **模型选择不是性能瓶颈、数据与任务定义才是**。不必为投稿再换 Random Forest。

---

## 六、如何复现

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project/scripts
export PATH="$HOME/.local/bin:$PATH"
python3 08_urat1_model_comparison.py
```

结果：`results/urat1_model_comparison/comparison_summary.csv`
