# 审稿人视角：模型质量质疑与回应

本文档用于 Methods / Discussion / rebuttal，说明 **URAT1 与 NLRP3 两个 ML 模型质量不等** 时，漏斗设计如何自洽，以及为何不应被解读为「弱模型硬筛库」。

---

## 1. 审稿人可能提出的质疑

| 质疑 | 实质 |
|------|------|
| URAT1 模型 R²≈0.51、benchmark 仅 2/4 | 「URAT1 ML 太差，整个双靶流程不可靠」 |
| NLRP3 AUROC 0.89 是否过拟合 | 「临床库外推性不足」 |
| colchicine / verinurad 高 P(active) | 「ML 把间接机制当成直接抑制」 |
| 对接与 ML 几乎无关 (ρ≈−0.04) | 「两路证据互相矛盾」 |
| Pareto 六分子多为 Phase 2 | 「输出没有临床可转化性」 |

---

## 2. 我们的立场（应在正文中明确写出）

### 2.1 URAT1 ML：** deliberately not used for main-screen ranking**

- 内部报告 `MODEL_QUALITY_REPORT.md` 结论：**URAT1_NO_GO** — 禁止单独用 ML 排序临床库。
- Benchmark：lesinurad、dotinurad 在严格阈值下 **fail**；与 ChEMBL 标注冲突、支架新颖 benchmark 有关（见 Discussion）。
- **主漏斗 URAT1 轴 = 9DKB 结构对接百分位 $S_U$**，不是 URAT1 ML。
- **8973 回顾轨**（AUC 0.705，EF@5%=4.23）独立验证「应对接 lead URAT1」，并展示 dotinurad「对接高、ML 低」——这恰恰支持 **对接主导**，而非 ML 主导。

> 回应话术：*We do not claim URAT1 ML is a standalone virtual screen. It is retained only for retrospective enrichment analysis on the 8973 distill set; all repurposing ranking on the clinical library uses 9DKB docking percentiles.*

### 2.2 NLRP3 ML：** library reduction only, not sole ranking**

- AUROC 0.89 / AUPRC 0.91（scaffold GroupKFold）→ 适合 **缩库**（8319→1588），不适合单独宣称「NLRP3 抑制剂排序」。
- $S_N = \max(P_{ML}, P_{dock})$：当 ML 与对接解耦（ρ≈−0.04）时，**取较强证据**，避免单一路径偏倚。
- colchicine **高 ML、低双轴** → 漏斗按设计 **拒绝** 其进入 Pareto，说明 ML 局限已在结果中 **自我暴露**。

> 回应话术：*High NLRP3 ML scores for colchicine illustrate phenotypic assay confounding; the composite axis and Pareto layer explicitly down-rank such compounds.*

### 2.3 双模型「质量差」≠ 流程无效

本工作的贡献是 **不对称双证据漏斗 + 可复现 benchmark 行为**，不是「训练两个 SOTA 单靶模型」。

| 组件 | 质量 | 在漏斗中的角色 |
|------|------|----------------|
| URAT1 ML | 中等 / benchmark 不稳 | **不用于**主筛排序 |
| NLRP3 ML | 较好 | **缩库** P≥0.5 |
| URAT1 Vina @ 9DKB | 回顾 AUC~0.7 | **主 URAT1 轴** |
| NLRP3 Vina @ 7ALV | 假设性 pose | **与 ML 取 max** |

---

## 3. 开源对接迁移（回应许可与可重复性）

因无 Schrödinger 商业许可，全部对接改为 **AutoDock Vina 1.2.5**（`exhaustiveness=32`，与 Glide SP 级搜索深度相当；见 `config/docking_open_source.yaml`）。

- **禁止**在正文将 Vina kcal/mol 与历史 Glide XP 数值混为一谈。
- **允许**在同一引擎、同一蛋白、同一配体池内做 **百分位排序**（$S_U$, $S_N^{dock}$）。
- 8973 回顾轨需用 **同一 Vina 协议** 重对接后更新 AUC/EF（脚本：`run_vina_batch.py`）。

---

## 4. 建议在正文增加的 Limitations 段

1. URAT1 regression ML recovers only 2/4 scaffold-novel uricosurics under strict thresholds.  
2. NLRP3 ML conflates direct and indirect inflammasome modulation (colchicine).  
3. Docking scores are pose hypotheses, not affinities; open-source Vina scores are not calibrated to experimental $K_i$.  
4. Pareto output is computationally non-dominated within the pool, not clinical validation.  
5. MD and experimental target engagement are required before mechanism claims.

---

## 5. 不推荐的辩护方式

- ❌ 「URAT1 AUROC 也很高」—— 回归任务上 EF@10%@p≥6 在 57% 阳性率下上限约 1.75，易误导。  
- ❌ 「两个模型都很好」—— 与内部 benchmark 矛盾，审稿人会抓。  
- ❌ 把 Glide 与 Vina 分数放在同一张表比较绝对值。

---

*数值来源：`data/models/training_report.json`，`docs/MODEL_QUALITY_REPORT.md`，`pareto_benchmark_report.json`。*
