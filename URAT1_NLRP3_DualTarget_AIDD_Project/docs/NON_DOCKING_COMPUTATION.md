# 非对接计算模块（Non-Docking Computational Evidence）

> **定位**：这些模块是加在**现有 Glide XP 对接/Pareto 结果之上的下游注释与验证**。
> 它们 **不重新对接、不改动任何 Glide 分、不改变 Pareto 归属、不覆盖生产 ML 分**。
> 因此可与已完成的对接筛选数据 **共存并直接引用**。
>
> **对接引擎说明**：现有对接分为 **Glide XP**（开发跑）。本页所有分析基于该数据；
> 若日后切换对接引擎，只需重跑对接层，本页脚本无需改动即可复用。

---

## 0. 一览：模块 ↔ 脚本 ↔ 输出

| 模块 | 脚本 | 输出目录 | 是否改动对接结果 |
|------|------|----------|:---:|
| A ML 严谨性 | `scripts/12_ml_rigor_validation.py` | `results/model_validation/` | ❌ |
| B 结构警报过滤 | `scripts/09_cheminformatics_filters.py` | `results/cheminformatics/` | ❌ |
| C ADMET/类药性 | `scripts/10_admet_druglikeness.py` | `results/cheminformatics/` | ❌ |
| D 化学空间/新颖性 | `scripts/11_chemical_space_novelty.py` | `results/cheminformatics/` | ❌ |
| E Pareto 稳健性 | `scripts/13_pareto_robustness.py` | `results/pareto_robustness/` | ❌ |

复现：

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
python3 scripts/09_cheminformatics_filters.py
python3 scripts/10_admet_druglikeness.py
python3 scripts/11_chemical_space_novelty.py
python3 scripts/12_ml_rigor_validation.py --n-permutations 20
python3 scripts/13_pareto_robustness.py
```

---

## A. ML 严谨性（y-scrambling + 适用域 + 校准）

| 检验 | 结果 | 含义 |
|------|------|------|
| URAT1 y-scrambling | 真实 Spearman **0.732** vs 置换最大 **0.065**，经验 p≈**0.048** | 模型学到真实信号,非指纹伪影 |
| NLRP3 y-scrambling | 真实 AUROC **0.891** vs 置换最大 **0.559**，p≈**0.048** | 同上 |
| NLRP3 概率校准 | Brier **0.128**（见 `nlrp3_calibration.csv`） | 概率基本可信,可作缩库阈值 |
| 适用域(AD)阈值 | URAT1 训练集内 NN Tanimoto 5% 分位 = **0.578** | 低于此值即"域外" |

**关键发现（重要叙事）**：Pareto 六分子对 URAT1 训练集的最近邻 Tanimoto 全部 **≈0.15–0.26 ≪ 0.578**,即 **全部落在 URAT1 ML 适用域之外**。
→ 这从数据上证明:**对这些重定位命中,URAT1 ML 预测本就不可靠,必须靠对接**——与全文"URAT1 docking-led"的设计完全自洽,是对 2/4 benchmark 的正面回应而非辩解。

> 提高 permutation 次数(如 `--n-permutations 100`)可把经验 p 压到更低,投稿建议 ≥100。

---

## B. 结构警报（PAINS / Brenk / NIH + 聚集体启发式）

| 集合 | n | PAINS | Brenk | NIH | 任一警报 |
|------|---|-------|-------|-----|----------|
| Pareto 短名单 | 6 | **1** | 4 | 3 | 4 (66.7%) |
| 双靶对接池 | 1451 | 58 | 554 | 200 | 597 (41.1%) |

**关键发现**:短名单里 **唯一命中 PAINS 的就是 EGCG**（PAINS_B + Brenk）。
→ 直接支撑上一轮结论:EGCG 应从"最优双靶 lead"降级为"漏斗自我纠错的诚实案例"。
详见 `results/cheminformatics/filters_shortlist.csv`（含命中子结构描述）。

---

## C. ADMET / 类药性

| 集合 | 中位 QED | Lipinski | Veber | Egan | 三规则全过 |
|------|---------|----------|-------|------|-----------|
| 短名单 | 0.222 | 4/6 | 1/6 | 2/6 | 1/6 |
| 对接池 | 0.436 | 1181 | 1176 | 1073 | 925/1451 |

**EGCG 明细**:MW 458、cLogP 2.23、**TPSA 197、HBD 8、HBA 11**,Lipinski/Veber/口服吸收 **全不过**。
→ 定量印证多酚"高极性、氢键过多、口服生物利用度低"的转化短板。
详见 `results/cheminformatics/admet_shortlist.csv`。

---

## D. 化学空间与新颖性

- ECFP4 + PCA 二维坐标见 `results/cheminformatics/chemical_space_pca.csv`（4 组:URAT1 活性、NLRP3 活性、对接池、短名单）。前两主成分方差占比仅 0.050 / 0.027,说明化学空间高度分散(指纹稀疏,正常)。
- **新颖性(最近邻 Tanimoto)**:短名单对 URAT1 已知活性 0.15–0.26、对 NLRP3 已知活性 0.15–0.21。

**双刃解读**:命中在结构上"新颖"(远离已知活性),但结合 A 的适用域结论——**新颖 = ML 外推不可靠**,再次说明为何 URAT1 必须结构证据主导。详见 `novelty_shortlist.csv`。

---

## E. Pareto 稳健性（不改动原前沿）

| 分析 | 结果 |
|------|------|
| 重算前沿大小 | 6（与生产一致,验证可复现） |
| Top-1% 双交集 | 1 |
| Top-2% | 3 |
| **Top-5%** | **10** |
| **Top-10%** | **46** |
| 阈值门 τ=85 / 90 / 95 | 83 / 46 / 10 |
| Bootstrap(500)前沿频率 | 六分子均 ≈0.60–0.67 |

**结论**:原始 6 分子前沿虽薄,但 bootstrap 显示其成员**中等稳定**(频率 ~0.6);若正文需要更可讨论的短名单,可改用 **Top-5%(10 个)或 Top-10%(46 个)双百分位交集**。详见 `results/pareto_robustness/`。

---

## 对稿件的落点建议

| 稿件位置 | 用哪个模块 |
|----------|-----------|
| Methods 新增"计算过滤与模型验证"小节 | A–E 全部 |
| Results:模型可信度(y-scramble、校准、AD) | A |
| Results/Discussion:EGCG 降级论证 | B + C + A(域外) |
| Results:短名单化学空间与新颖性图 | D |
| Results:Pareto 敏感性(SI 或正文) | E |
| Discussion:局限性(域外、PAINS、薄前沿) | A + B + E |

> 这些模块共同把论文从"仅对接 + MD"升级为"对接 + 模型验证 + 化学过滤 + 多目标稳健性"的完整计算证据链,且**不触碰现有对接筛选结果**。
