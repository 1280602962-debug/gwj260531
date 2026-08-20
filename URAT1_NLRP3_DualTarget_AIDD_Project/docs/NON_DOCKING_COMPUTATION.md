# 非对接计算模块

> 在 **gnina P2** 双成功合并表上运行。不要把已删除的 Glide XP 表当作输入。  
> 生产协议 **P2 / gnina CNNaffinity**。提名：`scripts/14_candidate_nomination.py`（默认读 `results/repurposing/`）。

---

## 0. 一览：模块 ↔ 脚本 ↔ 输出

| 模块 | 脚本 | 输出目录 | 是否改动对接结果 |
|------|------|----------|:---:|
| A ML 严谨性 | `scripts/12_ml_rigor_validation.py` | `results/model_validation/` | ❌ |
| B 结构警报过滤 | `scripts/09_cheminformatics_filters.py` | `results/cheminformatics/` | ❌ |
| C ADMET/类药性 | `scripts/10_admet_druglikeness.py` | `results/cheminformatics/` | ❌ |
| D 化学空间/新颖性 | `scripts/11_chemical_space_novelty.py` | `results/cheminformatics/` | ❌ |
| E Pareto 稳健性 | `scripts/13_pareto_robustness.py` | `results/pareto_robustness/` | ❌ |

复现（须先有 P2 合并表 `results/repurposing/pareto_merged_scores.csv`）：

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

这些数字来自训练集交叉验证，不依赖对接引擎。

| 检验 | 结果 | 含义 |
|------|------|------|
| URAT1 y-scrambling | 真实 Spearman **0.732** vs 置换最大 **0.065**，经验 p≈**0.048** | 模型学到真实信号,非指纹伪影 |
| NLRP3 y-scrambling | 真实 AUROC **0.891** vs 置换最大 **0.559**，p≈**0.048** | 同上 |
| NLRP3 概率校准 | Brier **0.128**（见 `nlrp3_calibration.csv`） | 概率基本可信,可作缩库阈值 |
| 适用域(AD)阈值 | URAT1 训练集内 NN Tanimoto 5% 分位 = **0.578** | 低于此值即"域外" |

对照药相对训练集的最近邻见 `results/model_validation/applicability_domain.csv`（benchmark 行）。短名单适用域在 P2 提名表生成后再填。

B–E 的计数表须在 P2 漏斗完成后重跑。仓库不保留历史 Glide 池的 PAINS / 类药性 / Pareto 稳健性数字。
