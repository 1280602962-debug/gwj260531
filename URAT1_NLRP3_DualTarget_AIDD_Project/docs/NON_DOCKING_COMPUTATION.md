# 非对接计算模块

> 在 **gnina P2** 双成功合并表上运行。不要把已删除的 Glide XP 表当作输入。  
> 生产协议 **P2 / gnina CNNaffinity**。提名：`scripts/14_candidate_nomination.py`（默认读 `data/repurposing/p2/`）。

---

## 0. 一览：模块 ↔ 脚本 ↔ 输出

| 模块 | 脚本 | 输出目录 | 是否改动对接结果 |
|------|------|----------|:---:|
| A ML 严谨性 | `scripts/12_ml_rigor_validation.py` | `results/model_validation/` | ❌ |
| B 结构警报过滤 | `scripts/09_cheminformatics_filters.py` | `results/cheminformatics/` | ❌ |
| C ADMET/类药性 | `scripts/10_admet_druglikeness.py` | `results/cheminformatics/` | ❌ |
| D 化学空间/新颖性 | `scripts/11_chemical_space_novelty.py` | `results/cheminformatics/` | ❌ |
| E Pareto 稳健性 | `scripts/13_pareto_robustness.py` | `results/pareto_robustness/` | ❌ |

复现（输入为已归档 P2 合并表 `data/repurposing/p2/pareto_merged_scores.csv`）：

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
python3 scripts/11_chemical_space_novelty.py --pool data/repurposing/p2/pareto_merged_scores.csv --shortlist data/repurposing/p2/pareto_shortlist.csv --output-dir data/repurposing/p2
python3 scripts/13_pareto_robustness.py --pool data/repurposing/p2/pareto_merged_scores.csv --output-dir data/repurposing/p2/pareto_robustness
python3 scripts/14_candidate_nomination.py --tau 90
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

B–E 已在 P2 完整案例（n=1,580）上归档：`data/repurposing/p2/`。

| 模块 | 关键计数 |
|------|----------|
| B 结构警报 | PAINS 78 / Brenk 626 / 1,580 |
| C 类药性 | Lipinski 752；Veber 1,254；MW 200–550：1,199 |
| D 新颖性 | GSK-3008348 NN Tanimoto URAT1/NLRP3 ≈ 0.21 / 0.20；Vecabrutinib ≈ 0.25 / 0.25 |
| E 稳健性 | 审计轴 τ=90 为 77；双对接门控 51；裸前沿 4（大环） |
| F 提名 | 优选 7；跟进 GSK-3008348 + Vecabrutinib |
