# DualFourClass-Bench analysis pack v1 — CI paper packaging

> Generated 2026-07-29 · no new docking · exploration/evaluation pool  
> Script: `scripts/build_benchmark_analysis_v1.py` · figures: `scripts/plot_forest_ci_v1.py`  
> Claim ceiling unchanged: evaluation/benchmark Article, not a universal scorer.

---

## 0. 一句话

对接已经够了；这一轮把现有 K=4 结果收成**可引用基准表 + bootstrap 95% CI + 失败模式诊断**。  
主结论不变，但措辞必须更严：**只有 PIK3CA/mTOR 的点估计超过平凡基线，且其 Δ 的 95% CI 仍跨 0**——论文卖的是任务/指标/基线协议，不是“对接显著胜利”。

---

## 1. 本轮补了什么（相对 POST_DOCKING_VERDICT）

| 模块 | 表/图 | 论文用途 |
|------|-------|----------|
| 统一 assembled 长表 | `tables/assembled_*` | Zenodo 主数据 |
| Bootstrap CI（B=2000） | `bootstrap_directional_ci_v1.csv`, `forest_summary_min_ci_v1.csv` | Fig 主 forest |
| 基线门控 Δ+CI | `baseline_gate_bootstrap_v1.csv` | Fig / Table：对接是否显著优于基线 |
| 端不对称 + 池化欺骗 | `asymmetry_pooled_vs_directional_v1.csv` | 方法学结果 |
| Top10 硬负组成 + CI | `top10_hardneg_bootstrap_v1.csv` | 决策失败案例 |
| 连续 Spearman | `continuous_spearman_v1.csv` | 标签/分数相关性附录 |
| θ 敏感度 | `threshold_sensitivity_v1.csv` | 标签稳健性 |
| AChE TPSA 捷径 | `ache_feature_dual_vs_hardneg_v1.csv` | 阴性对解释 |
| Prep 敏感 | `pm48_directional_by_prep_v1.csv` | 协议冻结证据 |
| 失败模式 JSON | `analysis_meta_v1.json` | 写作提纲 |

---

## 2. 主结果（带 CI）

`summary_min = min(AUROC dual–A_only, dual–B_only)`；配体 bootstrap 95% CI。

### 2.1 Vina vs 最佳平凡基线

| pair | vina summary_min [CI] | best baseline | Δ [CI] | 门控读法 |
|------|----------------------|---------------|--------|----------|
| **PIK3CA/mTOR** | **0.671 [0.43, 0.81]** | heavy 0.463 | **+0.208 [−0.02, +0.43]** | 点估计 PASS；CI 跨 0 → **不能宣称显著优于体积** |
| AChE/BChE | 0.530 [0.37, 0.68] | **TPSA 0.733** | −0.203 [−0.40, +0.01] | 点估计 FAIL；CI 刚跨 0 |
| PIK3CA/PIK3CB | 0.412 [0.26, 0.56] | heavy 0.622 | **−0.210 [−0.40, −0.01]** | **显著输给体积基线** |
| EGFR/HER2 | 0.282 [0.16, 0.42] | cLogP 0.482 | **−0.200 [−0.34, −0.03]** | **显著输给 cLogP**；弱臂 D/B |

RTM / GNINA：在 EGFR、AChE 上相对最佳基线多为显著负；在 PM 上点估计略正但 CI 更宽。**换引擎不改变评测文叙事。**

### 2.2 失败模式（写入 Discussion）

1. **PIK3CA/mTOR** — `beats_baseline_point_but_CI_includes_0`（N=44 有效四类；功率有限）  
2. **EGFR/HER2** — 弱 B 臂 + 池化掩盖 + 显著低于 cLogP  
3. **PIK3CA/PIK3CB** — 同模式：D/A≈0.70 vs D/B≈0.41；池化抬高；显著低于 heavy  
4. **AChE/BChE** — TPSA 捷径：dual 平均 TPSA≈75 vs 硬负≈51；TPSA dual-vs-hardneg AUROC≈0.77，远高于 vina≈0.56

---

## 3. 诊断摘要（论文可直接用）

### 3.1 AChE/BChE · TPSA chemotype confound

- dual vs 硬负（A∪B）：TPSA AUROC **0.769**；vina **0.559**；RTM/GNINA <0.5  
- 结论：该对上“极性/化学型”可当捷径；对接不能当作已学到双靶决策。必须在文中作为**阴性对照对**，不是方法失败的遮羞布。

### 3.2 Prep 敏感（PM48）

| prep | vina summary_min | rtm_min_z summary_min |
|------|------------------|----------------------|
| LigPrep (old) | 0.597 | 0.611 |
| RDKit (primary) | **0.671** | **0.520** |

主表必须统一 RDKit；RTM 对 prep 更敏感——与 EGFR M4 一致。

### 3.3 Top10 硬负（vina）

| pair | Top10 hardneg (point) | bootstrap CI |
|------|----------------------|--------------|
| EGFR/HER2 | 9/10 | [7, 10] |
| PIK3CA/PIK3CB | 7/10 | [3, 9] |
| AChE/BChE | 4/10 | [1, 8] |
| PIK3CA/mTOR | 4/10 | [1, 8] |

EGFR 顶部几乎全是硬负——方向失败的直观证据。

### 3.4 θ 敏感度（双侧 pChEMBL 重标）

- EGFR：θ=5.5→6.0→6.5 时 vina `summary_min` ≈ 0.37 / 0.28 / 0.30（弱臂问题稳健）  
- PM：θ=5.5 欠功率（B_only=5）；θ=6.0/6.5 稳定在 ~0.67  
- 严格面板（AChE、PIK3CB）按建造规则预筛选，简单 θ 重标变化小——报告时注明“面板已按 strict 配额冻结”

### 3.5 连续 Spearman（附录）

见 `continuous_spearman_v1.csv`：对接分数与 `min_pChEMBL` 相关普遍弱/不稳定；平凡描述符在部分对上更强——强化“任务难、基线必要”。

---

## 4. 对 JCIM 包装的判决（更新）

| 缺口 | 本轮前 | 本轮后 |
|------|--------|--------|
| 点估计 forest | 有 | 有 |
| Bootstrap CI | 缺 | **有** |
| 基线门控显著性 | 点比较 | **Δ + CI** |
| 失败模式分类 | 叙述 | **结构化 JSON + 表** |
| Zenodo 可复现表 | 半成品 | **assembled + 脚本** |
| 英文稿 / 图注 | 缺 | 下一步（图已出 PNG/PDF） |

**不再缺一轮巨大对接。** 缺的是：英文 Article 正文、Zenodo DOI、（可选）PM48 在严格配额下扩样本以收窄 CI——那是增强，不是门槛。

### 仍建议在文中写清的限制

- PM 的“超过体积”**不是** CI 显著  
- K=4 是冻结评测集，不是泛化证明  
- 严格硬负供给仍稀缺（J0/J1 审计保留为供给结果）

---

## 5. 复现

```bash
cd Dual_Target_Docking
python3 data/jcim_bench_v0/scripts/build_benchmark_analysis_v1.py
python3 data/jcim_bench_v0/scripts/plot_forest_ci_v1.py
```

依赖：`numpy`, `rdkit`, `scipy`（Spearman）, `matplotlib`（作图）。
