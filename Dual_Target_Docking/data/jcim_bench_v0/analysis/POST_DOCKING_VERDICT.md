# JCIM 对接阶段结果分析（K=4 完成后）

> 数据：`ache_bche_panel_v0`、`pik3ca_pik3cb_panel_v0`、`pik3ca_mtor_panel48_rdkit_v0`、EGFR 统一 prep EH110 + GNINA  
> 汇总包：`jcim_bench_v0/`  
> 本分析补表：`tables/directional_with_baselines_v1.csv`（含平凡基线；原 forest 缺 EGFR）  
> 日期：2026-07-29  
> 数据池：Exploration / evaluation（评测集，非方法选臂确认）

---

## 0. 一句话

**对接阶段按计划做完了（Vina + RTM + GNINA），而且强化了评测文叙事：四对靶里只有 PIK3CA/mTOR 上对接稳定超过体积基线；其余三对在方向分解 + 基线下都暴露失败。**  
这够支撑 **JCIM 评测/基准 Article（审慎 claim）**，不够支撑「新打分/决策臂胜利」文。

---

## 1. 完成度盘点

| 对 | 面板 | Prep | Vina | RTM | GNINA | 角色 |
|----|------|------|------|-----|-------|------|
| EGFR/HER2 | N=110（θ=6，供给受限） | 统一 RDKit | 已有 | 已有 | mode_01 有 | 案例 |
| PIK3CA/mTOR | N=48 | **RDKit 重跑 96/96** | 有 | 有 | 有 | 主开发 |
| AChE/BChE | N=100 strict | RDKit | 191/200（9 fail） | 有 | 96 配体 | 新主对 |
| PIK3CA/PIK3CB | N=100 strict | RDKit | 199/200 | 有 | 100 | 同工酶对照 |

授权文件：`protocols/PAIR_ROLES_APPROVED_JCIM.yaml`  
Claim 天花板：`jcim_bench_v0/CLAIM_CEILING.md`（正确：评测文，不卖通吃臂）

---

## 2. 主结果：方向分解 × 三通道 × 基线

`summary_min = min(AUROC dual-vs-A_only, dual-vs-B_only)`

| pair | vina | rtm_min_z | gnina_cnn_min | 最佳平凡基线 | 对接 vs 基线 |
|------|------|-----------|--------------|--------------|--------------|
| **PIK3CA/mTOR** | **0.671** (0.722/0.671) | 0.520 | 0.563 | heavy 0.463 | **全部 PASS** |
| AChE/BChE | 0.530 (0.530/0.585) | 0.409 | 0.372 | **TPSA 0.753** | 全部 FAIL |
| PIK3CA/PIK3CB | 0.412 (0.703/**0.412**) | 0.439 | 0.506 | heavy 0.599 | 全部 FAIL |
| EGFR/HER2 | 0.282 (0.680/**0.282**) | 0.253 | 0.263 | cLogP 0.482 | 全部 FAIL |

（括号内为 D/A / D/B。）

### 可写成论文的结论

1. **池化指标会骗人（跨对靶复现）**  
   EGFR：D/A≈0.68 与 D/B≈0.28 平均成假随机；  
   PIK3CA/PIK3CB：D/A≈0.70 与 D/B≈0.41 同样一端反转。  
   → 「端不对称」不是 EGFR 特例。

2. **对接增量高度对靶依赖**  
   仅 **PIK3CA/mTOR** 上 Vina/RTM/GNINA 的 `summary_min` 超过体积基线；  
   另外三对在该汇总指标下都输给平凡描述符（AChE 上 TPSA 尤其强，需在文中讨论 chemotype/极性捷径）。

3. **换引擎不等于救任务**  
   GNINA CNN rescore 在 AChE、EGFR 上并不优于 Vina；在 PIK3CA/PIK3CB 上 `summary_min` 提到 0.51，仍低于 heavy 0.60。  
   → 支持「问题在任务/指标/标签，不单是 Vina 太差」。

4. **Prep 协议有实证**  
   PM48 RDKit 后 Vina `summary_min` 0.597→**0.671**（变好），RTM 0.611→**0.520**（变差）。  
   → 主表必须统一 prep；RTM 对 prep 敏感（与 EGFR M4 一致）。

5. **Strict 面板可建，但供给仍紧**  
   AChE/BChE、PIK3CA/PIK3CB 均做成 dual/A/B/neither = 28/28/28/16；  
   结合 J0：公开数据里够厚的对仍然极少。

---

## 3. 对「能不能投 JCIM」的判决

### 可以冲的形态

**Evaluation / Benchmark Article**，主贡献是：

- 四类双靶决策任务形式化  
- 方向分解评测 + 必报平凡基线  
- 49 对供给审计（strict 硬负稀缺）  
- K=4 统一 prep × 三通道实证（含阴性）  
- 开放 DualFourClass-Bench  

### 投前缺口（2026-07-29 CI 包后更新）

| 缺口 | 状态 |
|------|------|
| 主文级 bootstrap CI / 森林图（含 EGFR + 基线） | **已补** → `BENCHMARK_ANALYSIS_V1.md` + `figures/` |
| AChE TPSA 化学型捷径诊断 | **已补** → `ache_feature_dual_vs_hardneg_v1.csv` |
| 基线门控 Δ±CI / 失败模式 | **已补** → `baseline_gate_bootstrap_v1.csv`；**PM Δ CI 仍跨 0** |
| Zenodo 正式 DOI 上传 | 表/脚本已齐；待人工 deposit |
| 英文稿 + JCIM 版式 | **下一步（写作）** |
| PM 扩到 ~110（收窄 CI） | 可选加分，非门槛 |

### 不要用现在的数据去卖

- 「找到了通用双靶决策臂」  
- 「GNINA/RTM 全面提升」  
- 「四对靶上都验证了协议」——实际是 **一对点估计有对接增量（CI 未显著）、三对暴露失败**  
- 「对接显著优于平凡基线」——仅 PM 点估计为正，且 95% CI 跨 0

### 综合判断

| 问题 | 答案 |
|------|------|
| 数据是否支持开始写 JCIM？ | **是（评测文）** |
| 分析包是否够写主图主表？ | **是**（见 `BENCHMARK_ANALYSIS_V1.md` / `FIGURE_PLAN_V1.md`） |
| 现在直接投稿是否稳？ | 分析侧可写；差英文稿 + Zenodo DOI |
| 命中率（评测定位） | **中等可冲**；若审稿人要「新方法+显著提升」会拒——那不是你们该卖的 |
| 保底刊 | 仍可用同一稿投 Molecular Informatics / JCAMD |

---

## 4. 下一步（白话）

### 第一步 — 已完成（零对接分析收束）
见 [`BENCHMARK_ANALYSIS_V1.md`](BENCHMARK_ANALYSIS_V1.md)：bootstrap CI、基线门控、TPSA、Top10、θ 敏感、prep、figure plan。

### 第二步 — 写作（当前主线）
- Intro：任务缺口（对标大库双靶发现文：他们有湿实验，你们补评测轴）  
- Results：供给审计 → 方向抵消 → **带 CI 的基线对照森林图** → prep 敏感性 → 三通道  
- 诚实写：PM 超过体积仅为点估计；EGFR/同工酶显著输基线  
- 不写 Track B 选臂胜利  

### 第三步（可选，加分）
- PM 按 strict 扩到 ~110 以收窄 Δ CI  
- Zenodo DOI + GitHub release  

**不建议再开：** 扩 EGFR、再找第 5 对赌方法赢、短板阈值调参、再开一轮巨大对接。

---

## 5. 与路线文件的关系

本轮结果 **符合** `JCIM_NEXT_ROUTE_PLAIN.md`；分析打包（CI/基线/诊断）已收官。  
Track B 仍不必复活：阴性结果本身就是评测文的证据。
