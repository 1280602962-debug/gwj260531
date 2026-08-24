# 已有结果清单 ↔ Discussion 口径

> 只收录仓库里已经跑完、可以写进稿的结果。讨论句必须不越过 `CLAIM_CEILING.md`。  
> 配套：`RESULTS_DRAFT_ZH_JCIM_V1.md`（3.1–3.6）、`DISCUSSION_DRAFT_ZH_JCIM_V1.md`（4.1–4.6）、`DISCUSSION_LIMITATIONS_DRAFT_ZH_JCIM_V1.md`。

全文定位：不是 “docking can/cannot identify dual-target ligands”，也不是新打分函数；而是评价对接分数在四状态硬负面板、双向 pairwise 主终点上的**可靠性与边界**。资源名 DualFourClass-Bench。

定量主结论（Abstract / 4.6 共用，不要写成更强）：

> Docking showed limited and strongly target-pair-dependent discrimination. Only PIK3CA/mTOR showed a point estimate above chance and its strongest trivial descriptor baseline, but the uncertainty interval remained compatible with chance.

PIK3CA/mTOR 的 boxed 说法：

> limited directional signal + ligand-panel persistence, **not** receptor invariance.

---

## 0. 现行骨架（Results 6 节 + Discussion 6 节）

| 段 | 功能 | 吃进哪些结果 |
|---|---|---|
| 3.1 | 公开数据能不能建成严格双靶基准 | 49 对审计、K=4、S12 |
| 3.2 | docking 能不能识别 | Table 2、Fig 3、θ 网格→S4、GNINA S15 |
| 3.3 | 表观信号是不是口袋特异 | 描述符 Δ、TPSA、协变量、ECFP4、匹配子集 |
| 3.4 | PM 信号在哪些条件下成立 | E8、PM110、holdout、换晶、Cα |
| 3.5 | specificity control 是否可靠 | 主面板 matched>wrong；holdout 反转；S11/S13 |
| 3.6 | 结构背景线索（探索性） | Table S7 序列一致性；T2/T5 姿态事实 |
| 4.1 | 评测的是什么任务 | 3.1 + DUD/LIT-PCBA/Ahmed |
| 4.2 | 为什么难于两个单靶任务 | 3.2–3.3 + Ballester 2023 |
| 4.3 | PM 是有限方向信号 | 3.4 + 3.6；Schaller 只作 receptor 变量 |
| 4.4 | 对 VS / 生成式的含义 | DualDiff/FuseDiff 用途句；Sindt 后处理类比 |
| 4.5 | 未解决的错口袋反转 | 3.5 |
| 4.6 | 五条局限 | 库存前五条；收束段已移至 Conclusions |

旧 Results 3.8（面板构成层重抽样）已移入 4.6 / Limitations，不再作为 Results 小节。

---

## 1. 供给与面板构成（Results 3.1）

| 已有结果 | 关键数字 | 源 |
|---|---|---|
| 49 对 ChEMBL 严格供给审计 | 两端硬负 ≥ 50 仅 4 对；去掉 HDAC1/HDAC6 后厚面板 3 对 | J0 / Table 1 / Fig 2 |
| EGFR/HER2 进 K=4 的理由 | 严格 B_only = **7**，供给受限案例，不是厚面板 | J0 |
| BindingDB / PubChem 计数核对（零对接） | `equal_only` 厚门槛不翻转（min HN 76/92/58 与 86/97/61）；EGFR 升至 ~30，仍 < 50；`as_is` 可过 50 是因为 `>` 截尾 | Table S12 |

**讨论应写：** 可平衡的四状态硬负面板在公开数据里很少；K=4 是 *constrained but experimentally grounded* 的供给冻结集，不是对全部双靶任务的抽样。

**不要写：** 我们建立了一个非常完整的 benchmark；公开数据硬负供给 “ChEMBL-invariant” 而不提 S12；ChEMBL 漏了约 80 个干净 HER2 选择性配体。

---

## 2. 主对接结果（Results 3.2）

**Table 2（θ = 6.0，口袋匹配 Vina）**

| 靶对 | n (D/A/B) | summary_min [95% CI] | 相对随机 | 相对最强描述符 |
|---|---|---|---|---|
| EGFR/HER2 | 28/38/32 | **0.430** [0.284, 0.576] | 点估计 < 0.5 | 不超过 |
| AChE/BChE | 27/25/28 | **0.606** [0.440, 0.740] | CI 含 0.5 | 低于 TPSA **0.733** |
| PIK3CA/PIK3CB | 28/27/28 | **0.500** [0.347, 0.648] | 随机 | 不超过 |
| PIK3CA/mTOR | 18/14/12 | **0.692** [0.464, 0.802] | CI 下界贴 0.5 | 点估计 > heavy **0.463**，但 Δ 的 CI **含 0** |

其它必须跟主表一起讨论的对照：

| 结果 | 数字 | 源 |
|---|---|---|
| 统一 θ = 6.0 | AChE、PIK3CB 上与严格 6.5/5.5 **分类和 AUROC 完全相同** | Table S4 |
| 阈值网格 | 整张网格排序不变：PM 最高，其余 ≤ 0.61 | Table S4 / Fig S1A |
| 池化会掩盖弱臂 | EGFR 分臂 B 端 0.430 | Table S6 / Fig 4A |
| EGFR Top-10（池化 Vina） | 9/10 为硬负；bootstrap 均值 ≈ 8.9，CI 7–10 | 3.2 |
| GNINA 真口袋匹配 best-of-9 | EGFR 0.290；AChE 0.413；PIK3CB 0.533；PM 0.655 | Table S15 |

**讨论应写：** docking 判别有限且高度依赖靶对。不要写 “PIK3CA/mTOR performed significantly better”（CI 跨 0.5）。不要写 “docking fails across all pairs”。

---

## 3. 混淆主导（Results 3.3）— 核心发现，不是附加分析

| 结果 | 数字 | 源 |
|---|---|---|
| 配对 Δ（matched Vina − best descriptor） | −0.052 / −0.128 / −0.122 / **+0.229**；四对 CI 均含 0 | Table S19 / Fig S3C |
| AChE TPSA | dual vs 硬负 ~75 vs 51；TPSA ~0.769 > Vina ~0.56；+heavy 后 D vs B 0.606 → 0.807；OR ≈ 1.18 | Fig 4C / Fig 7C |
| PM 控制尺寸/极性后 | AUROC +0.07～+0.11；OR ≈ 2.19 / 3.08 = **有限残余口袋信号**，须与含 0 的 Δ 一起读 | Fig 7C |
| ECFP4 GroupKFold | 多方向 0.78–0.91；EGFR D vs B 指纹 0.85 vs docking 0.43 | Fig 7A |
| 支架 vs 随机 | 平均 +0.011 | Table S20 / Fig S3D |

**讨论应写：** 表观 dual signal 大量可由 ligand properties / chemotype 解释。PM 的残余 OR 不是已确证的独立优势。

**不要写：** 描述符对照是补充实验；对接显著优于平凡描述符。

---

## 4. 配体层持续 vs 受体层崩溃（Results 3.4）

| 已有结果 | 关键数字 | 源 |
|---|---|---|
| PM exhaustiveness 16→8 | 0.692 → 0.660（Δ ≈ 0.03） | Fig S1D |
| PM110 稳定性核对 | Vina 0.648 [0.51, 0.76]，Δ ≈ −0.04 | Fig S1C |
| Holdout PM / AChE / PIK3CB | **0.765** [0.603, 0.891] / 0.618 [0.422, 0.759] / 0.425 [0.241, 0.618] | Table S8 / Fig S5 |
| 换 PIK3CA（PM） | 0.692 → **0.486** / **0.505**；4JSX 0.639；D/A 仍 0.714 | Table S9 / Fig 5A |
| 换 PIK3CA（PIK3CB） | 0.500 → **0.691** / **0.685**；弱臂切换 | Table S30 / Fig 5B |
| 换 mTOR | **0.639**；CI 含 0.5 | S9 |
| Cognate QC | 4JPS 0.607 Å、5DXT 0.624 Å、4JSX 0.515 Å | S9 |
| 口袋局域 Cα | 5DXT **0.343 Å** 仍使 PM 到 0.505；同一晶体使 PIK3CB 上升 | Table S10 |
| A4 max vs median | 翻转 7/110、1/95、1/99、0/48；pair-level 基本不动 | Table S29 |

**讨论应写：** holdout **supports persistence of the observed signal in an unused ligand pool**，不验证 benchmark。PM = limited directional signal + ligand-panel persistence，不是 receptor invariance。同一 PIK3CA 扰动可升可降。pose QC ≠ screening invariance。PAB_034 100/99/1 timeout。

**不要写：** holdout 是独立文献/跨库验证；PM 是结构不变的 positive case；Cα 定量解释了 AUROC；receptor swap 一律 collapse；max pChEMBL 仍是 fatal ground-truth 风险。

---

## 5. 错口袋反转（Results 3.5）

| 已有结果 | 关键数字 | 源 |
|---|---|---|
| 主面板 matched − wrong | 0.170 / 0.161 / 0.151 / 0.090；仅 EGFR 与 AChE 的 CI 不含 0 | Table S17 / Fig 6A / Fig S3A |
| Holdout wrong ≥ matched | 0.788 vs 0.765；0.643 vs 0.618；0.520 vs 0.425；三对 Δ 的 CI **均含 0** | Fig 6B / Fig S3B |
| 效价/尺寸匹配后仍不翻转 | Table S13 | Fig 6C |
| contact_count | B 臂 0.698–0.714；不能按幅度解释 Vina（PM 0.788 vs 0.552） | Table S11 / Fig 6D |

**讨论应写：** 未解决的 failure mode。不要解释掉。不要写成 wrong-pocket validation。

---

## 6. 探索性结构线索（Results 3.6）

| 已有结果 | 关键数字 | 源 |
|---|---|---|
| 全链序列一致性 | PM 18.1% / PIK3CB 40.5% / AChE 51.9% / EGFR 71.4% | Table S7 |
| 姿态分型 | T2 硬负两端干净 ATP/hinge-like pose；T5 dual 在 Vina 强、重打分偏离 | failure typology |

**Results 只写看到什么。Discussion 4.3 才写 ATP-site chemotype transferability 可能。** n = 4，禁止当相关。

---

## 7. 只出现在 Introduction/Discussion 的“用途句”

| 陈述 | 允许 | 禁止 |
|---|---|---|
| DualDiff / FuseDiff 用两端 Vina 相对**参考配体**定义 Dual High Affinity | 本基准可作下游诚实评测 | 未重打分那些生成物；Dual High Affinity ≠ 均值池化 |
| Ahmed / Ballester / Schaller / Sindt | 见 `DISCUSSION_REFS_JCIM_V1.md` 边界 | 写成 DualFourClass 使用了他们的数据或引擎 |

---

## 8. 明确不是本稿结果

| 未做 / 不可当结果 | 怎么处理 |
|---|---|
| A4 max vs median | Results 3.4；Table S29；Limitations 3（controlled） |
| 1000 个互不重叠独立 panel | Limitations 1 |
| 主面板 PLIF；rotamer | Limitations 4 |
| BindingDB 对接面板 | 只有计数 S12 |
| 湿实验 | Limitations 5 |

---

## 9. 主张强度总表

**可以主张：**

- 建立了带硬负、口袋匹配主指标、混淆对照的评测协议，并释放 DualFourClass-Bench。
- 在该冻结集上，对接方向信号总体有限、高度依赖靶对。
- 仅 PIK3CA/mTOR 点估计同时高于随机与最强描述符，幅度有限，CI 与随机相容，Δ 的 CI 含 0。
- 多数对上表观信号可被描述符或化学型解释。
- 同一配体协议在 unused-pool 上对 PM **同向**，但换晶体后 **受体依赖**。
- holdout 错口袋反转保持开放。

**Conclusions（两段，见 `CONCLUSIONS_SECTION_JCIM_EN_V1.md`）只保留：** 0.430–0.692；PM 最强点估计 + unused-pool 同向；不确定度与受体敏感性排除可迁移决策规则；配体/化学型/受体实现可实质影响表观信号；holdout 错口袋反转未解决且配对 CI 含 0；贡献 = reliability boundary。

**不可以主张：**

- 对接能/不能识别双靶配体（全称）。
- 对接显著优于平凡描述符。
- PM 是结构不变、可重复的成功。
- 错口袋悖论或换晶崩盘已被机制解释。
- RTM/GNINA 验证了通用决策臂。
- 结论外推到全部双靶对接任务。
