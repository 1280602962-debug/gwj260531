# 与 PLK1/NLRP3 不对称双靶论文的差异化说明

> **目的**：明确本项目的 URAT1/NLRP3 框架与 Gu et al.「Reliability-driven virtual screening for PLK1/NLRP3 dual-target candidate prioritization under asymmetric data conditions」在 **问题设定、算法选择、证据融合与候选来源** 上的本质区别，避免方法学雷同与审稿撞车。

---

## 一、两篇工作的共同背景（可承认，不必回避）

| 共性 | 说明 |
|------|------|
| 双靶组合 | 代谢/增殖层 + NLRP3 炎症层 |
| 数据不对称 | 一靶 ChEMBL 数据较丰富，NLRP3 测定异质性高 |
| 纯计算 | 无自建湿实验，依赖回顾性 benchmark |
| 多阶段漏斗 | ML → 对接 → ADMET → 候选排序 |

**以上属于领域通用范式，不构成创新点。** 本文档重点界定 **哪些模块不能照搬**，以及 **本项目独有的方法学贡献**。

---

## 二、PLK1/NLRP3 论文的核心方法（需刻意避开）

根据手稿，PLK1/NLRP3 工作的方法学指纹为：

```
PLK1（数据丰富）:  SVR + PCA(50) + Optuna + bootstrap/kNN 不确定性
NLRP3（数据贫乏）: 5-anchor 混合相似性（ECFP4 + 物化，max-pooling）+ ESM-2 口袋表征
融合:              0.5/0.5 等权归一化线性融合
候选来源:          仅商业化合物库虚拟筛选
结构验证:          Vina → Gnina → MM/GBSA → 100 ns MD
```

**若 URAT1/NLRP3 项目复用上述组合，将被视为增量换靶，创新性不足。**

---

## 三、模块级对照：避开 vs 采用

| 模块 | PLK1/NLRP3（避开） | URAT1/NLRP3 — TAPE-GATE（采用） |
|------|-------------------|--------------------------------|
| **数据丰富靶点建模** | SVR + PCA 降维 | **MiniMol/Chemprop 分子基础模型** + **Conformal prediction** 区间（非 kNN 邻域 UQ） |
| **URAT1 结构** | 激酶式单结构对接 | **转运体构象系综** + **$S_{\text{trap}}$ 构象捕获分**（inward/occluded vs outward） |
| **NLRP3 弱数据策略** | 多锚点指纹 max-pooling 相似性 | **Assay-conditioned 分类/排序**（CLAMP/TwinBooster 思路）或 **结构主导高精度分类** |
| **蛋白表征** | ESM-2 口袋 embedding 为主 | **晶体/冷冻电镜结合模式** + MM-GBSA/MD 稳定性；ESM-2 仅作可选消融 |
| **双靶融合** | 固定 0.5/0.5 线性加权 | **可靠性加权证据融合** + **Pareto 非支配排序**（按各臂 UQ 动态调权） |
| **候选来源** | 单一路径：库筛 | **双路径**：Path A 库筛 + Path B 生成式（CLM+RL）→ 统一候选池 |
| **迁移学习** | 无家族迁移 | **SLC22 家族 → URAT1** 分层迁移（转运体专属） |
| **标题关键词** | reliability-driven asymmetric | **transporter-aware conformation-ensemble** + **assay-conditioned** + **paired-path** |

---

## 四、NLRP3 处理：为何不用锚点相似性

PLK1/NLRP3 论文对 NLRP3 采用「5 个已知抑制剂作锚点、ECFP4+物化 max-pooling」—— 本质是 **活性数据不足时的相似性外推**。

本项目 ChEMBL 实测（用户数据）显示：

| 指标 | 数值 |
|------|------|
| 清洗后独特 SMILES（IL-1β + Assay B） | **513**（609 records；39 assays） |
| 不同 assay 间 >1 log 差异的化合物比例 | **~7.2%** |
| 与 URAT1 重叠 SMILES | **0** |

**结论**：
1. NLRP3 有足够数据支撑 **监督学习**，但须 **按 assay 条件化**，而非全局回归；
2. 锚点相似性在 assay 异质性下易放大噪声，且与 PLK1/NLRP3 高度雷同；
3. 推荐 **双轨 NLRP3 证据**：Assay-conditioned 分类概率 $P_{\text{active}}^{(a)}$ + 结构对接/MD 稳定性分 $S_{\text{NLRP3}}^{\text{struct}}$，由可靠性权重合并。

---

## 五、URAT1 处理：为何不用激酶式 QSAR+对接

| 维度 | PLK1 | URAT1 |
|------|------|-------|
| 蛋白类型 | 激酶，ATP 口袋 | **SLC22 转运体**，alternating access |
| 抑制机制 | 催化位点阻断 | **构象捕获**，阻断转运循环 |
| 结构资源 | 单构象为主 | **9DKB/9B1K/9B1L 三态 cryo-EM**（+ 9B1H 备用） |
| 数据量 | 数百至千级 | 清洗后 **822** 条（用户 ChEMBL） |

URAT1 侧创新核心：**$S_{\text{trap}}$ 构象系综评分** 与 **SLC22 家族迁移**，而非 SVR/PCA 管线。

---

## 六、双路径候选发现（本项目独有）

```
                    ┌─────────────────────────────────────┐
                    │         TAPE-GATE 双路径入口          │
                    └─────────────────────────────────────┘
                                      │
              ┌───────────────────────┴───────────────────────┐
              ▼                                               ▼
    ┌─────────────────────┐                     ┌─────────────────────┐
    │  Path A: 库筛        │                     │  Path B: 生成式      │
    │  Enamine/ChEMBL ~10⁶ │                     │  CLM fine-tune + RL  │
    │  ML/UQ 初筛          │                     │  双靶奖励函数优化     │
    │  → 构象系综对接      │                     │  → 同一对接漏斗      │
    └──────────┬──────────┘                     └──────────┬──────────┘
               │                                            │
               └────────────────────┬───────────────────────┘
                                    ▼
                    ┌─────────────────────────────────────┐
                    │  统一候选池 + 可靠性加权融合 + Pareto  │
                    │  回顾性 benchmark（lesinurad, MCC950） │
                    └─────────────────────────────────────┘
```

PLK1/NLRP3 论文 **无生成式路径**；双路径 + 合并排序是本项目的 **明确增量创新**。

---

## 七、证据融合公式对比

### PLK1/NLRP3（避开）

$$
S_{\text{dual}} = 0.5 \cdot \tilde{S}_{\text{PLK1}} + 0.5 \cdot \tilde{S}_{\text{NLRP3}}
$$

### TAPE-GATE（采用）

设 URAT1 预测区间宽度 $w_U$（conformal）、NLRP3 assay-conditioned 置信度 $c_N$、结构分 $S_U^{\text{struct}}, S_N^{\text{struct}}$：

$$
\omega_U = \frac{1/w_U}{\sum_k 1/w_k}, \quad \omega_N = \frac{c_N}{\sum_k c_k}
$$

$$
S_{\text{dual}}^{\text{fuse}} = \omega_U \cdot \hat{y}^U + \omega_N \cdot P_{\text{active}}^{N} + \gamma \cdot \sqrt{S_U^{\text{struct}} \cdot S_N^{\text{struct}}}
$$

最终排序：对 $(S_{\text{dual}}^{\text{fuse}}, \text{QED}, \text{SA}^{-1})$ 做 **Pareto rank**，取非支配前沿 Top-K。

---

## 八、论文叙事与关键词差异化

| 元素 | PLK1/NLRP3 | URAT1/NLRP3 (TAPE-GATE) |
|------|-----------|-------------------------|
| 疾病 | 肿瘤/炎症交叉 | **高尿酸血症/痛风** |
| 主靶生物学 | 细胞周期激酶 | **肾脏尿酸转运体** |
| 方法卖点 | reliability-driven asymmetric | **transporter-aware ensemble + assay-conditioned + paired-path generative** |
| 必做消融 | — | **vs PLK1-style baseline**（SVR+锚点相似性+0.5融合）作为阴性对照 |

---

## 九、审稿人「与 PLK1/NLRP3 有何不同」标准答复

1. **靶点生物学不同**：URAT1 是膜转运蛋白，必须构象系综评分，不能用激酶式单结构对接。
2. **NLRP3 策略不同**：我们用 assay-conditioned 监督模型 + 结构证据，而非锚点指纹相似性。
3. **候选空间不同**：库筛与生成式双路径扩展化学空间，非单库漏斗。
4. **融合机制不同**：可靠性加权 + Pareto，非固定等权。
5. **实证对照**：消融实验 **计划** 包含 PLK1-style pipeline（`06_retrospective_validation.py` 骨架），跑通后再定量对比 benchmark 回收。

---

## 十、参考文献（PLK1/NLRP3 对标）

- Gu et al., ECUST — PLK1/NLRP3 reliability-driven asymmetric screening（用户提供的内部手稿）
- 本项目框架：见 [`TAPE_GATE_FRAMEWORK.md`](TAPE_GATE_FRAMEWORK.md)、[`ALGORITHM_FRAMEWORK.md`](ALGORITHM_FRAMEWORK.md)
