# 自研双靶对接打分算法：实现指南

> 配套文献：`DUAL_MULTI_TARGET_DOCKING_SURVEY.md`  
> 现有挂钩：`config/docking_ensemble.yaml`（公式已写，无 runner）、`scripts/06_virtual_screening.py`（ML 复合分）、`config/targets.yaml`（权重）

---

## 1. 设计原则（先定“分什么”，再定“怎么合”）

自研双靶打分 = **单靶能量/几何组件** × **跨靶聚合算子** × **可选药化惩罚**。

不要一开始就写“端到端神经网络对接”。更稳的路径：

1. 用现成引擎（Vina / GNINA）只负责 **采样姿态**；
2. 你自己定义 **重打分（rescoring）组件**；
3. 用可学习或可调权重把两靶合成为 `S_dual`；
4. 用已知双靶/单靶分子校准权重。

这正是文献里 consensus scoring、CompScore、DualDiff 评估协议的工程化版本。

---

## 2. 推荐架构（四层）

```
SMILES
  │
  ├─[L0] 构象准备 (RDKit ETKDG → 可选 Meeko PDBQT)
  │
  ├─[L1] 姿态采样 (Vina/GNINA；每靶 top-K poses)   ← 可外包，不自研
  │
  ├─[L2] 单靶分量打分 φ_k(pose, pocket)            ← 你自研的核心
  │         例如: vina_aff, cnn_affinity, hbond,
  │               hydrophobic, clash, strain, pharmacophore_match
  │
  ├─[L3] 单靶聚合  S_A = g(φ(pose*)) ; S_B = g(φ(pose*))
  │
  └─[L4] 双靶合成  S_dual = F(S_A, S_B, props)      ← 双靶专用
```

**你该自研的是 L2–L4**；L1 用现成工具即可。

---

## 3. 打分组成：建议的分量菜单

每个姿态算一个向量 \(\phi \in \mathbb{R}^d\)，再线性或非线性合成。

| 分量 ID | 含义 | 典型来源 | 归一化建议 |
|---------|------|----------|------------|
| `vina` | 经验场对接分 | Vina / Smina | Z-score 或分位数 |
| `gnina_cnn` | CNN 亲和/姿态质量 | GNINA CNNaffinity / CNNscore | [0,1] 或 Z |
| `hbond` | 关键氢键计数/能量 | RDKit / PLIP / 自定义几何 | /max_expected |
| `hydrophobic` | 疏水接触面积 | 几何计数 | Z |
| `clash` | 立体冲突 | VDW overlap | 越大越差 → 取负 |
| `strain` | 配体应变能 | 对接构象 vs 真空最低能 | 越大越差 → 取负 |
| `pharm` | 药效团匹配分 | 相对共晶或共识药效团 | [0,1] |
| `buried` | 埋藏表面积 / 形状互补 | MSMS / RDKit | Z |

单靶分：

\[
S_{\text{target}} = \sum_k w_k\,\tilde\phi_k
\quad\text{或}\quad
S_{\text{target}} = \mathrm{MLP}(\tilde\phi)
\]

其中 \(\tilde\phi_k\) 是归一化后的分量。

---

## 4. 双靶合成算子 \(F(S_A, S_B)\)（关键设计选择）

设对接分“越负越好”（Vina 风格）时，先统一符号约定：  
**内部一律转成“越大越好”**：`score_up = -vina`。

| 算子 | 公式（越大越好） | 适用 |
|------|------------------|------|
| **Arithmetic mean** | \((S_A+S_B)/2\) | 两靶同等重要 |
| **Geometric mean** | \(\sqrt{S_A S_B}\)（需先平移到正） | 惩罚偏科 |
| **Max-penalty (类 Max Vina)** | \(\min(S_A,S_B)\) | 强制两边都不差（推荐默认） |
| **Softmin** | \(-\tau\log(e^{-S_A/\tau}+e^{-S_B/\tau})\) | 可微、可调尖锐度 |
| **Weighted sum** | \(\alpha S_A+(1-\alpha)S_B\) | 主靶/辅靶权重不同 |
| **Selectivity-style** | \(S_A - \max(S_{off})\) | 你仓库 JNK1 选择性已用此型 |
| **Pareto rank** | 非支配排序 + crowding | 多目标（亲和A/B + QED + SA） |
| **Rank fusion** | 对两靶排名做 ECR / 平均秩 | 尺度不可比时最稳 |
| **Learned fusion** | GBDT/logistic on \([S_A,S_B,\phi_A,\phi_B]\) | 有标注双靶数据时 |

**实践建议：**

- 数据少、两靶都要强 → 先用 **`min(S_A,S_B)` + 药化惩罚**；
- 只要“别太偏科” → softmin 或几何平均；
- 有一批已知双靶阳性 + 单靶阴性 → 学一个小融合模型。

与你现有 JNK 配置的关系：

```yaml
# config/docking_ensemble.yaml 已有（选择性，不是双激活）
scoring:
  selectivity: Score_JNK1 - max(Score_JNK2, Score_JNK3)
```

双**激活**（两靶都要抑制）应改成例如：

```yaml
scoring:
  dual_activation: min(Score_TargetA, Score_TargetB)   # 越大越好约定下
  # 或: softmin(Score_A, Score_B, tau=0.5)
```

---

## 5. 最小可运行实现骨架（可直接落进 `scripts/`）

### 5.1 配置扩展示例

```yaml
# config/dual_target_scoring.yaml
targets:
  A:
    name: JNK1   # 或 NLRP3
    receptors: [3ELJ, 4L7F]
    aggregate: mean          # 多晶体 → 单靶分
  B:
    name: JNK2               # 或另一疾病靶
    receptors: [3E7O]
    aggregate: identity

pose_engine: vina            # vina | gnina
poses_per_ligand: 9

components:
  - id: vina
    weight: 0.45
    direction: higher_better   # 若原始 vina 为负，pipeline 内取反
  - id: hbond
    weight: 0.20
  - id: clash
    weight: 0.15
    direction: higher_better   # 已取负后
  - id: strain
    weight: 0.10
  - id: pharm
    weight: 0.10

normalization: robust_z       # zscore | robust_z | quantile | none

dual_fusion:
  operator: softmin           # min | mean | softmin | weighted | rank_ecr
  tau: 0.5
  alpha: 0.5                  # for weighted
  property_penalty:
    qed_weight: 0.15
    sa_weight: 0.10           # 对 SA 高惩罚

gates:
  min_score_A: 0.0
  min_score_B: 0.0
  require_both_pass: true
```

### 5.2 核心伪代码

```python
# scripts/dual_target_rescoring.py（逻辑示意）

def normalize(df, cols, method="robust_z"):
    ...

def component_scores(pose, pocket, ligand) -> dict:
    return {
        "vina": -pose.vina_affinity,          # → higher better
        "hbond": count_hbonds(pose, pocket),
        "clash": -clash_energy(pose, pocket),
        "strain": -ligand_strain(pose, ligand),
        "pharm": pharmacophore_match(pose, pocket.ref_pharm),
    }

def score_target(poses, weights) -> float:
    # 选 top pose：可用 vina 粗排，再对 top-3 重打分取 max
    best = max(poses, key=lambda p: dot(weights, component_scores(p)))
    return dot(weights, component_scores(best))

def fuse_dual(s_a, s_b, op="softmin", tau=0.5, alpha=0.5):
    if op == "min":
        return min(s_a, s_b)
    if op == "mean":
        return 0.5 * (s_a + s_b)
    if op == "weighted":
        return alpha * s_a + (1 - alpha) * s_b
    if op == "softmin":
        # differentiable “差不多取较差侧”
        import numpy as np
        return -tau * np.logaddexp(-s_a / tau, -s_b / tau)
    raise ValueError(op)

def score_molecule(smiles, cfg):
    poses_a = dock(smiles, cfg.targets.A)   # 外包 Vina
    poses_b = dock(smiles, cfg.targets.B)
    s_a = score_target(poses_a, cfg.weights)
    s_b = score_target(poses_b, cfg.weights)
    s_dual = fuse_dual(s_a, s_b, **cfg.dual_fusion)
    s_dual += cfg.qed_weight * qed(smiles)
    s_dual -= cfg.sa_weight * sa_score(smiles)
    return {"S_A": s_a, "S_B": s_b, "S_dual": s_dual}
```

### 5.3 与现有漏斗衔接

```
06_virtual_screening.py  → topN SMILES
        ↓
08_dual_target_dock_score.py  (新)
        ↓
  all_rescored.csv  (S_A, S_B, S_dual, components...)
        ↓
  校准 / 多样性 / 采购列表
```

JNK **选择性**可继续用已有公式；若做 **NLRP3–JNK 双抑制**，L4 换成 `min`/`softmin`，不要用 `S1 - max(S2,S3)`。

---

## 6. 权重怎么定（三种难度）

### Level A — 手工 + 网格（无训练）

1. 收集：若干已知双靶阳性、仅 A 活性、仅 B 活性、双阴 decoy。  
2. 网格搜索 \(w_k\) 与 fusion 超参，优化：
   - 早期富集 EF@1% / BEDROC；或
   - 双阳性排在仅单阳性之前的 pair-accuracy。  
3. 固定权重上线。

### Level B — 线性/树模型融合（有中等标注）

特征：\([S_A, S_B, \phi_A, \phi_B, \mathrm{QED}, \mathrm{SA}]\)  
标签：`is_dual_active` 或 \(\min(pIC50_A, pIC50_B)\)  
模型：Logistic / XGBoost（你仓库已有 XGB 栈，可复用 `utils_ml.py` 风格）。

### Level C — 可学习打分函数（进阶）

- 姿态级：原子对距离指纹 → 预测 ΔG（RF-Score / NNScore / GNINA 微调思路）。  
- 双靶：两塔编码器 + fusion head，loss =  
  \(\mathcal{L}_A + \mathcal{L}_B + \lambda \mathcal{L}_{\text{dual}}\)  
  其中 \(\mathcal{L}_{\text{dual}}\) 可用 softmin 亲和或 pairwise ranking。

**建议从 Level A 起步**，验证协议通了再上 B。

---

## 7. 必须做的校准与验收（否则权重无意义）

| 阶段 | 做什么 | 通过标准 |
|------|--------|----------|
| 姿态协议 | 共晶自对接 | RMSD ≤ 2 Å（与 `docking_ensemble.yaml` validation_gate 一致） |
| 单靶筛选力 | 已知抑制剂 vs decoy | EF@1% 显著 > 随机 |
| 双靶分辨 | 双阳 vs 单阳/双阴 | Dual High Affinity 或 PR-AUC 提升 |
| 消融 | 去掉某一分量 | 证明该分量贡献 |
| 稳健性 | 换 fusion：min/mean/softmin | 头部 hit 重叠率报告 |

输出列建议：

```text
smiles, vina_A, vina_B, S_A, S_B, S_dual,
hbond_A, hbond_B, clash_A, clash_B, strain, qed, sa,
gate_pass_A, gate_pass_B, gate_pass_dual
```

---

## 8. Linked 分子的额外项（若做 linker 型双靶）

在 L2 增加：

- `bridge_ok`：两端锚点距离是否可被 linker 覆盖（0/1 或连续惩罚）；  
- `linker_strain`：linker 扭转应变；  
- `clash_protein_protein`：若模拟三元/双结构域。

合成时：

\[
S_{\text{dual}} = \min(S_A,S_B) - \lambda_1\mathrm{linker\_strain} - \lambda_2(1-\mathrm{bridge\_ok})
\]

这是 TwistDock 思想的打分化，不必一开始就上完整构象系综搜索。

---

## 9. 落在本仓库的具体步骤清单

1. **冻结接口**：新增 `config/dual_target_scoring.yaml`（分量 + fusion）。  
2. **实现 rescoring 库**：`scripts/lib/scoring_components.py`（hbond/clash/strain/pharm）。  
3. **实现 runner**：`scripts/08_dual_target_dock_score.py`  
   - 输入：`results/screening_v2/top*_diverse.csv`  
   - 对接：先 subprocess 调 Vina/GNINA（或读已有 pose 文件）  
   - 输出：带全部分量的 CSV + 简短 JSON 报告。  
4. **校准脚本**：`scripts/09_calibrate_dual_fusion.py`（网格/随机搜权重，benchmark 集）。  
5. **验收**：用 `Dual_Target_Docking/data/` 扩展“双阳/单阳”标签列，报告 EF 与消融。  
6. **再接入漏斗**：可选把 `S_dual` 写回类似 `06` 的 `final_score` 加权项。

当前缺口：仓库 **已有选择性公式配置，没有对接 runner 与分量打分代码**——从第 2–3 步开工即可。

---

## 10. 一个可发表的最小算法表述（便于写进方法学）

> We dock each candidate to pockets A and B with engine E, retaining the top-K poses.  
> Each pose is rescored by a weighted sum of physicochemically motivated terms  
> \(\phi=(\phi_{\mathrm{vina}},\phi_{\mathrm{HB}},\phi_{\mathrm{clash}},\phi_{\mathrm{strain}},\phi_{\mathrm{pharm}})\).  
> Target scores \(S_A,S_B\) use the best rescored pose (or ensemble mean over receptors).  
> Dual-target score is \(S_{\mathrm{dual}}=\mathrm{softmin}_\tau(S_A,S_B)\) plus drug-likeness penalties.  
> Weights and \(\tau\) are selected by maximizing enrichment of known dual-actives against single-target actives and decoys.

这就是“自研双靶打分算法”的标准叙事；创新点通常落在：**分量设计、融合算子、校准目标**，而不是重写采样引擎。

---

## 11. 常见坑

1. **两靶原始分尺度不同** → 先分靶 Z-score / 分位数，再 fusion；或改用 **排名融合**。  
2. **用平均值掩盖偏科** → 双激活场景优先 min/softmin。  
3. **只优化对接分忽略成药性** → linked 分子会系统性刷高“可对接性”。  
4. **没有单阳对照集** → 分不清“真双靶”和“口袋泛粘”。  
5. **分量强相关**（vina 与 hydrophobic）→ 权重不稳定；先看相关矩阵再降权。

---

## 12. 下一步（若要我直接写代码）

可按下列优先级直接实现到仓库：

1. `scoring_components.py` + 单元测试（不依赖 Vina，用假 pose）；  
2. `08_dual_target_dock_score.py` 读 YAML，对已有 CSV 做 **rescoring-only** 演示；  
3. `09_calibrate_dual_fusion.py` 接 benchmark。

指定靶对（例如 NLRP3–JNK1 或 JNK1 双位点）后即可开工。
