# 化合物 / 文献 / PDB 全量核对报告

核对日期：2026-07-10  
范围：Phase 0 种子表、主表、REFERENCES.md、`docs/report/`

---

## 一、Phase 0 种子分子（7 条）

| ID | 分子 | 结构 | 文献 | PDB | 结论 |
|----|------|------|------|-----|------|
| ACT001 | YL5084 | ✅ ChEMBL5398852 = PubChem 153610906；MW 600.73 | **Lu 2023** [R1] | 无独立 PDB；pose 来自 8ELC | ✅ |
| ACT002 | YL2056 | ✅ ChEMBL5413843；MW 586.70；丙烯酰胺 | **Lu 2023** [R1] | **8ELC** 共晶配体（PDB Y56 为 post-covalent butanamide） | ✅ |
| ACT003 | JNK-IN-8 | ✅ ChEMBL2216824 = PubChem 57340686；MW 507.60；**丙烯酰胺**（非 chloroacetamide） | **Zhang 2012** [R2] | 无独立 PDB；叠合 **3V6S** / 8ELC | ✅ |
| ACT004 | 56d | ✅ 用户/Wydra SMILES；MW 516.56；CHEMBL6167401 同结构 | **Wydra 2025** [R20] | 无 PDB；叠合 4WHZ+3V6S | ✅ |
| NEG002 | YL5084R | ✅ ChEMBL5436581；butanamide | Lu 2023 [R1] | — | ✅ |
| NEG001 | JNK-IN-6 | ✅ PubChem 57340684；propionamide；MW 438.49 | Zhang 2012 [R2] | — | ✅ |
| NEG003 | 56a | ✅ **已补** CHEMBL6151222；propionamide；MW 518.58 | Wydra 2025 [R20] | — | ✅ 新补 |

### 已修正的历史错误

| 问题 | 正确信息 |
|------|----------|
| 56d 曾绑 CHEMBL5947460 | 嘌呤结构，Tc≈0.17，**非 56d** |
| REPORT 写 YL5084/YL2056 → Wydra 2014 | **Lu 2023** J Med Chem |
| REPORT 写 JNK-IN-8/56d → chloroacetamide | 均为 **丙烯酰胺 Michael 受体** |
| REPORT 写 JNK-IN-6 为阳性 hit | **丙酰胺阴性对照** |
| REPORT 写 JNK1 对应 Met146 | 同位点为 **Ile106**（Met146 为 JNK3 Leu144 同源） |
| REPORT [R7b] 重复指 COValid | COValid 应为 **[R15]** |
| fig02 中 26k 2D 结构错误 | 已换为 4WHZ 配体 **3NL（26k）** RCSB SMILES |
| 8ELC  caption 写 YL5084 共晶 | 共晶配体为 **YL2056** |

---

## 二、文献编号对照（REFERENCES.md）

| 编号 | 内容 | 关联分子/结构 |
|------|------|---------------|
| **[R1]** | Lu W et al., *J Med Chem* **2023** — YL5084/YL2056 | YL5084, YL2056, YL5084R |
| **[R2]** | Zhang T et al., *Chem Biol* **2012** — JNK-IN 系列 | JNK-IN-2/6/7/8/11/12 |
| **[R3]** | PDB **8ELC** — JNK2–YL2056 | 主共价模板 |
| **[R4]** | PDB **7N8T** — JNK2–AMP | 对照结构 |
| **[R5]** | PDB **3V6R** — JNK3–JNK-IN-2 | 共价叠合参照 |
| **[R7]** | PDB **3NPC** — BIRB796 DFG-out | 共价初筛禁用 |
| **[R7b]** | PDB **4WHZ** — JNK3–**26k**（配体 3NL） | 可逆 DFG-in 选择性 |
| **[R15]** | Shamir et al., COValid / AF3 共价门控 | Phase 0 AF3 EF@1% |
| **[R20]** | Wydra et al., *J Med Chem* **2025** — 56d/56a 系列 | 56d, 56a, 21b, 51d |

> **注意：** [R7b] 在 REFERENCES.md 中标题为 4WHZ/26k；勿与 [R15] COValid 混用。

---

## 三、PDB 面板核对

| PDB | 蛋白 | 配体 | DFG | 文献 | 项目角色 |
|-----|------|------|-----|------|----------|
| 8ELC | JNK2 | YL2056 (Y56) | in | [R3]/[R1] | **共价主模板** |
| 3NPC | JNK2 | BIRB796 | out | [R7] | 共价初筛排除 |
| 4WHZ | JNK3 | 26k (3NL) | in | [R7b] | 可逆选择性参照 |
| 3V6S | JNK3 | JNK-IN-7 | in | [R2] | 共价叠合参照 |
| 7N8T | JNK2 | AMP | in | [R4] | 活化环参照 |
| 3V6R | JNK3 | JNK-IN-2 | in | [R2]/[R5] | IN-2 共晶（非 3V6S） |

---

## 四、Warhead 类型（勿混淆）

| 分子 | Warhead | 说明 |
|------|---------|------|
| YL5084, YL2056, JNK-IN-8, 56d | **丙烯酰胺** | Michael 受体；需 bondedAtomPairs |
| YL5084R, 56a, JNK-IN-6 | **饱和酰胺 / 丙酰胺** | 共价阴性对照 |
| BIRB796, 26k | **无共价弹头** | 可逆 ATP 竞争 |

JNK 系列 **无 chloroacetamide** 文献报道；勿在报告中写错。

---

## 五、ChEMBL 交叉验证摘要

```
CHEMBL5398852  YL5084   MW 600.73  ✅
CHEMBL5413843  YL2056   MW 586.70  ✅
CHEMBL2216824  JNK-IN-8 MW 507.60  ✅
CHEMBL5436581  YL5084R  MW 602.74  ✅
CHEMBL6167401  56d      MW 516.56  ✅
CHEMBL6151222  56a      MW 518.58  ✅
PubChem 57340684 JNK-IN-6 MW 438.49 ✅
PubChem 57340686 JNK-IN-8 InChIKey 与 seed 一致（E 型丙烯酰胺）✅
```

---

## 六、仍建议人工复核项

1. **56a/56d**：ChEMBL6151222 与 Wydra 2025 正文编号 56a 的对应关系建议在 SI Table 中再对一次化合物 ID。
2. **AF3 输入表**（用户本地 CSV）：YL2056 行的 `bonded_atom_pairs` / `use_bonded_atom_pairs` 需在跑 AF3 前修正。
3. **26k 2D 图**：已改用 4WHZ/3NL SMILES；若需无立体标绘，使用 `[C@H]5CCNC5` 形式。

---

*维护：结构变更后同步更新 `phase0_compounds_seed.csv`、`STRUCTURE_AUDIT.md` 与本文件。*
