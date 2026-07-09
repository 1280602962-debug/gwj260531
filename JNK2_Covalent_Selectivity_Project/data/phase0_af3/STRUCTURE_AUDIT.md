# Phase 0 分子结构核对报告

核对日期：2026-07-09  
方法：ChEMBL / PubChem / PDB 8ELC Y56 / Lu 2023 & Zhang 2012 SI / RDKit InChIKey 比对

## 汇总

| ID | 分子 | 结论 | 行动 |
|----|------|------|------|
| ACT001 | YL5084 | ✅ 正确 | 保持 |
| ACT002 | YL2056 | ✅ 正确（游离丙烯酰胺） | 补充 8ELC 共价加合说明 |
| ACT003 | JNK-IN-8 | ✅ 正确 | 保持 |
| ACT004 | 56d | ❌ **原库错误** | 已换为用户/Wydra 结构 |
| NEG002 | YL5084R | ✅ 正确 | 保持 |
| NEG001 | JNK-IN-6 | ⚠️ 原为空 | 已补 PubChem 57340684 |
| NEG003 | 56a | ❌ 仍缺 | 待 Wydra 2025 SI |

---

## 逐条详情

### ACT001 YL5084 — ✅ 正确

- **库中 SMILES** 与 ChEMBL5398852 一致
- PubChem 153610906 IsomericSMILES 与 ChEMBL **InChIKey 相同**（写法不同）
- MW 600.73；`(3S,4S)`-3-methyl 吡咯烷 + `(E)`-4-(dimethylamino)but-2-enamide，与 Lu 2023 SI 一致

### ACT002 YL2056 — ✅ 正确（注意与 8ELC 配体形式）

- **ChEMBL5413843** 与 Lu SI 命名一致：  
  `(S,E)-4-(dimethylamino)-N-[4-(3-{[4-(2-phenylpyrazolo[1,5-a]pyridin-3-yl)pyrimidin-2-yl]amino}pyrrolidine-1-carbonyl)phenyl]but-2-enamide`
- MW 586.70；含 **丙烯酰胺** warhead

**8ELC 注意：** PDB 配体 Y56 在 CCD 中登记为 **butanamide**（饱和），化学上接近 **YL2056R** 或 **Michael 加成后的 warhead 区**，并非游离 acrylamide SMILES：

```
8ELC Y56 (PDB): CN(C)CCCC(=O)N...  → MW≈588, 无 C=C
YL2056 (ChEMBL): CN(C)C/C=C/C(=O)N... → MW≈587, 有 acrylamide
```

Lu 2023 正文：8ELC 为 YL2056 浸泡共晶，**Cys116 连续密度证实共价键**。因此：

- **AF3 / Glide 输入**：用 ChEMBL **acrylamide** YL2056 + `bondedAtomPairs`
- **RMSD vs 8ELC**：比对 **骨架**（吡咯烷–苯甲酰胺–杂环），或从 8ELC 提取配体坐标并理解 warhead 为 **post-covalent** 构象

### ACT003 JNK-IN-8 — ✅ 正确

- ChEMBL2216824 = PubChem 57340686 = Chemical Probes Portal
- SMILES: `Cc1cc(NC(=O)c2cccc(NC(=O)/C=C/CN(C)C)c2)ccc1Nc1nccc(-c2cccnc2)n1`
- MW 507.60；相对 JNK-IN-7 多 **flag methyl**（3-甲基苯胺）

### ACT004 56d — ❌ 原库错误 → 已修正

| | 原库 (CHEMBL5947460) | 正确 (Wydra 2025 / 用户提供) |
|--|---------------------|------------------------------|
| SMILES | `Nc1ncnc2c1ccn2CC(=O)N(CC(=O)Nc1cccc(Br)n1)C1CC1` | `C=CC(=O)Nc1cccc(C(=O)Nc2cccc(-n3cc(NC(=O)Nc4cccc5ccccc45)cn3)c2)c1` |
| MW | 444 | **516.56** |
| 特征 | 嘌呤+溴吡啶 | 氨基吡唑+萘基酰胺+meta-丙烯酰胺 |
| Tc vs 正确结构 | 0.17 | — |

**CHEMBL5947460 为错误关联**，非 Wydra 56d。IPMS Mr≈516 与正确 MW 一致。

### NEG002 YL5084R — ✅ 正确

- ChEMBL5436581 与 Lu SI 一致：**butanamide**（饱和 warhead）
- 命名：`4-(dimethylamino)-N-[4-((3S,4S)-3-methyl-4-...pyrrolidine...)phenyl]butanamide`
- 与 YL5084 配对；AF3 **不设** bondedAtomPairs

### NEG001 JNK-IN-6 — ✅ 已补（原为空）

- PubChem **57340684**（Zhang 2012 系列）
- SMILES: `CCC(=O)Nc1cccc(C(=O)Nc2ccc(Nc3nccc(-c4cccnc4)n3)cc2)c1`
- **propionamide** 替换 JNK-IN-7 的 acrylamide；~100× 生化失活
- 注意：基于 **JNK-IN-5/7 scaffold**，**非** JNK-IN-8（无 flag methyl）

### NEG003 56a — ❌ 仍缺

- Wydra 2025 56 系列 **丙酰胺**阴性；需 SI 补 SMILES

---

## 使用建议

1. **56d decoy 锚点**：MW≈517（非 601），需单独算 6 性质窗口，不可与 YL5084 共用 decoy 规则。
2. **8ELC RMSD**：勿直接用 acrylamide YL2056 与 PDB Y56 做 warhead 原子 RMSD；比骨架或 post-covalent 模型。
3. **ChEMBL ID**：56d 勿再引用 CHEMBL5947460。
