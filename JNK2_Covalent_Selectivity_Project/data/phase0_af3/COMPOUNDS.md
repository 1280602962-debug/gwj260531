# Phase 0 分子清单与用途说明

完整结构核对见 → [`STRUCTURE_AUDIT.md`](STRUCTURE_AUDIT.md)

## 已验证分子

| ID | 分子 | 状态 | 用途 |
|----|------|------|------|
| ACT001 | YL5084 | ✅ | Phase 0 主阳性 + decoy 锚点 |
| ACT002 | YL2056 | ✅ | 8ELC 共晶系列；**AF3 用 acrylamide**；PDB 中为 post-covalent |
| ACT003 | JNK-IN-8 | ✅ | Pan-JNK 共价阳性 |
| ACT004 | 56d | ✅ 已修正 | Ligand-first JNK2/3>>JNK1 阳性 |
| NEG002 | YL5084R | ✅ | 共价阴性（butanamide） |
| NEG001 | JNK-IN-6 | ✅ 已补 | 丙酰胺阴性（JNK-IN-7 系列） |

## 仍待补

| 分子 | 用途 |
|------|------|
| **56a** | 56 系列丙酰胺阴性 |

## 重要更正（2026-07-09）

1. **56d**：原 CHEMBL5947460 完全错误，已替换为 Wydra 2025 结构（MW 516.56）。
2. **JNK-IN-6**：已补 PubChem 57340684；**不是** JNK-IN-8 的 SMILES。
3. **YL2056 vs 8ELC**：ChEMBL 丙烯酰胺正确；PDB 配体 Y56 登记为 butanamide（共价加成后 warhead 区），RMSD 请比骨架。

## Warhead 与 AF3

| 分子 | Warhead | bondedAtomPairs |
|------|---------|-----------------|
| YL5084, YL2056, JNK-IN-8, 56d | 丙烯酰胺 | 需要 |
| YL5084R, JNK-IN-6 | 饱和/丙酰胺 | 不需要 |
| 56a | 丙酰胺 | 不需要（待补） |

## 参考

- Lu 2023: [10.1021/acs.jmedchem.2c01834](https://doi.org/10.1021/acs.jmedchem.2c01834)
- Zhang 2012: [10.1016/j.chembiol.2011.11.010](https://doi.org/10.1016/j.chembiol.2011.11.010)
- Wydra 2025: [10.1021/acs.jmedchem.5c00884](https://doi.org/10.1021/acs.jmedchem.5c00884)
