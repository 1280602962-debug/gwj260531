# Phase 0 分子清单与用途说明

完整结构核对见 → [`STRUCTURE_AUDIT.md`](STRUCTURE_AUDIT.md)  
文献 / PDB 全量审计 → [`COMPOUND_LITERATURE_AUDIT.md`](COMPOUND_LITERATURE_AUDIT.md)

## 已验证分子（7/7）

| ID | 分子 | 状态 | 用途 |
|----|------|------|------|
| ACT001 | YL5084 | ✅ | Phase 0 主阳性 + decoy 锚点 |
| ACT002 | YL2056 | ✅ | 8ELC 共晶系列；AF3 用 acrylamide |
| ACT003 | JNK-IN-8 | ✅ | Pan-JNK 共价阳性（丙烯酰胺） |
| ACT004 | 56d | ✅ | Ligand-first JNK2/3>>JNK1 阳性 |
| NEG002 | YL5084R | ✅ | 共价阴性（butanamide） |
| NEG001 | JNK-IN-6 | ✅ | 丙酰胺阴性（JNK-IN-7 系列） |
| NEG003 | 56a | ✅ | 丙酰胺阴性（56d 配对） |

## 重要更正（2026-07-10）

1. **56a**：已补 CHEMBL6151222（propionamide，MW 518.58）。
2. **文献**：YL5084/YL2056 → Lu 2023 [R1]；JNK-IN 系列 → Zhang 2012 [R2]；56d/56a → Wydra 2025 [R20]。
3. **Warhead**：JNK 共价 hit 均为 **丙烯酰胺**，非 chloroacetamide。
4. **8ELC**：共晶配体为 **YL2056**（PDB Y56 为 post-covalent butanamide）。

## Warhead 与 AF3

| 分子 | Warhead | bondedAtomPairs |
|------|---------|-----------------|
| YL5084, YL2056, JNK-IN-8, 56d | 丙烯酰胺 | 需要 |
| YL5084R, JNK-IN-6, 56a | 饱和/丙酰胺 | 不需要 |

## 参考

- Lu 2023 [R1]: [10.1021/acs.jmedchem.2c01834](https://doi.org/10.1021/acs.jmedchem.2c01834)
- Zhang 2012 [R2]: [10.1016/j.chembiol.2011.11.010](https://doi.org/10.1016/j.chembiol.2011.11.010)
- Wydra 2025 [R20]: [10.1021/acs.jmedchem.5c00884](https://doi.org/10.1021/acs.jmedchem.5c00884)
