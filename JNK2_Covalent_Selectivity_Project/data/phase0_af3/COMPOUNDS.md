# Phase 0 分子清单与用途说明

本文件说明 `phase0_compounds_seed.csv` 中**已获取 SMILES** 的分子在项目中的角色。  
用途与 `JNK2选择性共价抑制剂筛选方案.md` §5.3（AF3 COValid 式 gate）及 Tier-0 对照一致。

## 已获取分子（可直接用于 AF3 / Glide 输入）

| ID | 分子 | 标签 | 用途（在本项目中的作用） | 数据来源 | 备注 |
|----|------|------|--------------------------|----------|------|
| **ACT001** | **YL5084** | 阳性 · decoy 锚点 | **Phase 0 主阳性对照**：严格 JNK2> JNK1 共价 lead；AF3 gate 中 mPAE 须优于 decoy 中位数；decoy 6 性质匹配的参照分子 | ChEMBL5398852 / PubChem 153610906 | kinact/KI(JNK2/JNK1)≈21×；无独立 PDB，pose 参照 8ELC |
| **ACT002** | **YL2056** | 阳性 | **8ELC 共晶配体**：JNK2 Cys116 共价几何与 Leu106 选择性结构证据；用于 AF3 pose vs 8ELC RMSD 对照（<2 Å） | ChEMBL5413843 | 与 YL5084 同系，无 flag methyl |
| **ACT003** | **JNK-IN-8** | 阳性 | **Pan-JNK 共价阳性**：验证 AF3 对经典丙烯酰胺共价系列的识别；可选第二 decoy 锚点（50 decoys） | ChEMBL2216824 / PubChem 57340686 | 三亚型近等，**不提供** JNK2 选择性信息 |
| **ACT004** | **56d** | 阳性 | **Ligand-first 共价阳性**：JNK2/3>>JNK1 第二独立 scaffold；验证 AF3 对氨基吡唑+meta-acrylamide 系列的富集 | ChEMBL5947460 | SMILES 待 Wydra 2025 SI 复核 |
| **NEG002** | **YL5084R** | 阴性 | **共价机制阴性对照**：饱和/还原 warhead，无法 Michael 加成 Cys116；用于分离共价 vs 非共价贡献（对应 Lu 2023 YL5084R） | ChEMBL5436581 | 应与 YL5084 配对使用；AF3 中**不设** bondedAtomPairs |

## 待补分子（占位，暂无可靠 SMILES）

| ID | 分子 | 预期用途 | 计划来源 |
|----|------|----------|----------|
| **NEG001** | **JNK-IN-6** | 丙酰胺阴性：JNK-IN-5 系列共价必要性对照（~100× 生化失活） | Zhang 2012 Chem Biol SI |
| **—** | **56a** | 56 系列丙酰胺阴性（与 56d 配对） | Wydra 2025 SI |

## 按实验环节的使用方式

```
AF3 Phase 0 gate
├── 必跑阳性：ACT001 (YL5084)          ← gate 排名核心
├── 结构对照：ACT002 (YL2056)          ← 8ELC RMSD
├── 扩展阳性：ACT003, ACT004           ← 跨 scaffold 验证
├── 阴性：NEG002 (YL5084R)             ← 无共价键约束
└── Decoy 库（待生成）：50× ACT001 性质匹配丙烯酰胺

Glide 共价 baseline（并行）
└── 受体：8ELC.pdb + 同上分子列表
```

## Warhead 与 AF3 `bondedAtomPairs`

| 分子 | Warhead | AF3 共价约束 |
|------|---------|--------------|
| YL5084, YL2056, JNK-IN-8, 56d | 丙烯酰胺 | **需要** Cys113-SG ↔ acrylamide β-C |
| YL5084R | 饱和酰胺（无 Michael 受体） | **不需要** |
| JNK-IN-6, 56a | 丙酰胺 | **不需要**（待补 SMILES） |

## 参考文献

- YL5084 / YL2056 / YL5084R：Lu et al., *J. Med. Chem.* 2023 ([R1](https://doi.org/10.1021/acs.jmedchem.2c01883))
- JNK-IN-8 / JNK-IN-6：Zhang et al., *Chem. Biol.* 2012 ([R2](https://doi.org/10.1016/j.chembiol.2011.11.010))
- 56d / 56a：Wydra et al., *J. Med. Chem.* 2025 ([R20](https://doi.org/10.1021/acs.jmedchem.5c00884))
