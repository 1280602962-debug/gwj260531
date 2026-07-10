# Phase 0 AF3 共价富集测试 — 数据包

本目录存放 **AF3 @ JNK2 Cys116** 回顾性 gate（COValid 式）可直接使用的种子数据。

**分子用途详见 → [`COMPOUNDS.md`](COMPOUNDS.md)**

## 已推送文件

| 文件 | 内容 |
|------|------|
| `jnk2_kinase_domain_4-364.fasta` | AF3 蛋白输入（筛选方案 §2.1）；**Cys116（UniProt）= 第 113 位** |
| `jnk2_uniprot_P45984_full.fasta` | 全长 UniProt P45984（备用） |
| `8ELC.pdb` | Glide 共价对接 / AF3 pose RMSD 对照受体 |
| `phase0_compounds_seed.csv` | 分子种子表（含 SMILES、用途 `role`、AF3 是否需共价键） |
| `COMPOUNDS.md` | **各分子用途中文说明**（阳性/阴性/decoy 锚点） |
| `STRUCTURE_AUDIT.md` | **结构核对报告**（ChEMBL/PubChem/PDB 交叉验证） |
| `COMPOUND_LITERATURE_AUDIT.md` | **化合物 + 文献 + PDB 全量审计** |
| `pdb_references.csv` | 结构模板索引 |
| `af3_input_template.json` | Cys 编号、`bondedAtomPairs`、gate 阈值 |
| `decoy_property_targets.yaml` | COValid 6 性质 + Tc 过滤规则 |

## 已获取 SMILES 的分子（6 个，2026-07-09 核对）

| 分子 | 角色 | 状态 |
|------|------|------|
| YL5084 | Phase 0 主阳性 + decoy 锚点 | ✅ confirmed |
| YL2056 | 8ELC 共晶系列（见 STRUCTURE_AUDIT） | ✅ confirmed |
| JNK-IN-8 | Pan-JNK 共价阳性 | ✅ confirmed |
| 56d | Ligand-first 共价阳性 | ✅ corrected (was wrong ChEMBL) |
| YL5084R | 共价机制阴性 | ✅ confirmed |
| JNK-IN-6 | 丙酰胺阴性 | ✅ added (PubChem 57340684) |

## 仍需准备

| 类别 | 数量 | 状态 |
|------|------|------|
| Property-matched decoy | 50 × YL5084（Phase 0 最小） | ❌ 待生成 |
| **56a** 阴性 SMILES | 1 | ❌ 待 Wydra SI |
| AF3 运行 | ~57 jobs（5 分子 + 50 decoy） | ❌ 需算力 |
| Glide baseline | 同上 | ❌ 需 Schrödinger |

## Gate 通过标准

- YL5084 **mPAE** < decoy 中位数
- **EF@1% ≥ 2**
- YL5084 AF3 pose vs **8ELC** 重原子 RMSD **< 2 Å**

## Cys 编号

| 构建方式 | Cys（SG）残基编号 |
|----------|------------------|
| UniProt / PDB 8ELC 全链 | **116** |
| 激酶域 4–364 FASTA（推荐 AF3 输入） | **113** |

`bondedAtomPairs` 必须与所选序列编号一致。
