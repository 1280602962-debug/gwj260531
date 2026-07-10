# JNK2 选择性共价抑制剂筛选项目

**项目代号：** JNK2-CovSel  
**靶点：** 人源 JNK2（MAPK9, UniProt P45984）  
**策略：** 靶向 **Cys116** 的共价抑制剂，优先 **JNK2 > JNK1**（可兼顾 JNK3 评估）  
**创建日期：** 2026-07

---

## 目录结构

```
JNK2_Covalent_Selectivity_Project/
├── README.md                              # 本文件
├── JNK2选择性共价抑制剂筛选方案.md          # 主方案（计算 + 实验 + 里程碑）
├── JNK2项目决策五阶段.md                    # 五阶段决策脉络（含 DFG-out 不作共价主筛理由）
├── JNK2共价抑制剂文献调研综述.md            # JNK2 抑制剂文献复盘（共价 + 可逆，v2.0, 38 refs）
├── JNK2_covalent_inhibitors_master_table.csv # JNK2 共价抑制剂主表（对齐 NLRP3 格式，12 条）
├── JNK2_covalent_inhibitors_README.md        # 主表字段说明 / Cys116 vs Cys154 编号
├── REFERENCES.md                          # 参考文献（含 DOI / PDB）
├── config/
│   └── targets.yaml                       # 结构模板与 benchmark 配置
└── data/phase0_af3/                       # AF3 Phase 0 种子数据（见 COMPOUNDS.md 分子用途）
    ├── COMPOUNDS.md
    ├── phase0_compounds_seed.csv
    └── 8ELC.pdb
```

## 科学定位

本项目与仓库内 **JNK1 非共价选择性筛选**（`JNK1_Selectivity_Project/`）互补：

| 维度 | JNK1 非共价项目 | 本课题（JNK2 共价） |
|------|----------------|---------------------|
| 结合模式 | 可逆 ATP 竞争 | **不可逆 Cys116 共价** |
| 选择性来源 | Gly87 位阻等（计算 Δsel 已证伪） | **kinact/KI + 非共价预定位**（Leu106 等） |
| 结构模板 | 3ELJ / 4L7F / 3E7O 等 | **8ELC（DFG-in 共晶）为主** |
| 先导路径 | ML + Glide 漏斗 | **JNK-IN-8 → YL5084 骨架 MedChem + 可选 AF3** |

## 关键文献入口

- **共价 JNK2 选择性先导：** Lu et al., *J. Med. Chem.* 2023 — YL5084（doi:[10.1021/acs.jmedchem.2c01834](https://doi.org/10.1021/acs.jmedchem.2c01834)）
- **Pan-JNK 共价起点：** Zhang et al., *Chem. Biol.* 2012 — JNK-IN-8（doi:[10.1016/j.chembiol.2011.11.010](https://doi.org/10.1016/j.chembiol.2011.11.010)）
- **共晶结构：** PDB [8ELC](https://www.rcsb.org/structure/8ELC)（JNK2–YL2056，DFG-in）

## 快速开始

0. 阅读 `JNK2项目决策五阶段.md` 了解从文献到 **DFG-in 共价 + 禁用 DFG-out 主筛** 的完整推理  
1. 阅读 `JNK2选择性共价抑制剂筛选方案.md` §一–§三 明确科学目标与 gate 标准  
2. 阅读 `JNK2共价抑制剂文献调研综述.md` 了解 **共价 + 可逆** JNK2 抑制剂的发现路径、活性数据与湿实验范式  
3. 使用 `data/phase0_af3/` 中的 **8ELC**、FASTA 与分子种子（各分子用途见 `COMPOUNDS.md`）  
4. Phase 0：AF3 gate（YL5084 mPAE vs decoy；EF@1%≥2）通过后，再扩大库筛选  

## 关联分支

- 非共价 IFP/MM-GBSA 工作流：`origin/cursor/jnk-step1-step2-workflow-cd8a`（**不含**共价配体，仅供交叉验证）
