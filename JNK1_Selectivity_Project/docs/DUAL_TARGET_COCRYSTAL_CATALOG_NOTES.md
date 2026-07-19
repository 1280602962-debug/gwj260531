# 双靶共晶结构数据收集与分类（v0.1）

数据路径：`data/benchmarks/dual_target_structures/`

## 结论

可以、也应该收集已发表双靶共晶来评估对接方法；与自有 PROTAC 活性集互补：

- **共晶集** → 姿态是否对（尤其 Tier A 两端 RMSD）
- **PROTAC 集** → linked 双功能有效/无效分辨

## 分类结果摘要

### 按结构完整度

| Tier | 含义 | 本版条目 |
|------|------|----------|
| **A_both_ends** | 同一双靶配体两端都有 PDB | DT-A-001–005（Mcl-1/Bcl-xL；LpxA/LpxD；PknA/PknB；EGFR/HER2 TAK-285） |
| **B_single_end** | 双靶分子但目前只结晶一端 | BET–HDAC 系列；ER/CA 调节剂；MurD/MurE；PD-L1/VISTA 等 |
| **C_series_related** | 同系列不同配体分别结晶在两靶 | GyrB/ParE 吡咯并嘧啶系列 |
| **D_claim_only** | 仅设计谱系，不作双靶姿态金标准 | Mcl-1 前体 3WIX |

### 按设计类型（Morphy）

| 类型 | 本版代表 |
|------|----------|
| **linked** | Mcl-1/Bcl-xL compound 10；部分 BET–HDAC hydroxamate hybrid |
| **merged** | LpxA/LpxD；PknA/PknB；TAK-285；多数 BET–HDAC merged 系列 |
| **fused** | 本种子集暂少，待扩充 |

## 评估时怎么用

1. **协议验证（必做）**：Tier A 两端自对接 RMSD  
2. **压力测试**：Tier B 单端 RMSD + Tier C 交叉对接  
3. **功能验证**：另接双靶活性回顾集 + 你们 PROTAC 标签  
4. **分类报告**：按 linked/merged 与 Tier 分别统计成功率  

## 重要提醒

- PDB 标题含 “dual” 噪声极大（双特异性磷酸酶、单靶 dual binding mode 等），必须人工核对。  
- 许多双靶药化文**没有**两端共晶；Tier A 天然偏少，属正常。  
- 本表为 **种子集 v0.1**，可按 README 中的 expansion protocol 继续追加。
