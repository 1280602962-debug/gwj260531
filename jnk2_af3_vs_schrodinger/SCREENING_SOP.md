# JNK2 Cys116 共价筛选 SOP（冻结版 v0.1）

**状态：** Phase 0 校准通过后冻结主规则  
**模板：** PDB 8ELC，Cys116（激酶域 FASTA 输入时为残基 113）  
**依据：** `AF3_vs_薛定谔共价对接对比总结.md`

## 1. 角色分工

| 工具 | 角色 | 主指标 |
|------|------|--------|
| **AF3 共价共折叠** | **主门控 / 主排序** | **mPAE**（protein–ligand PAE min） |
| **薛定谔 Covalent Docking** | **辅助 QC** | 共价几何 + 相互作用目视/规则；GScore **不单独决策** |
| WSL/RDKit | 建库、过滤、汇总 | SMARTS / 物化 / 可采购标记 |

## 2. Gate（Phase 0 已验证）

回顾性标准（已达到）：
- mPAE AUC = 1.0；EF@1% = 36.8；decoy 压过 active = 0%

前瞻性使用时建议：
1. 按 mPAE **升序**排序（越小越好）  
2. 初筛阈值：优于历史 decoy 中位数，或试点 **mPAE ≤ 1.05–1.20 Å**（需在 Layer1 上再校准）  
3. **禁止**单独用 ipTM / ranking_score 作 gate  
4. Glide：仅对 AF3 TopN 做 pose QC；**56d 邻域禁止 Glide-only 决策**

## 3. 漏斗

```
Layer1 库（丙烯酰胺）
  → AF3 @ Cys116（bondedAtomPairs 正确）
  → mPAE 排序 + gate
  → Top N → Glide covalent QC（可选交集）
  → 人工检查 Leu106/warhead 几何
  → 可采购/合成清单 → 湿实验
```

## 4. 阴性与对照

- 饱和酰胺阴性（YL5084R、JNK-IN-6、56a）不得设共价键约束  
- Pan 对照：JNK-IN-8  
- 选择性锚点：YL5084 / YL2056；ligand-first：56d（AF3 优先）

## 5. 变更控制

改主指标、改受体（非 8ELC）、启用 DFG-out，必须重做 Phase 0 级校准，不得口头修改。
