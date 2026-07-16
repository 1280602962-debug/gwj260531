# JNK2 Cys116 共价筛选 SOP（冻结版 v0.2）

**状态：** Phase 0 校准通过后冻结主指标；建库与算力分层见完善漏斗  
**模板：** PDB 8ELC，Cys116（激酶域 FASTA 输入时为残基 113）  
**依据：** `AF3_vs_薛定谔共价对接对比总结.md` · `LIBRARY_FUNNEL_COMPLETE.md`

## 1. 角色分工

| 工具 | 角色 | 主指标 |
|------|------|--------|
| **AF3 共价共折叠** | **最终排序 / 门控** | **mPAE** |
| **薛定谔 Covalent Docking** | **高通量粗筛 + pose QC** | 几何合格优先；GScore **不单独定采购** |
| WSL/RDKit | 建库、胺→丙烯酰胺枚举、合成性/双轨过滤 | SMARTS / SA / Tc / 可采购标记 |

## 2. Gate（Phase 0 已验证）

回顾性标准（已达到）：
- mPAE AUC = 1.0；EF@1% = 36.8；decoy 压过 active = 0%

前瞻性使用时：
1. 按 mPAE **升序**排序  
2. 试点阈值：**mPAE ≤ 1.05–1.20 Å** 或优于校准 decoy 分位（首轮后回调）  
3. **禁止**单独用 ipTM / ranking_score 作 gate  
4. Glide 粗筛阈值必须 **放宽**；**56d / Track-Novel 提高 AF3 配额**

## 3. 完整漏斗（v0.2）

```
Stage1  双源建库
        A. 商业丙烯酰胺（千–万，可采购快车道）
        B. ChEMBL/ZINC 类药 → 胺 → 装丙烯酰胺（虚拟扩展）
          ↓
Stage2  统一 QC：单弹头 · PAINS/过反应 · MW/cLogP · 合成性/可买
          ↓
Stage3  双轨分仓：Track-Sim（YL5084/56d 邻域）| Track-Novel（低相似新骨架）
          ↓
Stage4  薛定谔共价对接粗筛：5k–50k → 500–2000（松阈值）
          ↓
Stage5  AF3 mPAE 精排与门控（服务器）
          ↓
Stage6  AF3∩几何QC · 人工面板
          ↓
Stage7  purchase_list + synthesize_list → kinact/KI · C116S
```

细节与风险表：`LIBRARY_FUNNEL_COMPLETE.md`。

## 4. 阴性与对照

- 饱和酰胺阴性不得设共价键约束  
- Pan：JNK-IN-8；选择性锚点：YL5084 / YL2056；ligand-first：56d（AF3 优先）

## 5. 变更控制

改主指标、改受体（非 8ELC）、启用 DFG-out，必须重做 Phase 0 级校准。
