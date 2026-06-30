# URAT1 三态对接方法说明（事实校正版）

> **用途**：Teacher M-CPDL、$S_{\text{trap}}$、MASFL Stage 2 Glide 对接的 **唯一权威** PDB 映射。  
> **配置源**：`config/docking_ensemble.yaml` → `urat1_ensemble.three_state_primary`

---

## 1. 核心结论

| 问题 | 答案 |
|------|------|
| 用几个蛋白文件？ | **3 个独立 PDB**（同一条人 URAT1，三种构象） |
| 薛定谔能从一个 9JDZ 自动分态吗？ | **不能** — 需手动导入 3 个 PDB、建 3 个 Grid |
| 9JDZ 是 occluded/outward 吗？ | **否** — RCSB 标注为 *lesinurad-bound inward-open* |
| occluded / outward 用谁？ | **9B1K** / **9B1L**（Dai *Cell Res* 2024，RCSB 标题明确） |
| inward 抑制剂对接用谁？ | **9DKB**（lesinurad，2.55 Å，Fedor/Suo 2025） |

---

## 2. 三态主映射（对接执行表）

| 符号 | 态 | PDB | 共晶配体 | EMDB | 文献 |
|------|-----|-----|----------|------|------|
| $\pi_{\text{in}}$ | inward-open | **9DKB** | lesinurad | EMD-46950 | Suo *Nat Commun* 2025 |
| $\pi_{\text{occ}}$ | occluded | **9B1K** | urate | EMD-44082 | Dai *Cell Res* 2024 |
| $\pi_{\text{out}}$ | outward-open | **9B1L** | urate | EMD-44083 | Dai *Cell Res* 2024 |

备用 inward：**9B1H**（lesinurad）、**9B1J**（urate inward）— 敏感性分析，非默认 grid。

---

## 3. 为何修正 9JDZ

### 3.1 RCSB 事实

- **9JDZ**：*Human URAT1 bound to lesinurad* — inward-open，EMDB-61402  
- **9B1K**：*Urate bound human URAT1 in the **occluded state***  
- **9B1L**：*Urate bound human URAT1 in the **outward-facing state***

### 3.2 Wu *Cell Discov* 2025

该文 cryo-EM 解析了 urate 的 inward / outward / occluded，但 PDB 公开条目仅包括药物 inward（9JDZ/9JDY/9JE0/9JE1）和单一 urate 条目 **9JDV**。**Outward（4.1 Å）与 occluded（4.7 Å）未作为独立 PDB 编号公开入库**（截至 2025-10 RCSB 检索）。

因此不能从 9JDZ「分离 EM map」代替已标注的 9B1K/9B1L，除非自行从作者 supplementary 重建（非默认流程）。

### 3.3 机制一致性

URAT1 为 **alternating access** 转运体（Dai 2024 Fig. 2–3；Wu 2025 Fig. 6）：

1. **Outward-open**：胞外腔开放，底物进入  
2. **Occluded**：胞外门关闭（TM7 旋向 TM1），腔两侧不通  
3. **Inward-open**：胞内释放；**抑制剂**稳定 inward/occluded，阻断循环  

对接假说：强效 URAT1 抑制剂应在 **inward + occluded** 得分相对高，**outward** 相对低 → $S_\pi = \pi_{\text{in}}+\pi_{\text{occ}}-\pi_{\text{out}}$。

---

## 4. Schrödinger / Glide 标准流程

```
1. 下载 9DKB.pdb, 9B1K.pdb, 9B1L.pdb
2. 对每个 PDB：
   Protein Prep Wizard → 去配体 → pH 7.4 → 优化
   Grid：以共晶配体（lesinurad 或 urate）为中心，22 Å 立方体
3. 对 distill_manifest.csv（8973 条）：
   Glide SP → XP × 3 grids
4. 每分子记录：GlideScore_in, GlideScore_occ, GlideScore_out
5. pi_T(s|x) = softmax_s(score_s)；outward 在 S_pi 中取负
```

**Grid 中心**：occluded/outward 用 **urate** 位点（9B1K/9B1L）；与 inward 抑制剂口袋空间重叠（Dai 2024：三态 urate 结合位点相似）。

---

## 5. Teacher Gate 质控（对接后必做）

| 检验 | 标准 | 结构 |
|------|------|------|
| lesinurad redock | RMSD ≤ 2.0 Å | 对 **9DKB** grid |
| 四药方向 | inward $S_\pi$ > outward | lesinurad, benzbromarone, verinurad, dotinurad |
| 子集 D 负样本 | median($\pi_{\text{in}}+\pi_{\text{occ}}$) < 活性集 | distill subset D (8000) |

任一失败 → 停用 CPDL Teacher 标签，回退 v2 手工 $S_{\text{trap}}$。

---

## 6. 结构叠合验证（可选，PyMOL）

```python
# 伪代码
align 9B1L, 9B1K, object=9DKB and chain A and resi 1-290  # NTD 参考
```

预期：9B1L 胞外开口最大；9B1K TM7 靠近 TM1 腔变窄；9DKB 胞内开口、ECD 可见。

---

## 7. 与项目其他模块关系

| 模块 | 使用的三态 PDB |
|------|----------------|
| Teacher M-CPDL (Stage 2) | 9DKB + 9B1K + 9B1L |
| PC-Student 口袋图 $E_s$ | 同上三结构预计算 graph |
| $S_{\text{trap}}$ / Path A L2 | 同上 |
| Benchmark redock | 药物特异性 PDB（9DKA benzbromarone 等） |

---

## 8. 参考文献（PDB 出处）

1. Dai Y, Lee CH. Transport mechanism and structural pharmacology of human urate transporter URAT1. *Cell Res* 2024. doi:[10.1038/s41422-024-01023-1](https://doi.org/10.1038/s41422-024-01023-1) — **9B1K, 9B1L, 9B1H, 9B1J**  
2. Suo Y et al. Structural pharmacology of human urate transporter URAT1. *Nat Commun* 2025. doi:[10.1038/s41467-025-60480-3](https://doi.org/10.1038/s41467-025-60480-3) — **9DKB**  
3. Wu C et al. Molecular mechanisms of urate transport by native human URAT1. *Cell Discov* 2025. doi:[10.1038/s41421-025-00779-z](https://doi.org/10.1038/s41421-025-00779-z) — **9JDZ** (inward lesinurad only)
