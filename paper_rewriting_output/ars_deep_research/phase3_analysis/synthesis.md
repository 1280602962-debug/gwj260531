# Phase 3 — Synthesis + Devil’s Advocate

## Thematic findings

### Theme 1 — 对接前过滤是标准，但“标准”≠你们的 Step0–7
文献一致支持：超大库必须先做 **物化/合成性/反应枚举过滤**（RosettaAMRLD 显式允许对接前 MW/RB/logP；Deep Docking/HASTEN 用代理减对接）。  
**缺口：** 几乎没有论文把「ECFP 对四个 JNK 锚点 + 固定 Tc 窗 + Murcko 双轨」写成可复现且经校准的 JNK2 专用 SOP。

### Theme 2 — 共价筛选的真正轴是「识别 × 反应性 × 几何」
NRDD 2022、DOCKovalent、Nat Commun 2024、extended warheads JCIM 2024 共同指向：  
丙烯酰胺枚举只解决 warhead **存在**；还需 **反应性窗口** 与 **Cys 矢量几何**。ECFP 相似不编码后两者。

### Theme 3 — 类似物富集 vs 新骨架必须用不同表征
SwissSimilarity / ErG / SHAFTS / HybridSim：ECFP 擅类似物，ErG/药效团/shape 擅 hop。  
把 Murcko 不同直接叫 Novel，会高估新颖、低估噪声。

### Theme 4 — JNK 专属证据反对“只找像 YL 的分子就够”
Lu：选择性来自非共价预组织 + 动力学，不只 warhead。  
Wydra：ligand-first 再装弹头。  
AF3/JACS：最终排序应用 mPAE；传统共价对接弱。  
→ 初筛应服务 **AF3 配额**，并保留 **ligand-first/正交骨架** 通道。

### Theme 5 — Glide 在你们项目中的正确位置
Phase0：AF3 mPAE ≫ Glide；56d 邻域 Glide 假阳高。  
文献也警告共价对接误用。  
→ Glide = **松阈值几何粗筛**，不是活性初筛主裁判。

## Contradictions resolved
- “大库要对接” vs “AF3 更好” → 分层：ligand triage → 松 Glide → AF3。  
- “要新骨架” vs “要活性” → 分仓分配额，禁止混榜。  
- “丙烯酰胺最常见” vs “warhead 精度前沿” → Phase1 主用标准丙烯酰胺；反应性分桶；后续可选变体。

## Devil’s advocate checkpoint（Checkpoint 2）
**Verdict: REVISE**（不能 PASS as-is）

必须在宣称“理性初筛”前补齐：  
校准集/阈值敏感性、warhead-stripped 相似、反应性与选择性控制、可证伪成功标准、明确 Glide 非主排序。

详见并行 `devils_advocate` 意见记录于规划正文。
