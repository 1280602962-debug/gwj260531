# Research Dossier — Scene Analyst

**Scene:** `report_review`（技术报告 / 开题–立项评估，非投稿全文）  
**Target:** JNK2 Cys116 选择性共价抑制剂立项与推进评估  
**Date:** 2026-07-16  
**Channel:** local_first + WebSearch（无 PubMed/Scholar MCP）

## Venue Requirements（本任务类型）

本任务不是单一期刊投稿稿，而是 **贡献导向的立项/推进评估报告**。审稿人视角应模拟：

1. **Medicinal Chemistry / Chemical Biology 审稿人**：claim 是否超过证据（尤其疾病模型与虚筛发现）。
2. **结构生物学 / CADD 审稿人**：模板选择、共价几何、校准 vs 发现是否混淆。
3. **生物学审稿人**：JNK2 适应症叙事是否可辩护。

可接受交付形态：技术报告、开题报告核心章节、未来 J. Med. Chem. / ACS Med. Chem. Lett. / RSC MedChem 类文章的 **贡献骨架**。

## Review Criteria（硬标准）

| 标准 | 通过条件 |
|------|----------|
| Contribution-first | 有可证伪的主贡献句，而非“做了很多计算” |
| Results-as-validation | 每个结果小节对应一条贡献承诺 |
| Claim boundary | AF3 校准、MedChem 类似物、疾病模型分层声明 |
| Reproducibility | 8ELC、Cys 编号、种子库、decoy 规则可复现 |
| Honesty on biology | 不把 pan-JNK / CC-90001 故事写成 JNK2 独有适应症 |

## Accepted Paper Patterns（近三年对标）

1. **Lu 2023 (JMC)**：从 pan-JNK 共价起点 → 动力学选择性 → 共晶 → **明确承认细胞抗增殖不完全依赖 JNK2**。贡献类型 = new chemical matter + isoform kinetic selectivity。
2. **Wydra 2025 (JMC)**：ligand-first 可逆选择性 → 再装 warhead → kinact/KI + kinome。贡献类型 = new scaffold path to JNK2/3 covalent probes。
3. **Shamir/London 2025–2026 (JACS / bioRxiv COValid)**：共价虚筛基准 + AF3 mPAE。贡献类型 = new analysis/benchmark/method。

本课题若投方法学或探针文章，应 **显式选择** 上述三类之一为主贡献，避免三者平均用力却都站不住。

## Constraints for This Paper / Project

1. **不得**以 DSS/IBD 或纤维化临床故事作为主贡献。
2. **不得**把 Phase 0 AF3 回顾性 AUC/EF 写成前瞻性 hit discovery。
3. **必须**保留 8ELC DFG-in / 禁用 3NPC 共价主筛的决策。
4. WSL 路径必须与“可选 AF3/Schrödinger”分层，避免不可复现的算力依赖。
5. 与仓库内 JNK1 非共价项目结论衔接：**Δsel 不作采购 gate**。
