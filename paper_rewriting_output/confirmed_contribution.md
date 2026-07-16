# Confirmed Contribution

**Confirmed:** 2026-07-16  
**User choice:** B + 可复现框架必达  
**Scene:** report_review / 未来期刊文骨架

## Core Contribution

| Field | Content |
|---|---|
| Main contribution statement | 建立一个以 JNK2–YL2056 共晶（PDB 8ELC）与 Cys116 丙烯酰胺化学为锚、经 COValid 式回顾性校准、可在常规算力（含 WSL 开源漏斗）复现的 JNK2>JNK1 选择性共价筛选与候选推进框架，并在该框架上扩展 ligand-first / 邻域化学空间，以产出可检验的新共价候选（对标但不等同于 YL5084/56d）。 |
| Contribution type | **new system**（可复现筛选—验证框架，必达） + **new chemical matter**（新共价候选，升级达） / 兼 **new analysis-or-benchmark**（JNK2-Cys116 校准集） |
| One-sentence reviewer payoff | 审稿人得到的是「可复现的 JNK2 共价筛选框架 + 在框架内可证伪的新化学推进路径」，而不是又一次不可复核的对接截图或未归因的疾病药效。 |

## Why This Contribution Is Needed

| Field | Content |
|---|---|
| Field problem | JNK1/JNK2 功能可能相反，但 ATP 口袋 ~98% 相同，亚型选择性化学工具长期稀缺。 |
| Specific gap | （1）缺少以 8ELC/Cys116 为锚、公开可复现的校准筛选闭环；（2）YL5084/56d 之后仍需新骨架或更优属性候选；（3）疾病模型无法归因是因为缺探针，不是因为已证明 JNK2 成药。 |
| Concrete challenge | 共价几何与 Leu106 预定位耦合；kinact/KI 与 C116S 实验门槛高；AF3/商业对接算力与许可不均；易把回顾性富集误写成发现。 |
| Why prior work leaves it unresolved | Lu 2023 提供先导与共晶但非开放筛选框架；Wydra 2025 示范 ligand-first 但化学空间仍有限；COValid 证明 AF3/mPAE 通法，但非 JNK2 位点专用可复现工程包；JNK1 非共价 Δsel 已证伪，不能直接搬用。 |

## How This Paper Responds

| Field | Content |
|---|---|
| Design response | **双层贡献：** Layer-F（Framework，必达）= 种子/阴性对照 + property-matched decoy + 8ELC 模板规则 + WSL 开源漏斗（RDKit→共价对接）+ 可选 AF3/Schrödinger 精修分层；Layer-C（Chemistry，升级达）= ligand-first 邻域库与 MedChem 扩展，用同一 gate（几何/校准/可采购）推进新候选。 |
| Evidence required | （F）EF@1%/mPAE 或开源对接富集报告；Cys 编号与 bonded pairs 可复现；排除 3NPC 的论证；（C）至少一批可采购/可合成候选的对接或 AF3 排序 + 计划中的 kinact/KI、C116S 验证设计；与 YL5084/56d 的差异说明。 |
| Evidence available | 决策五阶段与主表审计；Phase 0 种子 7/7；8ELC 主模板；AF3 回顾性校准结果（用户侧）；排除 DFG-out；与 Lu/Wydra/COValid 文献对齐的方案文档。 |
| Evidence missing | 完整 decoy 库与规范化 EF 报告；Layer1 邻域库实物；自有合成与 kinact/KI；细胞 engagement；疾病模型归因数据（故意不作为主证据）。 |

## Claim Boundary

| Field | Content |
|---|---|
| Strong claims allowed | 8ELC/Cys116/丙烯酰胺路线有结构与文献依据；可复现框架的设计与校准标准可公开；Δsel 不作采购 gate；疾病叙事应改为探针驱动的开放问题。 |
| Claims to soften or avoid | 已完成前瞻性大库发现；已证明 DSS/纤维化/肿瘤中 JNK2 治疗优势；AF3 AUC=1 等于新药发现；新化学实体在湿实验前写“优效先导”。 |
| Novelty risk | “已有 YL5084/56d/COValid”——回答：贡献是 **JNK2 位点工程化可复现闭环 + 框架约束下的新化学扩展**，不是重复发表同一分子。 |
| Significance risk | “只是流程没有新分子”——回答：框架是必达贡献；新化学是升级贡献，文章/项目里程碑按 Evidence missing 分阶段披露，避免空头化学 claim。 |

## Dual-Track Milestone Mapping

| Track | Must-have deliverable | Upgrade deliverable |
|---|---|---|
| **Framework (必达)** | Layer0 校准集 + decoy；模板/编号 SOP；WSL 漏斗脚本与报告；claim boundary 写入正文 | 可选 AF3/Schrödinger 精修模块 |
| **Chemistry (升级)** | 邻域库设计与 TopN 可采购清单 | 合成 + kinact/KI + C116S + 细胞 engagement |

## Disease-Model Narrative Rule（绑定本贡献）

正文允许写：本框架产出的探针将用于回答上皮应激/炎症等场景中 JNK1 vs JNK2 的未决问题。  
正文禁止写：本工作已证明 JNK2 在某疾病模型中具有治疗优势（除非未来出现 JNK2 依赖的硬证据并单独立项）。
