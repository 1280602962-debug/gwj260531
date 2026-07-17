# Phase 2 — Annotated core corpus（精选，非 OpenAlex 全量）

证据等级：A=直接决定初筛设计；B=强方法迁移；C=背景/警示。

| ID | Citation | Venue/Year | Grade | Why included |
|----|----------|------------|-------|--------------|
| A1 | Lu et al. — YL5084 / YL2056 / 8ELC | J Med Chem 2023 | A | 选择性锚点 + 共晶模板；kinact/KI 选择性逻辑 |
| A2 | Wydra et al. — 56d ligand-first | J Med Chem 2025 | A | 正交骨架；装 warhead 前先优化可逆识别 |
| A3 | Zhang et al. — JNK-IN-8 | Chem Biol 2012 | A | Pan 阳性 / Cys116 范式校准 |
| A4 | Shamir/London — AF3 covalent ligands | JACS 2025/26 | A | mPAE 主排；大丙烯酰胺库前瞻；对接非充分 |
| A5 | Nat Commun — precision-guided JNK Cys116 warheads | 2024 | A | 仅标准丙烯酰胺已非唯一前沿；反应性/几何精度重要 |
| B1 | Boike et al. — Advances in covalent drug discovery | Nat Rev Drug Discov 2022 | B | 共价发现总图：识别 vs 反应性平衡 |
| B2 | Zhu et al. — CovDock-VS | JCIM 2014 | B | 共价对接 VS + decoy 富集；相互作用后过滤 |
| B3 | London et al. — DOCKovalent | Nat Chem Biol / PMC | B | 大库共价对接；反应性不进打分 → 库设计要控 warhead |
| B4 | Extended warheads for CKIs | JCIM 2024 | B | acrylamide 仍最常见；extended warhead 化学可参考 |
| B5 | CovalentLab | JACS Au 2025 | B | warhead 安装/片段生成平台；强调几何与反应位点 |
| B6 | RosettaAMRLD | JCIM 2025 | B | 反应驱动组合库；**对接前可配体过滤**（MW/RB/logP） |
| B7 | SwissSimilarity 2021 | IJMS 2022 | B | ECFP + ErG + 2D pharmacophore 并用 |
| B8 | SHAFTS / HybridSim | JCIM / Bioinformatics | B | 3D shape+pharmacophore 补 scaffold hop |
| B9 | MolPAL | Chem Sci 2021 | B | 主动学习减昂贵 oracle 调用 |
| B10 | Deep Docking / ultra-large VS protocols | Nat Protoc 2022 | B | 超大库分层：代理模型→精筛 |
| B11 | CovalentInDB 2.0 | NAR 2024/25 | B | 外部共价库与已知共价药相似对照 |
| B12 | Enrichment confidence bands | J Cheminform 2022 | B | 富集不确定性量化 |
| C1 | Covalent docking pitfalls | Mol Divers 2022 | C | 警示错误实践 |
| C2 | Electrophile fragment screening | JACS 2019 | C | 片段/反应性测绘传统 |
| C3 | PAINS utility limits | ACS Chem Biol 2018 | C | PAINS 勿黑箱硬杀 |
| C4 | Project Phase0 AF3 vs Glide | internal | A/local | mPAE 主门控；Glide 粗筛；56d Glide 假阳高 |

完整 OpenAlex 原始检索：`openalex_search_raw.json`（Q1–Q7）。
