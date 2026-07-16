# PaperSpine 立项评估纪要（JNK2 共价选择性）

**生成：** 2026-07-16  
**Skill：** PaperSpine v4.0.0（已镜像安装到 Cursor skills）  
**阶段：** Research 完成 → **等待用户确认 motivation/contribution**（按 V4 硬门禁，未开写全文）

---

## 0. Cursor 安装说明

上游 `install.sh` 官方支持：Claude Code / Codex / OpenClaw / Hermes，**不含 Cursor**。  
本环境已做镜像安装：

| 位置 | 路径 |
|------|------|
| 用户级 | `~/.cursor/skills/paper-spine/` |
| 项目级 | `/workspace/.cursor/skills/paper-spine/` |
| 说明 | `/workspace/.cursor/paperspine.md` |
| 官方副本 | `~/.claude/skills/paper-spine/`、`~/.codex/skills/paper-spine/` |

本地桌面 Cursor 若要使用：把该 `paper-spine` 目录复制到本机 `~/.cursor/skills/`，或在对话中 `@` / 说明启用 PaperSpine。云端会话已可读项目内 skill。

---

## 1. 课题建立是否合理？（Contribution-first 结论）

**合理，但必须换主贡献。**

| 维度 | 分数 | 理由 |
|------|------|------|
| 化学/结构可成性 | 高 | Cys116、8ELC、丙烯酰胺、Leu106 轴文献扎实 |
| 与近三年方向一致性 | 高 | 对齐 Lu 2023、Wydra 2025 |
| 计算框架可辩护性 | 中高 | Phase 0 校准有意义；不可夸大 |
| 疾病成药叙事 | 低 | JNK2 优势适应症证据不足 |
| **总体立项（探针/框架）** | **7.5/10** | 建议继续 |
| **总体立项（疾病驱动）** | **3.5/10** | 不建议作为主线 |

**审稿人一句话：**  
可以做；请证明你交付的是 **探针/框架**，不是又一篇“我们对接了 JNK2”或“DSS 有效所以 JNK2 成药”。

---

## 2. 疾病模型中 JNK2 优势不突出 —— 叙事如何完善

### 错误叙事
> JNK2 是 IBD/纤维化/肿瘤的优势靶点 → 我们做共价药 → 模型有效。

### PaperSpine 推荐叙事（Introduction ladder）

1. **Problem：** JNK1/JNK2 功能可能拮抗，但 ATP 口袋 ~98% 相同。  
2. **Progress：** pan-JNK 共价（JNK-IN-8）、JNK1-bias 临床（CC-90001）、JNK2 共价先导（YL5084、56d）。  
3. **Gap：** 疾病模型（含 DSS）**缺乏 isoform-resolved 工具**，故“JNK2 是否更优”**无法回答**（这是缺口，不是你的失败）。  
4. **RQ：** 能否建立以 8ELC/Cys116 为锚、可校准、可复现的 JNK2>JNK1 共价筛选与验证闭环？  
5. **Contribution promise：** 框架 ± 新候选（待证据）。  
6. **Boundary：** 不宣称已解决适应症；疾病模型仅作探针应用场景展望。

### 可写 / 勿写

| 可写 | 勿写 |
|------|------|
| 探针使上皮应激/炎症/存活场景可拆分 JNK1 vs JNK2 | JNK2 是 DSS 结肠炎主靶 |
| Lu 也未能证明 MM 表型完全依赖 JNK2 | 我们的化合物已证明 JNK2 治疗优势 |
| CC-90001 = JNK1 轴对照，不是 JNK2 金标准 | 用 CC90001 临床成功外推本课题 |

---

## 3. 后续推进的具体方向（按 PaperSpine Evidence required）

### P0 — 关贡献门禁（本周可完成）
1. 用户确认 Motivation Option（见 `motivation_options_after_research.md`）  
2. 固化 `confirmed_contribution.md`（贡献句 + claim boundary）  
3. 完成 COValid 式 decoy 与 EF@1% **方法学报告**（你们已有回顾性数字，差规范交付）  
4. 修正 AF3 输入（YL2056 bonded pairs；Cys116 vs 113 编号）

### P1 — WSL 可闭环的主路径
1. **Layer0 校准集**（已有种子）+ 50× decoy/阳性  
2. **Layer1 邻域库** 500–2000：YL5084/56d/4WHZ 骨架 + 丙烯酰胺 SMARTS  
3. RDKit 过滤 → AutoDock4/ADFR/AMECovDock 或 DOCKovalent（小库）  
4. Top100 人工/规则检查 Cys 几何与 Leu106 占据  
5. 可采购清单（Enamine acrylamide / REAL 子空间）

### P2 — 文章上限取决于湿实验
1. 对标 YL5084/56d：IC50、**kinact/KI**、C116S、washout  
2. 细胞 engagement；可选 NanoBRET  
3. **一个** 生物学拆分实验（遗传背景），不要以 DSS 疗效作主图

### 明确后置
- 百万级 AF3 发现叙事  
- 3NPC 共价主筛  
- Δsel/Gly87 复活  
- 疾病药效主贡献

---

## 4. 共价筛选库：如何选、如何建

### 近三年文献实际做法

| 来源 | 库策略 |
|------|--------|
| Lu 2023 | **不做大库盲筛**；JNK-IN-8 衍生数十类似物 |
| Wydra 2025 | **Ligand-first**：可逆 JNK2/3 系列 → 叠合装 acrylamide |
| COValid/AF3 | 活性丙烯酰胺 + **性质匹配、拓扑不相似 decoy**（ZINC 骨架装 warhead） |
| DOCKovalent 传统 | 分 warhead 商业电泳库大规模对接 |
| Enamine | 现成 Cys/acrylamide plated 库（千～万）供实验或预筛 |

### 推荐三层库（与你们项目匹配）

```
Layer0  校准（必须）     ~7 actives/negatives + 50–150 decoys
Layer1  发现主库         200–2000 邻域丙烯酰胺（ligand-first）
Layer2  多样性扩展       2k–10k Enamine sACR / ZINC SMARTS（可选）
```

### 建库硬规则
1. Warhead 固定丙烯酰胺（主筛）  
2. Linker 几何对准 Cys116（meta 优先于盲目 para）  
3. 优先 Leu106 后口袋药效团  
4. 排除过反应弹头；保留饱和酰胺阴性  
5. Decoy：性质匹配 + 低 Tanimoto  
6. 可采购 > 纯虚拟

---

## 5. 普通电脑 WSL 能实现什么？

### 能做（建议写成可复现方法）
- RDKit：清洗、SMARTS、decoy、构象、物化过滤  
- OpenBabel / RDKit 3D  
- AutoDock Vina / Gnina（预筛，非共价终局）  
- AutoDock4 共价 / ADFR + AMECovDock（数百～数千）  
- DOCKovalent 小库本地或官方 web  
- GROMACS 短 MD（少量 pose）  
- 全部主表、审计、报告脚本

### 不能/不宜依赖
- 大规模 AlphaFold3 共折叠虚筛（需官方 server/大 GPU）  
- Schrödinger CovDock 大规模（许可+算力）  
- 百万级共价对接作为唯一故事

### WSL-first 漏斗（可写进方法）

```
WSL: 建库 → decoy校准 → ADFR/DOCKovalent(≤2k) → TopN
可选云/机房: AF3 mPAE 或 Schrödinger 精修 Top100–500
实验: 采购 → GSH/生化 → C116S → kinact/KI
```

---

## 6. 建议的 Core Contribution 草案（待你确认后写入正式文件）

> 本文/本项目建立一个以 JNK2–YL2056 共晶（8ELC）与 Cys116 丙烯酰胺化学为锚、经 COValid 式回顾性校准的选择性共价筛选与候选推进框架，并明确将疾病适应症问题重新表述为「有待 isoform-selective 探针回答的开放生物学问题」，而非已解决的疗效主张。

**Strong claims allowed:** 结构模板选择合理；校准流程可复现；与 Lu/Wydra 路线兼容。  
**Claims to soften/avoid:** 已发现优效药；DSS/纤维化中 JNK2 治疗优势；AF3 已完成前瞻性发现。

---

## 7. 用户确认结果（2026-07-16）

**已确认：B + 可复现框架必达**

- 已写入 `confirmed_motivation.md`
- 已写入 `confirmed_contribution.md`
- 双轨：
  - **Framework（必达）：** 8ELC/Cys116 校准筛选—验证闭环（含 WSL 开源漏斗）
  - **Chemistry（升级）：** ligand-first / 邻域扩展 → 可检验新共价候选

下一阶段（若继续 PaperSpine）：section blueprints / 开题或技术报告写作，仍受 claim boundary 约束。
