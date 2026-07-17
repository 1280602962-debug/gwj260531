# SOTA Gap Map

**User motivation:** 评估立项合理性、推进方向、疾病叙事、共价库、近三年做法、WSL 可行性。

## Gap Table

| Candidate Contribution | What SOTA Already Does | User Evidence | Real Gap | Claim Strength | Risk |
|---|---|---|---|---|---|
| A. 新的 JNK2>JNK1 共价先导分子 | Lu 2023 YL5084；Wydra 2025 56d | 主表/审计/8ELC 路线；尚无自有合成与 kinact/KI | **有空间，但必须做出新化学或更优属性** | 高（若有湿实验） | 仅计算重述文献 → 拒稿 |
| B. 可复现的 Cys116 共价筛选框架（8ELC + decoy 校准 + 开源/WSL 漏斗） | COValid 给了通用方法；JNK2 位点专用校准与开源漏斗少见 | Phase 0 种子 7/7；AF3 回顾性 AUC/EF；决策五阶段；排除 3NPC | **JNK2-specific calibrated open workflow** 是真实缺口 | 中–高（方法学/报告） | 若夸大 AF3 发现力则弱 |
| C. Ligand-first 邻域库 → 共价扩展（对齐 Wydra） | Wydra 已示范该范式 | 有 4WHZ/26k 与 56d 逻辑文档；库未建成 | **可执行建库+对接闭环** 仍缺 | 中 | 与 Wydra 过近需强调差异 |
| D. JNK2 疾病模型药效突破 | 领域整体缺 isoform-resolved in vivo | DSS/CC90001 讨论材料；无自有 JNK2 依赖证明 | **不是当前可防守贡献** | 低 | 硬写会伤害全文 |
| E. Gly87/Δsel 计算选择性 | 用户侧已证伪 | 负结果已记录 | 负结果可作边界，不可作贡献 | 低 | 勿复活 |
| F. 探针驱动的细胞亚型拆分平台 | Lu/Wydra 有细胞 readout，但公开“可复用实验包”少 | 方案已写 kinact/KI、C116S、biotin 竞争 | 可作次贡献 | 中 | 需至少一套自有数据 |

## Gap Summary（最值得做的 2 个）

1. **主贡献候选：B（+部分 C）** — 「以 8ELC/Cys116 为锚的、经过 COValid 式校准、可在常规算力复现的 JNK2 共价筛选与候选推进框架」，化学发现作为框架的验证而非唯一卖点。  
2. **升级贡献候选：A** — 仅当 Layer1 库筛或 MedChem 产出 **区别于 YL5084/56d** 的新骨架/更优选择性或 developability 时，再升为文章主贡献。

**明确放弃作为主贡献：D（疾病药效）**——应改写为「探针使争议场景可被提问」。

---

## Addendum 2026-07-17 — 对接前 Step0–7 初筛漏斗

详见：`prefilter_pipeline_deep_research.md`

| Candidate | SOTA 已有 | 本漏斗 | Real gap | Claim |
|---|---|---|---|---|
| G. ECFP+Murcko 胺→丙烯酰胺双轨 | 模块全是标配 | Step0–7 | **无方法缺口** | 低（勿作主创新） |
| H. 经校准的 JNK2/8ELC prefilter（core-Tc+ErG/3Dpharm+反应性+UQ） | COValid/ErG/Kin-Cov 分散存在 | 尚未做 | **位点专用校准 triage** | 中（支撑 Framework） |
| I. 选择性差分（JNK2 vs JNK1）接入枚举后分仓 | Lu MD/kinetics；少见大规模预筛 | 缺失 | **与贡献句对齐的关键缺口** | 中–高（若有数据） |

**判决：** Step0–7 是工程 SOP；创新在 H/I + 湿实验，不在再写一套 Tc 阈值。
