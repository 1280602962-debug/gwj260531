# 已归档路线（请勿再按此执行）

> 以下文档/脚本描述 **已废弃叙事**（TAPE-GATE / MASFL / **默认 Glide XP 主文** 等）。  
> **当前论文主文档**：[`MANUSCRIPT_OUTLINE_V2.md`](MANUSCRIPT_OUTLINE_V2.md)  
> **中文 Methods**：[`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)

---

## 为何归档

| 旧叙事 | 归档原因 |
|--------|----------|
| TAPE-GATE 库筛 + 生成式双路径 | 未跑通；与当前重定位稿无关 |
| 8973 双靶 Pareto / NLRP3 ML on distill | **科学错误**：8973 为 URAT1 偏置集 |
| Teacher M-CPDL / 8973×三态全库 | Gate 未过；算力与叙事均已切换 |
| OAT 迁移作为主创新 | Δρ≈0.004，已从论文删除 |
| $S_{\text{trap}}$ 三态主筛 | B1K/B1L 刚性失败；留给远期论文 B |
| geometric mean / 对称融合主图 | 数据不支持双靶监督 |

---

## 归档文档列表（只读参考）

| 文件 | 原用途 |
|------|--------|
| **`MANUSCRIPT_DRAFT_CN.md`** | **中文全文旧稿（Glide XP 主叙事）——易与 V2 混淆，已加过时横幅** |
| `MANUSCRIPT_OUTLINE_REVISED.md` | 修订大纲（已被 V2 取代） |
| `MANUSCRIPT_OUTLINE_CURRENT.md` | 旧「定稿思路」；现仅保留重定向页 |
| `LOCAL_AGENT_POST_DOCKING_PROMPT.md` | Maestro/Glide XP 后处理 prompt |
| `LOCAL_AGENT_8973_DOCKING_PROMPT.md` | 8973 Glide XP 回顾 prompt |
| `TAPE_GATE_FRAMEWORK.md` | v2.0 双路径框架 |
| `MASFL_V3_WORKFLOW.md` | Teacher 蒸馏 |
| `ALGORITHM_FRAMEWORK.md` | 全阶段公式 |
| `PAPER_A_PRIME_PLUS_LOGIC.md` | 8973 双证据（已修正为仅 URAT1 回顾） |
| `COMPLETE_WORKFLOW_AND_FILES.md` | 旧端到端索引 |
| `TWO_PAPER_STRATEGY.md` | 双轨 JCIM + 快速线（部分仍可参考期刊表） |
| `LOCAL_AGENT_P1P2_PROMPT.md` | OAT 迁移离线 prompt |
| `INNOVATION_POINTS.md` | 旧创新点表述 |

---

## 仍有效的子模块

| 模块 | 文档 |
|------|------|
| PDB / 9DKB 对接参数 | `URAT1_THREE_STATE_DOCKING.md`, `config/docking_ensemble.yaml` |
| Benchmark 定义 | `BENCHMARK_SELECTION_CRITERIA.md`, `data/benchmarks/` |
| URAT1 ML 评估 | `URAT1_ML_MODEL_ASSESSMENT.md` |
| 数据事实 | `DATA_FACT_CHECK.md` |

---

## 脚本状态

| 脚本 | 状态 |
|------|------|
| `run_tape_gate_pipeline.py` | 骨架，**不用** |
| `03_library_screening.py`, `04_generative_optimization.py` | 骨架 |
| `05_fusion_and_ranking.py` | 骨架 |
| `run_oat_transfer_ablation.py` | 可跑但 **不写入论文** |
