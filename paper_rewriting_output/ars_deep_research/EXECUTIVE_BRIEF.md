# ARS Deep-Research Executive Brief — JNK2 对接前初筛

**Skill:** academic-research-suite (Imbad0202/academic-research-skills-codex)  
**Mode:** deep-research → experiment plan  
**Date:** 2026-07-17

## Answer in one paragraph
OpenAlex 系统检索（7 条查询、约 160 篇去重）与核心非 OA/顶刊文献综合表明：胺→丙烯酰胺→物化 QC→相似性缩库是常见工程实践，但 **ECFP 四锚点 + 固定 Tc + Murcko 双轨 as-is 不足以称为经证据支撑的“活性初筛”**。对 JNK2/Cys116/8ELC 项目，理性方案应改为：全库多表征特征 → 小集校准定阈值 → core-Tc 驱动 Track-Sim、ErG/药效团辅助 Track-Novel → 反应性分桶 → ~10k 分仓交付 → 松 Glide → AF3 mPAE。Devil’s advocate 结论为 **REVISE**。

## Deliverables
| File | Content |
|------|---------|
| `phase1_scoping/RQ_brief.md` | RQ / FINER / scope |
| `phase1_scoping/methodology_blueprint.md` | 方法蓝图 |
| `phase2_investigation/search_strategy.md` | OpenAlex 策略 |
| `phase2_investigation/openalex_search_raw.json` | 原始检索 |
| `phase2_investigation/annotated_core_corpus.md` | 精选语料 |
| `phase3_analysis/synthesis.md` | 主题综合 + DA |
| `phase4_plan/prefilter_replan_v1.md` | **重规划 SOP + WSL 指令** |

## Next action for user
把 `phase4_plan/prefilter_replan_v1.md` §4 整段贴给本地 WSL agent；先要求它只跑 L0–L5 并回报校准表，再授权 L6–L7。
