# STAGE_M_VERDICT — DualFourClass measurement audit

**Date:** 2026-07-28  
**Pack:** `data/stage_m_v0/`  
**Authority:** `docs/PLAN_V2_REDTEAM_AND_REDESIGN.md`, `docs/AGENT_COMMAND_STAGE_M.md`

## Gate table

| 项 | 状态 | 关键数字 | 含义 |
|----|------|----------|------|
| M1 | **Go** | EGFR vina D/A=0.689, D/B=0.311; PM D/A=0.698, D/B=0.597 | 方向分解已固化为主指标；池化仅附录 |
| M2 | **Weak** | gray EH=46%; strict B_only n=7; oracle@σ0.5 仍可辨（D/A≈0.94） | margin 子集 underpowered / gray 大；阈值未翻转，但 Track B 高风险 |
| M3 | **EGFR No-Go / PM Go** | EH heavy_atoms summary_min=0.369 > vina 0.311; PM vina 0.597 > heavy 0.463 | 对接增量不成对通用；勿平均两对靶 |
| M4 | **Go** (M4-min) | EH40 RDKit: rtm summary_min 0.467 vs LigPrep 0.607 | 统一 prep 对照完成；RTM 绝对值对 prep 敏感 |
| M5 | **Go** | 主竞赛臂≤4（`CANDIDATE_ARMS_V1_STAGE_M.yaml`） | rank_consensus / rtm_min_z 移出主竞赛 |

## Track B 总门控

**判定：Track B = Weak（仅文档规划，不启动大批对接）**

规则对照：
- Full-Go 需要 M1∧M2∧(M3≥1对Go)∧M4∧M5 → **M2=Weak 阻断 Full-Go**
- Weak：M2=Weak（且 M1–M3 非全面 No-Go；M4 已完成）
- 非 No-Go：M2 非 No-Go；并非两对靶 M3 皆 No-Go

**不允许：** 未经再次门控批准就开第三/四对靶大批对接或 Track B B0 选臂竞赛。

## Track A 写作头条（必须采用）

1. **方向抵消**：EGFR D/A≈0.69 与 D/B≈0.31 池化成 ~0.52（假“接近随机”）  
2. **体积基线**：EGFR 上 heavy_atoms ≥ 对接臂（`min(D/A,D/B)`）  
3. **标签边界**：硬负大量贴 θ=6±0.5；strict margin 灰度高  
4. **prep 混杂**：M4-min 显示 RTM 绝对值强依赖 LigPrep vs RDKit；不得把混 prep 的 panel120 RTM 分裂写成已验证方法结论  

## 下一步

1. **默认：** 跑 Track A STEP0（`AGENT_COMMAND_STEP0_TRACK_A_STARTER.md`）——诊断论文启动包  
2. **可选加强（非必须）：** M4-full（全 110 统一 RDKit）仅当需要彻底消除 panel120 混 prep 叙述  
3. **禁止：** 同协议 EGFR 扩样赌显著；未门控开 Track B

## Reproduce

```bash
python Dual_Target_Docking/data/stage_m_v0/scripts/run_m1_directional.py
python Dual_Target_Docking/data/stage_m_v0/scripts/run_m2_labels.py
python Dual_Target_Docking/data/stage_m_v0/scripts/run_m3_baselines.py
# M4-min (local docking) already run; scripts in stage_m_v0/scripts/m4_*.py
```
