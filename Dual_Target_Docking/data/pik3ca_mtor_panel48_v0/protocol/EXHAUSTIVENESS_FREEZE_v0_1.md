# EXHAUSTIVENESS FREEZE v0.1 — pik3ca_mtor_panel48_v0

## 结论

**`exhaustiveness_v0_1 = 16`（仅 DTPAIR-01 / PIK3CA–mTOR panel48）。**

`egfr_her2_panel40_v0` 的 `exhaustiveness_v0_1 = 8` **不改**。

## 升 E 合法理由（cognate 采样 QC）

| 条件 | E=8（seed 20260727） | E=16 诊断（seed 20260727） |
|------|----------------------|----------------------------|
| 4L23 best_of_9 | 0.624 Å（过） | 0.624 Å（过） |
| 4JT6 best_of_9 | **5.003 Å（不过）** | **0.445 Å（过，mode3）** |
| 4JT6 mode1 | 7.118 Å | 7.118 Å（仍失败；打分排序问题） |

来源：`analysis/cognate_redock_v0/COGNATE_QC_VERDICT.md` 与 `tables/pm48_01_rmsd_E16_diag.csv`。

## 明确边界

- 禁止因刷分升到 E=32 做全面板；E=32 仅允许 `PM48_01@4JT6` 诊断。
- 全面板在 `COGNATE_QC_VERDICT_E16.md` 写明 **Verdict: Go** 之前禁止启动（现已 Go 并完成）。
- EGFR/HER2 panel40 保持 E=8。

## 主文报告臂（定稿）

并列报告 **`vina_mean` 与 `rtm_min_z`**；禁止只报 RTM。

## Limitations（协议层）

- T2：PM48_26 / 20 / 21（化学型同源假双靶）
- T5：Torin1 / Omipalisib 重打分误伤
- PM48_34@4JT6 仅 8 个有效 mode
- clash 门控阴性；shortfall / consensus 冻结消融无法同时压硬负并保护 T5
- 不宣称 C4 已成功外推

## 签字等价

- agent 运行时间戳：2026-07-27T17:30:50+08:00
- hostname：LAPTOP-3GOC1J6E
- 工作根：`/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_v0`
- P1 协议回写：2026-07-27（decision_ablation_v0 + warning_flags）
