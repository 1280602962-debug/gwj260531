# 论文图表（JCAMD / Springer large-journal）

生产图由 `scripts/plot_jcamd_publication_figures.py` 生成。  
**每一个绘制数字都从归档 CSV/JSON 读取，并在脚本内断言**；图面不写图注（图注见 `generated/CAPTIONS.md`）。

## 版式（与近三年 Springer / JCAMD 一致）

| 项 | 设定 |
|----|------|
| 栏宽 | 单栏 84 mm，双栏 174 mm |
| 字体 | 8 pt 无衬线（Liberation Sans，Arial 度量兼容；构建机无 Arial 许可） |
| 颜色 | Okabe–Ito；无网格 |
| 矢量 | PDF，`pdf.fonttype = 42`（TrueType 嵌入） |
| 位图 | PNG 300 dpi（预览与排版核对） |
| 图注 | 单独 markdown，不叠在图上 |

七个化学过滤名单在图上标为 **chemistry-filtered audit set / list**，不是 dual-node candidates 或 follow-up hits。不画 MD 轨迹。没有 OOF 预测表则不画 ROC 曲线。

## 主图

| 文件 | 内容 | 数据来源 |
|------|------|----------|
| `generated/main/fig01_protocol_enrichment_selfdock.*` | P0–P5 EF@1%（True vs Random，bootstrap 95% CI）+ lesinurad 自对接 RMSD | `data/si/protocol_enrichment_ci/protocol_ef_ci.csv`；`data/redock_smoke/redock_results_lesinurad_9DKB.csv` |
| `generated/main/fig02_funnel_dual_percentiles.*` | 漏斗 8319→1588→1580→51→7；\(S_U\) vs \(S_{N,\mathrm{dock}}\) | `nlrp3_screening_summary_clinical_all.json`；`pareto_merged_scores.csv`；`nominated_shortlist_diverse.csv` |
| `generated/main/fig03_pose_qc.*` | 羧酸–Arg477；URAT1 / NLRP3 质心位移 | `data/si/pose_qc/pose_qc_table.csv`；`pose_qc_dual.csv` |
| `generated/main/fig04_active_set_decoy_leakage.*` | 阳性集命名药缺席、骨架偏置、诱饵 Tanimoto | `actives.csv`；`decoy_leakage_audit/` |

## 附图

| 文件 | 内容 |
|------|------|
| `generated/si/figS01_gate_sensitivity.*` | τ 敏感性（不替换 τ = 90） |
| `generated/si/figS02_protocol_auc.*` | P0–P5 ROC-AUC |
| `generated/si/figS03_nlrp3_fold_metrics.*` | NLRP3 骨架 CV 折间 AUROC/AUPRC（无 OOF 曲线文件） |
| `generated/si/figS04_known_gout_percentiles.*` | lesinurad / verinurad / 秋水仙碱百分位 |
| `generated/si/figS05_nlrp3_qN_histogram.*` | 临床库 8319 的 \(q_N\) |

图注：`generated/CAPTIONS.md`（英）、`generated/CAPTIONS_CN.md`（中）。  
数值锁：`generated/DATA_LOCK.json`。  
表：`generated/tables/`。

## 重新生成

```bash
python3 scripts/plot_jcamd_publication_figures.py
```

断言失败即停止：例如 P2 True EF@1% 必须为 2.5867、完整案例必须为 1,580、优选名单必须为 7。  
旧脚本 `scripts/plot_available_figures.py` 仍含“follow-up”用语，**不要**再用于投稿主图。
