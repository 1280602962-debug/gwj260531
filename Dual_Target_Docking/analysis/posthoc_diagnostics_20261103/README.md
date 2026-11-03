# 投稿前支持性诊断（2026-11-03）

本目录是独立于正式结果的只读支持性分析。正式八对靶点、867 条合格观察、评分和模型分折均未修改。

## 基线

- 仓库：`https://github.com/1280602962-debug/gwj260531`
- 正式基线分支：`cursor/v4-2-final-audited-results-20260930`
- 分析时核实的基线提交：`d1ca17512eabe5c9d102ebe75c3f9c59668452dc`
- 工作分支：`cursor/jcim-posthoc-diagnostics-20261103`
- 评分映射以 `submission/SCORING_SCHEMES.csv` 为准：S1=M0，S2=M1，S3=M2，S4=M1b，S5=M3。

远端 HEAD 在执行时与上述提交一致，因此没有版本混用。

## 执行

在 `Dual_Target_Docking/` 下运行：

```bash
python3 analysis/posthoc_diagnostics_20261103/scripts/run_posthoc.py
python3 -m pytest analysis/posthoc_diagnostics_20261103/tests -q
```

脚本会连续计算两次并比较数值输出的 SHA-256。随机种子为 271828。输入文件哈希写在 `qa/INPUT_SHA256_BEFORE.csv` 和 `qa/INPUT_SHA256_AFTER.csv`。

## 分析状态

- 分析 A：已完成。七对靶点没有非空 Murcko 骨架跨折，骨架分组重训按预设规则未执行，因此没有 `A_SCAFFOLD_GROUP_CV_*.csv`。
- 分析 B：已完成。EGFR/HER2 的 D_vs_B AUROC 已按正式成员复现。
- 分析 C：已完成。Top 10% 复现通过后才计算 Top 5% 和 Top 20%。

未执行：ECFP4 加 S1 增量的 bootstrap；35 条跨阈值重复测定的排除敏感性。

## 输出索引

- `outputs/`：三项分析的 CSV。
- `figures/`：支持性图，标记为 Supporting analysis，不使用正文 Figure 编号。
- `reports/A_SCAFFOLD_AUDIT_CN.md`：化学相似性审计。同文另存为 `A_ANALYSIS_REPORT_CN.md`。
- `reports/B_EGFR_INVERSION_DIAGNOSTICS_CN.md`：EGFR/HER2 反向排序。
- `reports/C_TOPK_SENSITIVITY_CN.md`：截断比例敏感性。同文另存为 `C_TOPK_SENSITIVITY_REPORT_CN.md`。
- `reports/FINAL_SCIENTIFIC_ASSESSMENT_CN.md`：汇总评估。
- `reports/CLAIM_EVIDENCE_ALIGNMENT.csv`：结论与证据对照。
- `qa/`：输入哈希、连接审计、环境版本和数值输出哈希。

冻结参数见 `config/analysis_config.json`。
