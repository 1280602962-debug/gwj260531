# 中英文润色修改摘要（2026-09-12，审稿边界轮）

本轮读取 PR #32 的 14 受体统一 CalcRMS 补充，并按投稿前 A/B/C 清单同时改中英文。数字未改。未改官方论文 PR #32。

## 科学边界

1. **2.3 靶对筛选可复现。** 按 ChEMBL 37 无序对宇宙写出门槛：2,164,618 → 63,790 → 5,253 → 86 → 26 → 19 → 17/12 → G5 8/6；排除共价 CTSK/CTSS；EGFR/HER2 作为供给受限例外。主文不写 PIK3CA/PIK3CB 项目史。
2. **2.2 主数据版本。** Table 2：三对 2026-07-23 REST；五对 ChEMBL 37 dump（2026-05-01）。API 与 dump 只作敏感性。
3. **标签异质性。** 4.5 写明 endpoint types remain heterogeneous despite pChEMBL；max pChEMBL 只部分缓解；类别是 database-derived retrospective states。
4. **summary_min。** 两个方向 AUROC 为一级终点；summary_min 为较弱臂描述，不是 overall performance。SI 说明非分层 pooled bootstrap。
5. **Vina 字段。** 已核 replay 脚本：Table 2 为第一模型 `REMARK VINA RESULT` / 生产 `*_affinity`。
6. **RMSD。** 并入 PR #32 的 `all14_cognate_rmsd_calcrrms_v1.csv`。14 受体同一 CalcRMS；搜索覆盖数字不变。
7. **Figure 3。** competing-explanation，不是 head-to-head benchmark。
8. **3.4–3.6。** 压缩种子审计、标签一致率堆砌；“不是外部验证”只留一处。

## 其他

- dual-measured → ligand measured at both targets
- 摘要补异质性，弱化 summary_min 技术句
- Figure 1 图注改为中性分开陈述
- mismatched：imperfect specificity control
- Conclusion：both-end experimental measurements；holdout 未一致再现
