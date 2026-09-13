# 未解决问题（投稿冻结）

只列仍须作者或投稿阶段处理的事项。当前工作分支：`cursor/jcim-language-polish-0b1a`（PR #37）。

## 仍开放

1. **固定成员脚本只排除空字符串。** 未另检 NaN/inf 字面量；未发现改变已采用的五对结果。完整种子表留在仓库，不进入正式 SI。
2. **未做高可比性 / 同测定子集分析。** 讨论中仅作为公开活性记录异质性陈述，不是本轮缺口修补对象。
3. **位图未做像素级目视审核。** 正式 SI 图为 Figure S1–S4。`FigS6_detectable_effect` 已从投稿 SI 撤出，仅留仓库。
4. **英文投稿格式。** 内容已与中文对齐，投稿前仍需期刊模板转换和常规 copy-edit。
5. **主文已删除 PIK3CA/PIK3CB（2WXF）撤出史。** 仓库 QC / 历史目录仍保留该记录，不进入主评价。

## 本冻结已关闭

- **14 个主受体共晶 RMSD 已用同一套化学对应 CalcRMS 重算。** 来源：`all14_cognate_rmsd_calcrrms_v1.csv`。14/14 最低保存姿态 < 2 Å；EGFR 3POZ 为重建 QC。正文、Table S2 与 Figure S4 使用该统一表。早期坐标匈牙利匹配不再作为正式 RMSD。
- **新增五对生产姿态已入 git。** 路径：`data/jcim_chembl_universe_v0/local_track_b_v0/poses/`（9742 个 pdbqt，含 1094 个 `mode_01`）。这是原保存树，不是新对接。`logs/vina/` 仍不入库。未申请 Zenodo DOI。

上述开放项均不推翻当前主结论，也不要求本轮新增实验。

## 投稿冻结

2026-09-13 在本分支运行 `python3 scripts/audit/freeze_submission_v1.py`，随后打 annotated tag `dualfourclass-jcim-pr37-freeze-20260913`。冻结后默认不再改数据或正式稿数字。本提交只补原生产姿态文件，不改分数或主表数字。
