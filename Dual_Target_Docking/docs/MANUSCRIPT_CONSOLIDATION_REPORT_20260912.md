# 稿件收尾修订与核查报告

日期：2026-09-12  
仓库：`1280602962-debug/gwj260531`  
项目：`Dual_Target_Docking`

## 本轮基准

- PR #32 实际 head 分支：`cursor/jcim-submission-pack-0b1a`
- 完整 SHA：`381bfe1692e0a1b93d8105194bcdbb7401ea75a8`
- 独立修订分支：`cursor/jcim-manuscript-consolidation-0b1a`（由该 SHA 检出，不覆盖原分支）
- 同名 `codex/jcim-manuscript-consolidation` 未使用；本环境按既有 `cursor/…-0b1a` 命名。
- 工作区在检出时干净。未推送、未强制推送、未合并、未更新 PR。
- 本轮本地提交：`e57002b2c060692bc4b6cf810f52e48b5d5fe97c`（本文件随后如仅补记该 SHA，以分支最新提交为准）。

事实来源仅限该提交中的有效文件，不以旧对话数字为准。

## 主要修改及理由

1. **拆开 max / median / 高置信 / 转储四项检查。** 正文曾把 API 快照上的聚合检查与高置信字段筛查写成“结果未改变”，并漏写新增五对已完成的 dump max/median。现按数据版本分别报告。
2. **收窄 GNINA 与 QC 用语。** 英文摘要不再写 independent GNINA “reproduced” 固定通道效应；共晶 QC 去掉“最终认可 / 槽位 / finally accepted / slots”，并写明 8 姿态体系、重建 QC 与最低 RMSD 的含义。
3. **外部数据与可用性。** BindingDB 限定为“在所用来源、过滤和门槛下未形成外部评价集”；可用性区分已公开分数表、本地回放、checksum 与仍待发布的生产姿态 / DOI。
4. **图表说明。** 删除“图例置于绘图区外 / legend is placed outside the plotting area”。未重绘、未改绘图风格。
5. **唯一主张—证据清单。** 重写 `docs/PRIMARY_RESULT_INDEX_V1.md`，不再另建一批重复审计文件。

中文稿先稳定，再同步英文稿、英文 SI、图注和投稿包副本。

## 已修正的明确事实错误

| 错误 | 现行依据 | 修正 |
|---|---|---|
| 将 max/median 与高置信写成“标签和结果均未改变” | `assay_max_vs_median_agreement_v1.csv`；`high_confidence_summary_v1.csv` | 分项报告：EGFR 冻结 0.430 / API-max 0.417 / median 0.424（一致率 93.6%，含 EH120_060 快照差异）；AChE median 0.629 ≠ 主结果 0.606；高置信 253/253 不变 |
| SI / 正文只写“仅覆盖三对”，未区分仓库已完成的五对转储分析 | `five_pair_dump_gated_v1/max_vs_median_auroc_v1.csv` | SI Table S3 增补 dump 范围：仅 PPARA/PPARD CHEMBL121 翻转；五对 `summary_min` 点估计三位小数不变；其 CI 不进入正文 |
| 英文摘要将独立 GNINA 写成 reproduced 固定通道效应 | Table S9a；独立姿态生成，非同姿态重评分 | 改为另一姿态生成流程下仍观察到类似任务差异 |
| Results 使用“最终认可 / 槽位 / finally accepted / slots” | 内部工作用语 | 改为正常论文用语 |
| 英文 Results 3.6 写成“公开数据不足以构建外部集” | Table S11 | 改为所用规则下未形成外部评价集 |
| 图注含绘图操作指令 | `MANUSCRIPT_FIGURE_CAPTIONS.md`；Figure 5 中文图注 | 删除 |

未改任何科学数据表、标签或对接输出。

## 本轮实际执行的核查

以下命令在仓库内运行；统计类检查只读已提交 CSV，不覆盖正式数值产物。组装与打包会更新正式文稿和投稿包副本，这是本轮修订本身。

```text
python3 docs/assemble_manuscript_zh.py
python3 docs/assemble_manuscript_en.py
python3 data/jcim_novelty_v0/scripts/validate_revision_v1.py
python3 data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py   # 文稿稳定后重写
python3 data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py --check
python3 scripts/audit/audit_submission_five_rounds_v1.py
python3 scripts/audit/pack_submission_v1.py
```

定向阅读并核对（未重跑对接或新统计）：

- Table 2 / 3 锁定：`STATISTICAL_LOCK_V1.md` 与对应 CSV
- API max/median、高置信、五对 dump max/median
- Table S9 GNINA：EGFR `summary_min` 0.220 无 min-of-two CI；[0.109, 0.343] 为 dual–B-only
- 2018 时间切分：仅 JAK1/TYK2、JAK1/JAK2 达报告门槛
- 固定成员脚本：`multiseed_fixed_membership_v1.py` 仅排除空字符串，未另检 NaN/inf；未把该限制写成现有数据错误
- BindingDB：`external_slice_summary_v1.csv` 八对 `packaged_as_external_evaluation=0`

未目视检查位图本身的像素数据；图号与磁盘文件名映射沿用现稿（Figure S5 对应 `FigS6_detectable_effect.png`）。不宣称已完成视觉审核。

## 尚未解决、是否阻止进入英文完善与投稿准备

| 问题 | 是否阻止本轮进入投稿准备 |
|---|---|
| 新增五对生产姿态未入 git，仅有分数表、本机回放和哈希清单 | 否。须在可用性声明中如实写出；投稿前若期刊要求原始姿态，再单独归档。本轮不上传、不申请 DOI。 |
| 公开活性记录存在实验异质性；未做高可比性确认 | 否。已在讨论中作为数据限制陈述，不列为本轮完成前提。 |
| 固定成员脚本未显式过滤 NaN/inf 字面量 | 否。未发现其改变本轮采用的五对结果。 |
| 共晶 QC 方法在原六受体与新增八受体之间仍非完全同一套 CalcRMS 表 | 否。正文已按“各受体对应计算方法”陈述；S2b 为历史对照。不在本轮重对接或重绘 Figure S4。 |
| 多数 `summary_min` 区间含 0.5；对应口袋优势未在留出集稳定出现 | 否。属可如实报告的研究结果与限制，不是稿件错误。 |
| BindingDB 未形成外部评价集 | 否。属研究限制。 |
| 英文尚需投稿前语言润色 | 否。内容已与中文对齐，可进入最终英文完善。 |

无本轮范围内无法掩盖的、会推翻主要结论的数据错误。

## 明确未开展的工作

未开展高可比性确认或实验可比性子集分析。未新增分子对接、新随机种子、新对接软件、新评分模型、MD/MMGBSA/机器学习、新外部评价集、GNINA 固定评分通道新比较，或新的模拟 / 统计检验 / 置信区间方法。重跑旧检查脚本不记为新实验。

## 最终判断

- **稿件内容是否基本稳定：** 是。核心定位、三问结构和主表数字未改；已知事实错位已按现行文件纠正。
- **能否进入最终英文完善与投稿准备：** 能。前提是投稿材料按可用性声明区分已公开与仍待发布内容。
- **提交前仍须处理：** 若目标期刊要求原始生产姿态，需另做归档；英文需常规润色和投稿格式转换。
- **可如实报告的研究限制：** 混合测定异质性、无独立外部评价集、部分区间含 0.5、对应口袋优势未在留出集稳定出现、独立 GNINA 仅三对、五对姿态未公开。

不承诺期刊录用。未新增实验不视为本轮任务未完成。除明确数据错误和必要投稿文件问题外，本轮不再提出新的实验任务。
