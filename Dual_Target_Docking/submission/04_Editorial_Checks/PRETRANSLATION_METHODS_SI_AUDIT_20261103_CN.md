# 英译前方法与 SI 文字核查报告

核查日期：2026-11-03。仓库：1280602962-debug/gwj260531。输入分支：`cursor/v4-2-final-audited-results-20260930`；输入提交：`d1ca17512eabe5c9d102ebe75c3f9c59668452dc`。

本报告先于正文及 SI 修改写入。仅进行已有文件、文本、引用、字段和文件字节核查；不执行对接、相似性计算、模型拟合、AUROC、bootstrap 或其他新的科学分析。以下是修订前结论及限定修改清单；修订后的逐处差异和验证另附。

## 一、逐项结论与处理

| 请求 | 已核实事实 | 限定处理 | 状态 |
| --- | --- | --- | --- |
| 评分方案编号 | `submission/Supplementary_Methods_CN.md` 原第 25 行写作“S2 为 Vina mode 1”，与 `submission/SCORING_SCHEMES.csv` 冲突 | 写明 S1 为 Vina mode 1；S2 为该姿态的 GNINA CNNscore；S3 为固定姿态 RTMScore；S4 为全部保存姿态的最高 CNNscore；S5 为独立 GNINA 对接首位姿态的 CNNscore | 可修订 |
| 跨折最近邻面板 | 当前提交及两个补充 ZIP 未找到阶段 A 跨折最近邻原始结果及对应计算脚本。现有相似性字段是同一 pair×arm 内的最大异类相似性 | 停止新增 Table S7 Panel B；保留原 Table S7 及其全部数值，不从历史文字总结编造统计 | 因证据缺失未新增 |
| CNNaffinity | 已确认保存于正式 GNINA 长表、独立对接 master 和 Q110 日志。正式评分读取 CNNscore | 在补充方法写明“原始记录中另保存了 CNNaffinity，但未将其纳入本研究的正式评价” | 可修订 |
| 共晶二维相似性 | 已取得的仓库 Results、SI、图注及 2026-10-29 中文稿未发现相关结论；未找到可核实的共晶二维相似性结果 | 不新增计算，不补写结论；保留重对接 RMSD 文件、表格和图。Discussion 正文未取得，不能声称已完成清理 | 现有文本未见待删除句；Discussion 未完成 |
| 全文引用及方法口径 | 投稿图表已采用新编号，但仓库 Results 和 Supporting Information 工作稿仍使用旧编号及 M 名称 | 按当前 TABLE_INDEX、FIGURE_INDEX、RENUMBERING 和 FIGURE_RENUMBERING 更新活跃 Markdown 文稿引用及评分显示名；保留源文件内部字段与原始文件名 | 可修订 |

## 二、评分方案及 CNNaffinity 证据

| 投稿编号 | 来源方法 | 正式评价读数 |
| --- | --- | --- |
| S1 | M0 | Vina seed 42 mode 1，取 affinity 相反数 |
| S2 | M1 | 同一 Vina mode 1 姿态的 GNINA CNNscore |
| S3 | M2 | 同一 Vina mode 1 姿态的 RTMScore model 1 |
| S4 | M1b | 同次 Vina 搜索全部保存姿态中的最高 GNINA CNNscore |
| S5 | M3 | GNINA 独立对接、按 CNNscore 排序的首位姿态 CNNscore |

上述映射与 `submission/SCORING_SCHEMES.csv`、`submission/README.md`、`submission/TABLE_CAPTIONS.md`、`submission/FIGURE_CAPTIONS.md`、`CURRENT_RESULT_AUTHORITIES.json` 以及冻结协议一致。来源脚本和 CSV 中的 M 名称是既有内部标识，不作重命名。

CNNaffinity 的直接证据（以下均相对项目根目录 `Dual_Target_Docking/`）：

1. `reruns/UNIFORM_RERUN_V4_2_20260921/09_secondary_scoring/GNINA_VINA_POSE_RESCORE_POSE_LONG_VERIFIED.csv`：表头包含 `CNNscore,CNNaffinity`，首个实际姿态同时有两个非空数值；不是仅凭字段名推定保存。
2. `reruns/UNIFORM_RERUN_V4_2_20260921/10_gnina_independent_docking/GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv`：同样包含实际 CNNaffinity 读数。
3. `extension_jcim_validation/pik3ca_mtor_q110_v42/gnina_rescore/4JT6__QXDJAGKQYJREHU-UHFFFAOYSA-N/vina_mode1.gnina.log`：第 17–18 行分别记录 CNNscore 和 CNNaffinity。该任务的 `status.json` 只提取 M1/M1b CNNscore，不能把 CNNaffinity 未进入状态汇总误写成原始日志未保存。
4. `scripts/rerun_v4_2/formal_metrics_lib.py` 的 `method_valid` 和 `method_score`：M1、M1b、M3 分别读取 `M1_CNNscore`、`M1b_CNNscore`、`CNNscore`。
5. `extension_jcim_validation/pik3ca_mtor_q110_v42/q110_continue.py` 中的正式人群构建同样读取上述 CNNscore 字段；`q110_stats.py` 读取固定的 M 方法评分列。没有新增 CNNaffinity 评价通道。

补充说明只确认“原始记录中保存”，不声称每个失败任务或每个汇总状态文件都有 CNNaffinity，也不将 CNNscore 称为实验亲和力。

## 三、Table S7 Panel B 的证据门槛及实际缺口

检查范围：当前提交的完整受版本管理项目树；`manuscript/Supplementary_Source_Data.zip` 全部 91 个条目；`submission/03_Supplementary_Data/Supplementary_Data.zip` 全部 143 个条目。检索覆盖路径及文本中的 interpretation audit、stage A、cross-fold、nearest neighbor、跨折、最近邻及相应 CSV 表头。另查找了已有 GitHub 分支线索及相关文件记录，未取得该分析的可复核文件。

现有可用文件：

- `results/jcim_stage2/module_D_chemistry_oof_ligands.csv`：正式折外预测、fold_id 和成员。
- `results/canonical/model_fold_assignments.csv`：固定分折。
- `results/jcim_stage2/module_D1_composition.csv`：字段 `max_ecfp4_tanimoto_other_class`，属于标签层的跨类别相似性。
- `analysis_plan_jcim/CHEMISTRY_BASELINE_FREEZE.md`：明确上述相似性在同一 pair×arm 内对另一类别计算。
- `submission/07_Figure_Source_Data/Figure_S6_Cross_Class_Similarity.csv`：对应现有跨类别图，不能替代训练折最近邻分析。

缺少的必要材料：

1. 阶段 A 的逐 OOF 成员结果：pair、arm、配体 ID、测试折号、训练折最近邻 ID 和 Tanimoto 数值，且可与正式 OOF 人群对应。
2. 与该逐成员表一致的已保存分布汇总及分母定义（如是否按方向重复计数 dual）。
3. 生成上述结果的原始脚本、指纹参数、训练集成员定义，以及输入路径和版本记录。

因此本轮不计算相似性，不计算分位数或阈值比例，不新增 Panel B，不将此前对话中的汇总数字写入 SI。这个结论是“当前可访问材料不足”，不推断历史上从未执行过该分析。

## 四、Discussion、正文 Methods 和共晶分析的可访问范围

当前提交未包含 `Discussion_CN.md` 或正文 `Methods_CN.md`。已读取的 `双靶对接第二版中文3(1).docx`（2026-10-29）包含引言、Methods 2.1–2.9 和 Results 3.1–3.7，末尾为联合排序结果，没有 Discussion 正文。本轮不拿旧分支讨论或讨论大纲冒充当前 Discussion。

该中文稿的评分方案定义与正式映射一致，Methods 对跨类别相似性的描述与 D1 文件定义一致，七对固定折与 PIK3CA/mTOR 骨架循环分折分别表述；当前引用采用投稿图表号。此文件仅作为正文方法核查参考，本轮不替换或重排 Word 文件。

取得的文本没有需要删除的共晶配体二维相似性实证句。当前共晶重对接证据是 `13_qa/redocking_master.csv`、Table S3 和 Figure S1，描述三维姿态 RMSD；它们不提供二维化学相似性或模板相似性偏差的结论。所有这些既有结果均保持不变。

Discussion 后续合稿时，跨折相似性及共晶二维相似性的实证结论仍不能引用不存在的 Panel B 或仅凭文字汇总保留。本轮未取得正文，故该部分列为未完成核查，不写成“全文全部通过”。

## 五、限定文件修订清单

| 文件 | 修改内容 | 保持不变 |
| --- | --- | --- |
| `submission/Supplementary_Methods_CN.md` | 修正并完整列明 S1–S5；补 CNNaffinity 未纳入评价；明确方向记号、summary_min 和跨类别相似性定义 | 软件参数、正式成员和全部科学数值 |
| `manuscript/Results_CN.md` | 更新旧图表引用及 M 显示名；明确点估计措辞 | 每个科学数值、方向性结论及其证据范围 |
| `manuscript/Supporting_Information_CN.md` | 将旧表号转为当前表号或补充数据文件名；统一评分显示名；补相似性定义和 CNNaffinity 说明 | 正式分析、历史支持分析事实及内部字段 |
| `manuscript/README.md`、`submission/README.md` | 说明活跃 Markdown 已对齐，历史工作簿/PDF/旧图注仍为来源版本；记录本轮证据缺口和报告入口 | 旧来源文件及投稿图表 |
| `submission/04_Editorial_Checks/` | 保存本报告、逐处差异及非科学计算的验证记录 | 不改旧审计记录 |

Table S7 的 Markdown、Excel、编号 CSV 和数据包不新增面板，不改任何单元格。图像、RMSD、人群、SMILES、标签、fold、评分、OOF 预测、PRIMARY/Stage 2/CAP5 指标、AUROC 及置信区间均不改。

旧 `manuscript/Figure_Captions_CN.md`、旧工作簿和审阅 PDF 是来源版本，不能与新的投稿图表混用；当前图注和表注分别以 `submission/FIGURE_CAPTIONS.md` 与 `submission/TABLE_CAPTIONS.md` 为准。它们保留历史编号是版本身份，不是本轮将评分方案重新编号。源文件中的内部键（例如 `pairwise_common_vs_M0`）保持原样。

## 六、限定修订执行结果

以上清单已执行到五个既有 Markdown 文件，共记录 62 次文字或引用替换操作。具体原文、修改后文字及理由见 `PRETRANSLATION_TEXT_CHANGES_20261103_CN.md`；结构化记录见同名 JSON。没有改动任何表格单元格、生成脚本、数据文件、模型结果、图或 ZIP。

另发现一项引用说明错误：现行投稿数据包中的 `Data_Index.csv` 表头为 `path_in_archive,bytes,data_rows`，并不提供旧 SI 所称的“包内 SHA256 清单”。已按实际字段改写文字；未创建或改动数据包。报告中的 SHA-256 验证是本轮对已有文件的字节核查，与包内索引分别说明。

`PRETRANSLATION_TEXT_VALIDATION_20261103.json` 中 53 项限定检查全部通过，涵盖：

- 修改仅涉及五个既有 Markdown 文件，科学来源路径没有变动。
- 17 个关键文件修订前后 SHA-256 相同，包括人群、标签映射、分折、OOF、PRIMARY 指标与差值、D1、重对接 RMSD、Table S7、整个 SI 工作簿及两个数据包。
- Results 与 SI 中原有科学小数值均保留；Results 仅补出原句已表达的 0.5 参考线，不产生新结果。
- 当前 Table/Panel 引用均能在 TABLE_INDEX 中定位，正文九张图及对应补充图已有引用；删除旧 S11–S31 长表引用，改为现行面板或数据文件。
- 明示的数据 CSV 名称均存在于当前补充数据包；真实内部字段 `pairwise_common_vs_M0` 保留。
- S1–S5 权威定义未改；错误语句已消除；CNNaffinity 和 summary_min 说明已补充。

“通过”只表示上述限定核查通过，不表示已取得缺失的跨折分析或已审校当前 Discussion。未新增 Table S7 Panel B；Discussion 仍属于未完成核查项。本轮未进行任何新的科学计算。

## 七、合并时的版本补记

本报告第 1–6 节对应输入提交 d1ca17512eabe5c9d102ebe75c3f9c59668452dc 的核查及文字修订。合并到原分支时，该分支已前进到 485ff1bf8e046ac3d8fb5ffdc7b75a9725d6476b，新增 `analysis/posthoc_diagnostics_20261103/` 支持性诊断目录。

已从该提交直接读取 `outputs/A_TEST_NEAREST_TRAIN.csv` 表头及数据行，确认存在测试折成员、训练折最近邻 ID 和相似性字段；同目录也包含 `scripts/analysis_a.py`、配置及报告。因此，前文“未取得跨折最近邻原始文件和脚本”的描述是旧基线的历史核查结论，不再代表合并后文件的可用性。此次合并保留全部 65 个新增支持性诊断文件，不运行它们，也不据此声称已独立核实该分析的全部数值。Table S7 Panel B 的编号整合尚未执行。

本次合并没有新科学计算，原分支已有的正式结果和支持性诊断文件均按原字节保留。Discussion 正文仍未取得。
