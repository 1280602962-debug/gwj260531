# 2026-10-26 当前证据整合记录

整合输入提交：7ef3b25fddb32b3d7089db90974fdfccb405847f。活性来源核对来自 196639e9，CAP5 来自 7ef3b25f。主结果继续使用 867 条合格观察，PIK3CA/mTOR 为 Q110；主分数和统计的来源提交仍为 dc028b3350e804cbf8f40acea08261ed4f89c4fa。

本次仅整合已完成的来源核对和支持分析，没有重选主评价成员，没有改变标签、分数、fold、PRIMARY、Stage 2 或主图，也没有重新对接、训练模型、重算 AUROC 或 bootstrap。

删除并替换了当前说明中“Q110 无逐 assay 来源导出”“完整来源仅覆盖旧 805 条”的过时结论。全部 867 条已有两端来源；855 条指定文献标识完整，12 条缺失。原 805 条导出仍是复用记录的原始来源证据，保留。

CAP5 保留 861 条、排除 6 条，仅作支持分析。更新 Results R1/R8、SI、表注、来源索引和权威入口。补充表新增 S19–S31，并把来源与 CAP5 文件纳入当前数据包。主表和主图保持原字节。旧审阅 PDF、补充工作簿、数据包和展示检查报告在原路径替换，避免并列旧版输出。

逐项核验见 results/current/CURRENT_VERIFICATION_REPORT.json。七对历史选样仍为 DOCUMENTED_NOT_REPLAYED，两项完整重放限制没有改写为已解决。

## 投稿表格分类与编号整理

以本次证据整合后的提交 1090ed86 为来源，新增独立 `submission/` 分类：正文 Tables 1–2、补充 Tables S1–S10 和完整补充数据。逐分子长表以数据文件提供，原 S1–S31 全部建立去向对照。新增受体及共晶重对接汇总使用既有配置和记录，未重新计算姿态。

投稿展示使用 S1=M0、S2=M1、S3=M2、S4=M1b、S5=M3。既有正文工作稿、图版和审阅 PDF 暂保留原编号，最终全文引用按新目录对齐。全部科学来源及既有工作文件保留原字节，投稿目录由单独的生成程序和核对报告管理。

## 2026-10-27 Results 核对与投稿图版

以 a7b9d752 为来源，对本轮 Results 文稿完成 205 项数值与汇总核对。正文结果图重新组织为 Figures 1–9，补充图为 Figures S1–S9，评分流程为补充 Scheme S1。所有新图由脚本读取既有数据绘制，提供矢量 PDF、SVG 及 600 dpi PNG。

Table S9 增加 D 面板，直接列出 CAP5 化学模型及 Vina 评分增量；其余投稿表号不变。补充工作簿与数据包同步更新。正式成员、标签、原始评分、分折、PRIMARY/Stage 2 统计及原始图版不变。本轮未重算对接、重评分、RMSD、bootstrap 或化学模型。

原 Results 文稿仍由作者随后修改；具体句子、共晶重对接结果的补充建议及图表引用映射保存在 `submission/04_Editorial_Checks/RESULTS_AND_FIGURE_AUDIT_CN.md` 和 `submission/RESULTS_CITATION_MAP_CN.md`。
