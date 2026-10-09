# ChEMBL 37 活性来源附件

本目录是现行来源附件，全仓库只保留这一处。它不改写冻结类别、mapping、alias、population、fold、PRIMARY、评分主表、权威 YAML 或第二阶段正式结果。

## 版本与日期

官方压缩包：

`https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/chembl_37_sqlite.tar.gz`

SHA-256：`33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`

库内 `version.name` 为 `ChEMBL_37`，创建日期 2026-05-01。`chembl_release` 的最后一条是 `CHEMBL_37`。这个日期是库内发布记录，不是本次下载日期。

压缩包于 2026-10-08 下载、拼接并完成 SHA-256 核对，同日解压并导出。工作库在 `/home/gwj/chembl_37_release/chembl_37/chembl_37_sqlite/chembl_37.db`。完整 SQLite、压缩包、下载分片和解压文件不进入 Git。

指定位置的既有副本检索没有找到可用数据库。一次更深的文件名查找被中断，该范围记为未完成，不能写成整机不存在。恢复使用的是上述官方 ChEMBL 37，不是更新的发行版，也不是在线 REST API。

## 提取条件

`export_chembl37_activity_rows.py` 只读打开上述数据库。805 条由 pair + panel_id 索引。基因按人类、SINGLE PROTEIN、GENE_SYMBOL、单一组分解析，每个基因只对应一个靶点。

活性行使用历史 STANDARD_OK 条件：`pchembl_value` 非空，且 `standard_type` 属于 IC50、Ki、Kd、EC50、Potency、IC50app、Ki app。没有追加 confidence、assay type 或 keep=1。最大值是未取整的 `MAX(pchembl_value)`。文献字段缺失时保留活性行。

前三个面板的历史面板值来自 2026-07-23 ChEMBL REST。本附件是 SQLite 上的 STANDARD_OK 对照，不是重放当时的 REST 请求。Track B 五个面板的历史规则就是该 dump 最大值。

## 覆盖

| 项目 | 计数 | 口径 |
| --- | --- | --- |
| 合格且有类别的 pair×ligand 观察 | 805 | pair + panel_id 唯一 |
| 不重复的分子 | 783 | 805 不是 805 个互不共享的分子 |
| 两端都有 STANDARD_OK 行 | 805 | 任一端为 0 的记录是 0 |
| 原始活性行 | 3555 | A 端 2076，B 端 1479 |
| 不同 activity_id | 3506 | 49 个 activity_id 因配体跨 pair 被重复引用，行未删除 |
| 未取整 MAX 与历史面板一致 | 805/805 | 两端绝对差最大为 0 |
| θ=6.0 回放与冻结类别一致 | 805/805 | |
| 两端全部并列最大值都有指定文献标识 | 794/805 | 见下一节 |
| 原始活性存在、文档无 DOI/PubMed/专利号 | 11 | 见数据集文档 |
| 后续审核差异实质理由待核实 | 3 | 冻结类别仍一致 |

该核对把每条观察对应到 ChEMBL 37 活性、测定和文档记录。它不是逐篇人工阅读原始论文，也不表示所有记录来自同一测定平台。

## 794 的口径

指定文献标识指 DOI、PubMed 或专利号至少一类非空。

794 的定义是：每一端的全部并列最大值记录都有指定文献标识。不是“每一端至少一条”。导出脚本用同一规则：并列集合中只要有一条缺少这三类标识，该端就记为 missing。

在本附件里，按“全部并列记录”和按“至少一条”计数都是 794。没有出现同一端并列记录中一部分有标识、一部分没有的情况。数字相同并不改变口径；报告使用全部并列记录。

## 11 条数据集文档

这 11 条两端都有原始活性。缺的是 DOI、PubMed 和专利号，不是活性行。两个文档在 ChEMBL 37 中都是 DATASET：

- CHEMBL1201862，标题 PubChem BioAssay data set，来源 PubChem BioAssays（src_id 7，PUBCHEM_BIOASSAY）。链接 `https://www.ebi.ac.uk/chembl/explore/document/CHEMBL1201862`。涉及 EH40_02、EH40_04、EH40_18、EH40_25、EH40_30、EH40_38、J1J2_039、J1TYK2_035、J1TYK2_040。
- CHEMBL1909046，标题 DrugMatrix in vitro pharmacology data，来源 DrugMatrix（src_id 15，DRUGMATRIX）。链接 `https://www.ebi.ac.uk/chembl/explore/document/CHEMBL1909046`。涉及 EH120_077 的 B 端和 EH120_102 的两端。

逐条判定在 `document_identifier_coverage_805.csv`。

## 三条历史审核差异

只读核对分支 `cursor/methods-sentence-audit-c7cc` 的提交 `17435413410290960da3437de640509705f00b51`。没有运行 `adjudicate_activity_records.py`。

该脚本先只接受 `high_confidence_activity_audit_v1.csv` 的 keep=1，再按 `primary_only_excluded_activity_rows_v1.csv` 排除。匹配键是 ligand、靶点、document_chembl_id 和四位小数 pChEMBL。下面三条的原库最大值本身是 keep=1，审核表的 `exclusion_reasons` 为空；数值变化发生在后一步排除表。因此不能把差异全部归因于 keep=1。

排除表给出的理由代码只有 `source_adjudication`，没有具体依据。排除步骤可追溯；实质审核理由待核实；冻结类别一致。这些原始记录没有被认定错误或不可靠。

| 记录 | 原库及历史面板 A/B | 后续审核聚合 A/B | 冻结类别 | 改变最大值的排除 |
| --- | --- | --- | --- | --- |
| PM48_04 | 8.85 / 9.35 | 8.40 / 9.35 | dual | A 端 activity 25847801，assay CHEMBL5505638，文档 CHEMBL5500428，IC50 = 1.4 nM。排除后 keep=1 最大值为 8.40（activity 26224982）。B 端另有 8.37 被排除，但最大值 9.35（activity 26139402）仍保留。 |
| PM48_05 | 10.00 / 8.52 | 9.12 / 7.32 | dual | A 端 activity 24958888，assay CHEMBL5215112，文档 CHEMBL5214883，IC50 = 0.1 nM。排除后两个 keep=1 最大值都是 9.12（activity 18088188、24954253）。同端 7.48 也被排除，但它不是原库最大值。B 端 activity 19139806，assay CHEMBL4376361，文档 CHEMBL4373732，IC50 = 3 nM。排除后两个 keep=1 最大值都是 7.32（activity 18088212、24958369）。 |
| PM48_22 | 8.77 / 5.83 | 8.77 / 5.52 | A_only | A 端 activity 26030793 未进入排除表。B 端 activity 22965963，assay CHEMBL4769547，文档 CHEMBL4765307，IC50 = 1461 nM。排除后 keep=1 最大值为 5.52（activity 22806053）。 |

活性、测定、文档和后继记录见 `historical_adjudication_differences.csv`。这三条例外不否定 805 条原始来源覆盖，也不表示全部问题已经解决。

## 本目录文件

- `ChEMBL37来源恢复说明.md`：本说明
- `chembl37_dual_endpoint_crosswalk_805.csv`：805 条两端最大值、冻结类别和回放类别
- `chembl37_activity_rows_long.csv`：3555 行原始活性
- `document_identifier_coverage_805.csv`：794 与 11 的判定
- `historical_adjudication_differences.csv`：三条审核差异
- `export_chembl37_activity_rows.py`：本次实际使用的只读导出脚本
