# Q110 新增分子的 ChEMBL 37 活性来源审计

审计提交是 `5d992b95ea121b482a9144f6ea7543f71422e975`。冻结名单 SHA-256 仍是 `1e85e8c25af4dbe14e756f7042993ae86ac7f89ae06f810379372d739246b550`。名单仍是 110 行、110 个 ChEMBL ID、110 个母体 InChIKey，类别 32/32/32/14，其中 9 个与 PM48 重叠，101 个是新增。本轮没有修改这些值或类别。

## 数据库

官方压缩包 `chembl_37_sqlite.tar.gz` 的 SHA-256 是 `33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`，与既定指纹一致。只读 SQLite 的 `version.name` 含 `ChEMBL_37`，`chembl_release` 最后一条是 `CHEMBL_37`。本轮没有重新解压这 29 GB 文件。

用同一套 `STANDARD_OK` 条件重查了已保存的 805 条观察：pChEMBL 非空，类型为 IC50、Ki、Kd、EC50、Potency、IC50app 或 Ki app，人源、`SINGLE PROTEIN`、单一蛋白组分。没有加入 confidence、assay type 或 keep=1。重查得到 3555 行，A 端 2076、B 端 1479、不同 activity_id 3506，与 2026-10-08 文件的行、activity_id 和两端最大值全部一致。因此数据库一致性记为 `CONFIRMED`，不是只凭版本名作出的判断。已有 757 条正式来源表没有重写。

PIK3CA 对应 CHEMBL4005（P42336），mTOR 对应 CHEMBL2842（P42345）。每个基因只有一个符合上述条件的靶点编号。PIK3CB、PIK3CD、PIK3CG 没有并入。

## 新增 101 人

101 个分子都能在 ChEMBL 37 中按 ChEMBL ID 找到，没有按 SMILES 替换。两端都有 `STANDARD_OK` 行，共 238 条原始活性，238 个不同 activity_id。

两端未取整最大 pChEMBL 与冻结值的绝对差都小于 1e-9，101/101 精确一致。按 θ=6.0 回放的四状态与冻结类别 101/101 一致。输出里个别数值显示为 `9.699999999999999` 这一类二进制浮点写法；比较使用原始差值，没有用四舍五入或 0.015 容差把差异判成一致。0.015 只作为与旧审计对照的参考列，本轮没有分子需要靠它通过。

文献覆盖沿用旧口径：每一端的全部并列最大值都要有 DOI、PubMed 或专利号之一。101 人中 100 人满足。不满足的是 Q110_A21：PIK3CA 最大值来自 CHEMBL3217752，有 DOI `10.1039/C2MD20044A`，IC50 6.6；mTOR 最大值是 activity 3919485，文档 CHEMBL1201862，来源 PUBCHEM_BIOASSAY，类型 DATASET，Potency 5.43，没有 DOI、PubMed 或专利号。该分子仍保留在冻结集中。

在 101 人中，6 人的两端最大值来自不同文档，14 人存在 IC50、Ki、Kd、EC50 或 Potency 的终点混合，3 人（Q110_A20、Q110_B09、Q110_B29）在同一端同时有不低于 6.0 和低于 6.0 的记录。没有分子的单端极差达到 2 个 pChEMBL 单位。这些标记可以同时存在，不表示最大值规则失效。

## 重叠的 9 人

9 个分子都映射回 2026-10-08 的 PM48 来源行。冻结的两端最大值和四状态与当时的 ChEMBL 37 最大值、θ=6.0 类别一致。

Q110_A24 / CHEMBL2396661 / PM48_22 的 ChEMBL 37 mTOR 最大值仍是 5.83，来自 activity 22965963、assay CHEMBL4769547、文档 CHEMBL4765307，IC50 1461 nM。后续审核聚合把这条记录以 `source_adjudication` 排除后，mTOR 最大值变为 5.52（activity 22806053）。排除步骤可以追溯，实质审核理由仍然没有具体依据。冻结类别在两种最大值下都是 A-only。本轮把这件事记入 `Q110_SOURCE_EXCEPTIONS.csv`，没有把 5.83 改成 5.52。
