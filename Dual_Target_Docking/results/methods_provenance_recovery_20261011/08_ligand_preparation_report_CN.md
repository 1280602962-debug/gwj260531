# 785 个生产配体的准备证据

解决程度：PARTIAL。785 个配体的文件身份、SDF 与 PDBQT 的重原子计数、以及历史方法记录都已逐个写出。把坐标重新生成出来的那一步没有恢复。文件存在没有被当成准备流程已经完整复现。

逐行表：`07_ligand_preparation_provenance_785.csv`。结构文件位于 `reruns/UNIFORM_RERUN_V4_2_20260921/02_ligands/`。身份对照使用 `00_protocol/global_ligand_registry.csv`。历史方法行来自 `a0ef2726` 的 `data/provenance/ligand_preparation_registry.csv`。

## 四项逐配体结果

结构身份：785/785 的 SDF 去氢后，连接层 InChIKey 与登记 canonical SMILES 的连接层一致。其中 618 个与登记 isomeric SMILES 的完整 InChIKey 一致；167 个连接图一致，而 SDF 的立体层与未写立体的 canonical SMILES 不同。这 167 个是立体信息所在位置不同，不是另一套分子。0 个连接图失败。

SDF 与 PDBQT：785 个 SDF 和 785 个 PDBQT 都在。756 个的非氢原子数与 SDF 重原子数相同。29 个含有 Meeko 的 `CG0` 和 `G0`；去掉 `G0` glue、保留 `CG0` 作为真实碳之后，重原子数与 SDF 相同。0 个计数不一致。这一步是计数核对。29 个大环配体没有逐原子重建 PDBQT 序号到冻结 SDF 的原子映射，因此原子映射仍未完成。

准备方法证据：历史登记里，110 行状态为 `CONFIRMED`，方法字符串是 ETKDGv3、种子 20260727，且这 110 行的靶点字段都是 EGFR/HER2，其中包括仍留在 785 个文件中的 EH120_059。其余 675 行状态为 `METADATA_ONLY`，种子字段写作 `20260727_metadata_only`。登记中的 RDKit 版本字段写明是分析用 2026.03.5，并注明不能把它当成当时对接使用的 RDKit。`ENV_PIN.md` 记录的准备环境是 2026-07-29 的 RDKit 2026.3.1。`CONFIRMED` 表示当时留下了方法快照，不表示该行坐标已经与现在的 SDF 比对过。

原始三维结构能否重建：785 行都是未重建。历史单构象脚本是 ETKDGv3、`randomSeed=20260727`、一次嵌入、MMFF 最多 200 步、Meeko 0.7.1。正式 SDF 是 V4.2 复用的 V4.1 资产，没有逐文件运行日志把坐标连到那次运行。共晶再对接脚本里的 20 构象、种子 42、MMFF94s 然后 UFF，只属于那条再对接路径，不用于这 785 个生产配体。

## 对正式结果的影响

没有替换 SDF、PDBQT，也没有重新对接。连接图一致支持“这 785 个文件对应当前登记分子”。立体层有 167 个只写在 SDF 里，方法快照与坐标之间仍有缺口。大环的 `G0` 计数核对不能代替 Meeko 原子映射的完整复现。
