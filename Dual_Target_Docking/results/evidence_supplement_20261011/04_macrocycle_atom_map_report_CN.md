# 29 个大环配体的原子映射

依据是正式冻结的 SDF、PDBQT，以及 PDBQT 里的 Meeko `REMARK SMILES` / `REMARK SMILES IDX`。映射不用“重原子个数相同”，也不用“元素序号碰巧一致”。对每个配体，在备注 SMILES 的全部同构中选取坐标偏差最小的那一个。逐行结果在 `03_macrocycle_atom_map_29.csv`。既有计数审计是 `13_qa/GNINA_LIGAND_REPRESENTATION_AUDIT.csv`，闭环的历史判定是 `13_qa/M3_MACROCYCLE_TOPOLOGY_DECISION.csv`。

29 个配体都有 2 个 `CG0` 和 2 个 `G0`。`CG0` 映射到 SDF 的碳，`G0` 不在 SMILES IDX 里。`G0` 落在互补 `CG0` 的坐标上。备注 SMILES 与 SDF 的连接图一致。两个 `CG0` 所对应的 SDF 原子之间有化学键。元素类型与映射到的 SDF 原子一致，不一致数为 0。

扭转树里，这两个真实碳原子都没有直接 `BRANCH`，也不在同一个刚性块中。这是 Meeko 用伪原子表示闭环的方式：化学键在 SDF 里，对接树用 `G0` 维持闭环。不能把“扭转树没有这条边”写成原子丢失。

已有的 M3 判定是另一件事。把 `G0` 删掉并把 `CG0` 改成 `C` 之后，转换后的 GNINA 输入不再保留这条闭环边。该结论已经写在 `M3_MACROCYCLE_TOPOLOGY_DECISION.csv` 的 29 行里，状态都是 `M3_INPUT_TOPOLOGY_FAIL`。它限制的是这条转换后的 GNINA 输入，不是冻结 PDBQT 的原子身份。

唯一的坐标例外是 AB_082（`SDHJEPJJARZDKU-UHFFFAOYSA-N`）。同构只有一个，闭环键仍在，元素也一致，但 PDBQT 序号 6、7、8、9 四个碳与 SDF 对应原子的距离约为 2.54、2.54、2.51、2.51 Å。其余 28 个配体的最大原子偏差小于 0.001 Å。AB_082 的图没有换，变的是这四个重原子的坐标。本次没有替换 PDBQT 或 SDF。
