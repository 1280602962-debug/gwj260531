# 受体与配体准备

命令和文件见同目录的三张表。这里只写可比性。

## 受体

正式生产受体是 Open Babel pH 7.4 准备的 `9DKB_receptor.pdbqt`（620900 字节，7759 行）和 `7ALV_receptor.pdbqt`（1480113 字节，4519 行）。水、晶体配体和其它 HETATM 被去掉。PDBQT 里的部分电荷全是 0。组氨酸互变异构没有记录，不能用软件默认值冒充。

phase-1 的 9839 分子基准用的是另一份 Maestro 准备的 9DKB，375040 字节、4688 行，带非零电荷。Maestro 的具体命令无法恢复。两份受体的分数不能放在一起比。

9DKB 是内向开口的 lesinurad 结合构象，分辨率 2.74 Å。`docking_ensemble_pdb.csv` 里写成 2.55 Å，与 CIF 和 Maestro 头不一致，以 2.74 Å 为准。文献 Arg477 在准备后的文件里是 A 链 476。7ALV 准备后被重新分段，文献残基必须用 `data/manuscript/tables/residue_mapping.csv`。8ETR 没有生产受体。

## 配体

三套准备不能混用。

1. 酸性战役和第三阶段 B：Dimorphite-DL pH 7.4，优先去质子羧酸，ETKDGv3 种子 `0xC0FFEE`，MMFF 最多 200 步，Meeko 默认写出 PDBQT。156 个和 32 个都是这条规则，失败数都是 0。历史那次运行的 RDKit、Dimorphite 和 Meeko 版本没有写入 `prep_summary.json`。
2. Rank 路线的 1583 个分子：不加 Dimorphite，直接给输入 SMILES 加氢。5 个嵌入失败。
3. phase-1：LigPrep 三维构象再交 Meeko。9849 进、9839 出。pH 和互变异构设置无法恢复。

第三阶段 B 的 32 个与历史 156 个因此可以在规则层面比较。不能和 Rank 路线或 phase-1 的分数比较。

32 个里 6 个选定微物种仍有未指定立体中心：RAZUPENEM、ISOPROPICILLIN、BENFURODIL HEMISUCCINATE、FENBENICILLIN、CEDIROGANT、EVATEGREL。5 个净电荷为 0，是羧酸根加铵的两性离子，羧酸仍然去质子。BENFURODIL HEMISUCCINATE 的名字像盐，SMILES 是单一组分。这些都只记录，没有因为对接结果改选微物种。156 个里 76 个有多于一个微物种，同样只对接了优先的羧酸形式。清单在 `PROTONATION_AND_STEREO_AUDIT.csv`。
