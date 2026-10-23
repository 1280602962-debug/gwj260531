# 其他七对靶点的选样证据审计

原始沉积面板和候选 JSON 来自 `a0ef2726036e35f5bf2b78412f9a3ba43506e8e9`。抽样脚本位于更早提交 `d824eb6b2b06fd56a9da327f8ab67a29b9eb1b78`。当前合格人数和身份以 V4.2 的 mapping 为准。后来的脚本说明只用来解释沉积名单，不把脚本里事后加上的排序说成原始抽样已经使用。

七对全部是 `DOCUMENTED_NOT_REPLAYED`。没有一对被写成 `EXACT_REPLAY`，也没有一对达到 `CONFIRMED_SELECTION_MISMATCH`。

## Track B 五对

F2/F10、JAK1/JAK2、JAK1/TYK2、PPARG/PPARA、PPARA/PPARD 的沉积面板各 110 人。当前合格观察各 110 人，ChEMBL ID 与沉积面板一致，沉积表中的 `theta6_class` 与当前 θ=6.0 类别一致。

`track_b_panel_summary_v1.csv` 记录种子 `20260729` 和配额 32/32/32/14。构建规则是严格 6.5/5.5，灰区不进入主面板；正式标签仍是 θ=6.0。这五对沉积成员的 strict 类与 θ=6.0 类相同。`extract_track_b_panels_v1.py` 写明使用 `random.Random(20260729).shuffle`，并在 2026-09-05 注明：交集排序是在五份已提交面板生成之后才固定的，不能重跑脚本来替换这些面板。分子级候选池缓存不在 Git 中。因此本轮没有重跑 Track B。小分子条件以该脚本记载的 MOL、Small molecule、分子量 150–750、重原子 10–60、无金属为准；脚本写明没有 Bemis–Murcko 上限。

## EGFR/HER2

沉积 `panel_v0_120.csv` 有 110 人，类别 28/38/32/12，阈值是 θ=6.0。其中 40 人来自 panel40，70 人的 role note 是 `stage1 expand hardneg-priority`。EH40_01 是具名的 TAK-285 构象锚点，EH40_02 是 lapatinib。另有一条 dual 记录写明替换了 imatinib 标签噪声候选。这是面板上已经写明的优先保留，不是本轮根据对接分数新发现的删减。

当前合格观察 109 人，类别与去掉 EH120_059 之后的沉积类别一致。EH120_059 / CHEMBL305194 / InChIKey `QQRLVNGJKMMEMW-POHAHGRESA-N` 在沉积面板中是 A_only（7.54/4.31）。`LABEL_PROVENANCE.md` 将它标为 `unresolved_missing_arm`，当前 `activity_eligible=0`。它不是 neither，也不是因为对接失败而被换成别人。

`build_egfr_her2_panel120.py` 的种子是 20260727，但分子信息可以来自当时的 ChEMBL API。本轮没有调用 API，所以没有做独立重放。

## AChE/BChE

沉积 `panel_v0_strict.csv` 有 100 人，规则是严格 6.5/5.5，配额 28/28/28/16。`build_strict_panels.py` 把 `gray_excluded=True` 写到每一行输出上，所以这一列不是逐分子排除标记。多样性键是 ChEMBL ID 前 8 个字符；脚本注释写明当时没有 SMILES，因此它不是 Bemis–Murcko。

本轮用 `a0ef2726` 中的 `mols_ACHE.json` 和 `mols_BCHE.json` 按该逻辑只读重放了一次，种子仍是 20260729，没有第二次尝试。重放得到 100 人，与沉积面板有 99 人相同。重放多出 CHEMBL173607，沉积面板多出 AB_056 / CHEMBL3960861，二者都是 strict A_only。AB_056 仍在当前合格集中。

当前合格观察是 98 人。沉积 100 人里的 AB_040 / CHEMBL222904 是 AB_046 / CHEMBL254016 的草酸盐记录，当前 mapping 把 AB_040 标成 AB_046 的 alias，不另算一条观察。AB_087 / CHEMBL344593 是 `unresolved_missing_arm`，`activity_eligible=0`。其余 98 人的 strict 类与当前 θ=6.0 类一致，因为满足 6.5/5.5 的四状态在 θ=6.0 下不会改变。人数变化由别名和技术排除闭合，不构成无法解释的成员替换。

这七对都没有发现类似 PM48 那种“先放入指定 ChEMBL 种子、种子还可以突破骨架上限”的清单。EGFR/HER2 有已经写在面板上的具名锚点和 hard-negative 优先扩充，需要在方法中说明，但不改变当前 109 人的身份。
