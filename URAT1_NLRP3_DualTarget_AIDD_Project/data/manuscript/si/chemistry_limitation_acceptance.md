# 条目2：接受为限制（不改 12/21）

决策：投稿前接受下列历史规则为方法学限制，不版本化重算名单。

1. 历史 40 个药化合格集合含按药名子串排除（FLOXACIN、PIROMIDIC、MK-5108 等，见 `config/chemistry_final.yaml`）。
2. 仅去掉药名排除、保留 PAINS/Brenk/β-lactam 与结构一致性时，会新增 MK-5108、Piromidic acid、Pradofloxacin。见 [对照表](structure_only_chemistry_counterfactual.csv)。
3. Arg 几何匹配器只识别 COOH/carboxylate 氧；156 中无该药效团的分子不通过不能解释为无 URAT1 活性。

12 primary 与 21 reserve 保持冻结。对照表只说明“若当时不用药名排除会怎样”，不是新的候选名单。
