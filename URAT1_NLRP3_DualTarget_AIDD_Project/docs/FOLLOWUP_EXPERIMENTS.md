# 后续需要补齐的内容与实验

状态：计划。`md_authorized` 仍为 false。没有轨迹、占有率或 MM/GBSA 数字。不改 12/21。

详细 MD 参数见 [`config/md_protocol.yaml`](../config/md_protocol.yaml)。PF 选择理由见 [`NOMINATION_RATIONALE.md`](NOMINATION_RATIONALE.md)。

## A. 投稿前已闭合（计算稿）

| 条 | 处理 |
|---|---|
| 2 药化 | 接受限制；SI 反事实表保留；不改 12/21 |
| 4 文献 | 三篇原始论文核对；compound 32 SMILES 改引 Nature SI；删除顶膜暴露与新靶点确认 |
| 3 GNINA | 已有 log 汇总表；无 log 标 unknown；未见 1.3.1 |
| 1 模型 | 冻结 1588；1377 仅 SI；不编多线程因果 |

## B. MD（审计第 5 条）：现有数据不能闭合

授权前不要跑。跑之前按 yaml 执行，跑完之前正文不写任何 ns、RMSD、占有率或结合能。

四个体系，各 3×200 ns（seed 42/43/44）：

1. 9DKB + lesinurad **晶体酸根姿** + POPC 膜（对照）
2. 9DKB + PF-03882845 **seed42 CNNscore Top-1（mode 1）** + POPC 膜
3. 7ALV + NP3-146 **共晶姿**（对照，无膜）
4. 7ALV + PF-03882845 **seed42 CNNscore Top-1（mode 1）**（无膜）

读出：配体 RMSD（受体坐标系）；URAT1 酸氧–Arg 氮距离分布与 Phe 笼接触；NLRP3 关键残基接触与口袋滞留。MD 中盐桥定义为酸氧–Arg NE/NH1/NH2 **< 4.0 Å**，与静态门 7.7027 Å 分开。分析窗口为生产段最后 150 ns。

对照失败则该靶侧不解释候选：

- URAT1：lesinurad 分析窗口内重原子 RMSD > 5 Å 的时间超过 50%，或酸–Arg < 4.0 Å 占有率 < 30%
- NLRP3：NP3-146 同样 RMSD 规则，或 7 个关键残基中持续接触少于 3 个

MM/GBSA 只作池内支持分析，不改名单、不当亲和力。

## C. 实验（MD 不能替代）

三个问题必须分开做，不能用 IL-1β 或肾病模型互相替代。

### C1. URAT1 功能（直接转运）

- HEK293（或同等）过表达人 URAT1，14C- 或荧光尿酸摄取，PF-03882845 剂量–反应。
- 对照：lesinurad 和溶剂；报告 IC50 或在最高可溶浓度下的抑制百分率。
- 阴性结果写成“本测定未检出抑制”，不写成无活性于一切 URAT1 测定。

### C2. NLRP3 直接结合 / ATPase（相对通路）

- 无细胞 NACHT 结合（SPR/MST 等）或 ATP 酶；和/或
- MR 敲低或 MR 拮抗剂背景中的巨噬细胞 IL-1β（LPS+ATP 或 MSU）。
- 仅有 IL-1β 下降且未做结合/ATPase，只能写通路水平，不能写 NACHT 配体。
- 已知混杂：Orena 等 2013 显示 PF 抑制醛固酮模型肾脏 *Il-6* 等基因。

### C3. 不在本项目用适应症代替暴露

- 糖尿病肾病或 Dahl 大鼠肾脏保护 **不是** 近端小管顶膜游离药物浓度证据。
- 若将来需要暴露论证，应单独设计肾小管/尿液或原位灌注实验；当前投稿不写顶膜暴露。

### C4. 任选、不作主线

- OAT1 对抗筛选（排除注释，不上标题）。
- 12 个 primary 的购样与核磁/LCMS 确认。
- 专利检索（PubMed 0 ≠ 自由实施）。
