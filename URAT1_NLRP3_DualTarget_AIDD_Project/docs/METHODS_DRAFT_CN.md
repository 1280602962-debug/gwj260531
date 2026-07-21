# Methods 正文草稿（中文）

> 接引言 [`INTRO_DRAFT_CN.md`](INTRO_DRAFT_CN.md)。只写设计与操作；协议胜出结果、富集终值与候选名单写入 Results。  
> 作者备忘（勿进投稿正文）：先验三阶段、勿写换方法日记；旧 Glide 稿见 `MANUSCRIPT_DRAFT_CN.md`（已过时）。

---

## 2.1 总体设计

本研究将计算流程预先划分为三个衔接阶段，并在分析开始前固定各阶段的判定规则。

第一阶段为对接协议筛选。以 URAT1 的 inward-open 结构（PDB 9DKB）为受体，分别构建性质匹配的 TrueDecoy 基准与等量随机的 RandomDecoy 对照；在预先定义的开源对接与重打分协议集合上计算富集及已知尿酸药回收，按下文“协议选优规则”选定用于后续生产排序的协议，记为 \(\Pi^{\ast}\)。

第二阶段为不对称临床库漏斗。对 ChEMBL 临床阶段化合物库，先由 NLRP3 分类模型给出活性概率并按预设阈值缩库；对进入对接池的分子，在 9DKB 与 NLRP3 NACHT 结构（PDB 7ALV）上按 \(\Pi^{\ast}\) 计算双靶对接分，再转换为池内百分位得分 \(S_{U}\) 与 \(S_{N}\)，据此求 Pareto 非支配前沿。

第三阶段为成药性审计与提名。对合并结果施加结构警报、类药性、化学适用域及模型稳健性检查。Pareto 前沿命中不自动等于最终提名：仅通过过滤器的分子进入假说短名单。

URAT1 回归模型仅用于基准回收对照。若其对已知尿酸药的回收未达到预设门槛，则不参与第二阶段的临床库主排序。全文对接相关分数只用于同一配体池内的相对比较，不换算为实验亲和力。

---

## 2.2 数据来源与分子预处理

URAT1 活性记录取自 ChEMBL。清洗关系符、统一浓度单位并按分子聚合后，以负对数形式定义活性：

\[
\mathrm{pActivity}=9-\log_{10}(\mathrm{IC}_{50,\mathrm{nM}}).
\]

若原始单位为微摩尔，则等价采用 \(\mathrm{pActivity}=6-\log_{10}(\mathrm{IC}_{50,\mu\mathrm{M}})\)。TrueDecoy 活性集默认取 \(\mathrm{pActivity}\ge 6\)。

NLRP3 训练表由多 assay 条件下的分类标签构建。临床重定位库由 ChEMBL 临床阶段小分子整理得到（\(n=8319\)）。诱饵候选池取自蒸馏子集中的未标记多样性分子。全部结构经 RDKit 解析、消毒并转为规范 SMILES；化学骨架按 Bemis–Murcko 定义。

---

## 2.3 TrueDecoy 与 RandomDecoy 基准

以 URAT1 活性分子集合 \(\mathcal{A}\) 为锚，从诱饵池中筛选性质匹配负样本。描述符集合取

\[
\mathcal{D}=\{\mathrm{MW},\,\log P,\,\mathrm{TPSA},\,\mathrm{HBD},\,\mathrm{HBA},\,\mathrm{NRot}\},
\]

各描述符给定窗口半宽 \(w_{d}\)（默认：分子量 \(\pm 40\)，\(\log P\pm 1\)，拓扑极性表面积 \(\pm 25\)，氢键供体 \(\pm 1\)，氢键受体 \(\pm 2\)，可旋转键 \(\pm 2\)）。候选诱饵 \(x\) 须对全部 \(d\in\mathcal{D}\) 满足 \(|a_{d}-x_{d}|\le w_{d}\)。性质接近程度用归一化欧氏距离衡量：

\[
d(a,x)=\sqrt{\frac{1}{|\mathcal{D}|}\sum_{d\in\mathcal{D}}\left(\frac{a_{d}-x_{d}}{w_{d}}\right)^{2}}.
\]

同时要求 Morgan 指纹（半径 2）相对活性集的最大 Tanimoto 相似度

\[
\mathrm{TC}_{\max}(x)=\max_{a\in\mathcal{A}}\mathrm{Tanimoto}\!\left(fp(x),fp(a)\right)\le 0.5.
\]

匹配按 round-robin 分配，避免少数活性分子垄断诱饵；若严格窗口不足，再以 \(1.5\,w_{d}\) 放宽一轮补齐。由此得到 TrueDecoy 集后，从同一诱饵池无放回抽取等量分子构成 RandomDecoy。两套基准共用同一批活性分子，仅负样本生成方式不同。

---

## 2.4 受体与配体准备及自对接报告

URAT1 对接采用 9DKB；NLRP3 对接采用 7ALV（共晶配体为 MCC950 类类似物 NP3-146）。受体去除结晶水、加氢并导出对接格式；配体经统一的三维嵌入与准备流程。搜索盒以共晶或参考配体几何中心为原点，边长约 \(22\,\text{Å}\)（以配置文件为准）。

自对接以 lesinurad–9DKB 为主。对重原子坐标，均方根偏差定义为

\[
\mathrm{RMSD}=\sqrt{\frac{1}{n}\sum_{k=1}^{n}\left\|\mathbf{r}_{k}^{\mathrm{dock}}-\mathbf{r}_{k}^{\mathrm{xtal}}\right\|^{2}}.
\]

预先规定报告三类数值：（1）协议原生排序第一构象（top-1）的 RMSD；（2）多构象集合中相对晶体最优构象的 RMSD；（3）RTMScore 所选构象的 RMSD。几何门控阈值取 top-1 RMSD \(\le 2\,\text{Å}\)。若 top-1 未达标而集合内存在近晶体构象，则富集结果仅解释为排序协议表现；结合模式讨论改用集合最优构象、RTMScore 选姿或晶体坐标。

---

## 2.5 候选对接与重打分协议

在同一受体、同一搜索盒与同一配体准备条件下，预先定义下列协议：

| 编号 | 构象来源 | 排序读出 | 角色 |
|------|----------|----------|------|
| P1 | AutoDock Vina | Vina affinity | 物理基线 |
| P2 | gnina | CNNaffinity | 主候选之一 |
| P3 | gnina | gnina affinity（kcal/mol） | 读出对照 |
| P4 | Vina 多构象 | RTMScore | 搜索–打分解耦 |
| P5 | gnina 多构象 | RTMScore | 搜索–打分解耦 |
| P0 | gnina | CNNscore | 负对照（不进入生产排序） |

每个分子保留不少于 9 个结合模式。富集评估与临床库排序使用各协议定义的排序构象分数；需要讨论几何时，另行报告 RTMScore 最高构象或晶体坐标。二者在设计上分离，不以单一 top-1 同时承担“排序最优”与“几何可信”两重主张。生产路径限定为上表开源协议；不将商业对接软件分数与开源分数做绝对值混比。

对接能 \(E_{\mathrm{dock}}\) 以越低越好时，统一转换为“越高越好”的排序分

\[
s=-E_{\mathrm{dock}};
\]

对 CNNaffinity 等本身越高越好的读出，直接取 \(s\) 等于该读出。

---

## 2.6 富集指标与协议选优规则

设标签 \(y\in\{0,1\}\)（1 表示活性），分数 \(s\) 越高越倾向活性。对库规模 \(N\)，取排名前比例 \(f\) 的分子（\(n_{f}=\lfloor Nf\rfloor\)），富集因子为

\[
\mathrm{EF}@f=\frac{\displaystyle\frac{1}{n_{f}}\sum_{i\in\mathrm{Top}\,f}y_{i}}{\displaystyle\frac{1}{N}\sum_{i=1}^{N}y_{i}}.
\]

主文预设报告 \(f=0.01\) 与 \(f=0.05\)，并同时给出 ROC-AUC。

选优规则在对接完成前锁定：（1）以 TrueDecoy 上的 \(\mathrm{EF}@1\%\) 为主判据，并列时依次比较 \(\mathrm{EF}@5\%\) 与 AUC；（2）若某协议在 RandomDecoy 上相对其 TrueDecoy 表现明显变差，则否决该协议；（3）仍并列时，比较 lesinurad、benzbromarone、verinurad、dotinurad 四药在排序中的百分位回收。胜出协议记为 \(\Pi^{\ast}\)，固定用于第二阶段的 URAT1 轴；NLRP3 轴采用同一开源引擎与准备流程，以保持漏斗内相对尺度一致。

---

## 2.7 NLRP3 机器学习缩库与 URAT1 对照模型

构建 assay-conditioned 的 NLRP3 二分类模型。以 Murcko 骨架为分组单位，进行 5 折 GroupKFold，报告受试者工作特征曲线下面积（AUROC）与精确率–召回率曲线下面积（AUPRC）。对临床库输出

\[
P(\mathrm{active})\in[0,1],
\]

预设阈值 \(P(\mathrm{active})\ge 0.5\) 定义对接池。

同步训练 URAT1 回归模型，报告均方根误差、决定系数 \(R^{2}\) 及 Spearman 相关系数，并检查已知尿酸药基准回收。按总体设计，该回归模型不参与临床库主排序，其作用限于检验“URAT1 轴应由对接主导”的数据前提。

---

## 2.8 双靶百分位得分与 Pareto 整合

对接池分子在 9DKB 与 7ALV 上按 \(\Pi^{\ast}\) 对接。对池 \(\mathcal{P}\) 内对接能 \(E\)（越低越好），百分位得分定义为

\[
S=100\times\mathrm{rankpct}(-E),
\]

即能量越低，百分位越高。URAT1 轴取对接百分位

\[
S_{U}=S_{\mathrm{URAT1\text{-}dock}}.
\]

NLRP3 轴取机器学习百分位与对接百分位的较大值

\[
S_{N}=\max\bigl(S_{\mathrm{NLRP3\text{-}ML}},\,S_{\mathrm{NLRP3\text{-}dock}}\bigr),
\]

其中 \(S_{\mathrm{NLRP3\text{-}ML}}=100\times\mathrm{rankpct}\!\bigl(P(\mathrm{active})\bigr)\)。取最大值是为在两类近正交证据并存时避免单轴埋没；二者相关性的数值分析列入 Results。

分子 \(i\) 位于 Pareto 非支配前沿，当且仅当不存在分子 \(j\) 同时满足

\[
S_{U}^{(j)}\ge S_{U}^{(i)}\quad\text{且}\quad S_{N}^{(j)}\ge S_{N}^{(i)},
\]

并且至少在一个坐标上严格更优。lesinurad、colchicine 等对照药进入同一百分位池，用于检查漏斗行为是否与机制预期一致。

---

## 2.9 成药性审计与假说提名

对合并结果标注 PAINS、Brenk 等结构警报，并计算类药性相关描述符。化学适用域采用相对训练集的最大 Morgan Tanimoto 相似度作简化判定：

\[
\mathrm{AD}=\mathbf{1}\!\left[\mathrm{TC}_{\max}\ge\theta_{\mathrm{AD}}\right],
\]

默认 \(\theta_{\mathrm{AD}}=0.35\)，敏感性分析置于补充材料。另进行标签置换（y-scrambling）及短名单稳健性检查。

预先规定：Pareto 前沿不等于最终提名。命中结构警报或成药性明显不佳的分子降级；仅在通过过滤器的子集上，按综合提名分给出跟进假说。提名分的示意形式为

\[
\tau=\alpha S_{U}+\beta S_{N}+\gamma\,\mathbf{1}_{\mathrm{clean}},
\]

其中 \(\mathbf{1}_{\mathrm{clean}}\) 表示通过警报与适用域等过滤器；权重 \(\alpha,\beta,\gamma\) 以最终提名脚本设定为准。

---

## 2.10 分子动力学（若纳入主文）

对 URAT1 基准药、NLRP3 工具药及审计后提名候选，构建相应蛋白–配体复合物并运行分子动力学。分析配体与口袋骨架的 RMSD、关键残基接触；必要时估算结合自由能近似

\[
\Delta G_{\mathrm{bind}}\approx\langle E_{\mathrm{MM}}\rangle+\langle G_{\mathrm{PB/GB}}\rangle+\langle G_{\mathrm{SA}}\rangle-T\langle S\rangle.
\]

初始坐标优先采用晶体构象或 RTMScore 选姿。模拟时长、力场与软件版本在定稿时据实填写。动力学结果用于构象讨论，不作为实验活性证明。

---

## 2.11 软件与可复现性

化学信息学处理使用 RDKit；机器学习使用 scikit-learn；对接使用 AutoDock Vina 与 gnina；构象重打分使用 RTMScore；动力学软件版本写入正文或补充材料。关键随机种子与对接盒参数见 `config/docking_open_source.yaml` 及仓库脚本。凡未特别声明的阈值，以相应脚本的默认设置为准。

---

## 附：投稿时公式取舍建议

| 建议放入主文 Methods | 可放补充材料 |
|----------------------|--------------|
| pActivity；EF@f；\(S_{U}\)/\(S_{N}\)；Pareto 定义；RMSD | TrueDecoy 距离 \(d(a,x)\) 细节；MM-PBSA 展开；提名分权重敏感性 |
| \(\mathrm{TC}_{\max}\le 0.5\)；适用域阈值 | 描述符窗口全表 |
