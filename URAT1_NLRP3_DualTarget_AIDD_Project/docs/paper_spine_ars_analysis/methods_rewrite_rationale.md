# Methods 重写依据（PaperSpine）

> 用途：记录 `METHODS_DRAFT_CN.md` 的写作决策；本文件不是投稿正文。

## 目标期刊范文中学到的组织方式

本轮参考了两类正式论文的 Methods 结构，而不复用其文字。Zengin 等发表于 *Journal of Computer-Aided Molecular Design*（2024；doi:10.1007/s10822-024-00554-4）的计算研究依次设置 benchmark data、docking protocol、rescoring、MD 和 free-energy calculation 等操作性小节，每一节直接说明输入、参数和输出，并未在开头以“阶段 A/阶段 B/阶段 C”重复整篇流程。Caba 等发表于 *Journal of Cheminformatics*（2024；doi:10.1186/s13321-024-00832-1）的虚拟筛选研究采用 data、protein-structure selection、molecular docking、featurization、model construction 和 performance measurement 的顺序；公式只在活性换算和评价指标实际使用处出现。两篇论文共同体现的写法是：Methods 以可复现操作为主线，工作流概述只保留一个短段，方法比较用连续文字说明，结果和作者判断不提前混入。

## 写作思路矩阵

| 单元 | 当前问题 | 与论文贡献的关系 | 范文结构启示 | 仓库证据 | 重写动作 | 终检 |
|---|---|---|---|---|---|---|
| 全节框架 | 作者备忘、开发说明和正文混排 | Methods 应支撑“可复现协议 + 不对称漏斗” | 正文只保留材料、操作、判据 | `confirmed_contribution.md`; `reviewer_audit.md` | 删除作者备忘、表格式协议说明和投稿建议 | PASS：正文可直接进入论文 |
| 2.1 研究设计 | 三个“阶段”被拆成四段，缺少因果连接 | 解释为什么先独立选协议、再应用于临床库 | 一段式 workflow overview | V2 outline；临床库不参与协议选优 | 改为“数据不对称 → 独立 benchmark → 固定协议 → 临床应用 → 审计”的连续段落 | PASS：无 A/B/C 标签 |
| 2.2 数据 | 数据来源和清洗过简 | 证明模型与基准输入可复现 | Data 小节写清筛选、聚合和标签 | `utils_ml.py`; `data_summary.json`; `build_repurposing_library.py` | 补充关系符、冲突删除、活性阈值、去盐和去重 | PASS：数字与脚本一致 |
| 2.3 诱饵基准 | 公式和操作分离；窗口曾写错 | 协议论文的主要 benchmark | benchmark data 单独成节 | `build_urat1_true_decoy.py` | 在匹配步骤中嵌入距离与 Tanimoto 公式 | PASS：窗口、seed、ratio 与脚本一致 |
| 2.4 机器学习 | 仅写“assay-conditioned”，不可复现 | 支撑 NLRP3 缩库与 URAT1 排除 | 模型特征、划分、聚合和指标依次写 | `02_train_asymmetric_models.py`; `screen_repurposing_library.py` | 补充 ECFP4、12 描述符、assay one-hot、XGBoost、骨架 CV | PASS：不提前报告性能结果 |
| 2.5 结构准备 | 两个搜索盒被笼统写成 22 Å | 结构计算的复现底线 | 受体、配体和网格按操作顺序写 | `docking_open_source.yaml`; prep scripts | 写入 9DKB/7ALV 精确中心和 22/20 Å 尺寸 | PASS：与配置一致 |
| 2.6 协议比较 | 文章正文出现协议表；P0 角色表述生硬 | 支撑 protocol-selection 主贡献 | docking/rescoring 用段落叙述 | Vina/gnina scripts；V2 plan | 将 P0–P5 写入连续段落；定义排序方向和 EF | PASS：无表格 |
| 2.7 Pareto | `rankpct` 公式过于口语化 | 支撑不对称双轴整合 | 先定义单轴，再定义非支配关系 | `merge_docking_pareto.py` | 用经验百分位公式定义三个分量和 Pareto | PASS：与代码方向一致 |
| 2.8 审计提名 | 原稿虚构加权提名分；AD 阈值写成 0.35 | “Pareto ≠ 提名”是核心贡献 | 审计规则须与实际实现一致 | scripts 09–14 | 删除虚构加权式；改写为双 90 百分位门控、清洁过滤和字典序排名；AD 改为训练集 5% 分位定义 | PASS：规则与脚本一致 |
| MD | 尚无终版运行参数 | PaperSpine 禁止制造证据 | 缺证据单元不写成已完成 Methods | V2 status | 暂不在正文伪写力场和时长；完成后按实际日志增加小节 | PASS：无虚构参数 |

## 需要在服务器结果返回后补齐的唯一技术项

gnina 最终 `num_modes`、RTMScore 输入构象数及各软件的实际版本应从服务器运行日志回填。当前正文不编造这些数值。分子动力学小节须在模拟实际完成后再依据真实体系、力场、平衡流程、生产时长和轨迹采样设置撰写。
