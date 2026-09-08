# URAT1–NLRP3临床库双节点候选：论文复现入口

研究问题：临床阶段及已上市小分子中，哪些酸性分子具有值得实验验证的URAT1与NLRP3口袋结构相容性？NLRP3-related activity classifier用于缩库；URAT1采用几何门，候选优先级不代表预测活性强弱。

现行路线：**8,319 → 1,588 → 303 → 156**，随后多种子双靶结构交集、历史药化审计与最终排除得到**12个primary**，另存**21个reserve/admissible分子**。PF-03882845是后续动态验证优先对象；MD尚未授权且没有轨迹结果。

| 入口 | 用途 |
|---|---|
| [MANUSCRIPT_MASTER](docs/MANUSCRIPT_MASTER.md) | 唯一当前稿件入口；整合Methods、Results和Discussion |
| [事实锁](MANUSCRIPT_FACTS.yaml) | 核心数字、定义、状态；由证据生成器产生 |
| [现行科学配置](config/campaign_final.yaml) | 当前解释规则和专题配置入口 |
| [冻结候选](data/frozen/README.md) | 唯一当前primary/reserve视图及模型哈希 |
| [论文数据](data/manuscript/README.md) | 六个证据包、主文表、SI与溯源 |
| [整理核对报告](docs/CONSOLIDATION_AUDIT.md) | 条目1–4已闭合；MD未跑 |
| [后续实验](docs/FOLLOWUP_EXPERIMENTS.md) | MD参数与湿实验清单；无轨迹数字 |
| [URAT1能力对接](docs/URAT1_CAPABILITY_DOCK.md) | 本机D1/D2；20个新作业；不改12/21 |
| [复现工作流](scripts/manuscript_pipeline/README.md) | 重建、核对、生成图表 |
| [开发归档](docs/archive/README.md) | 原文快照及兼容重定向；非当前主线 |

从项目根目录执行：

```bash
python scripts/manuscript_pipeline/build_evidence.py
python scripts/manuscript_pipeline/verify_evidence.py
python scripts/manuscript_pipeline/make_figures.py
```

依赖：Python、pandas、numpy、scipy、PyYAML、RDKit；制图另需matplotlib。实际整理版本见data/manuscript/audit/environment.json。此入口读取已完成结果，不重训模型、不重新对接、不覆盖原始数据。

当前证据可用于写作，但并非投稿就绪：历史药化排除含药名规则；酸性等价池比COOH姿势匹配器适用范围更宽；模型与历史生产池缺少同时生成的哈希证明；W4有失败诱饵和化学家族偏倚。详见审计报告。

原始data/raw、data/campaigns、P2及对接导出保留原址，生成器通过明确的来源清单读取；旧脚本和配置按登记表退出现行执行入口，保留导入兼容性。JNK1独立项目不在本次修改范围内。
