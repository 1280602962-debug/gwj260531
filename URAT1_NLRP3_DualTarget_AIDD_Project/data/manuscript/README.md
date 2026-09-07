# 论文证据层

只从该层的表和[事实锁](../../MANUSCRIPT_FACTS.yaml)取主文数字。表由[生成器](../../scripts/manuscript_pipeline/build_evidence.py)从冻结历史输入产生；[provenance.json](provenance.json)包含输入与输出SHA-256。生成器不会改写原始数据或重开候选排序。

| 证据包 | 主文数据 | SI/溯源 |
|---|---|---|
| 1. 数据与NLRP3模型 | [模型验证](tables/model_validation.json)、[漏斗](tables/screening_funnel.csv) | [39 assay来源](si/nlrp3_assay_provenance.csv)、[阈值敏感性](si/nlrp3_threshold_sensitivity.json)、[模型登记](../frozen/model_manifest.json) |
| 2. 酸性化学空间 | [描述符](tables/acid_descriptors.csv)、[总体分布](tables/chemical_space_summary.csv) | [微观状态](si/ligand_microstates.csv)、原303/156成员清单哈希见provenance |
| 3. URAT1结构方法 | [重对接](tables/urat1_redocking_summary.csv)、[A1/A2/W2验证](tables/urat1_gate_validation.csv) | [受体来源](tables/receptor_manifest.csv)、[残基映射](tables/residue_mapping.csv)、[规则](../../config/urat1_gate_definition.yaml) |
| 4. NLRP3结构方法 | [自对接](tables/nlrp3_selfdocking_summary.csv)、[分seed门验证](tables/nlrp3_gate_validation.csv) | [阳性骨架](si/w4_positive_scaffolds.csv)、[规则](../../config/nlrp3_gate_definition.yaml) |
| 5. 双臂交集与药化 | [156行决策台账](tables/screening_decision_ledger.csv)、[分seed候选指标](tables/candidate_seed_metrics.csv) | [primary](../frozen/primary_candidates_frozen.csv)、[reserve](../frozen/backup_candidates_frozen.csv)、[药化反事实](si/structure_only_chemistry_counterfactual.csv) |
| 6. 后续验证 | [PF选择理由](../../docs/NOMINATION_RATIONALE.md) | [MD计划](../../config/md_protocol.yaml)；没有MD结果 |

主图从tables生成至figures。敏感性和新颖性表位于si；历史P2/URAT1 ML证据见[data/archive](../archive/README.md)，不接入当前候选排序。关键限制见[核对报告](../../docs/CONSOLIDATION_AUDIT.md)。

统计口径：临床筛选以唯一repurposing_id计数；模型分子数以canonical_smiles计数；酸性官能团类别可重叠；Murcko计数包含空骨架类别。best-of-9 RMSD为诊断下界，不是可实际选择的Top-1。A1基于seed42，A2/NLRP3稳定性和W2注释分开计数。
