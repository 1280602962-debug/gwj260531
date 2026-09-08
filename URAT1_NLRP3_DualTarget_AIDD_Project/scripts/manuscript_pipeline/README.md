# 论文复现工作流

默认工作是从已存输入重建可投稿写作用的证据表，不重新启动模型训练、对接或MD。

```bash
python scripts/manuscript_pipeline/build_evidence.py
python scripts/manuscript_pipeline/verify_geometry.py
python scripts/manuscript_pipeline/make_figures.py
python scripts/manuscript_pipeline/verify_evidence.py
```

基础依赖：pandas、numpy、scipy、PyYAML、RDKit；几何复核另需gemmi，制图需matplotlib。使用独立环境安装依赖。实际核对环境写入data/manuscript/audit/environment.json；上游计算环境不同，重算注释不替代历史值。

build_evidence复核唯一ID、集合嵌套、609/513/39、40个历史eligible、12/21候选和A1/A2/W2/W4统计，输出事实锁、来源哈希、assay表、化学空间、受体来源/编号、冻结候选丰富视图及SI。--inventory额外对审计基点的所有Git blob生成SHA-256清单；只在需要更新基点清单时运行。

verify_geometry重新读取936份临床SDF，复核A1/A2/NLRP3门及33个冻结分子的W2，输出独立审计文件；这验证已存算法重现，不证明活性预测有效。

make_figures只读取data/manuscript/tables，输出静态论文图与来源哈希。verify_evidence为只读校验，检查来源/输出哈希、名单、归档字节和当前Markdown链接。任意不一致以非零退出。

## 从原始数据重新计算的依赖图

| 步骤 | 保留实现 | 重现边界 |
|---|---|---|
| ChEMBL清洗 | scripts/00_prepare_data.py | 当前生成器直接检查已存processed；重清洗另建输出 |
| NLRP3训练/筛库 | scripts/02_train_asymmetric_models.py、screen_repurposing_library.py | 重训是敏感性；不得覆盖1588池 |
| 酸池/配体准备 | build_c1_acid_clinical_pool.py、prepare_ligands_c1.py | 从冻结池出发；微观状态清单不可混用 |
| 双靶对接 | run_gnina_batch.py、run_c1_acid_dual_dock.py、run_c1_acid_dual_a2.py | 旧路径作为实现依赖；最终参数见config/docking_final.yaml |
| 门计算 | c1_acid_pose_selection.py、c1_nlrp3_pose_metrics.py、run_c5_w2_urat1_ifp_gate.py | 已存定义；A1/A2/W2作用分离 |
| 历史药化/冻结 | build_c5_tier_assignment.py、freeze_c5_shortlist.py | 新生成器独立复核历史eligible来源与集合 |
| 论文证据/图表 | 本目录脚本 | 唯一当前默认入口 |

retire_legacy.py仅用于本次一次性机械归档，不属于常规重建命令。不要执行旧手工提名脚本来更新当前候选。
