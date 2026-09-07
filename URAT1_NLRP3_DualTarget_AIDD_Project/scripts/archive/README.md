# 脚本退场登记

[script_registry.csv](script_registry.csv)登记原脚本的角色。为保持导入、PROJECT_ROOT和历史复现路径，原脚本保留原位。默认入口仅为[scripts/manuscript_pipeline](../manuscript_pipeline/README.md)。

build_c1_acid_shortlist_a2.py的硬编码primary/backup名单已退役；其40个eligible的历史生成逻辑只作为可追溯事实。原build_c5_tier_assignment.py仍依赖该40个集合，不能宣称完全没有历史药名规则。新生成器会独立重建并检查，而不会运行旧提名脚本。

编号开发脚本、P2流水线、placeholder或早期路线生成器均不是当前论文的快速开始命令。需要重新执行历史研究时，应使用登记的基点快照及原环境，输出到新的审计目录。
