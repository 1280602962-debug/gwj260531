# 历史数据分层登记

为保持脚本、原始SDF与审计表的路径引用，历史数据在原址逻辑归档，不批量移动或删除。

| 原目录 | 当前角色 |
|---|---|
| data/raw | 不可改写的原始输入 |
| data/processed | 标准化训练记录；论文provenance锁定输入哈希 |
| data/repurposing/p2 | P2 negative-transfer/回顾SI；51/7仅legacy audit set |
| data/campaigns/c1 | 已执行对接、A1/A2、自对接及原始药化审计；当前生成器的明确来源 |
| data/campaigns/c5 | 原始C5门验证和冻结成员；审计源，不是第二套当前候选入口 |
| data/si | 历史SI；只有被论文层列出的结果进入当前稿 |
| docking_export_20260820 | 原始benchmark姿势和协议导出 |

当前候选只由[data/frozen](../frozen/README.md)提供；正文/SI表由[data/manuscript](../manuscript/README.md)提供。原始数据生命周期与基点清单见data/manuscript/audit/base_file_inventory.csv。
