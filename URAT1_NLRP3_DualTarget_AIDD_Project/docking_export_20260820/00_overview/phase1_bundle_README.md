# 对接结果汇总包

路径：`/home/hww/gwj/NLRP3_URAT1/docking_results_bundle`

受体：9DKB（URAT1） | 库规模：9839（active 469 / true_decoy 4687 / random_decoy 4683）

## 目录说明

| 目录/文件 | 内容 |
|-----------|------|
| `scores/mol_protocol_scores.csv` | **主表**：每分子各协议 top 分（Vina / gnina CNN / RTM） |
| `scores/rtm_*_pose_scores.csv` | RTM 逐 pose 原始分 |
| `metrics/protocol_metrics.csv` | 各协议 EF@1%/EF@5%/AUC（True & Random） |
| `benchmarks/` | True/Random 基准与对接池 |
| `meta/mol_index_map.csv` | mol_id ↔ SMILES ↔ role |
| `poses/vina_pdbqt` | Vina 构象（**符号链接**到原目录） |
| `poses/gnina_sdf` | gnina 构象（**符号链接**） |
| `receptor/` | 受体 PDB/PDBQT 与 RTM 口袋 |
| `config/` | 对接配置与盒子 |
| `logs/` | 进度、超时列表；per-mol log 为链接 |
| `SUMMARY.txt` | 运行覆盖率与指标速览 |

## 协议列含义（`mol_protocol_scores.csv`）

| 列 | 含义 | 方向 |
|----|------|------|
| P1_vina_affinity | Vina 最优亲和力 | 越低越好 |
| P2_CNNaffinity | gnina CNN 亲和力 | 越高越好 |
| P3_gnina_affinity | gnina minimizedAffinity | 越低越好 |
| P0_CNNscore | gnina CNNscore | 越高越好 |
| P4_RTM_vina | RTM 打在 Vina pose 上 | 越高越好（覆盖不全） |
| P5_RTM_gnina | RTM 打在 gnina pose 上 | 越高越好 |

## 注意

- `poses/` 下是符号链接，删除本汇总包**不会**删原始构象；但若移动/删除原始 `work/vina` 或 `work/gnina`，链接会失效。
- 若需要可携带的完整拷贝（非链接），再说一声即可另做一份。
