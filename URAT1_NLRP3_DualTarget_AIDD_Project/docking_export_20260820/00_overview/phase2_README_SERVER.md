# 服务器 GPU 对接包：URAT1–NLRP3 双靶 P2

上传本文件夹到服务器后即可跑。协议锁定为 **Π\* = P2**（gnina，CNNaffinity，`cnn_scoring=rescore`，`exhaustiveness=32`）。

## 内容

| 路径 | 说明 |
|------|------|
| `data/structures/prepared/9DKB_receptor.pdbqt` | URAT1 受体（已准备） |
| `data/structures/prepared/7ALV_receptor.pdbqt` | NLRP3 受体（已准备） |
| `results/repurposing/ligands_p05/*.pdbqt` | 1583 个临床池配体（已 PDBQT） |
| `results/repurposing/ligands_p05/ligand_manifest.csv` | 配体清单 |
| `results/repurposing/docking_p2/9dkb/poses/` | 本机已完成的 9DKB 姿态（续跑会跳过） |
| `config/docking_production_p2_gpu.yaml` | GPU 版 P2 配置（`no_gpu: false`） |
| `run_server_gpu.sh` | 一键：9DKB → 7ALV → Pareto |

查看打包时进度：`PROGRESS.json`。

## 服务器依赖

- `gnina`（GPU 版，CUDA 匹配）
- `python3` + `pandas` + `pyyaml`
- NVIDIA 驱动 / `nvidia-smi` 可用

## 运行

```bash
cd server_p2_gpu_upload   # 或你解压后的目录名

# 单卡建议 JOBS=1 或 2（多卡再加大）
JOBS=2 bash run_server_gpu.sh
```

后台跑示例：

```bash
nohup env JOBS=2 bash run_server_gpu.sh > results/repurposing/logs/server_gpu.log 2>&1 &
tail -f results/repurposing/logs/server_gpu.log
```

## 完成后带回本地的文件

至少带回：

```
results/repurposing/docking_p2/9dkb/docking_9dkb_gnina.csv
results/repurposing/docking_p2/7alv/docking_7alv_gnina.csv
results/repurposing/docking_p2/*/poses/          # 姿态（可选但建议）
data/repurposing/pareto/pareto_merged_scores.csv
data/repurposing/pareto/pareto_shortlist.csv
```

## 注意

- **不要**降低 `exhaustiveness` 或关掉 `cnn_scoring=rescore`，否则与本地锁定的 P2 不一致。
- 中断后再次执行 `run_server_gpu.sh` 会跳过已有 pose，可安全续跑。
- 本机若还在 CPU 对接，上传后请以服务器结果为准，或停掉本机任务避免混写。
