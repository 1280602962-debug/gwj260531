# ENV_PIN — 可复现环境快照

> 采集：2026-07-29（本机 WSL2）

| 组件 | 版本 / 路径 |
|------|-------------|
| Python | `/home/gwj/miniconda3/bin/python3` (3.x) |
| RDKit | 2026.3.1 (`pip show rdkit`) |
| meeko | 0.7.1 |
| AutoDock Vina | v1.2.7 (`/home/gwj/miniconda3/bin/vina`) |
| GNINA | v1.3.2 master:f23dd2b (`/mnt/d/CADD paper exercise/gnina/bin/gnina`, CPU `--no_gpu`) |
| RTMScore | `/home/gwj/miniconda3/envs/rtmscore/bin/python` + `/home/gwj/software/RTMScore/trained_models/rtmscore_model1.pth` |
| Open Babel | 用于 GNINA 前 pdbqt→sdf（系统 obabel） |

## 当前零对接复分析环境（2026-08-26）

原始冻结对接仍由上表环境生成，未在本次修订中重跑。用于重算统计表、置信区间、化学描述符与姿态审计的 Python 环境固定在仓库根目录 `requirements-analysis.txt`：Python 3.12.13、NumPy 2.5.2、pandas 3.0.5、SciPy 1.18.1、scikit-learn 1.9.0、RDKit 2026.3.5、Meeko 0.7.1 与 gemmi 0.7.5。bootstrap 子种子使用 SHA-256 稳定偏移，不再依赖跨进程变化的 Python `hash()`。

## 协议常数

- ligand prep: RDKit ETKDGv3, seed **20260727**
- panel sampling seed: **20260729**
- Vina: n_modes=9, energy_range=3
- PM 主面板: E=16; exhaustiveness 对照: E=8
- GNINA: `--cnn_scoring rescore --minimize --seed 20260727`

## 激活命令

```bash
source "/mnt/d/CADD paper exercise/gnina/activate.sh"
export VINA=/home/gwj/miniconda3/bin/vina
export PYTHON=/home/gwj/miniconda3/bin/python
export RTM_PY=/home/gwj/miniconda3/envs/rtmscore/bin/python
```
