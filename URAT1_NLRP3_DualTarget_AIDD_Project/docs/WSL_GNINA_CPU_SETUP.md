# WSL2 + 无 NVIDIA 显卡：GNINA 对接配置指南

适用于 **Windows + WSL2 Ubuntu**，**没有独显** 时用 GNINA **CPU 模式**（`--no_gpu`）跑对接。

> ⚠️ CPU 很慢：1588×2 全库可能需要 **1–3 周**。建议先 `--limit 5` 烟雾测试，再挂机全库。

---

## 一、复制给 Agent 的完整命令（直接粘贴）

### A. 无 GitHub 连接（本机 WSL 常用）

项目代码需 **事先已在 WSL 本地**（U 盘拷贝、`/mnt/c/` 复制、或从能联网的机器打包传过来）。**不要** `git clone` / `git pull`。

把下面整段发给 **本机 WSL 终端里的 Cursor Agent**（只改「项目路径」）：

```text
你在本机 WSL2 Ubuntu 终端中工作。无法访问 GitHub，不要使用任何 git 命令。请仅用本地已有文件配置 GNINA CPU 对接（无 NVIDIA 显卡，必须 --no_gpu），并跑通烟雾测试。不要跑全库 1588×2。

【项目根目录】（改成你的实际路径）
/home/你的用户名/gwj260531/URAT1_NLRP3_DualTarget_AIDD_Project

【先检查本地是否有所需脚本，缺则停止并列出缺失文件】
必须存在：
  scripts/setup_gnina_wsl_cpu.sh
  scripts/run_gnina_batch.py
  scripts/prepare_receptor_vina.py
  scripts/prepare_ligands_vina.py
  config/docking_open_source.yaml
  data/structures/pdb/9DKB.cif 或 data/structures/pdb/7ALV.pdb
若不存在，告诉我缺什么，不要假装已从 GitHub 拉取。

【执行步骤】
1. cd 到项目根目录；pwd 确认路径
2. 确认在 /home/... 下（尽量不在 /mnt/c/ 跑大批量对接）
3. bash scripts/setup_gnina_wsl_cpu.sh
   - sudo apt 需要密码时提示我输入
   - GNINA 二进制下载地址（仅需访问 github.com/releases，若也失败则改用手动）：
     https://github.com/gnina/gnina/releases/download/v1.3.1/gnina-1.3.1-x86_64-ubuntu22.04
     保存为 tools/gnina && chmod +x tools/gnina
4. 验证：
   tools/gnina --help
   ls data/structures/prepared/9DKB_receptor.pdbqt
   cat results/repurposing/smoke_gnina/dock/docking_9dkb_gnina.csv
   （lesinurad 应为 docking_status=docked，dock_score 非空）
5. 仅测 10 个分子估时（不全库）：
   python3 scripts/prepare_ligands_vina.py \
     --input data/repurposing/screening/docking_pool_p05.csv \
     --output-dir results/repurposing/ligands_p05
   python3 scripts/run_gnina_batch.py --target urat1_9dkb \
     --manifest results/repurposing/ligands_p05/ligand_manifest.csv \
     --output-dir results/repurposing/docking_gnina/9dkb --jobs 2 --limit 10

【若 pip / apt 可用但 GitHub 不可用】
- Python 依赖：pip3 install -r requirements.txt（PyPI 一般仍可访问）
- 仅 GNINA 二进制需从 releases 下载；若 releases 也打不开，请我写出手动下载步骤（Windows 浏览器下载后 cp 到 WSL）

【成功标准】
烟雾测试通过 + 10 分子测速；最后只文字说明全库命令（不执行）：
  JOBS=2 bash scripts/run_gnina_docking_pipeline.sh
```

### B. 有 GitHub 时（可选）

```text
请在 WSL2 Ubuntu 里为项目 URAT1_NLRP3_DualTarget_AIDD_Project 配置 GNINA CPU 对接环境并跑通烟雾测试。

环境约束：
- Windows + WSL2，无 NVIDIA 显卡，必须用 gnina --no_gpu
- 项目路径：~/gwj260531/URAT1_NLRP3_DualTarget_AIDD_Project（若不同请先 cd 到实际路径）
- 分支：cursor/urat1-nlrp3-dualtarget-aidd-e43d

请按顺序执行并汇报每步结果：

1. 进入项目根目录，git pull 最新代码
2. 运行 bash scripts/setup_gnina_wsl_cpu.sh
   - 若 GNINA 二进制下载失败，从 https://github.com/gnina/gnina/releases 手动下载 gnina-1.3.1-x86_64-ubuntu22.04 到 tools/gnina 并 chmod +x
3. 确认 tools/gnina --help 可用；烟雾测试 CSV 中 lesinurad dock_score 非空
4. 不要立刻跑全库 1588×2；仅当烟雾测试通过后，告诉我如何用：
   bash scripts/run_gnina_docking_pipeline.sh
   以及如何用 JOBS=2 控制并行

5. 若 setup 失败，请诊断：python3 版本、pip install -r requirements.txt、受体 PDBQT 是否生成

成功标准：
- results/repurposing/smoke_gnina/dock/docking_9dkb_gnina.csv 有 1 行 docked
- data/structures/prepared/9DKB_receptor.pdbqt 存在

Methods 将写：AutoDock Vina sampling + GNINA CNN rescoring, CPU mode (--no_gpu), exhaustiveness=16。
```

### C. 代码还没进 WSL：先拷项目（不用 Git）

在 **Windows 浏览器** 打开 https://github.com/1280602962-debug/gwj260531 ，Download ZIP，解压后在 WSL 执行：

```bash
mkdir -p ~/gwj260531
cp -r /mnt/c/Users/你的Windows用户名/Downloads/gwj260531-main/URAT1_NLRP3_DualTarget_AIDD_Project ~/gwj260531/
cd ~/gwj260531/URAT1_NLRP3_DualTarget_AIDD_Project
```

或用 Cursor Cloud Agent 已改好的仓库：在能联网环境打包 `URAT1_NLRP3_DualTarget_AIDD_Project` 文件夹，拷进 WSL 再跑上面 **A** 段 Agent 指令。

---

## 二、你自己在 WSL 终端里手动跑（不通过 Agent）

```bash
# 1. 进入项目（按你的实际路径改）
cd ~/gwj260531/URAT1_NLRP3_DualTarget_AIDD_Project

# 2. 一键安装 + 烟雾测试
bash scripts/setup_gnina_wsl_cpu.sh

# 3. 全库（很慢，建议 screen/tmux 挂机）
# tmux new -s gnina
JOBS=2 bash scripts/run_gnina_docking_pipeline.sh
```

只跑前 10 个分子试时间：

```bash
python3 scripts/prepare_ligands_vina.py \
  --input data/repurposing/screening/docking_pool_p05.csv \
  --output-dir results/repurposing/ligands_p05
python3 scripts/run_gnina_batch.py --target urat1_9dkb \
  --manifest results/repurposing/ligands_p05/ligand_manifest.csv \
  --output-dir results/repurposing/docking_gnina/9dkb --jobs 2 --limit 10
```

---

## 三、配置说明

| 项 | 值 |
|----|-----|
| GNINA 二进制 | `tools/gnina`（setup 脚本下载） |
| CPU 模式 | `config/docking_open_source.yaml` → `gnina.no_gpu: true` |
| exhaustiveness | 16（CPU 建议；有 GPU 可改 32） |
| 并行 | `JOBS=2`（无 GPU 不宜超过 CPU 物理核数一半） |
| 打分 | `dock_score` = GNINA affinity (kcal/mol)，越低越好 |

---

## 四、常见问题

**Q: WSL 里没有 nvidia-smi？**  
A: 正常。使用 `--no_gpu`，不要装 CUDA。

**Q: gnina: error while loading shared libraries: libcudart**  
A: 换用带 CPU 后端的 release，或 Docker：`docker run ... gnina --no_gpu`（见 [GNINA releases](https://github.com/gnina/gnina/releases)）。

**Q: 太慢怎么办？**  
A: 租带 GPU 的云机只跑对接；或先用纯 Vina（`run_vina_batch.py`）缩时间，Top 10% 再用 GNINA rescore。

**Q: 和 Glide XP 数据能混用吗？**  
A: **不能。** Methods 写 GNINA 就必须用 GNINA 重跑的全库分数。

---

## 五、与论文章节对应

Methods 建议表述：

> Molecular docking used GNINA 1.3.1 in CPU mode (`--no_gpu`) with exhaustiveness 16 and CNN rescoring. Receptors were prepared with gemmi and Open Babel; ligands with RDKit and Meeko. Scores were converted to within-library percentiles for Pareto analysis.
