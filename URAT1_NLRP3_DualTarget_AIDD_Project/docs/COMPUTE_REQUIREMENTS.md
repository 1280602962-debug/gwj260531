# 算力需求评估

---

## 一、总览

| 阶段 | GPU | CPU | 内存 | 存储 | 墙钟时间（单卡 4090） |
|------|-----|-----|------|------|---------------------|
| 数据准备 | 否 | 4 核 | 8 GB | 5 GB | 数小时 |
| 数据集分析 | 否 | 4 核 | 16 GB | 2 GB | 1–2 小时 |
| MiniMol 指纹提取 | 可选 | 8 核 | 32 GB | 10 GB | 2–4 小时 |
| MTL 训练 | **是** | 8 核 | 32 GB | 5 GB | 4–8 小时 |
| Chemprop baseline | **是** | 8 核 | 32 GB | 5 GB | 2–6 小时 |
| Vina 系综对接 (10^6) | 否 | **32 核** | 64 GB | 20 GB | **2–5 天** |
| Glide SP (10^4) | 否 | 16 核 | 64 GB | 10 GB | 1–2 天 |
| MD 100 ns × 20 复合物 | 可选 GPU | 16 核 | 32 GB | 50 GB | 3–7 天 |
| CLM RL 生成 (5000 step) | **是** | 8 核 | 32 GB | 10 GB | 1–2 天 |
| **总计（完整 pipeline）** | 1–2 卡 | 32 核推荐 | 64 GB | **~100 GB** | **~2–3 周** |

---

## 二、分项详解

### 2.1 分子基础模型 / MTL 训练

| 配置 | 最低 | 推荐 | 理想 |
|------|------|------|------|
| GPU | RTX 3060 12GB | RTX 3090/4090 24GB | A100 40GB |
| 批大小 | 32 | 64–128 | 256 |
| 训练时间 | ~8 h | ~4 h | ~1 h |

MiniMol 指纹提取可 CPU 运行（~6M 分子/hour/8 核），GPU 加速约 3–5×。

**内存**：全库 10^6 分子指纹约 10^6 × 768 × 4 bytes ≈ **3 GB**（可批处理）。

### 2.2 虚拟对接（算力瓶颈）

**AutoDock Vina**（CPU 并行）：

```
10^6 分子 × 3 URAT1 构象 × 1 NLRP3 构象 ≈ 4×10^6 对接任务
单任务 ~30 s（含准备）→ 32 核并行约 40–80 小时
```

**优化策略**：
1. ML 预筛到 10^4 → 对接量降 100×
2. 仅对 Top 500 做 MD
3. 使用 GNINA CNN 重打分替代部分 Glide

**Glide SP**（需 Schrödinger 许可）：
- 10^4 分子 × 2 靶点 ≈ 2×10^4 任务
- 16 核：~24–48 小时

### 2.3 分子动力学

| 体系 | 原子数 | 100 ns 墙钟（GPU） | 存储 |
|------|--------|-------------------|------|
| URAT1-配体-膜 | ~80,000 | 2–4 天 | ~5 GB/体系 |
| NLRP3 NACHT-配体 | ~50,000 | 1–2 天 | ~3 GB/体系 |

**建议**：仅对 **Top 20–30** 候选做 MD，不做全库。

**软件**：GROMACS + CHARMM36m 力场；配体参数化用 CGenFF / GAFF2。

### 2.4 生成式 RL

| 组件 | 算力 |
|------|------|
| CLM 预训练权重加载 | 4–8 GB VRAM |
| RL fine-tune (5000 steps) | 8–16 GB VRAM, ~24 h |
| 每 100 步对接评估 | 额外 CPU 时间（可降频到每 500 步）|

**技巧**：结构奖励 **稀疏计算**（每 500 RL 步对接一次），可节省 80% 对接算力。

---

## 三、三档配置方案

### 方案 A：最低可行（预算有限）

- **硬件**：1× RTX 3060 12GB，16 GB RAM，8 核 CPU，100 GB SSD
- **策略**：
  - 跳过 MD，仅 Vina + MM-GBSA（静态）
  - 筛选库降至 10^5（Enamine 子集）
  - 跳过生成式模块
- **可发表**：JCIM / Pharmaceutics 级别
- **估计成本**：云 GPU ~$50–100（RunPod/Vast.ai）

### 方案 B：推荐配置（稳健发表）

- **硬件**：1× RTX 4090 24GB，64 GB RAM，32 核 CPU，200 GB SSD
- **策略**：
  - 完整 MTL + 系综对接
  - Top 20 MD 验证
  - CLM RL 生成
- **可发表**：J Cheminformatics / Briefings in Bioinformatics
- **估计成本**：云 GPU ~$200–400；或学校工作站

### 方案 C：理想配置（冲高）

- **硬件**：2× A100 40GB，128 GB RAM，64 核，500 GB NVMe
- **策略**：
  - 10^6 全库对接 + Glide XP 精修
  - 50 复合物 MD + FEP（可选）
  - 多 seed 生成 + 大规模消融
- **可发表**：冲击更高档期刊
- **估计成本**：云 ~$1000+

---

## 四、云算力平台参考（2026）

| 平台 | GPU 型号 | 价格约 |
|------|---------|--------|
| RunPod | RTX 4090 | $0.4–0.7/h |
| Vast.ai | RTX 3090 | $0.2–0.4/h |
| Google Colab Pro+ | A100 限量 | $50/月 |
| 阿里云/AutoDL | A100 | ¥3–8/h |

**对接任务**建议租 **高 CPU 实例**（32+ 核）而非 GPU。

---

## 五、存储规划

```
data/
  raw/           ~2 GB   (ChEMBL CSV, PDB)
  processed/     ~1 GB   (清洗数据, 指纹)
  structures/    ~5 GB   (预处理 PDB, 膜体系)
  libraries/     ~5 GB   (10^6 SMILES)
  docking/       ~30 GB  (对接 pose, log)
  md/            ~50 GB  (轨迹, 仅 top 候选)
results/
  models/        ~2 GB
  figures/       ~500 MB
  screening/     ~5 GB
```

**总计**：建议预留 **100–150 GB**。

---

## 六、并行化建议

```bash
# 对接：GNU parallel 或 SLURM array job
parallel -j 32 "vina --receptor {} --ligand {} " ::: receptors ::: ligands

# MTL：多 seed 并行
for seed in 0 1 2 3 4; do
  python train_mtl.py --seed $seed &
done

# MD：按候选分配不同 GPU
CUDA_VISIBLE_DEVICES=0 gmx mdrun -deffnm md_1 &
CUDA_VISIBLE_DEVICES=1 gmx mdrun -deffnm md_2 &
```

---

## 七、时间与算力权衡决策树

```
是否有 GPU (≥12GB)?
├─ 否 → 方案 A：XGBoost + Vina CPU，小库
└─ 是
    └─ 是否有 32 核 CPU / 3 周时间?
        ├─ 否 → ML 筛到 10^3 + 单构象对接
        └─ 是 → 完整 STAD-AIDD（方案 B）
            └─ 是否做生成式?
                ├─ 否 → 仍可发 JCIM（方法偏筛选）
                └─ 是 → 方案 B/C，冲 AI 方法期刊
```
