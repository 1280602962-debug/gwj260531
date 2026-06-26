# JNK 选择性筛选 — Windows 第 1/2 步工作流

在 **Windows PowerShell** 中运行（需安装 Schrödinger 2025）。

## 文件清单

| 文件 | 作用 |
|------|------|
| `jobs_step1.json` | 导出 + Prep + IFP 配置（排除 VSW，8 个非共价配体） |
| `export_complexes_batch.py` | 从 `benchmarks_*_pv.maegz` 导出复合物 |
| `prep_complexes_batch.py` | Protein Prep（固定 ligand，加氢） |
| `calc_ifp_batch.py` | SIFt IFP 批量计算 |
| `run_step1.ps1` / `run_step1.bat` | 一键执行第 1 步 |
| `jobs_mmgbsa.json` | MM-GBSA + ΔΔG 配置 |
| `run_mmgbsa_batch.py` | Prime MM-GBSA 批量提交（可选，已有结果时跳过） |
| `scan_mmgbsa_inventory.py` | 扫描已有 MM-GBSA 结果清单 |
| `calc_ddg_selectivity.py` | 解析已有 MM-GBSA，输出 ΔΔG(JNK1/2/3) |
| `run_step2.ps1` / `run_step2_parse_only.ps1` | 第 2 步（支持 `-ParseOnly` 仅解析） |

## 一键运行

```powershell
cd "D:\CADD paper exercise\JNK1_2_3\Docking\jnk_docking_export"
Set-ExecutionPolicy -Scope Process Bypass

# 第 1 步：导出 → Prep → SIFt IFP
.\run_step1.ps1 -DockingDir "D:\CADD paper exercise\JNK1_2_3\Docking"

# 第 2 步：解析已有 MM-GBSA（不重算）→ ΔΔG
.\run_step2_parse_only.ps1 -DockingDir "D:\CADD paper exercise\JNK1_2_3\Docking"

# 或
.\run_step2.ps1 -ParseOnly -DockingDir "D:\CADD paper exercise\JNK1_2_3\Docking"
```

## 分步运行

```powershell
cd "D:\CADD paper exercise\JNK1_2_3\Docking"
$env:SCHRODINGER = "D:\Schrodinger2025"
$run = "$env:SCHRODINGER\run.exe"

# 1a 导出 maegz
& $run python3 jnk_docking_export\export_complexes_batch.py --config jnk_docking_export\jobs_step1.json --out complexes_mae --format maegz

# 1b Protein Prep
& $run python3 jnk_docking_export\prep_complexes_batch.py --config jnk_docking_export\jobs_step1.json

# 1c SIFt IFP
& $run python3 jnk_docking_export\calc_ifp_batch.py --config jnk_docking_export\jobs_step1.json

# 2a MM-GBSA
& $run python3 jnk_docking_export\run_mmgbsa_batch.py --config jnk_docking_export\jobs_mmgbsa.json

# 2b ΔΔG
& $run python3 jnk_docking_export\calc_ddg_selectivity.py --config jnk_docking_export\jobs_mmgbsa.json
```

## 输出目录

```
Docking/
├── complexes_mae/          # 原始导出（maegz）
├── complexes_prepped/      # Prep 后复合物
├── ifp_results/            # SIFt IFP
│   ├── ifp_all.csv
│   ├── ifp_interactions.csv
│   └── ifp_summary.tsv
└── mmgbsa_results/         # MM-GBSA 作业 + ΔΔG
    ├── mmgbsa_jobs.tsv
    ├── ddg_selectivity.tsv
    └── ddg_selectivity_detail.tsv
```

## 配体范围

**包含（8 个非共价 benchmark）**：AS602801, CC-401, CC-90001, CC-930, E1, Q63, SP600125, TCS JNK 6O  
**排除**：JNK-IN-8（共价）

## 给 Cursor Agent 的命令模板

见下文「复制给 Agent 的 Prompt」。

---

## 复制给 Agent 的 Prompt

### Prompt A — 只写/检查第 1 步脚本

```
请为 JNK 非共价抑制剂选择性项目在 Windows + Schrödinger 2025 上实现第 1 步批处理脚本：

工作目录：D:\CADD paper exercise\JNK1_2_3\Docking
Schrödinger：D:\Schrodinger2025（用 run.exe python3）

任务：
1. 用 jobs_export.json / jobs_step1.json 从 benchmarks_*_pv.maegz 导出复合物到 complexes_mae/（maegz 格式，排除 vsw/top_5000/prime_mmgbsa）
2. prep_complexes_batch.py：Protein Prep 加氢，-fix ligand 保持对接 pose
3. calc_ifp_batch.py：对 complexes_prepped/ 批量算 SIFt IFP，输出 ifp_results/{ifp_all.csv, ifp_interactions.csv, ifp_summary.tsv}
4. 只处理 8 个非共价 benchmark，排除 JNK-IN-8
5. 提供 run_step1.ps1 一键运行

约束：不修改原始 maegz；脚本放 jnk_docking_export/；用 Schrödinger Python API（structutils.interactionfp）
```

### Prompt B — 只写/检查第 2 步脚本

```
请为 JNK 项目实现第 2 步 MM-GBSA 批处理（Windows + Schrödinger 2025）：

输入：complexes_prepped/（第 1 步输出）
配置：jobs_mmgbsa.json

任务：
1. run_mmgbsa_batch.py：对 8 个非共价 benchmark × 5 个 PDB（3ELJ,4L7F,3E7O,3TTI,4WHZ）提交 prime_mmgbsa
   - job_type REAL_MIN, -ligand ligand, -rflexdist 5.0
2. calc_ddg_selectivity.py：解析输出，按亚型平均（JNK1=3ELJ+4L7F, JNK2=3E7O, JNK3=3TTI+4WHZ）
   输出 ddg_selectivity.tsv（含 ddg_jnk1_minus_jnk2 等列）
3. run_step2.ps1 一键运行

排除 JNK-IN-8。脚本放 jnk_docking_export/。
```

### Prompt C — 两步一起（若仓库还没有脚本）

```
在 GitHub 仓库 jnk_docking_export/ 为 JNK1/2/3 非共价选择性筛选写 Windows 批处理工作流：

第 1 步：export → Protein Prep → SIFt IFP
第 2 步：Prime MM-GBSA → ΔΔG(JNK1/2/3) 表

路径：D:\CADD paper exercise\JNK1_2_3\Docking
Schrödinger：D:\Schrodinger2025
8 个非共价 benchmark，排除 JNK-IN-8 和 VSW 文件

请创建：jobs_step1.json, jobs_mmgbsa.json, prep_complexes_batch.py, calc_ifp_batch.py, run_mmgbsa_batch.py, calc_ddg_selectivity.py, run_step1.ps1, run_step2.ps1, WORKFLOW_STEP1_STEP2.md

写完后给出 Windows PowerShell 运行命令。
```

---

## 常见问题

1. **`complexes_prepped` 为空**：先确认 `benchmarks_*_pv.maegz` 存在，且 `jobs_step1.json` 的 `pose_glob` 能匹配到。
2. **IFP 报错 `ligand` ASL**：导出的是否为蛋白+配体合并结构（maegz）；配体需有 `ligand` ASL 识别。
3. **MM-GBSA 很慢**：40 个作业（8×5）可分批；先用 `run_step2.ps1 -DryRun` 检查任务列表。
4. **ΔΔG 表为空**：等 MM-GBSA 全部完成后重跑 `calc_ddg_selectivity.py`；检查 `mmgbsa_results/<PDB>/<ligand>/` 下是否有 `.csv` 或 `.out`。
