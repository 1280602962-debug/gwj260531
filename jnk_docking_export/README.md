# JNK 复合物 PDB 批量导出

## 文件来源

这些脚本由 Cursor Agent 生成，位于仓库 `jnk_docking_export/`。
您之前 IFP 用的 `calc_ifp_batch.py` 是另一套脚本；**导出 PDB 需要本目录下的文件**。

## 最快用法（Windows PowerShell）

1. 将整个 `jnk_docking_export` 文件夹复制到：

   ```
   D:\CADD paper exercise\JNK1_2_3\Docking\
   ```

2. 在 PowerShell 执行：

   ```powershell
   cd "D:\CADD paper exercise\JNK1_2_3\Docking\jnk_docking_export"
   Set-ExecutionPolicy -Scope Process Bypass
   .\RUN_EXPORT_PDB.ps1 -DockingDir "D:\CADD paper exercise\JNK1_2_3\Docking"
   ```

   或双击 `RUN_EXPORT_PDB.bat`（需把 bat 放在 Docking 根目录并改路径）。

3. 若只有 `RUN_EXPORT_PDB.ps1` 在 Docking 根目录，脚本会**自动下载**缺失的 py/json。

## 手动运行（不下载）

```powershell
cd "D:\CADD paper exercise\JNK1_2_3\Docking"
$env:SCHRODINGER = "D:\Schrodinger2025"
& "$env:SCHRODINGER\run.exe" python3 export_complexes_batch.py --config jobs_export.json --out complexes --format pdb
```

无 jobs.json 时自动发现 pv：

```powershell
& "$env:SCHRODINGER\run.exe" python3 export_complexes_batch.py --auto --out complexes --format pdb
```

## 输出

- `complexes/{PDB}/*.pdb` — 蛋白+配体复合物（默认每配体最佳 pose）
- `complexes/export_summary.tsv` — 汇总表

## 修复说明（相对初版）

- **递归搜索** `**/*_pv.maegz`，支持 `_XP_1\` 等子目录
- **`--auto`** 模式：无需 jobs.json，按文件名中的 PDB ID 分组
- **RUN_EXPORT_PDB.ps1**：缺文件时从 GitHub 自动下载
