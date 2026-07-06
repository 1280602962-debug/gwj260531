# =============================================================================
# JNK 对接复合物批量导出 — Windows PowerShell 完整脚本
# 从 Glide *_pv.maegz 导出蛋白+配体合并 PDB（5 PDB × 10 配体）
#
# 用法（在 PowerShell 中）：
#   cd "D:\CADD paper exercise\JNK1_2_3\Docking"
#   Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
#   .\export_jnk_complexes.ps1
#
# 导出全部 pose：
#   .\export_jnk_complexes.ps1 -AllPoses
# =============================================================================

param(
    [string]$DockingDir = "D:\CADD paper exercise\JNK1_2_3\Docking",
    [string]$SchrodingerHome = "D:\Schrodinger2025",
    [string]$OutDir = "complexes",
    [ValidateSet("pdb", "mae", "maegz")]
    [string]$Format = "pdb",
    [switch]$AllPoses
)

$ErrorActionPreference = "Stop"

# -----------------------------------------------------------------------------
# 路径检查
# -----------------------------------------------------------------------------
Write-Host "=== JNK 复合物批量导出 ===" -ForegroundColor Cyan
Write-Host ""

if (-not (Test-Path -LiteralPath $DockingDir)) {
    Write-Error "对接目录不存在: $DockingDir"
}

$RunExe = Join-Path $SchrodingerHome "run.exe"
if (-not (Test-Path -LiteralPath $RunExe)) {
    Write-Error "找不到薛定谔 run.exe: $RunExe`n请修改 `$SchrodingerHome 变量"
}

Set-Location -LiteralPath $DockingDir
$env:SCHRODINGER = $SchrodingerHome

$PyScript = Join-Path $DockingDir "export_complexes_batch.py"
$Config   = Join-Path $DockingDir "jobs_export.json"

if (-not (Test-Path -LiteralPath $PyScript)) {
    Write-Error "缺少 export_complexes_batch.py，请复制到: $DockingDir"
}
if (-not (Test-Path -LiteralPath $Config)) {
    Write-Error "缺少 jobs_export.json，请复制到: $DockingDir"
}

# -----------------------------------------------------------------------------
# 检查 pv 文件
# -----------------------------------------------------------------------------
Write-Host "工作目录   : $DockingDir"
Write-Host "薛定谔路径 : $SchrodingerHome"
Write-Host "输出目录   : $OutDir"
Write-Host "格式       : $Format"
Write-Host "仅最佳pose : $(-not $AllPoses)"
Write-Host ""

$pvFiles = Get-ChildItem -Path $DockingDir -Filter "*_pv.maegz" -File
Write-Host "发现 pv 文件 : $($pvFiles.Count) 个"
if ($pvFiles.Count -eq 0) {
    Write-Warning "当前目录下没有 *_pv.maegz，请确认路径正确"
}

# -----------------------------------------------------------------------------
# 构建命令
# -----------------------------------------------------------------------------
$argsList = @(
    "python3",
    "export_complexes_batch.py",
    "--config", "jobs_export.json",
    "--out", $OutDir,
    "--format", $Format
)
if ($AllPoses) {
    $argsList += "--all-poses"
}

Write-Host "执行命令:" -ForegroundColor Yellow
Write-Host "  & `"$RunExe`" $($argsList -join ' ')"
Write-Host ""

# -----------------------------------------------------------------------------
# 运行
# -----------------------------------------------------------------------------
& $RunExe @argsList
$exitCode = $LASTEXITCODE

Write-Host ""
if ($exitCode -eq 0) {
    Write-Host "=== 导出成功 ===" -ForegroundColor Green
} elseif ($exitCode -eq 2) {
    Write-Host "=== 导出完成（有部分错误，见 export_errors.log）===" -ForegroundColor Yellow
} else {
    Write-Host "=== 导出失败 (exit $exitCode) ===" -ForegroundColor Red
}

# -----------------------------------------------------------------------------
# 结果统计
# -----------------------------------------------------------------------------
$outPath = Join-Path $DockingDir $OutDir
if (Test-Path -LiteralPath $outPath) {
    $pdbCount = (Get-ChildItem -Path $outPath -Recurse -Include "*.pdb","*.mae","*.maegz" -File).Count
    Write-Host "导出文件数 : $pdbCount"
    Write-Host "汇总表     : $(Join-Path $outPath 'export_summary.tsv')"

    foreach ($pdb in @("3ELJ", "4L7F", "3E7O", "3TTI", "4WHZ")) {
        $sub = Join-Path $outPath $pdb
        if (Test-Path -LiteralPath $sub) {
            $n = (Get-ChildItem -Path $sub -File).Count
            Write-Host "  $pdb : $n 个文件"
        }
    }

    $summary = Join-Path $outPath "export_summary.tsv"
    if (Test-Path -LiteralPath $summary) {
        Write-Host ""
        Write-Host "前 5 行汇总:" -ForegroundColor Cyan
        Get-Content -LiteralPath $summary -TotalCount 6
    }
}

Write-Host ""
Write-Host "WSL 访问路径: /mnt/d/CADD paper exercise/JNK1_2_3/Docking/$OutDir" -ForegroundColor DarkGray

exit $exitCode
