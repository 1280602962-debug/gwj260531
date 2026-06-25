# =============================================================================
# JNK 对接复合物 → PDB 一键导出
#
# 功能：
#   1. 若缺少脚本，自动从 GitHub 下载 export_complexes_batch.py / jobs_export.json
#   2. 递归搜索子目录中的 *_pv.maegz（含 _XP_1\ 等）
#   3. 导出蛋白+配体合并 PDB 到 complexes\
#
# 用法（PowerShell）：
#   cd "D:\CADD paper exercise\JNK1_2_3\Docking"
#   Set-ExecutionPolicy -Scope Process Bypass
#   .\RUN_EXPORT_PDB.ps1
# =============================================================================

param(
    [string]$DockingDir = "D:\CADD paper exercise\JNK1_2_3\Docking",
    [string]$SchrodingerHome = "D:\Schrodinger2025",
    [string]$OutDir = "complexes",
    [switch]$AllPoses,
    [switch]$SkipDownload
)

$ErrorActionPreference = "Stop"

$RepoBase = "https://raw.githubusercontent.com/1280602962-debug/gwj260531/cursor-agent-nlrp3-report-89b1/jnk_docking_export"

function Write-Step($msg) { Write-Host $msg -ForegroundColor Cyan }
function Write-Ok($msg)   { Write-Host $msg -ForegroundColor Green }
function Write-Warn($msg) { Write-Host $msg -ForegroundColor Yellow }

Write-Step "=== JNK 复合物 PDB 批量导出 ==="
Write-Host ""

if (-not (Test-Path -LiteralPath $DockingDir)) {
    throw "目录不存在: $DockingDir"
}

Set-Location -LiteralPath $DockingDir
$env:SCHRODINGER = $SchrodingerHome

$RunExe   = Join-Path $SchrodingerHome "run.exe"
$PyScript = Join-Path $DockingDir "export_complexes_batch.py"
$Config   = Join-Path $DockingDir "jobs_export.json"

if (-not (Test-Path -LiteralPath $RunExe)) {
    throw "找不到薛定谔: $RunExe`n请修改 -SchrodingerHome 参数"
}

# --- 自动补全缺失脚本 ---
if (-not $SkipDownload) {
    if (-not (Test-Path -LiteralPath $PyScript)) {
        Write-Warn "缺少 export_complexes_batch.py，正在下载..."
        Invoke-WebRequest -Uri "$RepoBase/export_complexes_batch.py" -OutFile $PyScript -UseBasicParsing
        Write-Ok "已保存: $PyScript"
    }
    if (-not (Test-Path -LiteralPath $Config)) {
        Write-Warn "缺少 jobs_export.json，正在下载..."
        Invoke-WebRequest -Uri "$RepoBase/jobs_export.json" -OutFile $Config -UseBasicParsing
        Write-Ok "已保存: $Config"
    }
}

if (-not (Test-Path -LiteralPath $PyScript)) {
    throw @"
仍缺少 export_complexes_batch.py。
请手动复制 jnk_docking_export\export_complexes_batch.py 到:
  $DockingDir
或检查网络后重试（脚本会从 GitHub 自动下载）。
"@
}

# --- 查找 pv 文件 ---
$pvAll = Get-ChildItem -Path $DockingDir -Recurse -Filter "*_pv.maegz" -File -ErrorAction SilentlyContinue
Write-Host "工作目录 : $DockingDir"
Write-Host "薛定谔   : $SchrodingerHome"
Write-Host "pv 文件  : $($pvAll.Count) 个（含子目录）"
if ($pvAll.Count -gt 0) {
    $pvAll | Select-Object -First 5 | ForEach-Object { Write-Host "  - $($_.FullName.Replace($DockingDir, '.'))" }
    if ($pvAll.Count -gt 5) { Write-Host "  - ..." }
}
Write-Host ""

if ($pvAll.Count -eq 0) {
    throw "未找到任何 *_pv.maegz，请确认对接结果在此目录或子目录中。"
}

# --- 运行导出 ---
$pyArgs = @(
    "python3", "export_complexes_batch.py",
    "--out", $OutDir,
    "--format", "pdb"
)

if (Test-Path -LiteralPath $Config) {
    $pyArgs += @("--config", "jobs_export.json")
} else {
    Write-Warn "无 jobs_export.json，使用 --auto 自动发现 pv 文件"
    $pyArgs += "--auto"
}

if ($AllPoses) { $pyArgs += "--all-poses" }

Write-Step "执行: & `"$RunExe`" $($pyArgs -join ' ')"
Write-Host ""

& $RunExe @pyArgs
$code = $LASTEXITCODE

Write-Host ""
$outPath = Join-Path $DockingDir $OutDir

if ($code -eq 0 -and (Test-Path -LiteralPath $outPath)) {
    Write-Ok "=== 导出成功 ==="
    $pdbs = Get-ChildItem -Path $outPath -Recurse -Filter "*.pdb" -File
    Write-Host "PDB 总数: $($pdbs.Count)"
    foreach ($pdb in @("3ELJ","4L7F","3E7O","3TTI","4WHZ")) {
        $sub = Join-Path $outPath $pdb
        if (Test-Path -LiteralPath $sub) {
            $n = (Get-ChildItem -Path $sub -Filter "*.pdb" -File).Count
            Write-Host "  $pdb : $n 个 PDB"
        }
    }
    Write-Host ""
    Write-Host "汇总表: $(Join-Path $outPath 'export_summary.tsv')"
    Write-Host "WSL路径: /mnt/d/CADD paper exercise/JNK1_2_3/Docking/$OutDir"
} elseif ($code -eq 2) {
    Write-Warn "=== 部分 PDB 导出失败，见 $OutDir\export_errors.log ==="
} else {
    throw "=== 导出失败 (exit $code) ==="
}

exit $code
