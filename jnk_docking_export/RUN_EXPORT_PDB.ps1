# JNK docked complexes -> PDB batch export
# Usage:
#   cd "D:\CADD paper exercise\JNK1_2_3\Docking"
#   Set-ExecutionPolicy -Scope Process Bypass
#   .\RUN_EXPORT_PDB.ps1

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

Write-Step "=== JNK PDB export ==="
Write-Host ""

if (-not (Test-Path -LiteralPath $DockingDir)) {
    throw "Directory not found: $DockingDir"
}

Set-Location -LiteralPath $DockingDir
$env:SCHRODINGER = $SchrodingerHome

$RunExe   = Join-Path $SchrodingerHome "run.exe"
$PyScript = Join-Path $DockingDir "export_complexes_batch.py"
$Config   = Join-Path $DockingDir "jobs_export.json"

if (-not (Test-Path -LiteralPath $RunExe)) {
    throw "Schrodinger run.exe not found: $RunExe"
}

if (-not $SkipDownload) {
    if (-not (Test-Path -LiteralPath $PyScript)) {
        Write-Warn "Downloading export_complexes_batch.py ..."
        Invoke-WebRequest -Uri "$RepoBase/export_complexes_batch.py" -OutFile $PyScript -UseBasicParsing
        Write-Ok "Saved: $PyScript"
    }
    if (-not (Test-Path -LiteralPath $Config)) {
        Write-Warn "Downloading jobs_export.json ..."
        Invoke-WebRequest -Uri "$RepoBase/jobs_export.json" -OutFile $Config -UseBasicParsing
        Write-Ok "Saved: $Config"
    }
}

if (-not (Test-Path -LiteralPath $PyScript)) {
    throw "Missing export_complexes_batch.py in $DockingDir. Run without -SkipDownload or copy from GitHub jnk_docking_export."
}

$pvAll = Get-ChildItem -Path $DockingDir -Recurse -Filter "*_pv.maegz" -File -ErrorAction SilentlyContinue
Write-Host "Work dir : $DockingDir"
Write-Host "Schrodinger: $SchrodingerHome"
Write-Host "PV files : $($pvAll.Count)"
if ($pvAll.Count -gt 0) {
    $pvAll | Select-Object -First 5 | ForEach-Object {
        $rel = $_.FullName.Replace($DockingDir, ".")
        Write-Host "  - $rel"
    }
    if ($pvAll.Count -gt 5) { Write-Host "  - ..." }
}
Write-Host ""

if ($pvAll.Count -eq 0) {
    throw "No *_pv.maegz found under $DockingDir"
}

$pyArgs = @(
    "python3", "export_complexes_batch.py",
    "--out", $OutDir,
    "--format", "pdb"
)

if (Test-Path -LiteralPath $Config) {
    $pyArgs += @("--config", "jobs_export.json")
} else {
    Write-Warn "No jobs_export.json, using --auto"
    $pyArgs += "--auto"
}

if ($AllPoses) { $pyArgs += "--all-poses" }

Write-Step "Run: $RunExe $($pyArgs -join ' ')"
Write-Host ""

& $RunExe @pyArgs
$code = $LASTEXITCODE

Write-Host ""
$outPath = Join-Path $DockingDir $OutDir

if ($code -eq 0 -and (Test-Path -LiteralPath $outPath)) {
    Write-Ok "=== Export OK ==="
    $pdbs = Get-ChildItem -Path $outPath -Recurse -Filter "*.pdb" -File
    Write-Host "PDB count: $($pdbs.Count)"
    foreach ($pdb in @("3ELJ", "4L7F", "3E7O", "3TTI", "4WHZ")) {
        $sub = Join-Path $outPath $pdb
        if (Test-Path -LiteralPath $sub) {
            $n = (Get-ChildItem -Path $sub -Filter "*.pdb" -File).Count
            Write-Host "  $pdb : $n"
        }
    }
    Write-Host "Summary: $(Join-Path $outPath 'export_summary.tsv')"
} elseif ($code -eq 2) {
    Write-Warn "=== Partial errors, see $OutDir\export_errors.log ==="
} else {
    throw "Export failed, exit code $code"
}

exit $code
