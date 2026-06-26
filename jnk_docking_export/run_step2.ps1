param(
    [string]$DockingDir = "D:\CADD paper exercise\JNK1_2_3\Docking",
    [string]$Schrodinger = "D:\Schrodinger2025",
    [string]$WorkflowDir = "",
    [switch]$ParseOnly,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

if (-not $WorkflowDir) {
    $WorkflowDir = Join-Path $DockingDir "jnk_docking_export"
}

Set-Location $DockingDir

if ($ParseOnly) {
    & (Join-Path $WorkflowDir "run_step2_parse_only.ps1") -DockingDir $DockingDir -WorkflowDir $WorkflowDir
    exit $LASTEXITCODE
}

$env:SCHRODINGER = $Schrodinger
$run = Join-Path $env:SCHRODINGER "run.exe"

if (-not (Test-Path $run)) {
    throw "Schrödinger run.exe not found: $run"
}

if (-not (Test-Path "complexes_prepped")) {
    Write-Warning "complexes_prepped/ not found. If MM-GBSA is already done, use: .\run_step2.ps1 -ParseOnly"
}

Write-Host "=== Step 2a: Prime MM-GBSA (submit new jobs) ===" -ForegroundColor Cyan
$mmgbsaArgs = @(
    "python3", (Join-Path $WorkflowDir "run_mmgbsa_batch.py"),
    "--config", (Join-Path $WorkflowDir "jobs_mmgbsa.json")
)
if ($DryRun) { $mmgbsaArgs += "--dry-run" }
& $run @mmgbsaArgs

Write-Host "=== Step 2b: scan + ΔΔG summary ===" -ForegroundColor Cyan
python3 (Join-Path $WorkflowDir "scan_mmgbsa_inventory.py") `
    --config (Join-Path $WorkflowDir "jobs_mmgbsa.json")

python3 (Join-Path $WorkflowDir "calc_ddg_selectivity.py") `
    --config (Join-Path $WorkflowDir "jobs_mmgbsa.json") `
    --inventory mmgbsa_results/mmgbsa_inventory.tsv

Write-Host "=== Step 2 complete ===" -ForegroundColor Green
