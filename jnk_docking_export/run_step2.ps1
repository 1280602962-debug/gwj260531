param(
    [string]$DockingDir = "D:\CADD paper exercise\JNK1_2_3\Docking",
    [string]$Schrodinger = "D:\Schrodinger2025",
    [string]$WorkflowDir = "",
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

if (-not $WorkflowDir) {
    $WorkflowDir = Join-Path $DockingDir "jnk_docking_export"
}

$env:SCHRODINGER = $Schrodinger
$run = Join-Path $env:SCHRODINGER "run.exe"

if (-not (Test-Path $run)) {
    throw "Schrödinger run.exe not found: $run"
}

Set-Location $DockingDir

if (-not (Test-Path "complexes_prepped")) {
    throw "complexes_prepped/ not found. Run run_step1.ps1 first."
}

Write-Host "=== Step 2a: Prime MM-GBSA (8 non-covalent benchmarks) ===" -ForegroundColor Cyan
$mmgbsaArgs = @(
    "python3", (Join-Path $WorkflowDir "run_mmgbsa_batch.py"),
    "--config", (Join-Path $WorkflowDir "jobs_mmgbsa.json")
)
if ($DryRun) { $mmgbsaArgs += "--dry-run" }
& $run @mmgbsaArgs

Write-Host "=== Step 2b: ΔΔG summary ===" -ForegroundColor Cyan
& $run python3 (Join-Path $WorkflowDir "calc_ddg_selectivity.py") `
    --config (Join-Path $WorkflowDir "jobs_mmgbsa.json")

Write-Host "=== Step 2 complete ===" -ForegroundColor Green
Write-Host "Outputs:"
Write-Host "  mmgbsa_results/mmgbsa_jobs.tsv"
Write-Host "  mmgbsa_results/ddg_selectivity.tsv"
Write-Host "  mmgbsa_results/ddg_selectivity_detail.tsv"
