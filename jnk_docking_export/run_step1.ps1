param(
    [string]$DockingDir = "D:\CADD paper exercise\JNK1_2_3\Docking",
    [string]$Schrodinger = "D:\Schrodinger2025",
    [string]$WorkflowDir = ""
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

Write-Host "=== Step 1a: Export benchmark complexes (maegz, exclude VSW) ===" -ForegroundColor Cyan
& $run python3 (Join-Path $WorkflowDir "export_complexes_batch.py") `
    --config (Join-Path $WorkflowDir "jobs_step1.json") `
    --out complexes_mae `
    --format maegz

Write-Host "=== Step 1b: Protein Prep (fix ligand, add H) ===" -ForegroundColor Cyan
& $run python3 (Join-Path $WorkflowDir "prep_complexes_batch.py") `
    --config (Join-Path $WorkflowDir "jobs_step1.json") `
    --in-dir complexes_mae `
    --out-dir complexes_prepped

Write-Host "=== Step 1c: SIFt IFP ===" -ForegroundColor Cyan
& $run python3 (Join-Path $WorkflowDir "calc_ifp_batch.py") `
    --config (Join-Path $WorkflowDir "jobs_step1.json") `
    --in-dir complexes_prepped `
    --out-dir ifp_results

Write-Host "=== Step 1 complete ===" -ForegroundColor Green
Write-Host "Outputs:"
Write-Host "  complexes_mae/"
Write-Host "  complexes_prepped/"
Write-Host "  ifp_results/ifp_all.csv"
Write-Host "  ifp_results/ifp_interactions.csv"
Write-Host "  ifp_results/ifp_summary.tsv"
