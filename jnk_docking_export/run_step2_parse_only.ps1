param(
    [string]$DockingDir = "D:\CADD paper exercise\JNK1_2_3\Docking",
    [string]$WorkflowDir = ""
)

$ErrorActionPreference = "Stop"

if (-not $WorkflowDir) {
    $WorkflowDir = Join-Path $DockingDir "jnk_docking_export"
}

Set-Location $DockingDir

Write-Host "=== Step 2 (parse only): scan existing MM-GBSA ===" -ForegroundColor Cyan

# Inventory scan — plain Python is enough (no Schrödinger API needed)
python3 (Join-Path $WorkflowDir "scan_mmgbsa_inventory.py") `
    --config (Join-Path $WorkflowDir "jobs_mmgbsa.json") `
    --out mmgbsa_results/mmgbsa_inventory.tsv

Write-Host "=== Step 2 (parse only): ΔΔG summary ===" -ForegroundColor Cyan
python3 (Join-Path $WorkflowDir "calc_ddg_selectivity.py") `
    --config (Join-Path $WorkflowDir "jobs_mmgbsa.json") `
    --inventory mmgbsa_results/mmgbsa_inventory.tsv

Write-Host "=== Step 2 complete (existing MM-GBSA parsed) ===" -ForegroundColor Green
Write-Host "Outputs:"
Write-Host "  mmgbsa_results/mmgbsa_inventory.tsv"
Write-Host "  mmgbsa_results/ddg_selectivity.tsv"
Write-Host "  mmgbsa_results/ddg_selectivity_detail.tsv"
