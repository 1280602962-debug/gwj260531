@echo off
setlocal EnableExtensions

if not defined SCHRODINGER (
  set "SCHRODINGER=D:\Schrodinger2025"
)

set "SCRIPT_DIR=%~dp0"
cd /d "%SCRIPT_DIR%"

if "%~1"=="--all-poses" (
  "%SCHRODINGER%\run.exe" python3 export_complexes_batch.py --config jobs_export.json --out complexes --all-poses
) else (
  "%SCHRODINGER%\run.exe" python3 export_complexes_batch.py --config jobs_export.json --out complexes %*
)

endlocal
