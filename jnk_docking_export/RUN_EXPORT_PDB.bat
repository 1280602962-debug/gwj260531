@echo off
setlocal
if not defined SCHRODINGER set "SCHRODINGER=D:\Schrodinger2025"
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0RUN_EXPORT_PDB.ps1" %*
endlocal
