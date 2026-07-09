@echo off
setlocal
set DOCKING_DIR=%~dp0..
if "%~1" neq "" set DOCKING_DIR=%~1
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_step2_parse_only.ps1" -DockingDir "%DOCKING_DIR%"
endlocal
