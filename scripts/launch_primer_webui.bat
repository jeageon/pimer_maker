@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
set "PROJECT_DIR=%SCRIPT_DIR%.."
call "%PROJECT_DIR%run-PimerMaker-WebUI.bat"

endlocal
