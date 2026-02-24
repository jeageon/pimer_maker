@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
set "PROJECT_DIR=%SCRIPT_DIR%.."
set "RELEASE_DIR=%PROJECT_DIR%\release\pimer_maker_release"
set "APP_EXE=%RELEASE_DIR%\PrimerMaker.exe"
set "ICON_FILE=%RELEASE_DIR%\PrimerMaker.ico"
set "DESKTOP=%USERPROFILE%\Desktop"
set "SHORTCUT_NAME=Primer Maker.lnk"

if not exist "%APP_EXE%" (
  echo [ERROR] 실행파일을 찾을 수 없습니다. 먼저 prepare_release_windows.bat 를 먼저 실행하세요.
  echo 경로: %APP_EXE%
  pause
  exit /b 1
)

powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "$WshShell = New-Object -ComObject WScript.Shell; " ^
  "$Shortcut = $WshShell.CreateShortcut('%DESKTOP%\\%SHORTCUT_NAME%'); " ^
  "$Shortcut.TargetPath = '%APP_EXE%'; " ^
  "$Shortcut.WorkingDirectory = (Resolve-Path '%RELEASE_DIR%').Path; " ^
  "$Shortcut.WindowStyle = 1; " ^
  "$Shortcut.Description = 'Primer Maker'; " ^
  "$Shortcut.Hotkey = 'Ctrl+Alt+P'; " ^
  "if (Test-Path '%ICON_FILE%') { $Shortcut.IconLocation = '%ICON_FILE%'; } " ^
  "$Shortcut.Save()"

echo [OK] 바탕화면 바로가기가 생성되었습니다.
pause
endlocal
