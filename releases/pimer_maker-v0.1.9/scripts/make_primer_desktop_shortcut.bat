@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
set "PROJECT_DIR=%SCRIPT_DIR%.."
set "TARGET=%PROJECT_DIR%\run-PimerMaker-WebUI.bat"
set "SHORTCUT_NAME=Primer Maker.lnk"
set "DESKTOP=%USERPROFILE%\Desktop"

if not exist "%TARGET%" (
  echo [ERROR] 실행 파일을 찾을 수 없습니다.
  pause
  exit /b 1
)

powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "$WshShell = New-Object -ComObject WScript.Shell; " ^
  "$Shortcut = $WshShell.CreateShortcut('%DESKTOP%\\%SHORTCUT_NAME%'); " ^
  "$Shortcut.TargetPath = '%TARGET%'; " ^
  "$Shortcut.WorkingDirectory = (Resolve-Path '%PROJECT_DIR%').Path; " ^
  "$Shortcut.WindowStyle = 1; " ^
  "$Shortcut.Description = 'Primer Maker (Streamlit)'; " ^
  "$Shortcut.Save()"

if not exist "%DESKTOP%\%SHORTCUT_NAME%" (
  echo [ERROR] 바로가기가 생성되지 않았습니다.
  pause
  exit /b 1
)

echo [OK] 바탕화면에 "%DESKTOP%\%SHORTCUT_NAME%" 생성 완료
pause

endlocal
