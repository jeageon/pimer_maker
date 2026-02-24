@echo off
setlocal EnableExtensions

set "SCRIPT_DIR=%~dp0"
set "PROJECT_DIR=%SCRIPT_DIR%.."
set "VENV_DIR=%PROJECT_DIR%\.venv_release"
set "PYTHON_EXE=%VENV_DIR%\Scripts\python.exe"
set "RELEASE_BASE=%PROJECT_DIR%\release"
set "RELEASE_NAME=pimer_maker_release"
set "RELEASE_DIR=%RELEASE_BASE%\%RELEASE_NAME%"
set "ICON_DIR=%PROJECT_DIR%\resources"
set "ICON_FILE=%ICON_DIR%\primer_maker.ico"
set "APP_NAME=PrimerMaker"

if not exist "%PROJECT_DIR%\src\webui.py" (
  echo [ERROR] src\webui.py 를 찾을 수 없습니다.
  pause
  exit /b 1
)

if not exist "%VENV_DIR%" (
  echo [INFO] 릴리즈 전용 가상환경 생성 중: %VENV_DIR%
  py -3 -m venv "%VENV_DIR%"
  if errorlevel 1 (
    echo [ERROR] Python 3(py -3)을 찾을 수 없습니다.
    echo Python 3.10 이상 설치 후 다시 시도하세요.
    pause
    exit /b 1
  )
)

if not exist "%PYTHON_EXE%" (
  echo [ERROR] 가상환경 Python 실행 파일이 없습니다.
  pause
  exit /b 1
)

if not exist "%ICON_DIR%" mkdir "%ICON_DIR%"

"%PYTHON_EXE%" -m pip install --upgrade pip
"%PYTHON_EXE%" -m pip install -r "%PROJECT_DIR%\requirements.txt" pyinstaller
"%PYTHON_EXE%" -m pip install pyinstaller==6.12.0

if not exist "%ICON_FILE%" (
  echo [INFO] 실행 아이콘 생성 중: primer_maker.ico
  powershell -NoProfile -ExecutionPolicy Bypass -Command ^
    "Add-Type -AssemblyName System.Drawing; " ^
    "$shellIcon = [System.Drawing.Icon]::ExtractAssociatedIcon([System.IO.Path]::Combine($env:WINDIR,'system32','imageres.dll')); " ^
    "$stream = New-Object System.IO.FileStream('%ICON_FILE%', [System.IO.FileMode]::Create); " ^
    "$shellIcon.Save($stream); " ^
    "$stream.Close();"
)

if not exist "%ICON_FILE%" (
  echo [WARN] 아이콘 생성 실패. 기본 아이콘으로 진행합니다.
  set "ICON_ARG="
) else (
  set "ICON_ARG=--icon=%ICON_FILE%"
)

cd /d "%PROJECT_DIR%"
rmdir /s /q "%PROJECT_DIR%\build" >nul 2>&1
rmdir /s /q "%PROJECT_DIR%\dist" >nul 2>&1

"%PYTHON_EXE%" -m PyInstaller --noconfirm --onefile --name %APP_NAME% --collect-all streamlit %ICON_ARG% src\webui.py

if not exist "%PROJECT_DIR%\dist\%APP_NAME%.exe" (
  echo [ERROR] exe 생성 실패.
  pause
  exit /b 1
)

if exist "%RELEASE_DIR%" rmdir /s /q "%RELEASE_DIR%"
mkdir "%RELEASE_DIR%"

xcopy /E /I /Y "%PROJECT_DIR%\src" "%RELEASE_DIR%\source\src" >nul
xcopy /E /I /Y "%PROJECT_DIR%\scripts" "%RELEASE_DIR%\source\scripts" >nul
xcopy /E /I /Y "%PROJECT_DIR%\tests" "%RELEASE_DIR%\source\tests" >nul
copy /Y "%PROJECT_DIR%\requirements.txt" "%RELEASE_DIR%\source\requirements.txt" >nul
copy /Y "%PROJECT_DIR%\README.md" "%RELEASE_DIR%\source\README.md" >nul
if exist "%PROJECT_DIR%\data" xcopy /E /I /Y "%PROJECT_DIR%\data" "%RELEASE_DIR%\source\data" >nul

copy /Y "%PROJECT_DIR%\dist\%APP_NAME%.exe" "%RELEASE_DIR%\%APP_NAME%.exe" >nul
if exist "%ICON_FILE%" copy /Y "%ICON_FILE%" "%RELEASE_DIR%\%APP_NAME%.ico" >nul
copy /Y "%PROJECT_DIR%\scripts\install_primer_maker.bat" "%RELEASE_DIR%\install_primer_maker.bat" >nul
copy /Y "%PROJECT_DIR%\scripts\make_primer_desktop_shortcut.bat" "%RELEASE_DIR%\make_primer_desktop_shortcut.bat" >nul

if exist "%RELEASE_BASE%\pimer_maker_release.zip" del "%RELEASE_BASE%\pimer_maker_release.zip"
powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "Compress-Archive -Path '%RELEASE_DIR%\*' -DestinationPath '%RELEASE_BASE%\pimer_maker_release.zip' -Force"

echo [OK] 릴리즈 폴더 생성: %RELEASE_DIR%
echo [OK] 실행 파일: %RELEASE_DIR%\%APP_NAME%.exe
echo [OK] 릴리즈 압축본: %RELEASE_BASE%\pimer_maker_release.zip
pause
endlocal
