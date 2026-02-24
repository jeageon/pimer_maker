@echo off
setlocal EnableExtensions

set "SCRIPT_DIR=%~dp0"
set "PROJECT_DIR=%SCRIPT_DIR%.."
set "VENV_DIR=%PROJECT_DIR%\.venv"
set "PYTHON_EXE=%VENV_DIR%\Scripts\python.exe"

if not exist "%PROJECT_DIR%\src\webui.py" (
  echo [ERROR] src\webui.py 를 찾을 수 없습니다.
  echo [%PROJECT_DIR%]
  pause
  exit /b 1
)

if not exist "%VENV_DIR%" (
  echo [INFO] 가상환경이 없어서 생성합니다: %VENV_DIR%
  py -3 -m venv "%VENV_DIR%"
  if errorlevel 1 (
    echo [ERROR] Python 3 실행 파일(py -3)을 찾을 수 없습니다.
    echo Python 3.10 이상을 먼저 설치해 주세요.
    pause
    exit /b 1
  )
)

if not exist "%PYTHON_EXE%" (
  echo [ERROR] 가상환경 python 실행 파일을 찾을 수 없습니다.
  pause
  exit /b 1
)

"%PYTHON_EXE%" -m pip install --upgrade pip
"%PYTHON_EXE%" -m pip install -r "%PROJECT_DIR%\requirements.txt"
set "STREAMLIT_BROWSER_GATHER_USAGE_STATS=false"
cd /d "%PROJECT_DIR%"
"%PYTHON_EXE%" -m streamlit run src\webui.py

endlocal
