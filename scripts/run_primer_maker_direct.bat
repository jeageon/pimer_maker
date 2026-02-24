@echo off
setlocal EnableExtensions EnableDelayedExpansion
chcp 65001 >nul

set "SCRIPT_DIR=%~dp0"
if "%SCRIPT_DIR:~-1%"=="\\" set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
set "EXE_PATH=%SCRIPT_DIR%\PrimerMaker.exe"
if not exist "%EXE_PATH%" set "EXE_PATH=%SCRIPT_DIR%\pimer_maker_release\PrimerMaker.exe"

set "SOURCE_DIR="
if exist "%SCRIPT_DIR%\source\src\webui.py" (
    set "SOURCE_DIR=%SCRIPT_DIR%\source"
) else if exist "%SCRIPT_DIR%\..\source\src\webui.py" (
    set "SOURCE_DIR=%SCRIPT_DIR%\..\source"
) else if exist "%SCRIPT_DIR%\..\src\webui.py" (
    set "SOURCE_DIR=%SCRIPT_DIR%\.."
) else if exist "%SCRIPT_DIR%\src\webui.py" (
    set "SOURCE_DIR=%SCRIPT_DIR%"
)

if defined EXE_PATH (
    if exist "%EXE_PATH%" (
        echo [INFO] 실행 파일로 시작합니다: "%EXE_PATH%"
        "%EXE_PATH%"
        if not errorlevel 1 (
            echo [OK] 실행 파일이 시작되었습니다. 브라우저가 자동으로 뜨지 않으면
            echo     http://localhost:8501 으로 접속하세요.
            exit /b 0
        )
        echo [WARN] 실행 파일 시작 실패, Python 직접 실행 모드로 전환합니다.
    )
)

if "%SOURCE_DIR%"=="" (
    echo [ERROR] 실행에 필요한 소스 경로(src/webui.py)를 찾을 수 없습니다.
    echo 현재 경로: "%SCRIPT_DIR%"
    pause
    exit /b 1
)

set "VENV_DIR=%SOURCE_DIR%\.venv_portable"
set "PYTHON_EXE=%VENV_DIR%\Scripts\python.exe"

if exist "%VENV_DIR%" if not exist "%PYTHON_EXE%" rmdir /s /q "%VENV_DIR%" >nul 2>&1

if not exist "%PYTHON_EXE%" (
    set "PY_CMD="
    where py >nul 2>&1 && set "PY_CMD=py -3" || where python >nul 2>&1 && set "PY_CMD=python"
    if not defined PY_CMD (
        echo [ERROR] Python 3를 찾을 수 없습니다.
        echo py.exe 또는 python.exe가 PATH에 등록되어 있어야 합니다.
        pause
        exit /b 1
    )

    echo [INFO] 로컬 실행 전용 가상환경 생성: "%VENV_DIR%"
    %PY_CMD% -m venv "%VENV_DIR%"
    if errorlevel 1 (
        echo [ERROR] 가상환경 생성 실패.
        pause
        exit /b 1
    )
)

echo [INFO] pip 점검 및 의존성 설치 중...
if not exist "%SOURCE_DIR%\requirements.txt" (
    echo [ERROR] requirements.txt를 찾을 수 없습니다: "%SOURCE_DIR%\requirements.txt"
    pause
    exit /b 1
)

"%PYTHON_EXE%" -m pip install --upgrade pip -q
"%PYTHON_EXE%" -m pip install -r "%SOURCE_DIR%\requirements.txt" -q
if errorlevel 1 (
    echo [ERROR] 패키지 설치에 실패했습니다.
    echo "%SOURCE_DIR%\requirements.txt" 또는 네트워크 연결을 확인해 주세요.
    pause
    exit /b 1
)

echo [INFO] 브라우저에서 http://localhost:8501 로 접속하세요.
set "STREAMLIT_BROWSER_GATHER_USAGE_STATS=false"
cd /d "%SOURCE_DIR%"
"%PYTHON_EXE%" -m streamlit run src\webui.py
if errorlevel 1 (
    echo [ERROR] Streamlit 실행 실패.
    pause
    exit /b 1
)

exit /b 0
