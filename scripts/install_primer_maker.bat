@echo off
setlocal EnableExtensions

set "APP_NAME=Primer Maker"
set "APP_EXE=PrimerMaker.exe"
set "APP_ICON=PrimerMaker.ico"

set "SCRIPT_DIR=%~dp0"
if "%SCRIPT_DIR:~-1%"=="\\" set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
set "SOURCE_EXE=%SCRIPT_DIR%\%APP_EXE%"
set "SOURCE_ICON=%SCRIPT_DIR%\%APP_ICON%"

if "%1"=="--help" (
    echo Usage:
    echo   %~nx0 [install_directory]
    echo.
    echo Example:
    echo   %~nx0
    echo   %~nx0 "%%USERPROFILE%%\\AppData\\Local\\PrimerMaker"
    pause
    exit /b 0
)

if "%~1"=="" (
    set "INSTALL_DIR=%ProgramFiles%\%APP_NAME%"
) else (
    set "INSTALL_DIR=%~1"
)

set "INSTALL_ROOT=%INSTALL_DIR%"
call :ensure_dir "%INSTALL_DIR%"
if errorlevel 1 (
    if not defined LOCALAPPDATA set "LOCALAPPDATA=%USERPROFILE%\\AppData\\Local"
    echo [WARN] cannot create default path: %INSTALL_DIR%
    set "INSTALL_DIR=%LOCALAPPDATA%\\Programs\\%APP_NAME%"
    echo [INFO] fallback install path: %INSTALL_DIR%
    call :ensure_dir "%INSTALL_DIR%"
    if errorlevel 1 (
        echo [ERROR] cannot prepare install directory.
        echo Install path used: %INSTALL_DIR%
        pause
        exit /b 1
    )
)

if not exist "%SOURCE_EXE%" (
    echo [ERROR] executable not found: %SOURCE_EXE%
    pause
    exit /b 1
)

copy /Y "%SOURCE_EXE%" "%INSTALL_DIR%\%APP_EXE%" >nul
if exist "%SOURCE_ICON%" (
    copy /Y "%SOURCE_ICON%" "%INSTALL_DIR%\%APP_ICON%" >nul
)

set "INSTALL_EXE=%INSTALL_DIR%\%APP_EXE%"
set "INSTALL_ICON=%INSTALL_DIR%\%APP_ICON%"
if exist "%INSTALL_ICON%" (
    set "HAS_ICON=1"
) else (
    set "HAS_ICON="
)

set "SHORTCUT_SCRIPT=%TEMP%\\PrimerMakerShortcut.vbs"
set "UNINSTALLER=%INSTALL_DIR%\\uninstall_primer_maker.bat"

echo Set WshShell = CreateObject("WScript.Shell") >"%SHORTCUT_SCRIPT%"
echo DesktopPath = WshShell.SpecialFolders("Desktop") >>"%SHORTCUT_SCRIPT%"
echo StartMenuBase = WshShell.SpecialFolders("StartMenu") >>"%SHORTCUT_SCRIPT%"
echo StartMenuPath = StartMenuBase & "\\%APP_NAME%" >>"%SHORTCUT_SCRIPT%"
echo Set FSO = CreateObject("Scripting.FileSystemObject") >>"%SHORTCUT_SCRIPT%"
echo If Not FSO.FolderExists(StartMenuPath) Then >>"%SHORTCUT_SCRIPT%"
echo     FSO.CreateFolder(StartMenuPath) >>"%SHORTCUT_SCRIPT%"
echo End If >>"%SHORTCUT_SCRIPT%"

echo Set DesktopLink = WshShell.CreateShortcut(DesktopPath & "\\%APP_NAME%.lnk") >>"%SHORTCUT_SCRIPT%"
echo DesktopLink.TargetPath = "%INSTALL_EXE%" >>"%SHORTCUT_SCRIPT%"
echo DesktopLink.WorkingDirectory = "%INSTALL_DIR%" >>"%SHORTCUT_SCRIPT%"
echo DesktopLink.Description = "%APP_NAME%" >>"%SHORTCUT_SCRIPT%"
if defined HAS_ICON echo DesktopLink.IconLocation = "%INSTALL_ICON%,0" >>"%SHORTCUT_SCRIPT%"
echo DesktopLink.Save >>"%SHORTCUT_SCRIPT%"

echo Set MenuLink = WshShell.CreateShortcut(StartMenuPath & "\\%APP_NAME%.lnk") >>"%SHORTCUT_SCRIPT%"
echo MenuLink.TargetPath = "%INSTALL_EXE%" >>"%SHORTCUT_SCRIPT%"
echo MenuLink.WorkingDirectory = "%INSTALL_DIR%" >>"%SHORTCUT_SCRIPT%"
echo MenuLink.Description = "%APP_NAME%" >>"%SHORTCUT_SCRIPT%"
if defined HAS_ICON echo MenuLink.IconLocation = "%INSTALL_ICON%,0" >>"%SHORTCUT_SCRIPT%"
echo MenuLink.Save >>"%SHORTCUT_SCRIPT%"

cscript //nologo "%SHORTCUT_SCRIPT%"
del "%SHORTCUT_SCRIPT%" >nul 2>&1

(
    echo @echo off
    echo setlocal EnableExtensions
    echo set "INSTALL_DIR=%INSTALL_DIR%"
    echo set "APP_NAME=%APP_NAME%"
    echo set "APP_EXE=%APP_EXE%"
    echo set "START_MENU=%%APPDATA%%\\Microsoft\\Windows\\Start Menu\\Programs\\%%APP_NAME%%"
    echo set "DESKTOP=%%USERPROFILE%%\\Desktop"
    echo del /q "%%DESKTOP%%\\%%APP_NAME%%.lnk" 2^>nul
    echo del /q "%%START_MENU%%\\%%APP_NAME%%.lnk" 2^>nul
    echo if exist "%%START_MENU%%" rmdir /s /q "%%START_MENU%%" 2^>nul
    echo if exist "%%INSTALL_DIR%%\\%%APP_EXE%%" del /q "%%INSTALL_DIR%%\\%%APP_EXE%%"
    echo if exist "%%INSTALL_DIR%%\\PrimerMaker.ico" del /q "%%INSTALL_DIR%%\\PrimerMaker.ico"
    echo if exist "%%INSTALL_DIR%%" rmdir /q "%%INSTALL_DIR%%" 2^>nul
    echo echo Primer Maker has been removed.
    echo pause
) > "%UNINSTALLER%"

echo [OK] installation complete.
echo Install path: %INSTALL_DIR%
if not "%INSTALL_ROOT%"=="%INSTALL_DIR%" (
    echo [INFO] fallback install path was used.
)
echo Desktop and Start menu shortcuts were created.
echo Run: "%%USERPROFILE%%\\Desktop\\%APP_NAME%.lnk"
echo Uninstall: %UNINSTALLER%
pause
exit /b 0

:ensure_dir
setlocal
set "DIR=%~1"
if exist "%DIR%" exit /b 0
md "%DIR%" >nul 2>&1
if not exist "%DIR%" exit /b 1
exit /b 0
